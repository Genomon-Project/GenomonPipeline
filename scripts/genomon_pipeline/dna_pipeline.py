import os
import shutil
import yaml
import linecache
from ruffus import *
from genomon_pipeline.config.run_conf import *
from genomon_pipeline.config.genomon_conf import *
from genomon_pipeline.config.task_conf import *
from genomon_pipeline.config.sample_conf import *
from genomon_pipeline.dna_resource.bamtofastq import *
from genomon_pipeline.dna_resource.fastq_splitter import *
from genomon_pipeline.dna_resource.bwa_align import *
from genomon_pipeline.dna_resource.markduplicates import *
from genomon_pipeline.dna_resource.mutation_call import *
from genomon_pipeline.dna_resource.mutation_merge import *
from genomon_pipeline.dna_resource.sv_parse import *
from genomon_pipeline.dna_resource.sv_merge import *
from genomon_pipeline.dna_resource.sv_filt import *

# set task classes
bamtofastq = Bam2Fastq(task_conf.get("bam2fastq", "qsub_option"), run_conf.project_root + '/script')
fastq_splitter = Fastq_splitter(task_conf.get("split_fast", "qsub_option"), run_conf.project_root + '/script')
bwa_align = Bwa_align(task_conf.get("bwa_mem", "qsub_option"), run_conf.project_root + '/script')
markduplicates = Markduplicates(task_conf.get("markduplicates", "qsub_option"), run_conf.project_root + '/script')
mutation_call = Mutation_call(task_conf.get("mutation_call", "qsub_option"), run_conf.project_root + '/script')
mutation_merge = Mutation_merge(task_conf.get("mutation_merge", "qsub_option"), run_conf.project_root + '/script')
sv_parse = SV_parse(task_conf.get("sv_parse", "qsub_option"), run_conf.project_root + '/script')
sv_merge = SV_merge(task_conf.get("sv_merge", "qsub_option"), run_conf.project_root + '/script')
sv_filt = SV_filt(task_conf.get("sv_filt", "qsub_option"), run_conf.project_root + '/script')

# generate output list of 'linked fastq'
linked_fastq_list = []
for sample in sample_conf.fastq:
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/1.sorted.bam'): continue
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/' + sample + '.markdup.bam'): continue
    fastq_prefix, ext = os.path.splitext(sample_conf.fastq[sample][0][0])
    linked_fastq_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1' + ext,
                              run_conf.project_root + '/fastq/' + sample + '/1_2' + ext])

# generate output list of 'bam2fastq'
bam2fastq_output_list = []
for sample in sample_conf.bam_tofastq:
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/1.sorted.bam'): continue
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/' + sample + '.markdup.bam'): continue
    bam2fastq_output_list.append([run_conf.project_root + '/fastq/' + sample + '/1_1.fastq',
                                  run_conf.project_root + '/fastq/' + sample + '/1_2.fastq'])

# generate input list of 'mutation call'
markdup_bam_list = []
merge_mutation_list = []
for complist in sample_conf.mutation_call:
    if os.path.exists(run_conf.project_root + '/mutation/' + complist[0] + '/' + complist[0] + '_genomon_mutations.result.txt'): continue
    markdup_bam_list.append([run_conf.project_root + '/bam/' + complist[0] + '/' + complist[0] + '.markdup.bam',
                             run_conf.project_root + '/bam/' + complist[1] + '/' + complist[1] + '.markdup.bam',
                             run_conf.project_root + '/mutation/control_panel/' + complist[2] + ".control_panel.txt"])

# generate input list of 'SV parse'
parse_sv_bam_list = []
all_target_bams = []
unique_bams = []
for complist in sample_conf.sv_detection:
    tumor_sample = complist[0]
    if tumor_sample != None:
        all_target_bams.append(run_conf.project_root + '/bam/' + tumor_sample + '/' + tumor_sample + '.markdup.bam')
    normal_sample = complist[1]
    if normal_sample != None:
        all_target_bams.append(run_conf.project_root + '/bam/' + normal_sample + '/' + normal_sample + '.markdup.bam')
    panel_name = complist[2]
    if panel_name != None:
        for panel_sample in sample_conf.control_panel[panel_name]:
            all_target_bams.append(run_conf.project_root + '/bam/' + panel_sample + '/' + panel_sample + '.markdup.bam')
    unique_bams = list(set(all_target_bams))       
    
for bam in unique_bams:
    dir_name = os.path.dirname(bam)
    sample_name = os.path.basename(dir_name)
    if os.path.exists(run_conf.project_root + '/sv/' + sample_name + '/' + sample_name + '.junction.clustered.bedpe.gz'): continue
    parse_sv_bam_list.append(bam)

# generate input list of 'SV merge'
unique_complist = []
merge_bedpe_list = []
for complist in sample_conf.sv_detection:
    control_panel_name = complist[2]
    if control_panel_name != None and control_panel_name not in unique_complist:
        unique_complist.append(control_panel_name)

for control_panel_name in unique_complist:
    if os.path.exists(run_conf.project_root + '/sv/non_matched_control_panel/' + control_panel_name + '.merged.junction.control.bedpe.gz'): continue
    tmp_list = []
    tmp_list.append(run_conf.project_root + '/sv/config/' + control_panel_name + ".control.yaml")
    for sample in sample_conf.control_panel[control_panel_name]:
        tmp_list.append(run_conf.project_root+ "/sv/"+ sample +"/"+ sample +".junction.clustered.bedpe.gz")
    merge_bedpe_list.append(tmp_list)

# generate input list of 'SV filt'
filt_bedpe_list = []
for complist in sample_conf.sv_detection:
    if os.path.exists(run_conf.project_root + '/sv/' + complist[0] +'/'+ complist[0] +'.genomonSV.result.txt'): continue
    filt_bedpe_list.append(run_conf.project_root+ "/sv/"+ complist[0] +"/"+ complist[0] +".junction.clustered.bedpe.gz")

# prepare output directories
if not os.path.isdir(run_conf.project_root): os.mkdir(run_conf.project_root)
if not os.path.isdir(run_conf.project_root + '/script'): os.mkdir(run_conf.project_root + '/script')
if not os.path.isdir(run_conf.project_root + '/log'): os.mkdir(run_conf.project_root + '/log')
if not os.path.isdir(run_conf.project_root + '/fastq'): os.mkdir(run_conf.project_root + '/fastq')
if not os.path.isdir(run_conf.project_root + '/bam'): os.mkdir(run_conf.project_root + '/bam')
if not os.path.isdir(run_conf.project_root + '/mutation'): os.mkdir(run_conf.project_root + '/mutation')
if not os.path.isdir(run_conf.project_root + '/mutation/control_panel'): os.mkdir(run_conf.project_root + '/mutation/control_panel')
if not os.path.isdir(run_conf.project_root + '/sv'): os.mkdir(run_conf.project_root + '/sv')
if not os.path.isdir(run_conf.project_root + '/sv/non_matched_control_panel'): os.mkdir(run_conf.project_root + '/sv/non_matched_control_panel')
if not os.path.isdir(run_conf.project_root + '/sv/config'): os.mkdir(run_conf.project_root + '/sv/config')
for outputfiles in (bam2fastq_output_list + linked_fastq_list):
    sample = os.path.basename(os.path.dirname(outputfiles[0]))
    fastq_dir = run_conf.project_root + '/fastq/' + sample
    bam_dir = run_conf.project_root + '/bam/' + sample
    if not os.path.isdir(fastq_dir): os.mkdir(fastq_dir)
    if not os.path.isdir(bam_dir): os.mkdir(bam_dir)

# prepare output directory for each sample and make mutation control panel file
for complist in sample_conf.mutation_call:
    # make dir
    mutation_dir = run_conf.project_root + '/mutation/' + complist[0]
    if not os.path.isdir(mutation_dir): os.mkdir(mutation_dir)
    # make the control panel text 
    control_panel_name = complist[2]
    if control_panel_name != None:
        control_panel_file = run_conf.project_root + '/mutation/control_panel/' + control_panel_name + ".control_panel.txt"
        with open(control_panel_file,  "w") as out_handle:
            for panel_sample in sample_conf.control_panel[control_panel_name]:
                out_handle.write(run_conf.project_root + '/bam/' + panel_sample + '/' + panel_sample + '.markdup.bam' + "\n")

# make SV configuration file
for complist in sample_conf.sv_detection:
    # make the control yaml file
    control_panel_name = complist[2]
    if control_panel_name != None:
        control_conf = run_conf.project_root + '/sv/config/' + control_panel_name + ".control.yaml"
        with open(control_conf,  "w") as out_handle:
            for sample in sample_conf.control_panel[control_panel_name]:
                out_handle.write(sample +": "+run_conf.project_root+ "/sv/"+ sample +"/"+ sample +".junction.clustered.bedpe.gz\n")

# link the import bam to project directory
@originate(sample_conf.bam_import.keys())
def link_import_bam(sample):
    bam = sample_conf.bam_import[sample]
    link_dir = run_conf.project_root + '/bam/' + sample
    bam_prefix, ext = os.path.splitext(bam)
    
    if not os.path.isdir(link_dir): os.mkdir(link_dir)
    if (not os.path.exists(link_dir +'/'+ sample +'.markdup.bam')) and (not os.path.exists(link_dir +'/'+ sample +'.markdup.bam.bai')): 
        os.symlink(bam, link_dir +'/'+ sample +'.markdup.bam')
        if (os.path.exists(bam +'.bai')):
            os.symlink(bam +'.bai', link_dir +'/'+ sample +'.markdup.bam.bai')
        elif (os.path.exists(bam_prefix +'.bai')):
            os.symlink(bam_prefix +'.bai', link_dir +'/'+ sample +'.markdup.bam.bai')

# convert bam to fastq
@originate(bam2fastq_output_list)
def bam2fastq(outputfiles):
    sample = os.path.basename(os.path.dirname(outputfiles[0]))
    output_dir = run_conf.project_root + '/fastq/' + sample
            
    arguments = {"biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "input_bam": sample_conf.bam_tofastq[sample],
                 "f1_name": outputfiles[0],
                 "f2_name": outputfiles[1],
                 "o1_name": output_dir + '/unmatched_first_output.txt',
                 "o2_name": output_dir + '/unmatched_second_output.txt',
                 "t": output_dir + '/temp.txt',
                 "s": output_dir + '/single_end_output.txt',
                 "log": run_conf.project_root + '/log'}
    bamtofastq.task_exec(arguments)


# link the input fastq to project directory
@originate(linked_fastq_list, sample_conf.fastq)
def link_input_fastq(output_file, sample_list_fastq):
    sample = os.path.basename(os.path.dirname(output_file[0]))
    fastq_dir = run_conf.project_root + '/fastq/' + sample
    fastq_prefix, ext = os.path.splitext(sample_list_fastq[sample][0][0])
    # Todo
    # 1. should compare the timestamps between input and linked file
    # 2. check md5sum ?
    if not os.path.exists(fastq_dir + '/1_1' + ext): os.symlink(sample_list_fastq[sample][0][0], fastq_dir + '/1_1' + ext)
    if not os.path.exists(fastq_dir + '/1_2' + ext): os.symlink(sample_list_fastq[sample][1][0], fastq_dir + '/1_2' + ext)


# split fastq
@subdivide([bam2fastq, link_input_fastq], formatter(".+/(.+).fastq"), "{path[0]}/*_*.fastq_split", "{path[0]}")
def split_files(input_files, output_files, target_dir):

    for oo in output_files:
        os.unlink(oo)

    arguments = {"lines": task_conf.get("split_fast", "split_fastq_line_number"),
                 "target_dir": target_dir,
                 "log": run_conf.project_root + '/log'}
    fastq_splitter.task_exec(arguments, 2)
    
    all_line_num = sum(1 for line in open(input_files[0]))
    with open(target_dir + "/fastq_line_num.txt",  "w") as out_handle:
        out_handle.write(str(all_line_num))

    os.unlink(input_files[0])
    os.unlink(input_files[1])


#bwa
@subdivide(split_files, formatter(".+/(.+)/1_0000.fastq_split"), add_inputs("{subpath[0][2]}/fastq/{subdir[0][0]}/2_0000.fastq_split"), "{subpath[0][2]}/bam/{subdir[0][0]}/{subdir[0][0]}_*.sorted.bam", "{subpath[0][2]}/fastq/{subdir[0][0]}", "{subpath[0][2]}/bam/{subdir[0][0]}")
def map_dna_sequence(input_files, output_files, input_dir, output_dir):

    sample_name = os.path.basename(output_dir)

    all_line_num = 0
    with open(input_dir + "/fastq_line_num.txt") as in_handle:
        tmp_num = in_handle.read()
        all_line_num = int(tmp_num)
    split_lines = task_conf.get("split_fast", "split_fastq_line_number")

    ans_quotient = all_line_num / int(split_lines)
    ans_remainder = all_line_num % int(split_lines)
    max_task_id = ans_quotient if ans_remainder == 0 else ans_quotient + 1
    
    arguments = {"input_dir": input_dir,
                 "output_dir": output_dir,
                 "sample_name": sample_name,
                 "bwa": genomon_conf.get("SOFTWARE", "bwa"),
                 "bwa_params": task_conf.get("bwa_mem", "bwa_params"),
                 "ref_fa":genomon_conf.get("REFERENCE", "ref_fasta"),
                 "biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "log": run_conf.project_root + '/log'}

    bwa_align.task_exec(arguments, max_task_id) 

    for task_id in range(max_task_id):
        num = str(task_id).zfill(4)
        os.unlink(input_dir +'/1_'+str(num)+'.fastq_split')
        os.unlink(input_dir +'/2_'+str(num)+'.fastq_split')
        os.unlink(output_dir+'/'+sample_name+'_'+str(num)+'.bwa.sam')


# merge sorted bams into one and mark duplicate reads with biobambam
@collate(map_dna_sequence, formatter(), "{subpath[0][2]}/bam/{subdir[0][0]}/{subdir[0][0]}.markdup.bam")
def markdup(input_files, output_file):

    output_prefix, ext = os.path.splitext(output_file)

    input_bam_files = ""
    for input_file in input_files:
        input_bam_files = input_bam_files + " I=" + input_file

    arguments = {"biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "out_prefix": output_prefix,
                 "input_bam_files": input_bam_files,
                 "out_bam": output_file,
                 "log": run_conf.project_root + '/log'}

    markduplicates.task_exec(arguments)

    for input_file in input_files:
        os.unlink(input_file)
        os.unlink(input_file + ".bai")


# identify mutations
@follows( markdup )
@follows( link_import_bam )
@subdivide(markdup_bam_list, formatter(), "{subpath[0][2]}/mutation/{subdir[0][0]}/{subdir[0][0]}_mutations_candidate.*.hg19_multianno.txt", "{subpath[0][2]}/mutation/{subdir[0][0]}")
def identify_mutations(input_file, output_file, output_dir):

    sample_name = os.path.basename(output_dir)

    arguments = {
        # fisher mutation
        "fisher": genomon_conf.get("SOFTWARE", "fisher"),
        "map_quality": task_conf.get("fisher_mutation_call", "map_quality"),
        "base_quality": task_conf.get("fisher_mutation_call", "base_quality"),
        "min_allele_freq": task_conf.get("fisher_mutation_call", "disease_min_allele_frequency"),
        "max_allele_freq": task_conf.get("fisher_mutation_call", "control_max_allele_frequency"),
        "min_depth": task_conf.get("fisher_mutation_call", "min_depth"),
        "fisher_thres": task_conf.get("fisher_mutation_call", "fisher_thres_hold"),
        "post_10_q": task_conf.get("fisher_mutation_call", "post_10_q"),
        # realignment filter
        "mutfilter": genomon_conf.get("SOFTWARE", "mutfilter"),
        "realign_min_mismatch": task_conf.get("realignment_filter","disease_min_mismatch"),
        "realign_max_mismatch": task_conf.get("realignment_filter","control_max_mismatch"),
        "realign_score_diff": task_conf.get("realignment_filter","score_diff"),
        "realign_window_size": task_conf.get("realignment_filter","window_size"),
        "realign_max_depth": task_conf.get("realignment_filter","max_depth"),
        # indel filter
        "indel_search_length": task_conf.get("indel_filter","search_length"),
        "indel_neighbor": task_conf.get("indel_filter","neighbor"),
        "indel_base_quality": task_conf.get("indel_filter","base_quality"),
        "indel_min_depth": task_conf.get("indel_filter","min_depth"),
        "indel_min_mismatch": task_conf.get("indel_filter","max_mismatch"),
        "indel_min_allele_freq": task_conf.get("indel_filter","max_allele_freq"),
        # breakpoint filter
        "bp_max_depth": task_conf.get("breakpoint_filter","max_depth"),
        "bp_min_clip_size": task_conf.get("breakpoint_filter","min_clip_size"),
        "bp_junc_num_thres": task_conf.get("breakpoint_filter","junc_num_thres"),
        "bp_map_quality": task_conf.get("breakpoint_filter","map_quality"),
        # EB filter
        "EBFilter": genomon_conf.get("SOFTWARE", "ebfilter"),
        "eb_map_quality": task_conf.get("eb_filter","map_quality"),
        "eb_base_quality": task_conf.get("eb_filter","base_quality"),
        "control_bam_list": input_file[2],
        # annovar
        "active_annovar_flag": task_conf.get("annotation", "active_annovar_flag"),
        "annovar": genomon_conf.get("SOFTWARE", "annovar"),
        "table_annovar_params": task_conf.get("annotation", "table_annovar_params"),
        # commmon
        "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
        "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
        "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
        "ref_fa":genomon_conf.get("REFERENCE", "ref_fasta"),
        "interval_list": genomon_conf.get("REFERENCE", "interval_list"),
        "simple_repeat_db":genomon_conf.get("REFERENCE", "simple_repeat_tabix_db"),
        "disease_bam": input_file[0],
        "control_bam": input_file[1],
        "out_prefix": output_dir + '/' + sample_name,
        "samtools": genomon_conf.get("SOFTWARE", "samtools"),
        "blat": genomon_conf.get("SOFTWARE", "blat"),
        "log": run_conf.project_root + '/log'}

    interval_list = genomon_conf.get("REFERENCE", "interval_list")
    max_task_id = sum(1 for line in open(interval_list))

    mutation_call.task_exec(arguments, max_task_id)
     
    for task_id in range(1,(max_task_id + 1)):
        if os.path.exists(output_dir+'/'+sample_name+'.fisher_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.fisher_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.realignment_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.realignment_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.indel_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.indel_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.breakpoint_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.breakpoint_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.simplerepeat_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.simplerepeat_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.ebfilter_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.ebfilter_mutations.'+str(task_id)+'.txt')


# merge candidate mutations
@collate(identify_mutations, formatter(".+/(.+)_mutations_candidate.(.+).hg19_multianno.txt"), "{path[0]}/{subdir[0][0]}_genomon_mutations.result.txt", "{subpath[0][2]}/mutation/{subdir[0][0]}")
def merge_mutation(input_files, output_file, output_dir):

    sample_name = os.path.basename(output_dir)

    interval_list = genomon_conf.get("REFERENCE", "interval_list")
    max_task_id = sum(1 for line in open(interval_list))

    arguments = {
        "active_annovar_flag": task_conf.get("annotation", "active_annovar_flag"),
        "filecount": max_task_id,
        "out_prefix": output_dir + '/' + sample_name,
        "log": run_conf.project_root + '/log'}

    mutation_merge.task_exec(arguments)
     
    for task_id in range(1,(max_task_id + 1)):
        input_file = output_dir+'/'+sample_name+'_mutations_candidate.'+str(task_id)+'.hg19_multianno.txt'
        os.unlink(input_file)


# parse SV 
@follows( link_import_bam )
@follows( markdup )
@transform(parse_sv_bam_list, formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.junction.clustered.bedpe.gz")
def parse_sv(input_file, output_file):

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)

    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    yaml_flag = False
    sv_sampleConf = {"target": {}, "matched_control": {}, "non_matched_control_panel": {}}
    for complist in sample_conf.sv_detection:
        if sample_name == complist[0]:

            # tumor:exist, matched_normal:exist, non-matched-normal:exist
            if complist[0] != None and complist[1] != None and complist[2] != None:
                sv_sampleConf["target"]["label"] = sample_name
                sv_sampleConf["target"]["path_to_bam"] = input_file
                sv_sampleConf["target"]["path_to_output_dir"] = dir_name
                sv_sampleConf["matched_control"]["use"] = True
                sv_sampleConf["matched_control"]["path_to_bam"] = run_conf.project_root + '/bam/' + complist[1] + '/' + complist[1] + '.markdup.bam'
                sv_sampleConf["non_matched_control_panel"]["use"] = True
                sv_sampleConf["non_matched_control_panel"]["matched_control_label"] = complist[1]
                sv_sampleConf["non_matched_control_panel"]["data_path"] = run_conf.project_root +"/sv/non_matched_control_panel/"+ complist[2] +".merged.junction.control.bedpe.gz"
                yaml_flag = True
                break

            # tumor:exist, matched_normal:exist, non-matched-normal:none
            elif complist[0] != None and complist[1] != None and complist[2] == None:
                sv_sampleConf["target"]["label"] = sample_name
                sv_sampleConf["target"]["path_to_bam"] = input_file
                sv_sampleConf["target"]["path_to_output_dir"] = dir_name
                sv_sampleConf["matched_control"]["use"] = True
                sv_sampleConf["matched_control"]["path_to_bam"] = run_conf.project_root + '/bam/' + complist[1] + '/' + complist[1] + '.markdup.bam'
                sv_sampleConf["non_matched_control_panel"]["use"] = False
                yaml_flag = True
                break

            # tumor:exist, matched_normal:none, non-matched-normal:exist
            elif complist[0] != None and complist[1] == None and complist[2] != None:
                sv_sampleConf["target"]["label"] = sample_name
                sv_sampleConf["target"]["path_to_bam"] = input_file
                sv_sampleConf["target"]["path_to_output_dir"] = dir_name
                sv_sampleConf["matched_control"]["use"] = False
                sv_sampleConf["matched_control"]["path_to_bam"] = None
                sv_sampleConf["non_matched_control_panel"]["use"] = True
                sv_sampleConf["non_matched_control_panel"]["matched_control_label"] = None
                sv_sampleConf["non_matched_control_panel"]["data_path"] = run_conf.project_root +"/sv/non_matched_control_panel/"+ complist[2] +".merged.junction.control.bedpe.gz"
                yaml_flag = True
                break

    if not yaml_flag:
        sv_sampleConf["target"]["label"] = sample_name
        sv_sampleConf["target"]["path_to_bam"] = input_file
        sv_sampleConf["target"]["path_to_output_dir"] = dir_name
        sv_sampleConf["matched_control"]["use"] = False
        sv_sampleConf["matched_control"]["path_to_bam"] = None
        sv_sampleConf["non_matched_control_panel"]["use"] = False

    sample_yaml = run_conf.project_root + "/sv/config/" + sample_name + ".yaml"
    hOUT = open(sample_yaml, "w")
    print >> hOUT, yaml.dump(sv_sampleConf, default_flow_style = False)
    hOUT.close()

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "sample_conf": sample_yaml,
                 "param_conf": task_conf.get("genomon_sv", "param_file"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "log": run_conf.project_root + '/log'}
    sv_parse.task_exec(arguments)


# merge SV
@follows( parse_sv )
@transform(merge_bedpe_list, formatter(".+/(?P<NAME>.+).control.yaml"), "{subpath[0][2]}/sv/non_matched_control_panel/{NAME[0]}.merged.junction.control.bedpe.gz")
def merge_sv(input_files,  output_file):

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "control_conf": input_files[0],
                 "bedpe": output_file,
                 "param_conf": task_conf.get("genomon_sv", "param_file"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "log": run_conf.project_root + '/log'}
    sv_merge.task_exec(arguments)


# filt SV
@follows( merge_sv )
@transform(filt_bedpe_list, formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.genomonSV.result.txt")
def filt_sv(input_files,  output_file):

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    sample_yaml = run_conf.project_root + "/sv/config/" + sample_name + ".yaml"

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "sample_conf": sample_yaml,
                 "param_conf": task_conf.get("genomon_sv", "param_file"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "log": run_conf.project_root + '/log'}
    sv_filt.task_exec(arguments)

