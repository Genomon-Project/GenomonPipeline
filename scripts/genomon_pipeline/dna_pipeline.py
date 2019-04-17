import os
import shutil
import glob
from ruffus import *
from genomon_pipeline.config.run_conf import *
from genomon_pipeline.config.genomon_conf import *
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
from genomon_pipeline.dna_resource.qc_bamstats import *
from genomon_pipeline.dna_resource.qc_coverage import *
from genomon_pipeline.dna_resource.qc_merge import *
from genomon_pipeline.dna_resource.post_analysis import *
from genomon_pipeline.dna_resource.pre_pmsignature import *
from genomon_pipeline.dna_resource.pmsignature import *
from genomon_pipeline.dna_resource.paplot import *

# set task classes
bamtofastq = Bam2Fastq(genomon_conf.get("bam2fastq", "qsub_option"), run_conf.drmaa)
fastq_splitter = Fastq_splitter(genomon_conf.get("split_fastq", "qsub_option"), run_conf.drmaa)
bwa_align = Bwa_align(genomon_conf.get("bwa_mem", "qsub_option"), run_conf.drmaa)
markduplicates = Markduplicates(genomon_conf.get("markduplicates", "qsub_option"), run_conf.drmaa)
mutation_call = Mutation_call(genomon_conf.get("mutation_call", "qsub_option"), run_conf.drmaa)
mutation_merge = Mutation_merge(genomon_conf.get("mutation_merge", "qsub_option"), run_conf.drmaa)
sv_parse = SV_parse(genomon_conf.get("sv_parse", "qsub_option"), run_conf.drmaa)
sv_merge = SV_merge(genomon_conf.get("sv_merge", "qsub_option"), run_conf.drmaa)
sv_filt = SV_filt(genomon_conf.get("sv_filt", "qsub_option"), run_conf.drmaa)
r_qc_bamstats = Res_QC_Bamstats(genomon_conf.get("qc_bamstats", "qsub_option"), run_conf.drmaa)
r_qc_coverage = Res_QC_Coverage(genomon_conf.get("qc_coverage", "qsub_option"), run_conf.drmaa)
r_qc_merge = Res_QC_Merge(genomon_conf.get("qc_merge", "qsub_option"), run_conf.drmaa)
r_paplot = Res_PA_Plot(genomon_conf.get("paplot", "qsub_option"), run_conf.drmaa)
r_post_analysis = Res_PostAnalysis(genomon_conf.get("post_analysis", "qsub_option"), run_conf.drmaa)
r_pre_pmsignature = Res_PrePmsignature(genomon_conf.get("pre_pmsignature", "qsub_option"), run_conf.drmaa)
r_pmsignature_ind = Res_Pmsignature(genomon_conf.get("pmsignature_ind", "qsub_option"), run_conf.drmaa)
r_pmsignature_full = Res_Pmsignature(genomon_conf.get("pmsignature_full", "qsub_option"), run_conf.drmaa)

_debug = False
if genomon_conf.has_section("develop"):
    if genomon_conf.has_option("develop", "debug") == True:
        _debug = genomon_conf.getboolean("develop", "debug")

# generate output list of 'linked fastq'
linked_fastq_list = []
for sample in sample_conf.fastq:
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/1.sorted.bam'): continue
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/' + sample + '.markdup.bam'): continue

    link_fastq_arr1 = []
    link_fastq_arr2 = []
    for (count, fastq_file) in enumerate(sample_conf.fastq[sample][0]):
        fastq_prefix, ext = os.path.splitext(fastq_file)
        link_fastq_arr1.append(run_conf.project_root + '/fastq/' + sample + '/' + str(count+1) + '_1' + ext)
        link_fastq_arr2.append(run_conf.project_root + '/fastq/' + sample + '/' + str(count+1) + '_2' + ext)
    linked_fastq_list.append([link_fastq_arr1,link_fastq_arr2])

# generate output list of 'bam2fastq'
bam2fastq_output_list = []
for sample in sample_conf.bam_tofastq:
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/1.sorted.bam'): continue
    if os.path.exists(run_conf.project_root + '/bam/' + sample + '/' + sample + '.markdup.bam'): continue
    bam2fastq_arr1 = []
    bam2fastq_arr2 = []
    bam2fastq_arr1.append(run_conf.project_root + '/fastq/' + sample + '/1_1.fastq')
    bam2fastq_arr2.append(run_conf.project_root + '/fastq/' + sample + '/1_2.fastq')
    bam2fastq_output_list.append([bam2fastq_arr1,bam2fastq_arr2])

# generate input list of 'mutation call'
markdup_bam_list = []
merge_mutation_list = []
for complist in sample_conf.mutation_call:
     if os.path.exists(run_conf.project_root + '/mutation/' + complist[0] + '/' + complist[0] + '.genomon_mutation.result.filt.txt'): continue
     tumor_bam  = run_conf.project_root + '/bam/' + complist[0] + '/' + complist[0] + '.markdup.bam'
     normal_bam = run_conf.project_root + '/bam/' + complist[1] + '/' + complist[1] + '.markdup.bam' if complist[1] != None else None
     panel = run_conf.project_root + '/mutation/control_panel/' + complist[2] + ".control_panel.txt" if complist[2] != None else None
     markdup_bam_list.append([tumor_bam, normal_bam, panel])


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
    if os.path.exists(run_conf.project_root + '/sv/' + sample_name + '/' + sample_name + '.junction.clustered.bedpe.gz') and os.path.exists(run_conf.project_root + '/sv/' + sample_name + '/' + sample_name + '.junction.clustered.bedpe.gz.tbi'): continue
    parse_sv_bam_list.append(bam)

# generate input list of 'SV merge'
unique_complist = []
merge_bedpe_list = []
for complist in sample_conf.sv_detection:
    control_panel_name = complist[2]
    if control_panel_name != None and control_panel_name not in unique_complist:
        unique_complist.append(control_panel_name)

for control_panel_name in unique_complist:
    if os.path.exists(run_conf.project_root + '/sv/non_matched_control_panel/' + control_panel_name + '.merged.junction.control.bedpe.gz') and os.path.exists(run_conf.project_root + '/sv/non_matched_control_panel/' + control_panel_name + '.merged.junction.control.bedpe.gz.tbi'): continue
    tmp_list = []
    tmp_list.append(run_conf.project_root + '/sv/control_panel/' + control_panel_name + ".control_info.txt")
    for sample in sample_conf.control_panel[control_panel_name]:
        tmp_list.append(run_conf.project_root+ "/sv/"+ sample +"/"+ sample +".junction.clustered.bedpe.gz")
    merge_bedpe_list.append(tmp_list)

# generate input list of 'SV filt'
filt_bedpe_list = []
for complist in sample_conf.sv_detection:
    if os.path.exists(run_conf.project_root + '/sv/' + complist[0] +'/'+ complist[0] +'.genomonSV.result.filt.txt'): continue
    filt_bedpe_list.append(run_conf.project_root+ "/sv/"+ complist[0] +"/"+ complist[0] +".junction.clustered.bedpe.gz")

# generate input list of 'qc'
qc_bamstats_list = []
qc_coverage_list = []
qc_merge_list = []
for sample in sample_conf.qc:
    if os.path.exists(run_conf.project_root + '/qc/' + sample + '/' + sample + '.genomonQC.result.txt'): continue
    qc_merge_list.append(
        [run_conf.project_root + '/qc/' + sample + '/' + sample + '.bamstats',
         run_conf.project_root + '/qc/' + sample + '/' + sample + '.coverage'])
    if not os.path.exists(run_conf.project_root + '/qc/' + sample + '/' + sample + '.bamstats'):
        qc_bamstats_list.append(run_conf.project_root + '/bam/' + sample +'/'+ sample +'.markdup.bam')
    if not os.path.exists(run_conf.project_root + '/qc/' + sample + '/' + sample + '.coverage'):
        qc_coverage_list.append(run_conf.project_root + '/bam/' + sample +'/'+ sample +'.markdup.bam')

### 
# input/output lists for post-analysis
###
genomon_conf_name, genomon_conf_ext = os.path.splitext(os.path.basename(run_conf.genomon_conf_file))
sample_conf_name, sample_conf_ext = os.path.splitext(os.path.basename(run_conf.sample_conf_file))

# generate input list of 'post analysis for mutation'
pa_outputs_mutation = r_post_analysis.output_files("mutation", sample_conf.mutation_call, run_conf.project_root, sample_conf_name, genomon_conf)

pa_inputs_mutation = []
if pa_outputs_mutation["run_pa"] == True:
    for complist in sample_conf.mutation_call:
        pa_inputs_mutation.append(run_conf.project_root + '/mutation/' + complist[0] +'/'+ complist[0] +'.genomon_mutation.result.filt.txt')
        
# generate input list of 'post analysis for SV'
pa_outputs_sv = r_post_analysis.output_files("sv", sample_conf.sv_detection, run_conf.project_root, sample_conf_name, genomon_conf)

pa_inputs_sv = []
if pa_outputs_sv["run_pa"] == True:
    for complist in sample_conf.sv_detection:
        pa_inputs_sv.append(run_conf.project_root + '/sv/' + complist[0] +'/'+ complist[0] +'.genomonSV.result.filt.txt')
        
# generate input list of 'post analysis for qc'
pa_outputs_qc = r_post_analysis.output_files("qc", sample_conf.qc, run_conf.project_root, sample_conf_name, genomon_conf)

pa_inputs_qc = []
if pa_outputs_qc["run_pa"] == True:
    for sample in sample_conf.qc:
        pa_inputs_qc.append(run_conf.project_root + '/qc/' + sample + '/' + sample + '.genomonQC.result.txt')

### 
# input/output lists for paplot
###
paplot_output = run_conf.project_root + '/paplot/' + sample_conf_name + '/index.html'

## mutation
use_mutations = []
if pa_outputs_mutation["case1"]["output_filt"] != "":
    use_mutations.append(pa_outputs_mutation["case1"]["output_filt"])
if pa_outputs_mutation["case2"]["output_filt"] != "" and genomon_conf.getboolean("paplot", "include_unpanel"):
    use_mutations.append(pa_outputs_mutation["case2"]["output_filt"])
if pa_outputs_mutation["case3"]["output_filt"] != "" and genomon_conf.getboolean("paplot", "include_unpair"):
    use_mutations.append(pa_outputs_mutation["case3"]["output_filt"])
if pa_outputs_mutation["case4"]["output_filt"] != "" and genomon_conf.getboolean("paplot", "include_unpanel") and genomon_conf.getboolean("paplot", "include_unpair"):
    use_mutations.append(pa_outputs_mutation["case4"]["output_filt"])

paplot_inputs_mutation = []
if os.path.exists(paplot_output) == False or pa_outputs_mutation["run_pa"] == True:
    paplot_inputs_mutation.extend(use_mutations)

## pmsignature
# ind
ind_outputs = []
ind_exists = True
for i in range(genomon_conf.getint("pmsignature_ind", "signum_min"), genomon_conf.getint("pmsignature_ind", "signum_max") + 1):
    fname = run_conf.project_root + '/pmsignature/' + sample_conf_name + '/pmsignature.ind.result.%d.json' % i
    ind_outputs.append(fname)
    if not os.path.exists(fname): ind_exists = False
        
run_ind = False
paplot_inputs_ind = []
if len(sample_conf.mutation_call) > 0 and genomon_conf.getboolean("pmsignature_ind", "enable") and len(use_mutations) > 0:
    if ind_exists == False: run_ind = True
    elif pa_outputs_mutation["run_pa"] == True: run_ind = True
    elif not os.path.exists(run_conf.project_root + '/pmsignature/' + sample_conf_name + '/mutation.cut.txt'): run_ind = True
    if os.path.exists(paplot_output) == False or run_ind == True:    
        paplot_inputs_ind.extend(ind_outputs)
    
# full
full_outputs = []
full_exists = True
for i in range(genomon_conf.getint("pmsignature_full", "signum_min"), genomon_conf.getint("pmsignature_full", "signum_max") + 1):
    fname = run_conf.project_root + '/pmsignature/' + sample_conf_name + '/pmsignature.full.result.%d.json' % i
    full_outputs.append(fname)
    if not os.path.exists(fname): full_exists = False
        
run_full = False
paplot_inputs_full = []
if len(sample_conf.mutation_call) > 0 and genomon_conf.getboolean("pmsignature_full", "enable") and len(use_mutations) > 0:
    if full_exists == False: run_full = True
    elif pa_outputs_mutation["run_pa"] == True: run_full = True
    elif not os.path.exists(run_conf.project_root + '/pmsignature/' + sample_conf_name + '/mutation.cut.txt'): run_full = True
    if os.path.exists(paplot_output) == False or run_full == True:
        paplot_inputs_full.extend(full_outputs)
 
pmsignature_inputs = []
if run_ind == True or run_full == True: 
    pmsignature_inputs.extend(use_mutations)

## sv
paplot_inputs_sv = []
if os.path.exists(paplot_output) == False or pa_outputs_sv["run_pa"] == True:

    if pa_outputs_sv["case1"]["output_filt"] != "":
        paplot_inputs_sv.append(pa_outputs_sv["case1"]["output_filt"])
    if pa_outputs_sv["case2"]["output_filt"] != "" and genomon_conf.getboolean("paplot", "include_unpanel"):
        paplot_inputs_sv.append(pa_outputs_sv["case2"]["output_filt"])
    if pa_outputs_sv["case3"]["output_filt"] != "" and genomon_conf.getboolean("paplot", "include_unpair"):
        paplot_inputs_sv.append(pa_outputs_sv["case3"]["output_filt"])
    if pa_outputs_sv["case4"]["output_filt"] != "" and genomon_conf.getboolean("paplot", "include_unpanel") and genomon_conf.getboolean("paplot", "include_unpair"):
        paplot_inputs_sv.append(pa_outputs_sv["case4"]["output_filt"])

## qc
paplot_inputs_qc = []
if os.path.exists(paplot_output) == False or pa_outputs_qc["run_pa"] == True:
    paplot_inputs_qc.extend(pa_outputs_qc["outputs"])

paplot_inputs = []
paplot_inputs.extend(paplot_inputs_qc)
paplot_inputs.extend(paplot_inputs_sv)
paplot_inputs.extend(paplot_inputs_mutation)
paplot_inputs.extend(paplot_inputs_ind)
paplot_inputs.extend(paplot_inputs_full)

if _debug:
    from pprint import pprint
    print ("post-analysis-mutation");  pprint (pa_outputs_mutation); print ("post-analysis-sv");  pprint (pa_outputs_sv); print ("post-analysis-qc");  pprint (pa_outputs_qc)
    print ("paplot"); pprint (paplot_inputs)
    print ("pmsignature"); pprint (pmsignature_inputs)

# prepare output directories
if not os.path.isdir(run_conf.project_root): os.mkdir(run_conf.project_root)
if not os.path.isdir(run_conf.project_root + '/script'): os.mkdir(run_conf.project_root + '/script')
if not os.path.isdir(run_conf.project_root + '/script/sv_merge'): os.mkdir(run_conf.project_root + '/script/sv_merge')
if not os.path.isdir(run_conf.project_root + '/log'): os.mkdir(run_conf.project_root + '/log')
if not os.path.isdir(run_conf.project_root + '/log/sv_merge'): os.mkdir(run_conf.project_root + '/log/sv_merge')
if not os.path.isdir(run_conf.project_root + '/fastq'): os.mkdir(run_conf.project_root + '/fastq')
if not os.path.isdir(run_conf.project_root + '/bam'): os.mkdir(run_conf.project_root + '/bam')
if not os.path.isdir(run_conf.project_root + '/mutation'): os.mkdir(run_conf.project_root + '/mutation')
if not os.path.isdir(run_conf.project_root + '/mutation/control_panel'): os.mkdir(run_conf.project_root + '/mutation/control_panel')
if not os.path.isdir(run_conf.project_root + '/mutation/hotspot'): os.mkdir(run_conf.project_root + '/mutation/hotspot')
if not os.path.isdir(run_conf.project_root + '/sv'): os.mkdir(run_conf.project_root + '/sv')
if not os.path.isdir(run_conf.project_root + '/sv/non_matched_control_panel'): os.mkdir(run_conf.project_root + '/sv/non_matched_control_panel')
if not os.path.isdir(run_conf.project_root + '/sv/control_panel'): os.mkdir(run_conf.project_root + '/sv/control_panel')
if not os.path.isdir(run_conf.project_root + '/qc'): os.mkdir(run_conf.project_root + '/qc')
for sample in sample_conf.qc:
    if not os.path.isdir(run_conf.project_root + '/qc/' + sample): os.mkdir(run_conf.project_root + '/qc/' + sample)

if (genomon_conf.getboolean("post_analysis", "enable") == True):
    if not os.path.exists(run_conf.project_root + '/post_analysis'): os.mkdir(run_conf.project_root + '/post_analysis')
    if not os.path.exists(run_conf.project_root + '/post_analysis/' + sample_conf_name): os.mkdir(run_conf.project_root + '/post_analysis/' + sample_conf_name)
    if not os.path.isdir(run_conf.project_root + '/script/post_analysis'): os.mkdir(run_conf.project_root + '/script/post_analysis')
    if not os.path.isdir(run_conf.project_root + '/log/post_analysis'): os.mkdir(run_conf.project_root + '/log/post_analysis')
    
    if (genomon_conf.getboolean("paplot", "enable") == True):
        if not os.path.isdir(run_conf.project_root + '/paplot/'): os.mkdir(run_conf.project_root + '/paplot/')
        if not os.path.isdir(run_conf.project_root + '/paplot/' + sample_conf_name): os.mkdir(run_conf.project_root + '/paplot/' + sample_conf_name)
        if not os.path.isdir(run_conf.project_root + '/script/paplot'): os.mkdir(run_conf.project_root + '/script/paplot')
        if not os.path.isdir(run_conf.project_root + '/log/paplot'): os.mkdir(run_conf.project_root + '/log/paplot')

    if (genomon_conf.getboolean("pmsignature_ind", "enable") == True) or (genomon_conf.getboolean("pmsignature_full", "enable") == True):
        if not os.path.isdir(run_conf.project_root + '/pmsignature/'): os.mkdir(run_conf.project_root + '/pmsignature/')
        if not os.path.isdir(run_conf.project_root + '/pmsignature/' + sample_conf_name): os.mkdir(run_conf.project_root + '/pmsignature/' + sample_conf_name)
        if not os.path.isdir(run_conf.project_root + '/script/pmsignature'): os.mkdir(run_conf.project_root + '/script/pmsignature')
        if not os.path.isdir(run_conf.project_root + '/log/pmsignature'): os.mkdir(run_conf.project_root + '/log/pmsignature')

if not os.path.isdir(run_conf.project_root + '/config'): os.mkdir(run_conf.project_root + '/config')

for outputfiles in (bam2fastq_output_list, linked_fastq_list):
    for outputfile in outputfiles:
        sample = os.path.basename(os.path.dirname(outputfile[0][0]))
        fastq_dir = run_conf.project_root + '/fastq/' + sample
        bam_dir = run_conf.project_root + '/bam/' + sample
        if not os.path.isdir(fastq_dir): os.mkdir(fastq_dir)
        if not os.path.isdir(bam_dir): os.mkdir(bam_dir)

for target_sample_dict in (sample_conf.bam_import, sample_conf.fastq, sample_conf.bam_tofastq):
    for sample in target_sample_dict:
        script_dir = run_conf.project_root + '/script/' + sample
        log_dir = run_conf.project_root + '/log/' + sample
        if not os.path.isdir(script_dir): os.mkdir(script_dir)
        if not os.path.isdir(log_dir): os.mkdir(log_dir)

shutil.copyfile(run_conf.genomon_conf_file, run_conf.project_root + '/config/' + genomon_conf_name +'_'+ run_conf.analysis_timestamp + genomon_conf_ext)
shutil.copyfile(run_conf.sample_conf_file, run_conf.project_root + '/config/' + sample_conf_name +'_'+ run_conf.analysis_timestamp + sample_conf_ext)

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
        control_conf = run_conf.project_root + '/sv/control_panel/' + control_panel_name + ".control_info.txt"
        with open(control_conf,  "w") as out_handle:
            for sample in sample_conf.control_panel[control_panel_name]:
                out_handle.write(sample+ "\t"+ run_conf.project_root+ "/sv/"+ sample +"/"+ sample+ "\n")

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
    sample = os.path.basename(os.path.dirname(outputfiles[0][0]))
    output_dir = run_conf.project_root + '/fastq/' + sample
            
    arguments = {"biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "param": genomon_conf.get("bam2fastq", "params"),
                 "input_bam": sample_conf.bam_tofastq[sample],
                 "f1_name": outputfiles[0][0],
                 "f2_name": outputfiles[1][0],
                 "o1_name": output_dir + '/unmatched_first_output.txt',
                 "o2_name": output_dir + '/unmatched_second_output.txt',
                 "t": output_dir + '/temp.txt',
                 "s": output_dir + '/single_end_output.txt'}
    bamtofastq.task_exec(arguments, run_conf.project_root + '/log/' + sample, run_conf.project_root + '/script/'+ sample)


# link the input fastq to project directory
@originate(linked_fastq_list)
def link_input_fastq(output_file):
    sample = os.path.basename(os.path.dirname(output_file[0][0]))
    fastq_dir = run_conf.project_root + '/fastq/' + sample
    fastq_prefix, ext = os.path.splitext(sample_conf.fastq[sample][0][0])
    # Todo
    # 1. should compare the timestamps between input and linked file
    # 2. check md5sum ?
    for (count, fastq_files) in enumerate(sample_conf.fastq[sample][0]):
        fastq_prefix, ext = os.path.splitext(fastq_files)
        if not os.path.exists(fastq_dir + '/'+str(count+1)+'_1'+ ext): os.symlink(sample_conf.fastq[sample][0][count], fastq_dir + '/'+str(count+1)+'_1'+ ext)
        if not os.path.exists(fastq_dir + '/'+str(count+1)+'_2'+ ext): os.symlink(sample_conf.fastq[sample][1][count], fastq_dir + '/'+str(count+1)+'_2'+ ext)


# split fastq
@subdivide([bam2fastq, link_input_fastq], formatter(), "{path[0]}/*_*.fastq_split", "{path[0]}")
def split_files(input_files, output_files, target_dir):

    sample_name = os.path.basename(target_dir)

    for oo in output_files:
        os.unlink(oo)

    split_lines = genomon_conf.get("split_fastq", "split_fastq_line_number")

    input_prefix, ext = os.path.splitext(input_files[0][0])
    arguments = {"lines": split_lines,
                 "fastq_filter": genomon_conf.get("split_fastq", "fastq_filter"),
                 "target_dir": target_dir,
                 "ext": ext}
    
    fastq_splitter.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/'+ sample_name, 2)
   
    file_list = glob.glob(target_dir + '/1_*.fastq_split')
    file_list.sort()
    last_file_lines = sum(1 for line in open(file_list[-1]))
    all_line_num = ((len(file_list)-1)*int(split_lines)) + last_file_lines
    
    with open(target_dir + "/fastq_line_num.txt",  "w") as out_handle:
        out_handle.write(str(all_line_num)+"\n")
    
    for input_fastq in input_files[0]:
        os.unlink(input_fastq)
    for input_fastq in input_files[1]:
        os.unlink(input_fastq)


#bwa
@subdivide(split_files, formatter(".+/(.+)/1_0000.fastq_split"), add_inputs("{subpath[0][2]}/fastq/{subdir[0][0]}/2_0000.fastq_split"), "{subpath[0][2]}/bam/{subdir[0][0]}/{subdir[0][0]}_*.sorted.bam", "{subpath[0][2]}/fastq/{subdir[0][0]}", "{subpath[0][2]}/bam/{subdir[0][0]}")
def map_dna_sequence(input_files, output_files, input_dir, output_dir):

    sample_name = os.path.basename(output_dir)

    all_line_num = 0
    with open(input_dir + "/fastq_line_num.txt") as in_handle:
        tmp_num = in_handle.read()
        all_line_num = int(tmp_num)
    split_lines = genomon_conf.get("split_fastq", "split_fastq_line_number")

    ans_quotient = all_line_num / int(split_lines)
    ans_remainder = all_line_num % int(split_lines)
    max_task_id = ans_quotient if ans_remainder == 0 else ans_quotient + 1
    
    arguments = {"input_dir": input_dir,
                 "output_dir": output_dir,
                 "sample_name": sample_name,
                 "bwa": genomon_conf.get("SOFTWARE", "bwa"),
                 "bwa_params": genomon_conf.get("bwa_mem", "bwa_params"),
                 "ref_fa":genomon_conf.get("REFERENCE", "ref_fasta"),
                 "biobambam": genomon_conf.get("SOFTWARE", "biobambam")}

    bwa_align.task_exec(arguments, run_conf.project_root + '/log/' + sample_name , run_conf.project_root + '/script/' + sample_name, max_task_id) 

    for task_id in range(max_task_id):
        num = str(task_id).zfill(4)
        os.unlink(input_dir +'/1_'+str(num)+'.fastq_split')
        os.unlink(input_dir +'/2_'+str(num)+'.fastq_split')
        os.unlink(output_dir+'/'+sample_name+'_'+str(num)+'.bwa.sam')


# merge sorted bams into one and mark duplicate reads with biobambam
@collate(map_dna_sequence, formatter(), "{subpath[0][2]}/bam/{subdir[0][0]}/{subdir[0][0]}.markdup.bam", "{subpath[0][2]}/bam/{subdir[0][0]}")
def markdup(input_files, output_file, output_dir):

    sample_name = os.path.basename(output_dir)

    output_prefix, ext = os.path.splitext(output_file)

    input_bam_files = ""
    for input_file in input_files:
        input_bam_files = input_bam_files + " I=" + input_file

    arguments = {"biobambam": genomon_conf.get("SOFTWARE", "biobambam"),
                 "out_prefix": output_prefix,
                 "input_bam_files": input_bam_files,
                 "out_bam": output_file}

    markduplicates.task_exec(arguments, run_conf.project_root + '/log/' + sample_name , run_conf.project_root + '/script/'+ sample_name)

    for input_file in input_files:
        os.unlink(input_file)
        os.unlink(input_file + ".bai")


# identify mutations
@follows( markdup )
@follows( link_import_bam )
@subdivide(markdup_bam_list, formatter(), "{subpath[0][2]}/mutation/{subdir[0][0]}/{subdir[0][0]}.genomon_mutation.result.filt.txt", "{subpath[0][2]}/mutation/{subdir[0][0]}")
def identify_mutations(input_file, output_file, output_dir):

    sample_name = os.path.basename(output_dir)

    active_inhouse_normal_flag = False
    if genomon_conf.has_option("annotation", "active_inhouse_normal_flag"):
        active_inhouse_normal_flag = genomon_conf.get("annotation", "active_inhouse_normal_flag")

    inhouse_normal_tabix_db = ""
    if genomon_conf.has_option("REFERENCE", "inhouse_normal_tabix_db"):
        inhouse_normal_tabix_db = genomon_conf.get("REFERENCE", "inhouse_normal_tabix_db")

    active_inhouse_tumor_flag = False
    if genomon_conf.has_option("annotation", "active_inhouse_tumor_flag"):
        active_inhouse_tumor_flag = genomon_conf.get("annotation", "active_inhouse_tumor_flag")

    inhouse_tumor_tabix_db = ""
    if genomon_conf.has_option("REFERENCE", "inhouse_tumor_tabix_db"):
        inhouse_tumor_tabix_db = genomon_conf.get("REFERENCE", "inhouse_tumor_tabix_db")

    active_HGMD_flag = False
    if genomon_conf.has_option("annotation", "active_HGMD_flag"):
        active_HGMD_flag = genomon_conf.get("annotation", "active_HGMD_flag")
        
    HGMD_tabix_db = ""
    if genomon_conf.has_option("REFERENCE", "HGMD_tabix_db"):
        HGMD_tabix_db = genomon_conf.get("REFERENCE", "HGMD_tabix_db")

    arguments = {
        # fisher mutation
        "fisher": genomon_conf.get("SOFTWARE", "fisher"),
        "fisher_pair_params": genomon_conf.get("fisher_mutation_call", "pair_params"),
        "fisher_single_params": genomon_conf.get("fisher_mutation_call", "single_params"),
        # realignment filter
        "mutfilter": genomon_conf.get("SOFTWARE", "mutfilter"),
        "realignment_params": genomon_conf.get("realignment_filter","params"),
        # indel filter
        "indel_params": genomon_conf.get("indel_filter", "params"),
        # breakpoint filter
        "breakpoint_params": genomon_conf.get("breakpoint_filter","params"),
        # simplerepeat filter
        "simple_repeat_db":genomon_conf.get("REFERENCE", "simple_repeat_tabix_db"),
        # EB filter
        "EBFilter": genomon_conf.get("SOFTWARE", "ebfilter"),
        "eb_map_quality": genomon_conf.get("eb_filter","map_quality"),
        "eb_base_quality": genomon_conf.get("eb_filter","base_quality"),
        "filter_flags": genomon_conf.get("eb_filter","filter_flags"),
        "control_bam_list": input_file[2],
        # hotspot mutation caller
        "hotspot": genomon_conf.get("SOFTWARE","hotspot"),
        "hotspot_database":genomon_conf.get("REFERENCE","hotspot_db"),
        "active_hotspot_flag":genomon_conf.get("hotspot","active_hotspot_flag"),
        "hotspot_params": genomon_conf.get("hotspot","params"),
        "mutil": genomon_conf.get("SOFTWARE", "mutil"),
        # original_annotations
        "mutanno": genomon_conf.get("SOFTWARE", "mutanno"),
        "active_inhouse_normal_flag": active_inhouse_normal_flag,
        "inhouse_normal_database":inhouse_normal_tabix_db,
        "active_inhouse_tumor_flag": active_inhouse_tumor_flag,
        "inhouse_tumor_database":inhouse_tumor_tabix_db,
        "active_HGVD_2013_flag": genomon_conf.get("annotation", "active_HGVD_2013_flag"),
        "HGVD_2013_database":genomon_conf.get("REFERENCE", "HGVD_2013_tabix_db"),
        "active_HGVD_2016_flag": genomon_conf.get("annotation", "active_HGVD_2016_flag"),
        "HGVD_2016_database":genomon_conf.get("REFERENCE", "HGVD_2016_tabix_db"),
        "active_ExAC_flag": genomon_conf.get("annotation", "active_ExAC_flag"),
        "ExAC_database":genomon_conf.get("REFERENCE", "ExAC_tabix_db"),
        "active_HGMD_flag": active_HGMD_flag,
        "HGMD_database": HGMD_tabix_db,
        # annovar
        "active_annovar_flag": genomon_conf.get("annotation", "active_annovar_flag"),
        "annovar": genomon_conf.get("SOFTWARE", "annovar"),
        "annovar_database": genomon_conf.get("annotation", "annovar_database"),
        "table_annovar_params": genomon_conf.get("annotation", "table_annovar_params"),
        "annovar_buildver": genomon_conf.get("annotation", "annovar_buildver"),
        # commmon
        "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
        "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
        "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
        "ref_fa":genomon_conf.get("REFERENCE", "ref_fasta"),
        "interval_list": genomon_conf.get("REFERENCE", "interval_list"),
        "disease_bam": input_file[0],
        "control_bam": input_file[1],
        "out_prefix": output_dir + '/' + sample_name,
        "samtools": genomon_conf.get("SOFTWARE", "samtools"),
        "blat": genomon_conf.get("SOFTWARE", "blat")}

    interval_list = genomon_conf.get("REFERENCE", "interval_list")
    max_task_id = sum(1 for line in open(interval_list))

    mutation_call.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name, max_task_id)
    
    arguments = {
        "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
        "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
        "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
        "control_bam": input_file[1],
        "control_bam_list": input_file[2],
        "active_annovar_flag": genomon_conf.get("annotation", "active_annovar_flag"),
        "annovar_buildver": genomon_conf.get("annotation", "annovar_buildver"),
        "active_HGVD_2013_flag": genomon_conf.get("annotation", "active_HGVD_2013_flag"),
        "active_HGVD_2016_flag": genomon_conf.get("annotation", "active_HGVD_2016_flag"),
        "active_ExAC_flag": genomon_conf.get("annotation", "active_ExAC_flag"),
        "active_HGMD_flag": active_HGMD_flag,
        "active_inhouse_normal_flag": active_inhouse_normal_flag,
        "active_inhouse_tumor_flag": active_inhouse_tumor_flag,
        "filecount": max_task_id,
        "mutil": genomon_conf.get("SOFTWARE", "mutil"),
        "pair_params": genomon_conf.get("mutation_util","pair_params"),
        "single_params": genomon_conf.get("mutation_util","single_params"),
        "active_hotspot_flag":genomon_conf.get("hotspot","active_hotspot_flag"),
        "hotspot_database":genomon_conf.get("REFERENCE","hotspot_db"),
        "meta_info_em": get_meta_info(["fisher", "mutfilter", "ebfilter", "mutil", "mutanno"]),
        "meta_info_m": get_meta_info(["fisher", "mutfilter", "mutil", "mutanno"]),
        "meta_info_ema": get_meta_info(["fisher", "mutfilter", "ebfilter", "mutil", "mutanno", "hotspot"]),
        "meta_info_ma": get_meta_info(["fisher", "mutfilter", "mutil", "mutanno", "hotspot"]),
        "out_prefix": output_dir + '/' + sample_name}

    mutation_merge.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name)

    annovar_buildver = genomon_conf.get("annotation", "annovar_buildver"),
    for task_id in range(1,(max_task_id + 1)):
        input_file = output_dir+'/'+sample_name+'_mutations_candidate.'+str(task_id)+'.'+annovar_buildver[0]+'_multianno.txt'
        os.unlink(input_file)

    for task_id in range(1,(max_task_id + 1)):
        if os.path.exists(output_dir+'/'+sample_name+'.fisher_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.fisher_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.hotspot_mutations.'+str(task_id)+'.txt'):
            os.unlink(output_dir+'/'+sample_name+'.hotspot_mutations.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.fisher_hotspot_mutations.'+str(task_id)+'.txt'):
           os.unlink(output_dir+'/'+sample_name+'.fisher_hotspot_mutations.'+str(task_id)+'.txt')
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
        if os.path.exists(output_dir+'/'+sample_name+'.inhouse_normal.'+str(task_id)+'.txt'):
           os.unlink(output_dir+'/'+sample_name+'.inhouse_normal.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.inhouse_tumor.'+str(task_id)+'.txt'):
           os.unlink(output_dir+'/'+sample_name+'.inhouse_tumor.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.HGVD_2013.'+str(task_id)+'.txt'):
           os.unlink(output_dir+'/'+sample_name+'.HGVD_2013.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.HGVD_2016.'+str(task_id)+'.txt'):
           os.unlink(output_dir+'/'+sample_name+'.HGVD_2016.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.ExAC.'+str(task_id)+'.txt'):
           os.unlink(output_dir+'/'+sample_name+'.ExAC.'+str(task_id)+'.txt')
        if os.path.exists(output_dir+'/'+sample_name+'.HGMD.'+str(task_id)+'.txt'):
           os.unlink(output_dir+'/'+sample_name+'.HGMD.'+str(task_id)+'.txt')

# parse SV 
@follows( link_import_bam )
@follows( markdup )
@transform(parse_sv_bam_list, formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.junction.clustered.bedpe.gz")
def parse_sv(input_file, output_file):

    dir_name = os.path.dirname(output_file)
    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    sample_name = os.path.basename(dir_name)

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "input_bam": input_file,
                 "output_prefix": output_file.replace(".junction.clustered.bedpe.gz", ""),
                 "param": genomon_conf.get("sv_parse", "params"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "htslib": genomon_conf.get("SOFTWARE", "htslib")}

    sv_parse.task_exec(arguments, run_conf.project_root + '/log/' + sample_name , run_conf.project_root + '/script/' + sample_name)


# merge SV
@follows( parse_sv )
@transform(merge_bedpe_list, formatter(".+/(?P<NAME>.+).control_info.txt"), "{subpath[0][2]}/sv/non_matched_control_panel/{NAME[0]}.merged.junction.control.bedpe.gz")
def merge_sv(input_files,  output_file):

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "control_info": input_files[0],
                 "merge_output_file": output_file,
                 "param": genomon_conf.get("sv_merge", "params"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "htslib": genomon_conf.get("SOFTWARE", "htslib")}

    sv_merge.task_exec(arguments, run_conf.project_root + '/log/sv_merge', run_conf.project_root + '/script/sv_merge')


# filt SV
@follows( merge_sv )
@transform(filt_bedpe_list, formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.genomonSV.result.filt.txt")
def filt_sv(input_files,  output_file):

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    #sample_yaml = run_conf.project_root + "/sv/config/" + sample_name + ".yaml"

    filt_param = ""

    for complist in sample_conf.sv_detection:
        if sample_name == complist[0]:

            if complist[1] != None:
                filt_param = filt_param + " --matched_control_bam " + run_conf.project_root + "/bam/" + complist[1] + '/' + complist[1] + ".markdup.bam"

            if complist[2] != None:
                filt_param = filt_param + " --non_matched_control_junction " + run_conf.project_root +"/sv/non_matched_control_panel/"+ complist[2] +".merged.junction.control.bedpe.gz"
                if complist[1] != None:
                    filt_param = filt_param + " --matched_control_label " + complist[1]

            break

    filt_param = filt_param.lstrip(' ') + ' ' + genomon_conf.get("sv_filt", "params")

    arguments = {"genomon_sv": genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "input_bam": run_conf.project_root + "/bam/" + sample_name + '/' + sample_name + ".markdup.bam",
                 "output_prefix": run_conf.project_root + "/sv/" + sample_name + '/' + sample_name,
                 "reference_genome": genomon_conf.get("REFERENCE", "ref_fasta"),
                 "param": filt_param,
                 "meta_info": get_meta_info(["genomon_sv", "sv_utils"]),
                 "sv_utils": genomon_conf.get("SOFTWARE", "sv_utils"),
                 "sv_utils_param": genomon_conf.get("sv_filt", "sv_utils_params"),
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "htslib": genomon_conf.get("SOFTWARE", "htslib"),
                 "blat": genomon_conf.get("SOFTWARE", "blat")}

    sv_filt.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name)


# qc
@follows( link_import_bam )
@follows( markdup )
@follows( filt_sv )
@follows( identify_mutations )
@transform(qc_bamstats_list, formatter(), "{subpath[0][2]}/qc/{subdir[0][0]}/{subdir[0][0]}.bamstats")
def bam_stats(input_file, output_file):
    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    
    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_qc": genomon_conf.get("SOFTWARE", "genomon_qc"),
                 "bamstats": genomon_conf.get("SOFTWARE", "bamstats"),
                 "perl": genomon_conf.get("ENV", "PERL"),
                 "perl5lib": genomon_conf.get("ENV", "PERL5LIB"),
                 "input_file": input_file,
                 "output_file": output_file}
    
    r_qc_bamstats.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name)


@follows( link_import_bam )
@follows( markdup )
@follows( filt_sv )
@follows( identify_mutations )
@transform(qc_coverage_list, formatter(), "{subpath[0][2]}/qc/{subdir[0][0]}/{subdir[0][0]}.coverage")
def coverage(input_file, output_file):
    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    
    data_type = "exome"
    if genomon_conf.get("qc_coverage", "wgs_flag") == "True":
        data_type = "wgs"

    arguments = {"data_type": data_type,
                 "pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_qc": genomon_conf.get("SOFTWARE", "genomon_qc"),
                 "coverage_text": genomon_conf.get("qc_coverage", "coverage"),
                 "i_bed_lines": genomon_conf.get("qc_coverage", "wgs_i_bed_lines"),
                 "i_bed_width": genomon_conf.get("qc_coverage", "wgs_i_bed_width"),
                 "incl_bed_width":genomon_conf.get("qc_coverage", "wgs_incl_bed_width"),
                 "genome_size_file": genomon_conf.get("REFERENCE", "genome_size"),
                 "gaptxt": genomon_conf.get("REFERENCE", "gaptxt"),
                 "bait_file": genomon_conf.get("REFERENCE", "bait_file"),
                 "samtools_params": genomon_conf.get("qc_coverage", "samtools_params"),
                 "bedtools": genomon_conf.get("SOFTWARE", "bedtools"),
                 "samtools": genomon_conf.get("SOFTWARE", "samtools"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "input_file": input_file,
                 "output_file": output_file}

    r_qc_coverage.task_exec(arguments, run_conf.project_root + '/log/' + sample_name , run_conf.project_root + '/script/' + sample_name)

@follows( bam_stats )
@follows( coverage )
@collate(qc_merge_list, formatter(), "{subpath[0][2]}/qc/{subdir[0][0]}/{subdir[0][0]}.genomonQC.result.txt")
def merge_qc(input_files, output_file):

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    
    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_qc": genomon_conf.get("SOFTWARE", "genomon_qc"),
                 "bamstats_file": input_files[0][0],
                 "coverage_file": input_files[0][1],
                 "output_file": output_file,
                 "meta": get_meta_info(["genomon_pipeline"]),
                 "fastq_line_num_file": run_conf.project_root +'/fastq/'+ sample_name +'/fastq_line_num.txt'}
    
    r_qc_merge.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name)

#####################
# post analysis stage
@active_if(genomon_conf.getboolean("post_analysis", "enable"))
@active_if(len(pa_inputs_mutation) > 0)
@follows(filt_sv)
@follows(identify_mutations)
@collate(pa_inputs_mutation, formatter(), pa_outputs_mutation["outputs"])
def post_analysis_mutation(input_files, output_file):
        
    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_pa":  genomon_conf.get("SOFTWARE", "genomon_pa"),
                 "mode": "mutation",
                 "genomon_root": run_conf.project_root,
                 "output_dir": run_conf.project_root + "/post_analysis/" + sample_conf_name,
                 "sample_sheet": os.path.abspath(run_conf.sample_conf_file),
                 "config_file": genomon_conf.get("post_analysis", "config_file"),
                 "samtools": genomon_conf.get("SOFTWARE", "samtools"),
                 "bedtools": genomon_conf.get("SOFTWARE", "bedtools"),
                 "input_file_case1": ",".join(pa_outputs_mutation["case1"]["samples"]),
                 "input_file_case2": ",".join(pa_outputs_mutation["case2"]["samples"]),
                 "input_file_case3": ",".join(pa_outputs_mutation["case3"]["samples"]),
                 "input_file_case4": ",".join(pa_outputs_mutation["case4"]["samples"]),
                }

    r_post_analysis.task_exec(arguments, run_conf.project_root + '/log/post_analysis', run_conf.project_root + '/script/post_analysis')
    
@active_if(genomon_conf.getboolean("post_analysis", "enable"))
@active_if(len(pa_inputs_sv) > 0)
@follows(filt_sv)
@follows(identify_mutations)
@collate(pa_inputs_sv, formatter(), pa_outputs_sv["outputs"])
def post_analysis_sv(input_files, output_file):

    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_pa":  genomon_conf.get("SOFTWARE", "genomon_pa"),
                 "mode": "sv",
                 "genomon_root": run_conf.project_root,
                 "output_dir": run_conf.project_root + "/post_analysis/" + sample_conf_name,
                 "sample_sheet": os.path.abspath(run_conf.sample_conf_file),
                 "config_file": genomon_conf.get("post_analysis", "config_file"),
                 "samtools": genomon_conf.get("SOFTWARE", "samtools"),
                 "bedtools": genomon_conf.get("SOFTWARE", "bedtools"),
                 "input_file_case1": ",".join(pa_outputs_sv["case1"]["samples"]),
                 "input_file_case2": ",".join(pa_outputs_sv["case2"]["samples"]),
                 "input_file_case3": ",".join(pa_outputs_sv["case3"]["samples"]),
                 "input_file_case4": ",".join(pa_outputs_sv["case4"]["samples"]),
                }
                 
    r_post_analysis.task_exec(arguments, run_conf.project_root + '/log/post_analysis', run_conf.project_root + '/script/post_analysis')

@active_if(genomon_conf.getboolean("post_analysis", "enable"))
@active_if(len(pa_inputs_qc) > 0)
@follows(merge_qc)
@collate(pa_inputs_qc, formatter(), pa_outputs_qc["outputs"])
def post_analysis_qc(input_files, output_file):

    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_pa":  genomon_conf.get("SOFTWARE", "genomon_pa"),
                 "mode": "qc",
                 "genomon_root": run_conf.project_root,
                 "output_dir": run_conf.project_root + "/post_analysis/" + sample_conf_name,
                 "sample_sheet": os.path.abspath(run_conf.sample_conf_file),
                 "config_file": genomon_conf.get("post_analysis", "config_file"),
                 "samtools": genomon_conf.get("SOFTWARE", "samtools"),
                 "bedtools": genomon_conf.get("SOFTWARE", "bedtools"),
                 "input_file_case1": ",".join(sample_conf.qc),
                 "input_file_case2": "",
                 "input_file_case3": "",
                 "input_file_case4": "",
                }
                 
    r_post_analysis.task_exec(arguments, run_conf.project_root + '/log/post_analysis', run_conf.project_root + '/script/post_analysis')
    
@active_if(genomon_conf.getboolean("post_analysis", "enable"))
@active_if(genomon_conf.getboolean("pmsignature_ind", "enable") or genomon_conf.getboolean("pmsignature_full", "enable"))
@active_if(len(pmsignature_inputs) > 0)
@follows(post_analysis_mutation)
@collate(pmsignature_inputs, formatter(), run_conf.project_root + '/pmsignature/' + sample_conf_name + "/mutation.cut.txt")
def pre_pmsignature(input_files, output_file):
        
    arguments = {"input_files" : " ".join(input_files),
                 "output_file" : run_conf.project_root + '/pmsignature/' + sample_conf_name + "/mutation.cut.txt"
                }

    r_pre_pmsignature.task_exec(arguments, run_conf.project_root + '/log/pmsignature', run_conf.project_root + '/script/pmsignature')
    
@active_if(genomon_conf.getboolean("post_analysis", "enable"))
@active_if(genomon_conf.getboolean("pmsignature_ind", "enable"))
@active_if(run_ind)
@follows(pre_pmsignature)
@transform(run_conf.project_root + '/pmsignature/' + sample_conf_name + "/mutation.cut.txt", formatter(), ind_outputs[0])
def pmsignature_ind(input_file, output_file):
    
    command = r_pmsignature_ind.ind_template.format(
                inputfile =  input_file,
                outputdir = run_conf.project_root + '/pmsignature/' + sample_conf_name,
                trdirflag = genomon_conf.get("pmsignature_ind", "trdirflag").upper(),
                trialnum = genomon_conf.getint("pmsignature_ind", "trialnum"),
                bs_genome = genomon_conf.get("pmsignature_ind", "bs_genome"),
                bgflag = genomon_conf.get("pmsignature_ind", "bgflag"),
                txdb_transcript = genomon_conf.get("pmsignature_ind", "txdb_transcript"),
                script_path = genomon_conf.get("SOFTWARE", "r_scripts"))
    
    sig_nums = range(genomon_conf.getint("pmsignature_ind", "signum_min"), genomon_conf.getint("pmsignature_ind", "signum_max") + 1)
    sig_num_text = ""
    for i in sig_nums: sig_num_text += "%d " % i

    arguments = {"r_path": genomon_conf.get("ENV", "R_PATH"),
                 "r_ld_library_path": genomon_conf.get("ENV", "R_LD_LIBRARY_PATH"),
                 "r_libs": genomon_conf.get("ENV", "R_LIBS"),
                 "command": command,
                 "sig_list": sig_num_text
                }
    max_task_id = len(sig_nums)
    r_pmsignature_ind.task_exec(arguments, run_conf.project_root + '/log/pmsignature', run_conf.project_root + '/script/pmsignature', max_task_id)

@active_if(genomon_conf.getboolean("post_analysis", "enable"))
@active_if(genomon_conf.getboolean("pmsignature_full", "enable"))
@active_if(run_full)
@follows(pre_pmsignature)
@transform(run_conf.project_root + '/pmsignature/' + sample_conf_name + "/mutation.cut.txt", formatter(), full_outputs[0])
def pmsignature_full(input_file, output_file):
    
    command = r_pmsignature_full.full_template.format(
                inputfile = input_file,
                outputdir = run_conf.project_root + '/pmsignature/' + sample_conf_name,
                trdirflag = genomon_conf.get("pmsignature_full", "trdirflag").upper(),
                trialnum = genomon_conf.getint("pmsignature_full", "trialnum"),
                bgflag = genomon_conf.get("pmsignature_full", "bgflag"),
                bs_genome = genomon_conf.get("pmsignature_full", "bs_genome"),
                txdb_transcript = genomon_conf.get("pmsignature_full", "txdb_transcript"),
                script_path = genomon_conf.get("SOFTWARE", "r_scripts"))
    
    sig_nums = range(genomon_conf.getint("pmsignature_full", "signum_min"), genomon_conf.getint("pmsignature_full", "signum_max") + 1)
    sig_num_text = ""
    for i in sig_nums: sig_num_text += "%d " % i

    arguments = {"r_path": genomon_conf.get("ENV", "R_PATH"),
                 "r_ld_library_path": genomon_conf.get("ENV", "R_LD_LIBRARY_PATH"),
                 "r_libs": genomon_conf.get("ENV", "R_LIBS"),
                 "command": command,
                 "sig_list": sig_num_text
                }
    max_task_id = len(sig_nums)
    r_pmsignature_full.task_exec(arguments, run_conf.project_root + '/log/pmsignature', run_conf.project_root + '/script/pmsignature', max_task_id)

@active_if(genomon_conf.getboolean("post_analysis", "enable"))
@active_if(genomon_conf.getboolean("paplot", "enable"))
@active_if(len(paplot_inputs) > 0)
@follows(post_analysis_sv)
@follows(post_analysis_qc)
@follows(pmsignature_ind)
@follows(pmsignature_full)
@collate(paplot_inputs, formatter(), run_conf.project_root + '/paplot/' + sample_conf_name + '/index.html')
def paplot(input_file, output_file):
    if not os.path.exists(paplot_output) and os.path.exists(run_conf.project_root + '/paplot/' + sample_conf_name + '/.meta.json'):
        os.unlink(run_conf.project_root + '/paplot/' + sample_conf_name + '/.meta.json')

    command = ""
    if len(paplot_inputs_qc) > 0:
        command += r_paplot.qc_template.format(
                    paplot = genomon_conf.get("SOFTWARE", "paplot"),
                    inputs = ",".join(paplot_inputs_qc),
                    output_dir = run_conf.project_root + "/paplot/" + sample_conf_name,
                    title = genomon_conf.get("paplot", "title"),
                    config_file = genomon_conf.get("paplot", "config_file"))
                        
    if len(paplot_inputs_sv) > 0:
        command += r_paplot.sv_template.format(
                    paplot = genomon_conf.get("SOFTWARE", "paplot"),
                    inputs = ",".join(paplot_inputs_sv),
                    output_dir = run_conf.project_root + "/paplot/" + sample_conf_name,
                    title = genomon_conf.get("paplot", "title"),
                    config_file = genomon_conf.get("paplot", "config_file"))
                        
    if len(paplot_inputs_mutation) > 0:
        command += r_paplot.mutation_template.format(
                    paplot = genomon_conf.get("SOFTWARE", "paplot"),
                    inputs = ",".join(paplot_inputs_mutation),
                    output_dir = run_conf.project_root + "/paplot/" + sample_conf_name,
                    title = genomon_conf.get("paplot", "title"),
                    config_file = genomon_conf.get("paplot", "config_file"),
                    annovar = genomon_conf.getboolean("annotation", "active_annovar_flag"))
    
    if genomon_conf.getboolean("pmsignature_ind", "enable"):
        for i in range(len(paplot_inputs_ind)):
            command += r_paplot.ind_template.format(
                        paplot = genomon_conf.get("SOFTWARE", "paplot"),
                        input = paplot_inputs_ind[i],
                        output_dir = run_conf.project_root + "/paplot/" + sample_conf_name,
                        title = genomon_conf.get("paplot", "title"),
                        config_file = genomon_conf.get("paplot", "config_file"))
    
    if genomon_conf.getboolean("pmsignature_full", "enable"):
        for i in range(len(paplot_inputs_full)):
            command += r_paplot.full_template.format(
                        paplot =genomon_conf.get("SOFTWARE", "paplot"),
                        input = paplot_inputs_full[i],
                        output_dir = run_conf.project_root + "/paplot/" + sample_conf_name,
                        title = genomon_conf.get("paplot", "title"),
                        config_file = genomon_conf.get("paplot", "config_file"))
    
    remark = genomon_conf.get("paplot", "remarks")
    remark += "<ul>"
    
    for item in genomon_conf.get("paplot", "software").split(","):
        key = item.split(":")[0].strip(" ").rstrip(" ")
        name = item.split(":")[1].strip(" ").rstrip(" ")
        try:
            version = get_version(key).split("-")
        except Exception:
            print ("[WARNING] paplot: %s is not defined." % (key))
            continue
        
        remark += "<li>" + name + " " + version[-1] + "</li>"

    remark += "</ul>"
    
    command += r_paplot.index_template.format(
                        paplot = genomon_conf.get("SOFTWARE", "paplot"),
                        output_dir = run_conf.project_root + "/paplot/" + sample_conf_name,
                        remarks = remark,
                        config_file = genomon_conf.get("paplot", "config_file"))

    arguments = {"pythonhome": genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": genomon_conf.get("ENV", "PYTHONPATH"),
                 "paplot":  genomon_conf.get("SOFTWARE", "paplot"),
                 "command": command
                }
                 
    r_paplot.task_exec(arguments, run_conf.project_root + '/log/paplot', run_conf.project_root + '/script/paplot')

