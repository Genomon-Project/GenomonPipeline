import os
import shutil
import glob
import ruffus
import genomon_pipeline.config.run_conf              as rc
import genomon_pipeline.config.genomon_conf          as gc
import genomon_pipeline.config.sample_conf           as sc
import genomon_pipeline.dna_resource.bamtofastq      as dr_bamtofastq
import genomon_pipeline.dna_resource.fastq_splitter  as dr_fastq_splitter 
import genomon_pipeline.dna_resource.bwa_align       as dr_bwa_align
import genomon_pipeline.dna_resource.markduplicates  as dr_markduplicates 
import genomon_pipeline.dna_resource.mutation_call   as dr_mutation_call
import genomon_pipeline.dna_resource.mutation_merge  as dr_mutation_merge
import genomon_pipeline.dna_resource.sv_parse        as dr_sv_parse
import genomon_pipeline.dna_resource.sv_merge        as dr_sv_merge
import genomon_pipeline.dna_resource.sv_filt         as dr_sv_filt
import genomon_pipeline.dna_resource.qc_bamstats     as dr_qc_bamstats
import genomon_pipeline.dna_resource.qc_coverage     as dr_qc_coverage
import genomon_pipeline.dna_resource.qc_merge        as dr_qc_merge
import genomon_pipeline.dna_resource.post_analysis   as dr_post_analysis
import genomon_pipeline.dna_resource.pre_pmsignature as dr_pre_pmsignature
import genomon_pipeline.dna_resource.pmsignature     as dr_pmsignature
import genomon_pipeline.dna_resource.paplot          as dr_paplot

# set task classes
bamtofastq = dr_bamtofastq.Bam2Fastq(gc.genomon_conf.get("bam2fastq", "qsub_option"), rc.run_conf.drmaa)
fastq_splitter = dr_fastq_splitter.Fastq_splitter(gc.genomon_conf.get("split_fastq", "qsub_option"), rc.run_conf.drmaa)
bwa_align = dr_bwa_align.Bwa_align(gc.genomon_conf.get("bwa_mem", "qsub_option"), rc.run_conf.drmaa)
markduplicates = dr_markduplicates.Markduplicates(gc.genomon_conf.get("markduplicates", "qsub_option"), rc.run_conf.drmaa)
mutation_call = dr_mutation_call.Mutation_call(gc.genomon_conf.get("mutation_call", "qsub_option"), rc.run_conf.drmaa)
mutation_merge = dr_mutation_merge.Mutation_merge(gc.genomon_conf.get("mutation_merge", "qsub_option"), rc.run_conf.drmaa)
sv_parse = dr_sv_parse.SV_parse(gc.genomon_conf.get("sv_parse", "qsub_option"), rc.run_conf.drmaa)
sv_merge = dr_sv_merge.SV_merge(gc.genomon_conf.get("sv_merge", "qsub_option"), rc.run_conf.drmaa)
sv_filt = dr_sv_filt.SV_filt(gc.genomon_conf.get("sv_filt", "qsub_option"), rc.run_conf.drmaa)
r_qc_bamstats = dr_qc_bamstats.Res_QC_Bamstats(gc.genomon_conf.get("qc_bamstats", "qsub_option"), rc.run_conf.drmaa)
r_qc_coverage = dr_qc_coverage.Res_QC_Coverage(gc.genomon_conf.get("qc_coverage", "qsub_option"), rc.run_conf.drmaa)
r_qc_merge = dr_qc_merge.Res_QC_Merge(gc.genomon_conf.get("qc_merge", "qsub_option"), rc.run_conf.drmaa)
r_paplot = dr_paplot.Res_PA_Plot(gc.genomon_conf.get("paplot", "qsub_option"), rc.run_conf.drmaa)
r_post_analysis = dr_post_analysis.Res_PostAnalysis(gc.genomon_conf.get("post_analysis", "qsub_option"), rc.run_conf.drmaa)
r_pre_pmsignature = dr_pre_pmsignature.Res_PrePmsignature(gc.genomon_conf.get("pre_pmsignature", "qsub_option"), rc.run_conf.drmaa)
r_pmsignature_ind = dr_pmsignature.Res_Pmsignature(gc.genomon_conf.get("pmsignature_ind", "qsub_option"), rc.run_conf.drmaa)
r_pmsignature_full = dr_pmsignature.Res_Pmsignature(gc.genomon_conf.get("pmsignature_full", "qsub_option"), rc.run_conf.drmaa)

_debug = False
if gc.genomon_conf.has_section("develop"):
    if gc.genomon_conf.has_option("develop", "debug") == True:
        _debug = gc.genomon_conf.getboolean("develop", "debug")

# generate output list of 'linked fastq'
linked_fastq_list = []
for sample in sc.sample_conf.fastq:
    if os.path.exists(rc.run_conf.project_root + '/bam/' + sample + '/1.sorted.bam'): continue
    if os.path.exists(rc.run_conf.project_root + '/bam/' + sample + '/' + sample + '.markdup.bam'): continue

    link_fastq_arr1 = []
    link_fastq_arr2 = []
    for (count, fastq_file) in enumerate(sc.sample_conf.fastq[sample][0]):
        fastq_prefix, ext = os.path.splitext(fastq_file)
        link_fastq_arr1.append(rc.run_conf.project_root + '/fastq/' + sample + '/' + str(count+1) + '_1' + ext)
        link_fastq_arr2.append(rc.run_conf.project_root + '/fastq/' + sample + '/' + str(count+1) + '_2' + ext)
    linked_fastq_list.append([link_fastq_arr1,link_fastq_arr2])

# generate output list of 'bam2fastq'
bam2fastq_output_list = []
for sample in sc.sample_conf.bam_tofastq:
    if os.path.exists(rc.run_conf.project_root + '/bam/' + sample + '/1.sorted.bam'): continue
    if os.path.exists(rc.run_conf.project_root + '/bam/' + sample + '/' + sample + '.markdup.bam'): continue
    bam2fastq_arr1 = []
    bam2fastq_arr2 = []
    bam2fastq_arr1.append(rc.run_conf.project_root + '/fastq/' + sample + '/1_1.fastq')
    bam2fastq_arr2.append(rc.run_conf.project_root + '/fastq/' + sample + '/1_2.fastq')
    bam2fastq_output_list.append([bam2fastq_arr1,bam2fastq_arr2])

# generate input list of 'mutation call'
markdup_bam_list = []
merge_mutation_list = []
for complist in sc.sample_conf.mutation_call:
     if os.path.exists(rc.run_conf.project_root + '/mutation/' + complist[0] + '/' + complist[0] + '.genomon_mutation.result.filt.txt'): continue
     tumor_bam  = rc.run_conf.project_root + '/bam/' + complist[0] + '/' + complist[0] + '.markdup.bam'
     normal_bam = rc.run_conf.project_root + '/bam/' + complist[1] + '/' + complist[1] + '.markdup.bam' if complist[1] != None else None
     panel = rc.run_conf.project_root + '/mutation/control_panel/' + complist[2] + ".control_panel.txt" if complist[2] != None else None
     markdup_bam_list.append([tumor_bam, normal_bam, panel])


# generate input list of 'SV parse'
parse_sv_bam_list = []
all_target_bams = []
unique_bams = []
for complist in sc.sample_conf.sv_detection:
    tumor_sample = complist[0]
    if tumor_sample != None:
        all_target_bams.append(rc.run_conf.project_root + '/bam/' + tumor_sample + '/' + tumor_sample + '.markdup.bam')
    normal_sample = complist[1]
    if normal_sample != None:
        all_target_bams.append(rc.run_conf.project_root + '/bam/' + normal_sample + '/' + normal_sample + '.markdup.bam')
    panel_name = complist[2]
    if panel_name != None:
        for panel_sample in sc.sample_conf.control_panel[panel_name]:
            all_target_bams.append(rc.run_conf.project_root + '/bam/' + panel_sample + '/' + panel_sample + '.markdup.bam')
    unique_bams = list(set(all_target_bams))       
    
for bam in unique_bams:
    dir_name = os.path.dirname(bam)
    sample_name = os.path.basename(dir_name)
    if os.path.exists(rc.run_conf.project_root + '/sv/' + sample_name + '/' + sample_name + '.junction.clustered.bedpe.gz') and os.path.exists(rc.run_conf.project_root + '/sv/' + sample_name + '/' + sample_name + '.junction.clustered.bedpe.gz.tbi'): continue
    parse_sv_bam_list.append(bam)

# generate input list of 'SV merge'
unique_complist = []
merge_bedpe_list = []
for complist in sc.sample_conf.sv_detection:
    control_panel_name = complist[2]
    if control_panel_name != None and control_panel_name not in unique_complist:
        unique_complist.append(control_panel_name)

for control_panel_name in unique_complist:
    if os.path.exists(rc.run_conf.project_root + '/sv/non_matched_control_panel/' + control_panel_name + '.merged.junction.control.bedpe.gz') and os.path.exists(rc.run_conf.project_root + '/sv/non_matched_control_panel/' + control_panel_name + '.merged.junction.control.bedpe.gz.tbi'): continue
    tmp_list = []
    tmp_list.append(rc.run_conf.project_root + '/sv/control_panel/' + control_panel_name + ".control_info.txt")
    for sample in sc.sample_conf.control_panel[control_panel_name]:
        tmp_list.append(rc.run_conf.project_root+ "/sv/"+ sample +"/"+ sample +".junction.clustered.bedpe.gz")
    merge_bedpe_list.append(tmp_list)

# generate input list of 'SV filt'
filt_bedpe_list = []
for complist in sc.sample_conf.sv_detection:
    if os.path.exists(rc.run_conf.project_root + '/sv/' + complist[0] +'/'+ complist[0] +'.genomonSV.result.filt.txt'): continue
    filt_bedpe_list.append(rc.run_conf.project_root+ "/sv/"+ complist[0] +"/"+ complist[0] +".junction.clustered.bedpe.gz")

# generate input list of 'qc'
qc_bamstats_list = []
qc_coverage_list = []
qc_merge_list = []
for sample in sc.sample_conf.qc:
    if os.path.exists(rc.run_conf.project_root + '/qc/' + sample + '/' + sample + '.genomonQC.result.txt'): continue
    qc_merge_list.append(
        [rc.run_conf.project_root + '/qc/' + sample + '/' + sample + '.bamstats',
         rc.run_conf.project_root + '/qc/' + sample + '/' + sample + '.coverage'])
    if not os.path.exists(rc.run_conf.project_root + '/qc/' + sample + '/' + sample + '.bamstats'):
        qc_bamstats_list.append(rc.run_conf.project_root + '/bam/' + sample +'/'+ sample +'.markdup.bam')
    if not os.path.exists(rc.run_conf.project_root + '/qc/' + sample + '/' + sample + '.coverage'):
        qc_coverage_list.append(rc.run_conf.project_root + '/bam/' + sample +'/'+ sample +'.markdup.bam')

### 
# input/output lists for post-analysis
###
genomon_conf_name, genomon_conf_ext = os.path.splitext(os.path.basename(rc.run_conf.genomon_conf_file))
sample_conf_name, sample_conf_ext = os.path.splitext(os.path.basename(rc.run_conf.sample_conf_file))

# generate input list of 'post analysis for mutation'
pa_outputs_mutation = r_post_analysis.output_files("mutation", sc.sample_conf.mutation_call, rc.run_conf.project_root, sample_conf_name, gc.genomon_conf)

pa_inputs_mutation = []
if pa_outputs_mutation["run_pa"] == True:
    for complist in sc.sample_conf.mutation_call:
        pa_inputs_mutation.append(rc.run_conf.project_root + '/mutation/' + complist[0] +'/'+ complist[0] +'.genomon_mutation.result.filt.txt')
        
# generate input list of 'post analysis for SV'
pa_outputs_sv = r_post_analysis.output_files("sv", sc.sample_conf.sv_detection, rc.run_conf.project_root, sample_conf_name, gc.genomon_conf)

pa_inputs_sv = []
if pa_outputs_sv["run_pa"] == True:
    for complist in sc.sample_conf.sv_detection:
        pa_inputs_sv.append(rc.run_conf.project_root + '/sv/' + complist[0] +'/'+ complist[0] +'.genomonSV.result.filt.txt')
        
# generate input list of 'post analysis for qc'
pa_outputs_qc = r_post_analysis.output_files("qc", sc.sample_conf.qc, rc.run_conf.project_root, sample_conf_name, gc.genomon_conf)

pa_inputs_qc = []
if pa_outputs_qc["run_pa"] == True:
    for sample in sc.sample_conf.qc:
        pa_inputs_qc.append(rc.run_conf.project_root + '/qc/' + sample + '/' + sample + '.genomonQC.result.txt')

### 
# input/output lists for paplot
###
paplot_output = rc.run_conf.project_root + '/paplot/' + sample_conf_name + '/index.html'

## mutation
use_mutations = []
if pa_outputs_mutation["case1"]["output_filt"] != "":
    use_mutations.append(pa_outputs_mutation["case1"]["output_filt"])
if pa_outputs_mutation["case2"]["output_filt"] != "" and gc.genomon_conf.getboolean("paplot", "include_unpanel"):
    use_mutations.append(pa_outputs_mutation["case2"]["output_filt"])
if pa_outputs_mutation["case3"]["output_filt"] != "" and gc.genomon_conf.getboolean("paplot", "include_unpair"):
    use_mutations.append(pa_outputs_mutation["case3"]["output_filt"])
if pa_outputs_mutation["case4"]["output_filt"] != "" and gc.genomon_conf.getboolean("paplot", "include_unpanel") and gc.genomon_conf.getboolean("paplot", "include_unpair"):
    use_mutations.append(pa_outputs_mutation["case4"]["output_filt"])

paplot_inputs_mutation = []
if os.path.exists(paplot_output) == False or pa_outputs_mutation["run_pa"] == True:
    paplot_inputs_mutation.extend(use_mutations)

## pmsignature
# ind
ind_outputs = []
ind_exists = True
for i in range(gc.genomon_conf.getint("pmsignature_ind", "signum_min"), gc.genomon_conf.getint("pmsignature_ind", "signum_max") + 1):
    fname = rc.run_conf.project_root + '/pmsignature/' + sample_conf_name + '/pmsignature.ind.result.%d.json' % i
    ind_outputs.append(fname)
    if not os.path.exists(fname): ind_exists = False
        
run_ind = False
paplot_inputs_ind = []
if len(sc.sample_conf.mutation_call) > 0 and gc.genomon_conf.getboolean("pmsignature_ind", "enable") and len(use_mutations) > 0:
    if ind_exists == False: run_ind = True
    elif pa_outputs_mutation["run_pa"] == True: run_ind = True
    elif not os.path.exists(rc.run_conf.project_root + '/pmsignature/' + sample_conf_name + '/mutation.cut.txt'): run_ind = True
    if os.path.exists(paplot_output) == False or run_ind == True:    
        paplot_inputs_ind.extend(ind_outputs)
    
# full
full_outputs = []
full_exists = True
for i in range(gc.genomon_conf.getint("pmsignature_full", "signum_min"), gc.genomon_conf.getint("pmsignature_full", "signum_max") + 1):
    fname = rc.run_conf.project_root + '/pmsignature/' + sample_conf_name + '/pmsignature.full.result.%d.json' % i
    full_outputs.append(fname)
    if not os.path.exists(fname): full_exists = False
        
run_full = False
paplot_inputs_full = []
if len(sc.sample_conf.mutation_call) > 0 and gc.genomon_conf.getboolean("pmsignature_full", "enable") and len(use_mutations) > 0:
    if full_exists == False: run_full = True
    elif pa_outputs_mutation["run_pa"] == True: run_full = True
    elif not os.path.exists(rc.run_conf.project_root + '/pmsignature/' + sample_conf_name + '/mutation.cut.txt'): run_full = True
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
    if pa_outputs_sv["case2"]["output_filt"] != "" and gc.genomon_conf.getboolean("paplot", "include_unpanel"):
        paplot_inputs_sv.append(pa_outputs_sv["case2"]["output_filt"])
    if pa_outputs_sv["case3"]["output_filt"] != "" and gc.genomon_conf.getboolean("paplot", "include_unpair"):
        paplot_inputs_sv.append(pa_outputs_sv["case3"]["output_filt"])
    if pa_outputs_sv["case4"]["output_filt"] != "" and gc.genomon_conf.getboolean("paplot", "include_unpanel") and gc.genomon_conf.getboolean("paplot", "include_unpair"):
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
if not os.path.isdir(rc.run_conf.project_root): os.mkdir(rc.run_conf.project_root)
if not os.path.isdir(rc.run_conf.project_root + '/script'): os.mkdir(rc.run_conf.project_root + '/script')
if not os.path.isdir(rc.run_conf.project_root + '/script/sv_merge'): os.mkdir(rc.run_conf.project_root + '/script/sv_merge')
if not os.path.isdir(rc.run_conf.project_root + '/log'): os.mkdir(rc.run_conf.project_root + '/log')
if not os.path.isdir(rc.run_conf.project_root + '/log/sv_merge'): os.mkdir(rc.run_conf.project_root + '/log/sv_merge')
if not os.path.isdir(rc.run_conf.project_root + '/fastq'): os.mkdir(rc.run_conf.project_root + '/fastq')
if not os.path.isdir(rc.run_conf.project_root + '/bam'): os.mkdir(rc.run_conf.project_root + '/bam')
if not os.path.isdir(rc.run_conf.project_root + '/mutation'): os.mkdir(rc.run_conf.project_root + '/mutation')
if not os.path.isdir(rc.run_conf.project_root + '/mutation/control_panel'): os.mkdir(rc.run_conf.project_root + '/mutation/control_panel')
if not os.path.isdir(rc.run_conf.project_root + '/mutation/hotspot'): os.mkdir(rc.run_conf.project_root + '/mutation/hotspot')
if not os.path.isdir(rc.run_conf.project_root + '/sv'): os.mkdir(rc.run_conf.project_root + '/sv')
if not os.path.isdir(rc.run_conf.project_root + '/sv/non_matched_control_panel'): os.mkdir(rc.run_conf.project_root + '/sv/non_matched_control_panel')
if not os.path.isdir(rc.run_conf.project_root + '/sv/control_panel'): os.mkdir(rc.run_conf.project_root + '/sv/control_panel')
if not os.path.isdir(rc.run_conf.project_root + '/qc'): os.mkdir(rc.run_conf.project_root + '/qc')
for sample in sc.sample_conf.qc:
    if not os.path.isdir(rc.run_conf.project_root + '/qc/' + sample): os.mkdir(rc.run_conf.project_root + '/qc/' + sample)

if (gc.genomon_conf.getboolean("post_analysis", "enable") == True):
    if not os.path.exists(rc.run_conf.project_root + '/post_analysis'): os.mkdir(rc.run_conf.project_root + '/post_analysis')
    if not os.path.exists(rc.run_conf.project_root + '/post_analysis/' + sample_conf_name): os.mkdir(rc.run_conf.project_root + '/post_analysis/' + sample_conf_name)
    if not os.path.isdir(rc.run_conf.project_root + '/script/post_analysis'): os.mkdir(rc.run_conf.project_root + '/script/post_analysis')
    if not os.path.isdir(rc.run_conf.project_root + '/log/post_analysis'): os.mkdir(rc.run_conf.project_root + '/log/post_analysis')
    
    if (gc.genomon_conf.getboolean("paplot", "enable") == True):
        if not os.path.isdir(rc.run_conf.project_root + '/paplot/'): os.mkdir(rc.run_conf.project_root + '/paplot/')
        if not os.path.isdir(rc.run_conf.project_root + '/paplot/' + sample_conf_name): os.mkdir(rc.run_conf.project_root + '/paplot/' + sample_conf_name)
        if not os.path.isdir(rc.run_conf.project_root + '/script/paplot'): os.mkdir(rc.run_conf.project_root + '/script/paplot')
        if not os.path.isdir(rc.run_conf.project_root + '/log/paplot'): os.mkdir(rc.run_conf.project_root + '/log/paplot')

    if (gc.genomon_conf.getboolean("pmsignature_ind", "enable") == True) or (gc.genomon_conf.getboolean("pmsignature_full", "enable") == True):
        if not os.path.isdir(rc.run_conf.project_root + '/pmsignature/'): os.mkdir(rc.run_conf.project_root + '/pmsignature/')
        if not os.path.isdir(rc.run_conf.project_root + '/pmsignature/' + sample_conf_name): os.mkdir(rc.run_conf.project_root + '/pmsignature/' + sample_conf_name)
        if not os.path.isdir(rc.run_conf.project_root + '/script/pmsignature'): os.mkdir(rc.run_conf.project_root + '/script/pmsignature')
        if not os.path.isdir(rc.run_conf.project_root + '/log/pmsignature'): os.mkdir(rc.run_conf.project_root + '/log/pmsignature')

if not os.path.isdir(rc.run_conf.project_root + '/config'): os.mkdir(rc.run_conf.project_root + '/config')

for outputfiles in (bam2fastq_output_list, linked_fastq_list):
    for outputfile in outputfiles:
        sample = os.path.basename(os.path.dirname(outputfile[0][0]))
        fastq_dir = rc.run_conf.project_root + '/fastq/' + sample
        bam_dir = rc.run_conf.project_root + '/bam/' + sample
        if not os.path.isdir(fastq_dir): os.mkdir(fastq_dir)
        if not os.path.isdir(bam_dir): os.mkdir(bam_dir)

for target_sample_dict in (sc.sample_conf.bam_import, sc.sample_conf.fastq, sc.sample_conf.bam_tofastq):
    for sample in target_sample_dict:
        script_dir = rc.run_conf.project_root + '/script/' + sample
        log_dir = rc.run_conf.project_root + '/log/' + sample
        if not os.path.isdir(script_dir): os.mkdir(script_dir)
        if not os.path.isdir(log_dir): os.mkdir(log_dir)

shutil.copyfile(rc.run_conf.genomon_conf_file, rc.run_conf.project_root + '/config/' + genomon_conf_name +'_'+ rc.run_conf.analysis_timestamp + genomon_conf_ext)
shutil.copyfile(rc.run_conf.sample_conf_file, rc.run_conf.project_root + '/config/' + sample_conf_name +'_'+ rc.run_conf.analysis_timestamp + sample_conf_ext)

# prepare output directory for each sample and make mutation control panel file
for complist in sc.sample_conf.mutation_call:
    # make dir
    mutation_dir = rc.run_conf.project_root + '/mutation/' + complist[0]
    if not os.path.isdir(mutation_dir): os.mkdir(mutation_dir)
    # make the control panel text 
    control_panel_name = complist[2]
    if control_panel_name != None:
        control_panel_file = rc.run_conf.project_root + '/mutation/control_panel/' + control_panel_name + ".control_panel.txt"
        with open(control_panel_file,  "w") as out_handle:
            for panel_sample in sc.sample_conf.control_panel[control_panel_name]:
                out_handle.write(rc.run_conf.project_root + '/bam/' + panel_sample + '/' + panel_sample + '.markdup.bam' + "\n")

# make SV configuration file
for complist in sc.sample_conf.sv_detection:
    # make the control yaml file
    control_panel_name = complist[2]
    if control_panel_name != None:
        control_conf = rc.run_conf.project_root + '/sv/control_panel/' + control_panel_name + ".control_info.txt"
        with open(control_conf,  "w") as out_handle:
            for sample in sc.sample_conf.control_panel[control_panel_name]:
                out_handle.write(sample+ "\t"+ rc.run_conf.project_root+ "/sv/"+ sample +"/"+ sample+ "\n")

# link the import bam to project directory
@ruffus.originate(sc.sample_conf.bam_import.keys())
def link_import_bam(sample):
    bam = sc.sample_conf.bam_import[sample]
    link_dir = rc.run_conf.project_root + '/bam/' + sample
    bam_prefix, ext = os.path.splitext(bam)
    
    if not os.path.isdir(link_dir): os.mkdir(link_dir)
    if (not os.path.exists(link_dir +'/'+ sample +'.markdup.bam')) and (not os.path.exists(link_dir +'/'+ sample +'.markdup.bam.bai')): 
        os.symlink(bam, link_dir +'/'+ sample +'.markdup.bam')
        if (os.path.exists(bam +'.bai')):
            os.symlink(bam +'.bai', link_dir +'/'+ sample +'.markdup.bam.bai')
        elif (os.path.exists(bam_prefix +'.bai')):
            os.symlink(bam_prefix +'.bai', link_dir +'/'+ sample +'.markdup.bam.bai')

# convert bam to fastq
@ruffus.originate(bam2fastq_output_list)
def bam2fastq(outputfiles):
    sample = os.path.basename(os.path.dirname(outputfiles[0][0]))
    output_dir = rc.run_conf.project_root + '/fastq/' + sample
            
    arguments = {"biobambam": gc.genomon_conf.get("SOFTWARE", "biobambam"),
                 "param": gc.genomon_conf.get("bam2fastq", "params"),
                 "input_bam": sc.sample_conf.bam_tofastq[sample],
                 "f1_name": outputfiles[0][0],
                 "f2_name": outputfiles[1][0],
                 "o1_name": output_dir + '/unmatched_first_output.txt',
                 "o2_name": output_dir + '/unmatched_second_output.txt',
                 "t": output_dir + '/temp.txt',
                 "s": output_dir + '/single_end_output.txt'}
    bamtofastq.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample, rc.run_conf.project_root + '/script/'+ sample)


# link the input fastq to project directory
@ruffus.originate(linked_fastq_list)
def link_input_fastq(output_file):
    sample = os.path.basename(os.path.dirname(output_file[0][0]))
    fastq_dir = rc.run_conf.project_root + '/fastq/' + sample
    fastq_prefix, ext = os.path.splitext(sc.sample_conf.fastq[sample][0][0])
    # Todo
    # 1. should compare the timestamps between input and linked file
    # 2. check md5sum ?
    for (count, fastq_files) in enumerate(sc.sample_conf.fastq[sample][0]):
        fastq_prefix, ext = os.path.splitext(fastq_files)
        if not os.path.exists(fastq_dir + '/'+str(count+1)+'_1'+ ext): os.symlink(sc.sample_conf.fastq[sample][0][count], fastq_dir + '/'+str(count+1)+'_1'+ ext)
        if not os.path.exists(fastq_dir + '/'+str(count+1)+'_2'+ ext): os.symlink(sc.sample_conf.fastq[sample][1][count], fastq_dir + '/'+str(count+1)+'_2'+ ext)


# split fastq
@ruffus.subdivide([bam2fastq, link_input_fastq], ruffus.formatter(), "{path[0]}/*_*.fastq_split", "{path[0]}")
def split_files(input_files, output_files, target_dir):

    sample_name = os.path.basename(target_dir)

    for oo in output_files:
        os.unlink(oo)

    split_lines = gc.genomon_conf.get("split_fastq", "split_fastq_line_number")

    input_prefix, ext = os.path.splitext(input_files[0][0])
    arguments = {"lines": split_lines,
                 "fastq_filter": gc.genomon_conf.get("split_fastq", "fastq_filter"),
                 "target_dir": target_dir,
                 "ext": ext}
    
    fastq_splitter.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name, rc.run_conf.project_root + '/script/'+ sample_name, 2)
   
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
@ruffus.subdivide(split_files, ruffus.formatter(".+/(.+)/1_0000.fastq_split"), ruffus.add_inputs("{subpath[0][2]}/fastq/{subdir[0][0]}/2_0000.fastq_split"), "{subpath[0][2]}/bam/{subdir[0][0]}/{subdir[0][0]}_*.sorted.bam", "{subpath[0][2]}/fastq/{subdir[0][0]}", "{subpath[0][2]}/bam/{subdir[0][0]}")
def map_dna_sequence(input_files, output_files, input_dir, output_dir):

    sample_name = os.path.basename(output_dir)

    all_line_num = 0
    with open(input_dir + "/fastq_line_num.txt") as in_handle:
        tmp_num = in_handle.read()
        all_line_num = int(tmp_num)
    split_lines = gc.genomon_conf.get("split_fastq", "split_fastq_line_number")

    ans_quotient = all_line_num / int(split_lines)
    ans_remainder = all_line_num % int(split_lines)
    max_task_id = ans_quotient if ans_remainder == 0 else ans_quotient + 1
    
    arguments = {"input_dir": input_dir,
                 "output_dir": output_dir,
                 "sample_name": sample_name,
                 "bwa": gc.genomon_conf.get("SOFTWARE", "bwa"),
                 "bwa_params": gc.genomon_conf.get("bwa_mem", "bwa_params"),
                 "ref_fa": gc.genomon_conf.get("REFERENCE", "ref_fasta"),
                 "biobambam": gc.genomon_conf.get("SOFTWARE", "biobambam")}

    bwa_align.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name , rc.run_conf.project_root + '/script/' + sample_name, max_task_id) 

    for task_id in range(max_task_id):
        num = str(task_id).zfill(4)
        os.unlink(input_dir +'/1_'+str(num)+'.fastq_split')
        os.unlink(input_dir +'/2_'+str(num)+'.fastq_split')
        os.unlink(output_dir+'/'+sample_name+'_'+str(num)+'.bwa.sam')


# merge sorted bams into one and mark duplicate reads with biobambam
@ruffus.collate(map_dna_sequence, ruffus.formatter(), "{subpath[0][2]}/bam/{subdir[0][0]}/{subdir[0][0]}.markdup.bam", "{subpath[0][2]}/bam/{subdir[0][0]}")
def markdup(input_files, output_file, output_dir):

    sample_name = os.path.basename(output_dir)

    output_prefix, ext = os.path.splitext(output_file)

    input_bam_files = ""
    for input_file in input_files:
        input_bam_files = input_bam_files + " I=" + input_file

    arguments = {"biobambam": gc.genomon_conf.get("SOFTWARE", "biobambam"),
                 "out_prefix": output_prefix,
                 "input_bam_files": input_bam_files,
                 "out_bam": output_file}

    markduplicates.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name , rc.run_conf.project_root + '/script/'+ sample_name)

    for input_file in input_files:
        os.unlink(input_file)
        os.unlink(input_file + ".bai")


# identify mutations
@ruffus.follows( markdup )
@ruffus.follows( link_import_bam )
@ruffus.subdivide(markdup_bam_list, ruffus.formatter(), "{subpath[0][2]}/mutation/{subdir[0][0]}/{subdir[0][0]}.genomon_mutation.result.filt.txt", "{subpath[0][2]}/mutation/{subdir[0][0]}")
def identify_mutations(input_file, output_file, output_dir):

    sample_name = os.path.basename(output_dir)

    active_inhouse_normal_flag = False
    if gc.genomon_conf.has_option("annotation", "active_inhouse_normal_flag"):
        active_inhouse_normal_flag = gc.genomon_conf.get("annotation", "active_inhouse_normal_flag")

    inhouse_normal_tabix_db = ""
    if gc.genomon_conf.has_option("REFERENCE", "inhouse_normal_tabix_db"):
        inhouse_normal_tabix_db = gc.genomon_conf.get("REFERENCE", "inhouse_normal_tabix_db")

    active_inhouse_tumor_flag = False
    if gc.genomon_conf.has_option("annotation", "active_inhouse_tumor_flag"):
        active_inhouse_tumor_flag = gc.genomon_conf.get("annotation", "active_inhouse_tumor_flag")

    inhouse_tumor_tabix_db = ""
    if gc.genomon_conf.has_option("REFERENCE", "inhouse_tumor_tabix_db"):
        inhouse_tumor_tabix_db = gc.genomon_conf.get("REFERENCE", "inhouse_tumor_tabix_db")

    active_HGMD_flag = False
    if gc.genomon_conf.has_option("annotation", "active_HGMD_flag"):
        active_HGMD_flag = gc.genomon_conf.get("annotation", "active_HGMD_flag")
        
    HGMD_tabix_db = ""
    if gc.genomon_conf.has_option("REFERENCE", "HGMD_tabix_db"):
        HGMD_tabix_db = gc.genomon_conf.get("REFERENCE", "HGMD_tabix_db")

    arguments = {
        # fisher mutation
        "fisher": gc.genomon_conf.get("SOFTWARE", "fisher"),
        "fisher_pair_params": gc.genomon_conf.get("fisher_mutation_call", "pair_params"),
        "fisher_single_params": gc.genomon_conf.get("fisher_mutation_call", "single_params"),
        # realignment filter
        "mutfilter": gc.genomon_conf.get("SOFTWARE", "mutfilter"),
        "realignment_params": gc.genomon_conf.get("realignment_filter","params"),
        # indel filter
        "indel_params": gc.genomon_conf.get("indel_filter", "params"),
        # breakpoint filter
        "breakpoint_params": gc.genomon_conf.get("breakpoint_filter","params"),
        # simplerepeat filter
        "simple_repeat_db": gc.genomon_conf.get("REFERENCE", "simple_repeat_tabix_db"),
        # EB filter
        "EBFilter": gc.genomon_conf.get("SOFTWARE", "ebfilter"),
        "eb_map_quality": gc.genomon_conf.get("eb_filter","map_quality"),
        "eb_base_quality": gc.genomon_conf.get("eb_filter","base_quality"),
        "filter_flags": gc.genomon_conf.get("eb_filter","filter_flags"),
        "control_bam_list": input_file[2],
        # hotspot mutation caller
        "hotspot": gc.genomon_conf.get("SOFTWARE","hotspot"),
        "hotspot_database": gc.genomon_conf.get("REFERENCE","hotspot_db"),
        "active_hotspot_flag": gc.genomon_conf.get("hotspot","active_hotspot_flag"),
        "hotspot_params": gc.genomon_conf.get("hotspot","params"),
        "mutil": gc.genomon_conf.get("SOFTWARE", "mutil"),
        # original_annotations
        "mutanno": gc.genomon_conf.get("SOFTWARE", "mutanno"),
        "active_inhouse_normal_flag": active_inhouse_normal_flag,
        "inhouse_normal_database":inhouse_normal_tabix_db,
        "active_inhouse_tumor_flag": active_inhouse_tumor_flag,
        "inhouse_tumor_database":inhouse_tumor_tabix_db,
        "active_HGVD_2013_flag": gc.genomon_conf.get("annotation", "active_HGVD_2013_flag"),
        "HGVD_2013_database": gc.genomon_conf.get("REFERENCE", "HGVD_2013_tabix_db"),
        "active_HGVD_2016_flag": gc.genomon_conf.get("annotation", "active_HGVD_2016_flag"),
        "HGVD_2016_database": gc.genomon_conf.get("REFERENCE", "HGVD_2016_tabix_db"),
        "active_ExAC_flag": gc.genomon_conf.get("annotation", "active_ExAC_flag"),
        "ExAC_database": gc.genomon_conf.get("REFERENCE", "ExAC_tabix_db"),
        "active_HGMD_flag": active_HGMD_flag,
        "HGMD_database": HGMD_tabix_db,
        # annovar
        "active_annovar_flag": gc.genomon_conf.get("annotation", "active_annovar_flag"),
        "annovar": gc.genomon_conf.get("SOFTWARE", "annovar"),
        "annovar_database": gc.genomon_conf.get("annotation", "annovar_database"),
        "table_annovar_params": gc.genomon_conf.get("annotation", "table_annovar_params"),
        "annovar_buildver": gc.genomon_conf.get("annotation", "annovar_buildver"),
        # commmon
        "pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
        "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),   
        "ld_library_path": gc.genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
        "ref_fa": gc.genomon_conf.get("REFERENCE", "ref_fasta"),
        "interval_list": gc.genomon_conf.get("REFERENCE", "interval_list"),
        "disease_bam": input_file[0],
        "control_bam": input_file[1],
        "out_prefix": output_dir + '/' + sample_name,
        "samtools": gc.genomon_conf.get("SOFTWARE", "samtools"),
        "blat": gc.genomon_conf.get("SOFTWARE", "blat")}

    interval_list = gc.genomon_conf.get("REFERENCE", "interval_list")
    max_task_id = sum(1 for line in open(interval_list))

    mutation_call.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name, rc.run_conf.project_root + '/script/' + sample_name, max_task_id)
    
    arguments = {
        "pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
        "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),   
        "ld_library_path": gc.genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
        "control_bam": input_file[1],
        "control_bam_list": input_file[2],
        "active_annovar_flag": gc.genomon_conf.get("annotation", "active_annovar_flag"),
        "annovar_buildver": gc.genomon_conf.get("annotation", "annovar_buildver"),
        "active_HGVD_2013_flag": gc.genomon_conf.get("annotation", "active_HGVD_2013_flag"),
        "active_HGVD_2016_flag": gc.genomon_conf.get("annotation", "active_HGVD_2016_flag"),
        "active_ExAC_flag": gc.genomon_conf.get("annotation", "active_ExAC_flag"),
        "active_HGMD_flag": active_HGMD_flag,
        "active_inhouse_normal_flag": active_inhouse_normal_flag,
        "active_inhouse_tumor_flag": active_inhouse_tumor_flag,
        "filecount": max_task_id,
        "mutil": gc.genomon_conf.get("SOFTWARE", "mutil"),
        "pair_params": gc.genomon_conf.get("mutation_util","pair_params"),
        "single_params": gc.genomon_conf.get("mutation_util","single_params"),
        "active_hotspot_flag": gc.genomon_conf.get("hotspot","active_hotspot_flag"),
        "hotspot_database": gc.genomon_conf.get("REFERENCE","hotspot_db"),
        "meta_info_em": gc.get_meta_info(["fisher", "mutfilter", "ebfilter", "mutil", "mutanno"]),
        "meta_info_m": gc.get_meta_info(["fisher", "mutfilter", "mutil", "mutanno"]),
        "meta_info_ema": gc.get_meta_info(["fisher", "mutfilter", "ebfilter", "mutil", "mutanno", "hotspot"]),
        "meta_info_ma": gc.get_meta_info(["fisher", "mutfilter", "mutil", "mutanno", "hotspot"]),
        "out_prefix": output_dir + '/' + sample_name}

    mutation_merge.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name, rc.run_conf.project_root + '/script/' + sample_name)

    annovar_buildver = gc.genomon_conf.get("annotation", "annovar_buildver"),
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
@ruffus.follows( link_import_bam )
@ruffus.follows( markdup )
@ruffus.transform(parse_sv_bam_list, ruffus.formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.junction.clustered.bedpe.gz")
def parse_sv(input_file, output_file):

    dir_name = os.path.dirname(output_file)
    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    sample_name = os.path.basename(dir_name)

    arguments = {"genomon_sv": gc.genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "input_bam": input_file,
                 "output_prefix": output_file.replace(".junction.clustered.bedpe.gz", ""),
                 "param": gc.genomon_conf.get("sv_parse", "params"),
                 "pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": gc.genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "htslib": gc.genomon_conf.get("SOFTWARE", "htslib")}

    sv_parse.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name , rc.run_conf.project_root + '/script/' + sample_name)


# merge SV
@ruffus.follows( parse_sv )
@ruffus.transform(merge_bedpe_list, ruffus.formatter(".+/(?P<NAME>.+).control_info.txt"), "{subpath[0][2]}/sv/non_matched_control_panel/{NAME[0]}.merged.junction.control.bedpe.gz")
def merge_sv(input_files,  output_file):

    arguments = {"genomon_sv": gc.genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "control_info": input_files[0],
                 "merge_output_file": output_file,
                 "param": gc.genomon_conf.get("sv_merge", "params"),
                 "pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": gc.genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "htslib": gc.genomon_conf.get("SOFTWARE", "htslib")}

    sv_merge.task_exec(arguments, rc.run_conf.project_root + '/log/sv_merge', rc.run_conf.project_root + '/script/sv_merge')


# filt SV
@ruffus.follows( merge_sv )
@ruffus.transform(filt_bedpe_list, ruffus.formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.genomonSV.result.filt.txt")
def filt_sv(input_files,  output_file):

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    #sample_yaml = rc.run_conf.project_root + "/sv/config/" + sample_name + ".yaml"

    filt_param = ""

    for complist in sc.sample_conf.sv_detection:
        if sample_name == complist[0]:

            if complist[1] != None:
                filt_param = filt_param + " --matched_control_bam " + rc.run_conf.project_root + "/bam/" + complist[1] + '/' + complist[1] + ".markdup.bam"

            if complist[2] != None:
                filt_param = filt_param + " --non_matched_control_junction " + rc.run_conf.project_root +"/sv/non_matched_control_panel/"+ complist[2] +".merged.junction.control.bedpe.gz"
                if complist[1] != None:
                    filt_param = filt_param + " --matched_control_label " + complist[1]

            break

    filt_param = filt_param.lstrip(' ') + ' ' + gc.genomon_conf.get("sv_filt", "params")

    arguments = {"genomon_sv": gc.genomon_conf.get("SOFTWARE", "genomon_sv"),
                 "input_bam": rc.run_conf.project_root + "/bam/" + sample_name + '/' + sample_name + ".markdup.bam",
                 "output_prefix": rc.run_conf.project_root + "/sv/" + sample_name + '/' + sample_name,
                 "reference_genome": gc.genomon_conf.get("REFERENCE", "ref_fasta"),
                 "param": filt_param,
                 "meta_info": gc.get_meta_info(["genomon_sv", "sv_utils"]),
                 "sv_utils": gc.genomon_conf.get("SOFTWARE", "sv_utils"),
                 "sv_utils_param": gc.genomon_conf.get("sv_filt", "sv_utils_params"),
                 "pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),   
                 "ld_library_path": gc.genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "htslib": gc.genomon_conf.get("SOFTWARE", "htslib"),
                 "blat": gc.genomon_conf.get("SOFTWARE", "blat")}

    sv_filt.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name, rc.run_conf.project_root + '/script/' + sample_name)


# qc
@ruffus.follows( link_import_bam )
@ruffus.follows( markdup )
@ruffus.follows( filt_sv )
@ruffus.follows( identify_mutations )
@ruffus.transform(qc_bamstats_list, ruffus.formatter(), "{subpath[0][2]}/qc/{subdir[0][0]}/{subdir[0][0]}.bamstats")
def bam_stats(input_file, output_file):
    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    
    arguments = {"pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_qc": gc.genomon_conf.get("SOFTWARE", "genomon_qc"),
                 "bamstats": gc.genomon_conf.get("SOFTWARE", "bamstats"),
                 "perl": gc.genomon_conf.get("ENV", "PERL"),
                 "perl5lib": gc.genomon_conf.get("ENV", "PERL5LIB"),
                 "input_file": input_file,
                 "output_file": output_file}
    
    r_qc_bamstats.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name, rc.run_conf.project_root + '/script/' + sample_name)


@ruffus.follows( link_import_bam )
@ruffus.follows( markdup )
@ruffus.follows( filt_sv )
@ruffus.follows( identify_mutations )
@ruffus.transform(qc_coverage_list, ruffus.formatter(), "{subpath[0][2]}/qc/{subdir[0][0]}/{subdir[0][0]}.coverage")
def coverage(input_file, output_file):
    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    
    data_type = "exome"
    if gc.genomon_conf.get("qc_coverage", "wgs_flag") == "True":
        data_type = "wgs"

    arguments = {"data_type": data_type,
                 "pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_qc": gc.genomon_conf.get("SOFTWARE", "genomon_qc"),
                 "coverage_text": gc.genomon_conf.get("qc_coverage", "coverage"),
                 "i_bed_lines": gc.genomon_conf.get("qc_coverage", "wgs_i_bed_lines"),
                 "i_bed_width": gc.genomon_conf.get("qc_coverage", "wgs_i_bed_width"),
                 "incl_bed_width": gc.genomon_conf.get("qc_coverage", "wgs_incl_bed_width"),
                 "genome_size_file": gc.genomon_conf.get("REFERENCE", "genome_size"),
                 "gaptxt": gc.genomon_conf.get("REFERENCE", "gaptxt"),
                 "bait_file": gc.genomon_conf.get("REFERENCE", "bait_file"),
                 "samtools_params": gc.genomon_conf.get("qc_coverage", "samtools_params"),
                 "bedtools": gc.genomon_conf.get("SOFTWARE", "bedtools"),
                 "samtools": gc.genomon_conf.get("SOFTWARE", "samtools"),
                 "ld_library_path": gc.genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "input_file": input_file,
                 "output_file": output_file}

    r_qc_coverage.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name , rc.run_conf.project_root + '/script/' + sample_name)

@ruffus.follows( bam_stats )
@ruffus.follows( coverage )
@ruffus.collate(qc_merge_list, ruffus.formatter(), "{subpath[0][2]}/qc/{subdir[0][0]}/{subdir[0][0]}.genomonQC.result.txt")
def merge_qc(input_files, output_file):

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)
    
    arguments = {"pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
                 "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_qc": gc.genomon_conf.get("SOFTWARE", "genomon_qc"),
                 "bamstats_file": input_files[0][0],
                 "coverage_file": input_files[0][1],
                 "output_file": output_file,
                 "meta": gc.get_meta_info(["genomon_pipeline"]),
                 "fastq_line_num_file": rc.run_conf.project_root +'/fastq/'+ sample_name +'/fastq_line_num.txt'}
    
    r_qc_merge.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name, rc.run_conf.project_root + '/script/' + sample_name)

#####################
# post analysis stage
@ruffus.active_if(gc.genomon_conf.getboolean("post_analysis", "enable"))
@ruffus.active_if(len(pa_inputs_mutation) > 0)
@ruffus.follows(filt_sv)
@ruffus.follows(identify_mutations)
@ruffus.collate(pa_inputs_mutation, ruffus.formatter(), pa_outputs_mutation["outputs"])
def post_analysis_mutation(input_files, output_file):
        
    arguments = {"pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": gc.genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_pa":  gc.genomon_conf.get("SOFTWARE", "genomon_pa"),
                 "mode": "mutation",
                 "genomon_root": rc.run_conf.project_root,
                 "output_dir": rc.run_conf.project_root + "/post_analysis/" + sample_conf_name,
                 "sample_sheet": os.path.abspath(rc.run_conf.sample_conf_file),
                 "config_file": gc.genomon_conf.get("post_analysis", "config_file"),
                 "samtools": gc.genomon_conf.get("SOFTWARE", "samtools"),
                 "bedtools": gc.genomon_conf.get("SOFTWARE", "bedtools"),
                 "input_file_case1": ",".join(pa_outputs_mutation["case1"]["samples"]),
                 "input_file_case2": ",".join(pa_outputs_mutation["case2"]["samples"]),
                 "input_file_case3": ",".join(pa_outputs_mutation["case3"]["samples"]),
                 "input_file_case4": ",".join(pa_outputs_mutation["case4"]["samples"]),
                }

    r_post_analysis.task_exec(arguments, rc.run_conf.project_root + '/log/post_analysis', rc.run_conf.project_root + '/script/post_analysis')
    
@ruffus.active_if(gc.genomon_conf.getboolean("post_analysis", "enable"))
@ruffus.active_if(len(pa_inputs_sv) > 0)
@ruffus.follows(filt_sv)
@ruffus.follows(identify_mutations)
@ruffus.collate(pa_inputs_sv, ruffus.formatter(), pa_outputs_sv["outputs"])
def post_analysis_sv(input_files, output_file):

    arguments = {"pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": gc.genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_pa":  gc.genomon_conf.get("SOFTWARE", "genomon_pa"),
                 "mode": "sv",
                 "genomon_root": rc.run_conf.project_root,
                 "output_dir": rc.run_conf.project_root + "/post_analysis/" + sample_conf_name,
                 "sample_sheet": os.path.abspath(rc.run_conf.sample_conf_file),
                 "config_file": gc.genomon_conf.get("post_analysis", "config_file"),
                 "samtools": gc.genomon_conf.get("SOFTWARE", "samtools"),
                 "bedtools": gc.genomon_conf.get("SOFTWARE", "bedtools"),
                 "input_file_case1": ",".join(pa_outputs_sv["case1"]["samples"]),
                 "input_file_case2": ",".join(pa_outputs_sv["case2"]["samples"]),
                 "input_file_case3": ",".join(pa_outputs_sv["case3"]["samples"]),
                 "input_file_case4": ",".join(pa_outputs_sv["case4"]["samples"]),
                }
                 
    r_post_analysis.task_exec(arguments, rc.run_conf.project_root + '/log/post_analysis', rc.run_conf.project_root + '/script/post_analysis')

@ruffus.active_if(gc.genomon_conf.getboolean("post_analysis", "enable"))
@ruffus.active_if(len(pa_inputs_qc) > 0)
@ruffus.follows(merge_qc)
@ruffus.collate(pa_inputs_qc, ruffus.formatter(), pa_outputs_qc["outputs"])
def post_analysis_qc(input_files, output_file):

    arguments = {"pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": gc.genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),
                 "genomon_pa":  gc.genomon_conf.get("SOFTWARE", "genomon_pa"),
                 "mode": "qc",
                 "genomon_root": rc.run_conf.project_root,
                 "output_dir": rc.run_conf.project_root + "/post_analysis/" + sample_conf_name,
                 "sample_sheet": os.path.abspath(rc.run_conf.sample_conf_file),
                 "config_file": gc.genomon_conf.get("post_analysis", "config_file"),
                 "samtools": gc.genomon_conf.get("SOFTWARE", "samtools"),
                 "bedtools": gc.genomon_conf.get("SOFTWARE", "bedtools"),
                 "input_file_case1": ",".join(sc.sample_conf.qc),
                 "input_file_case2": "",
                 "input_file_case3": "",
                 "input_file_case4": "",
                }
                 
    r_post_analysis.task_exec(arguments, rc.run_conf.project_root + '/log/post_analysis', rc.run_conf.project_root + '/script/post_analysis')
    
@ruffus.active_if(gc.genomon_conf.getboolean("post_analysis", "enable"))
@ruffus.active_if(gc.genomon_conf.getboolean("pmsignature_ind", "enable") or gc.genomon_conf.getboolean("pmsignature_full", "enable"))
@ruffus.active_if(len(pmsignature_inputs) > 0)
@ruffus.follows(post_analysis_mutation)
@ruffus.collate(pmsignature_inputs, ruffus.formatter(), rc.run_conf.project_root + '/pmsignature/' + sample_conf_name + "/mutation.cut.txt")
def pre_pmsignature(input_files, output_file):
        
    arguments = {"input_files" : " ".join(input_files),
                 "output_file" : rc.run_conf.project_root + '/pmsignature/' + sample_conf_name + "/mutation.cut.txt"
                }

    r_pre_pmsignature.task_exec(arguments, rc.run_conf.project_root + '/log/pmsignature', rc.run_conf.project_root + '/script/pmsignature')
    
@ruffus.active_if(gc.genomon_conf.getboolean("post_analysis", "enable"))
@ruffus.active_if(gc.genomon_conf.getboolean("pmsignature_ind", "enable"))
@ruffus.active_if(run_ind)
@ruffus.follows(pre_pmsignature)
@ruffus.transform(rc.run_conf.project_root + '/pmsignature/' + sample_conf_name + "/mutation.cut.txt", ruffus.formatter(), ind_outputs[0])
def pmsignature_ind(input_file, output_file):
    
    command = r_pmsignature_ind.ind_template.format(
                inputfile =  input_file,
                outputdir = rc.run_conf.project_root + '/pmsignature/' + sample_conf_name,
                trdirflag = gc.genomon_conf.get("pmsignature_ind", "trdirflag").upper(),
                trialnum = gc.genomon_conf.getint("pmsignature_ind", "trialnum"),
                bs_genome = gc.genomon_conf.get("pmsignature_ind", "bs_genome"),
                bgflag = gc.genomon_conf.get("pmsignature_ind", "bgflag"),
                txdb_transcript = gc.genomon_conf.get("pmsignature_ind", "txdb_transcript"),
                script_path = gc.genomon_conf.get("SOFTWARE", "r_scripts"))
    
    sig_nums = range(gc.genomon_conf.getint("pmsignature_ind", "signum_min"), gc.genomon_conf.getint("pmsignature_ind", "signum_max") + 1)
    sig_num_text = ""
    for i in sig_nums: sig_num_text += "%d " % i

    arguments = {"r_path": gc.genomon_conf.get("ENV", "R_PATH"),
                 "r_ld_library_path": gc.genomon_conf.get("ENV", "R_LD_LIBRARY_PATH"),
                 "r_libs": gc.genomon_conf.get("ENV", "R_LIBS"),
                 "command": command,
                 "sig_list": sig_num_text
                }
    max_task_id = len(sig_nums)
    r_pmsignature_ind.task_exec(arguments, rc.run_conf.project_root + '/log/pmsignature', rc.run_conf.project_root + '/script/pmsignature', max_task_id)

@ruffus.active_if(gc.genomon_conf.getboolean("post_analysis", "enable"))
@ruffus.active_if(gc.genomon_conf.getboolean("pmsignature_full", "enable"))
@ruffus.active_if(run_full)
@ruffus.follows(pre_pmsignature)
@ruffus.transform(rc.run_conf.project_root + '/pmsignature/' + sample_conf_name + "/mutation.cut.txt", ruffus.formatter(), full_outputs[0])
def pmsignature_full(input_file, output_file):
    
    command = r_pmsignature_full.full_template.format(
                inputfile = input_file,
                outputdir = rc.run_conf.project_root + '/pmsignature/' + sample_conf_name,
                trdirflag = gc.genomon_conf.get("pmsignature_full", "trdirflag").upper(),
                trialnum = gc.genomon_conf.getint("pmsignature_full", "trialnum"),
                bgflag = gc.genomon_conf.get("pmsignature_full", "bgflag"),
                bs_genome = gc.genomon_conf.get("pmsignature_full", "bs_genome"),
                txdb_transcript = gc.genomon_conf.get("pmsignature_full", "txdb_transcript"),
                script_path = gc.genomon_conf.get("SOFTWARE", "r_scripts"))
    
    sig_nums = range(gc.genomon_conf.getint("pmsignature_full", "signum_min"), gc.genomon_conf.getint("pmsignature_full", "signum_max") + 1)
    sig_num_text = ""
    for i in sig_nums: sig_num_text += "%d " % i

    arguments = {"r_path": gc.genomon_conf.get("ENV", "R_PATH"),
                 "r_ld_library_path": gc.genomon_conf.get("ENV", "R_LD_LIBRARY_PATH"),
                 "r_libs": gc.genomon_conf.get("ENV", "R_LIBS"),
                 "command": command,
                 "sig_list": sig_num_text
                }
    max_task_id = len(sig_nums)
    r_pmsignature_full.task_exec(arguments, rc.run_conf.project_root + '/log/pmsignature', rc.run_conf.project_root + '/script/pmsignature', max_task_id)

@ruffus.active_if(gc.genomon_conf.getboolean("post_analysis", "enable"))
@ruffus.active_if(gc.genomon_conf.getboolean("paplot", "enable"))
@ruffus.active_if(len(paplot_inputs) > 0)
@ruffus.follows(post_analysis_sv)
@ruffus.follows(post_analysis_qc)
@ruffus.follows(pmsignature_ind)
@ruffus.follows(pmsignature_full)
@ruffus.collate(paplot_inputs, ruffus.formatter(), rc.run_conf.project_root + '/paplot/' + sample_conf_name + '/index.html')
def paplot(input_file, output_file):
    if not os.path.exists(paplot_output) and os.path.exists(rc.run_conf.project_root + '/paplot/' + sample_conf_name + '/.meta.json'):
        os.unlink(rc.run_conf.project_root + '/paplot/' + sample_conf_name + '/.meta.json')

    command = ""
    if len(paplot_inputs_qc) > 0:
        command += r_paplot.qc_template.format(
                    paplot = gc.genomon_conf.get("SOFTWARE", "paplot"),
                    inputs = ",".join(paplot_inputs_qc),
                    output_dir = rc.run_conf.project_root + "/paplot/" + sample_conf_name,
                    title = gc.genomon_conf.get("paplot", "title"),
                    config_file = gc.genomon_conf.get("paplot", "config_file"))
                        
    if len(paplot_inputs_sv) > 0:
        command += r_paplot.sv_template.format(
                    paplot = gc.genomon_conf.get("SOFTWARE", "paplot"),
                    inputs = ",".join(paplot_inputs_sv),
                    output_dir = rc.run_conf.project_root + "/paplot/" + sample_conf_name,
                    title = gc.genomon_conf.get("paplot", "title"),
                    config_file = gc.genomon_conf.get("paplot", "config_file"))
                        
    if len(paplot_inputs_mutation) > 0:
        command += r_paplot.mutation_template.format(
                    paplot = gc.genomon_conf.get("SOFTWARE", "paplot"),
                    inputs = ",".join(paplot_inputs_mutation),
                    output_dir = rc.run_conf.project_root + "/paplot/" + sample_conf_name,
                    title = gc.genomon_conf.get("paplot", "title"),
                    config_file = gc.genomon_conf.get("paplot", "config_file"),
                    annovar = gc.genomon_conf.getboolean("annotation", "active_annovar_flag"))
    
    if gc.genomon_conf.getboolean("pmsignature_ind", "enable"):
        for i in range(len(paplot_inputs_ind)):
            command += r_paplot.ind_template.format(
                        paplot = gc.genomon_conf.get("SOFTWARE", "paplot"),
                        input = paplot_inputs_ind[i],
                        output_dir = rc.run_conf.project_root + "/paplot/" + sample_conf_name,
                        title = gc.genomon_conf.get("paplot", "title"),
                        config_file = gc.genomon_conf.get("paplot", "config_file"))
    
    if gc.genomon_conf.getboolean("pmsignature_full", "enable"):
        for i in range(len(paplot_inputs_full)):
            command += r_paplot.full_template.format(
                        paplot = gc.genomon_conf.get("SOFTWARE", "paplot"),
                        input = paplot_inputs_full[i],
                        output_dir = rc.run_conf.project_root + "/paplot/" + sample_conf_name,
                        title = gc.genomon_conf.get("paplot", "title"),
                        config_file = gc.genomon_conf.get("paplot", "config_file"))
    
    remark = gc.genomon_conf.get("paplot", "remarks")
    remark += "<ul>"
    
    for item in gc.genomon_conf.get("paplot", "software").split(","):
        key = item.split(":")[0].strip(" ").rstrip(" ")
        name = item.split(":")[1].strip(" ").rstrip(" ")
        try:
            version = gc.get_version(key).split("-")
        except Exception:
            print ("[WARNING] paplot: %s is not defined." % (key))
            continue
        
        remark += "<li>" + name + " " + version[-1] + "</li>"

    remark += "</ul>"
    
    command += r_paplot.index_template.format(
                        paplot = gc.genomon_conf.get("SOFTWARE", "paplot"),
                        output_dir = rc.run_conf.project_root + "/paplot/" + sample_conf_name,
                        remarks = remark,
                        config_file = gc.genomon_conf.get("paplot", "config_file"))

    arguments = {"pythonhome": gc.genomon_conf.get("ENV", "PYTHONHOME"),
                 "ld_library_path": gc.genomon_conf.get("ENV", "LD_LIBRARY_PATH"),
                 "pythonpath": gc.genomon_conf.get("ENV", "PYTHONPATH"),
                 "paplot":  gc.genomon_conf.get("SOFTWARE", "paplot"),
                 "command": command
                }
                 
    r_paplot.task_exec(arguments, rc.run_conf.project_root + '/log/paplot', rc.run_conf.project_root + '/script/paplot')

