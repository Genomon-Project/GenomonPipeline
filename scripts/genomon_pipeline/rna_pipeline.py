import os
import shutil
import ruffus
import genomon_pipeline.config.run_conf as rc
import genomon_pipeline.config.genomon_conf as gc
import genomon_pipeline.config.sample_conf as sc
import genomon_pipeline.rna_resource.star_align         as rr_star_align
import genomon_pipeline.rna_resource.fusionfusion       as rr_fusionfusion
import genomon_pipeline.rna_resource.post_analysis      as rr_post_analysis
import genomon_pipeline.rna_resource.paplot             as rr_paplot
import genomon_pipeline.rna_resource.fusion_count       as rr_fusion_count
import genomon_pipeline.rna_resource.fusion_merge       as rr_fusion_merge
import genomon_pipeline.rna_resource.genomon_expression as rr_genomon_expression
import genomon_pipeline.rna_resource.intron_retention   as rr_intron_retention
import genomon_pipeline.dna_resource.bamtofastq         as rr_bamtofastq

# set task classes
bamtofastq = rr_bamtofastq.Bam2Fastq(gc.genomon_conf.get("bam2fastq", "qsub_option"), rc.run_conf.drmaa)
star_align = rr_star_align.Star_align(gc.genomon_conf.get("star_align", "qsub_option"), rc.run_conf.drmaa)
fusionfusion = rr_fusionfusion.Fusionfusion(gc.genomon_conf.get("fusionfusion", "qsub_option"), rc.run_conf.drmaa)
fusion_count = rr_fusion_count.Fusion_count(gc.genomon_conf.get("fusion_count_control", "qsub_option"), rc.run_conf.drmaa)
fusion_merge = rr_fusion_merge.Fusion_merge(gc.genomon_conf.get("fusion_merge_control", "qsub_option"), rc.run_conf.drmaa)
genomon_expression = rr_genomon_expression.Genomon_expression(gc.genomon_conf.get("genomon_expression", "qsub_option"), rc.run_conf.drmaa)
intron_retention = rr_intron_retention.Intron_retention(gc.genomon_conf.get("intron_retention", "qsub_option"), rc.run_conf.drmaa)
r_paplot = rr_paplot.Res_PA_Plot(gc.genomon_conf.get("paplot", "qsub_option"), rc.run_conf.drmaa)
r_post_analysis = rr_post_analysis.Res_PostAnalysis(gc.genomon_conf.get("post_analysis", "qsub_option"), rc.run_conf.drmaa)

_debug = False
if gc.genomon_conf.has_section("develop"):
    if gc.genomon_conf.has_option("develop", "debug") == True:
        _debug = gc.genomon_conf.getboolean("develop", "debug")

# generate list of linked_fastq file path
linked_fastq_list = []
for sample in sc.sample_conf.fastq:
    if os.path.exists(rc.run_conf.project_root + '/star/' + sample + '/' + sample + '.Aligned.sortedByCoord.out.bam'): continue

    linked_fastq_list.append([rc.run_conf.project_root + '/fastq/' + sample + '/1_1.fastq',
                              rc.run_conf.project_root + '/fastq/' + sample + '/1_2.fastq'])

# generate output list of 'bam2fastq'
bam2fastq_output_list = []
for sample in sc.sample_conf.bam_tofastq:
    if os.path.exists(rc.run_conf.project_root + '/star/' + sample + '/' + sample + '.Aligned.sortedByCoord.out.bam'): continue

    bam2fastq_output_list.append([rc.run_conf.project_root + '/fastq/' + sample + '/1_1.fastq',
                                  rc.run_conf.project_root + '/fastq/' + sample + '/1_2.fastq'])

sample_list_fastq = sc.sample_conf.fastq

sample_conf_name, ext = os.path.splitext(os.path.basename(rc.run_conf.sample_conf_file))

# generate input list of 'post analysis for fusionfusion'
pa_outputs_fusion = r_post_analysis.output_files("fusion", sc.sample_conf.fusion, rc.run_conf.project_root, sample_conf_name, gc.genomon_conf)

pa_inputs_fusion = []
if pa_outputs_fusion["run_pa"] == True:
    for complist in sc.sample_conf.fusion:
        pa_inputs_fusion.append(rc.run_conf.project_root + '/fusion/' + complist[0] + '/' + complist[0] + '.genomonFusion.result.filt.txt')
        
# generate input list of 'post analysis for qc'
pa_outputs_starqc = r_post_analysis.output_files("starqc", sc.sample_conf.qc, rc.run_conf.project_root, sample_conf_name, gc.genomon_conf)

pa_inputs_starqc = []
if pa_outputs_starqc["run_pa"] == True:
    for sample in sc.sample_conf.qc:
        pa_inputs_starqc.append(rc.run_conf.project_root + '/star/' + sample + '/' + sample + '.Log.final.out')

# generate input list of paplot
paplot_output = rc.run_conf.project_root + '/paplot/' + sample_conf_name + '/index.html'

## fusionfusion
paplot_inputs_fusion = []
if os.path.exists(paplot_output) == False or pa_outputs_fusion["run_pa"] == True:

    if pa_outputs_fusion["case1"]["output_filt"] != "":
        paplot_inputs_fusion.append(pa_outputs_fusion["case1"]["output_filt"])
    if pa_outputs_fusion["case2"]["output_filt"] != "" and gc.genomon_conf.getboolean("paplot", "include_unpanel"):
        paplot_inputs_fusion.append(pa_outputs_fusion["case2"]["output_filt"])

## star-qc
paplot_inputs_starqc = []
if os.path.exists(paplot_output) == False or pa_outputs_starqc["run_pa"] == True:
    paplot_inputs_starqc.extend(pa_outputs_starqc["outputs"])

paplot_inputs = []
paplot_inputs.extend(paplot_inputs_starqc)
paplot_inputs.extend(paplot_inputs_fusion)

if _debug:
    from pprint import pprint
    print ("post-analysis-fusion");  pprint (pa_outputs_fusion); print ("post-analysis-starqc");  pprint (pa_outputs_starqc)
    print ("paplot"); pprint (paplot_inputs)

# prepare output directories
if not os.path.isdir(rc.run_conf.project_root): os.makedirs(rc.run_conf.project_root)
if not os.path.isdir(rc.run_conf.project_root + '/script'): os.mkdir(rc.run_conf.project_root + '/script')
if not os.path.isdir(rc.run_conf.project_root + '/script/fusion_merge'): os.mkdir(rc.run_conf.project_root + '/script/fusion_merge')
if not os.path.isdir(rc.run_conf.project_root + '/log'): os.mkdir(rc.run_conf.project_root + '/log')
if not os.path.isdir(rc.run_conf.project_root + '/log/fusion_merge'): os.mkdir(rc.run_conf.project_root + '/log/fusion_merge')
if not os.path.isdir(rc.run_conf.project_root + '/fastq'): os.mkdir(rc.run_conf.project_root + '/fastq')
if not os.path.isdir(rc.run_conf.project_root + '/star'): os.mkdir(rc.run_conf.project_root + '/star')
if not os.path.isdir(rc.run_conf.project_root + '/fusion'): os.mkdir(rc.run_conf.project_root + '/fusion')
if not os.path.isdir(rc.run_conf.project_root + '/fusion/control_panel'): os.mkdir(rc.run_conf.project_root + '/fusion/control_panel')
if not os.path.isdir(rc.run_conf.project_root + '/config'): os.mkdir(rc.run_conf.project_root + '/config')
if not os.path.isdir(rc.run_conf.project_root + '/expression'): os.mkdir(rc.run_conf.project_root + '/expression')
if not os.path.isdir(rc.run_conf.project_root + '/intron_retention'): os.mkdir(rc.run_conf.project_root + '/intron_retention')

if (gc.genomon_conf.getboolean("post_analysis", "enable") == True):
    if not os.path.exists(rc.run_conf.project_root + '/post_analysis'): os.mkdir(rc.run_conf.project_root + '/post_analysis')
    if not os.path.exists(rc.run_conf.project_root + '/post_analysis/' + sample_conf_name): os.mkdir(rc.run_conf.project_root + '/post_analysis/' + sample_conf_name)
    if not os.path.isdir(rc.run_conf.project_root + '/script/post_analysis'): os.mkdir(rc.run_conf.project_root + '/script/post_analysis')
    if not os.path.isdir(rc.run_conf.project_root + '/log/post_analysis'): os.mkdir(rc.run_conf.project_root + '/log/post_analysis')

    if (gc.genomon_conf.getboolean("paplot", "enable") == True):
        if not os.path.exists(rc.run_conf.project_root + '/paplot'): os.mkdir(rc.run_conf.project_root + '/paplot')
        if not os.path.exists(rc.run_conf.project_root + '/paplot/' + sample_conf_name): os.mkdir(rc.run_conf.project_root + '/paplot/' + sample_conf_name)
        if not os.path.isdir(rc.run_conf.project_root + '/script/paplot'): os.mkdir(rc.run_conf.project_root + '/script/paplot')
        if not os.path.isdir(rc.run_conf.project_root + '/log/paplot'): os.mkdir(rc.run_conf.project_root + '/log/paplot')

for target_sample_dict in (sc.sample_conf.bam_import, sc.sample_conf.fastq, sc.sample_conf.bam_tofastq):
    for sample in target_sample_dict:
        script_dir = rc.run_conf.project_root + '/script/' + sample
        log_dir = rc.run_conf.project_root + '/log/' + sample
        if not os.path.isdir(script_dir): os.mkdir(script_dir)
        if not os.path.isdir(log_dir): os.mkdir(log_dir)

genomon_conf_name, genomon_conf_ext = os.path.splitext(os.path.basename(rc.run_conf.genomon_conf_file))
sample_conf_name, sample_conf_ext = os.path.splitext(os.path.basename(rc.run_conf.sample_conf_file))
shutil.copyfile(rc.run_conf.genomon_conf_file, rc.run_conf.project_root + '/config/' + genomon_conf_name +'_'+ rc.run_conf.analysis_timestamp + genomon_conf_ext)
shutil.copyfile(rc.run_conf.sample_conf_file, rc.run_conf.project_root + '/config/' + sample_conf_name +'_'+ rc.run_conf.analysis_timestamp + sample_conf_ext)

expression_bams = []
# generate input list of genomon expression
for tumor_sample in sc.sample_conf.expression:
    if os.path.exists(rc.run_conf.project_root + '/expression/' + tumor_sample + '/' + tumor_sample + '.genomonExpression.result.txt'): continue
    expression_bams.append(rc.run_conf.project_root + '/star/' + tumor_sample + '/' + tumor_sample + '.Aligned.sortedByCoord.out.bam')

intron_retention_bams = []
# generate input list of intron_retention
for tumor_sample in sc.sample_conf.intron_retention:
    if os.path.exists(rc.run_conf.project_root + '/intron_retention/' + tumor_sample + '/' + tumor_sample + '.genomonIR.result.txt'): continue
    intron_retention_bams.append(rc.run_conf.project_root + '/star/' + tumor_sample + '/' + tumor_sample + '.Aligned.sortedByCoord.out.bam')

fusionfusion_bams = []
fusion_control_panel = []
chimeric_ctrl_sams = []
# generate input list of fusionfusion
for complist in sc.sample_conf.fusion:

    tumor_sample = complist[0]
    control_panel_name = complist[1]
    if os.path.exists(rc.run_conf.project_root + '/fusion/' + tumor_sample + '/' + tumor_sample + '.genomonFusion.result.filt.txt'): continue

    # generate input list of 'fusionfusion'
    tumor_bam = rc.run_conf.project_root + '/star/' + tumor_sample + '/' + tumor_sample + '.Aligned.sortedByCoord.out.bam'
    merged_count = rc.run_conf.project_root + '/fusion/control_panel/' + control_panel_name + ".merged.Chimeric.count" if control_panel_name != None else None
    fusionfusion_bams.append([tumor_bam, merged_count])
   
    if control_panel_name != None:
        # generate input list of 'fusion count'
        for panel_sample in sc.sample_conf.control_panel[control_panel_name]:
            chimeric_ctrl_sams.append(rc.run_conf.project_root + '/star/' + panel_sample + '/' + panel_sample + '.Chimeric.out.sam')
        chimeric_ctrl_sams = list(set(chimeric_ctrl_sams))

        # generate input list of 'fusion merge control'
        control_panel_file = rc.run_conf.project_root + '/fusion/control_panel/' + control_panel_name + ".Chimeric_count.list"
        fusion_control_panel.append(control_panel_file)
        fusion_control_panel = list(set(fusion_control_panel))
        with open(control_panel_file,  "w") as out_handle:
            for panel_sample in sc.sample_conf.control_panel[control_panel_name]:
                out_handle.write(rc.run_conf.project_root + '/fusion/' + panel_sample + '/' + panel_sample + '.Chimeric.count' + "\n")


# link the import bam to project directory
@ruffus.originate(list(sc.sample_conf.bam_import.keys()))
def link_import_bam(sample):
    bam = sc.sample_conf.bam_import[sample]
    link_dir = rc.run_conf.project_root + '/star/' + sample
    bam_prefix, ext = os.path.splitext(bam)
    input_dir_name = os.path.dirname(bam)
    sample_name = os.path.basename(input_dir_name)
    input_chimeric_sam = input_dir_name + '/' + sample_name + ".Chimeric.out.sam"
    input_log_final = input_dir_name + '/' + sample_name + ".Log.final.out"
    
    if not os.path.isdir(link_dir): os.mkdir(link_dir)
    if (not os.path.exists(link_dir +'/'+ sample +'.Aligned.sortedByCoord.out.bam')) and (not os.path.exists(link_dir +'/'+ sample +'.Aligned.sortedByCoord.out.bam.bai')): 
        os.symlink(bam, link_dir +'/'+ sample +'.Aligned.sortedByCoord.out.bam')
        if (os.path.exists(bam +'.bai')):
            os.symlink(bam +'.bai', link_dir +'/'+ sample +'.Aligned.sortedByCoord.out.bam.bai')
        elif (os.path.exists(bam_prefix +'.bai')):
            os.symlink(bam_prefix +'.bai', link_dir +'/'+ sample +'.Aligned.sortedByCoord.out.bam.bai')
    if not os.path.exists(link_dir +'/'+ sample +'.Chimeric.out.sam') and os.path.exists(input_chimeric_sam): 
        os.symlink(input_chimeric_sam, link_dir +'/'+ sample +'.Chimeric.out.sam')
    if not os.path.exists(link_dir +'/'+ sample +'.Log.final.out') and os.path.exists(input_log_final): 
        os.symlink(input_log_final, link_dir +'/'+ sample +'.Log.final.out')

# convert bam to fastq
@ruffus.originate(bam2fastq_output_list)
def bam2fastq(outputfiles):
    sample = os.path.basename(os.path.dirname(outputfiles[0]))
    output_dir = rc.run_conf.project_root + '/fastq/' + sample
            
    arguments = {
        "param": gc.genomon_conf.get("bam2fastq", "params"),
        "input_bam": sc.sample_conf.bam_tofastq[sample],
        "f1_name": outputfiles[0],
        "f2_name": outputfiles[1],
        "o1_name": output_dir + '/unmatched_first_output.txt',
        "o2_name": output_dir + '/unmatched_second_output.txt',
        "t": output_dir + '/temp.txt',
        "s": output_dir + '/single_end_output.txt'
    }
    bind = [rc.run_conf.project_root]
    if sample in sc.sample_conf.bam_tofastq_src:
        bind.extend(sc.sample_conf.bam_tofastq_src[sample])
        
    singularity_params = {
        "image": gc.genomon_conf.get("bam2fastq", "image"),
        "option": gc.genomon_conf.get("bam2fastq", "singularity_option"),
        "bind": bind,
    }
    if not os.path.isdir(output_dir): os.mkdir(output_dir)
    bamtofastq.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample, rc.run_conf.project_root + '/script/'+ sample, singularity_params)

# link the input fastq files
@ruffus.originate(linked_fastq_list, sample_list_fastq)
def link_input_fastq(output_file, sample_list_fastq):
    sample = os.path.basename(os.path.dirname(output_file[0]))
    link_dir = rc.run_conf.project_root + '/fastq/' + sample

    # Todo
    # 1. should compare the timestamps between input and linked file
    # 2. check md5sum ?
    if not os.path.isdir(link_dir): os.mkdir(link_dir)
    if not os.path.exists(link_dir + '/1_1.fastq'): os.symlink(sample_list_fastq[sample][0][0], link_dir + '/1_1.fastq')
    if not os.path.exists(link_dir + '/1_2.fastq'): os.symlink(sample_list_fastq[sample][1][0], link_dir + '/1_2.fastq')


@ruffus.transform([bam2fastq,link_input_fastq], ruffus.formatter(), "{subpath[0][2]}/star/{subdir[0][0]}/{subdir[0][0]}.Aligned.sortedByCoord.out.bam")
def task_star_align(input_files, output_file):

    dir_name = os.path.dirname(output_file)
    sample_name = os.path.basename(dir_name)

    arguments = {
        "star_genome": gc.genomon_conf.get("REFERENCE", "star_genome"),
        "additional_params": gc.genomon_conf.get("star_align", "star_params"),
        "samtools_sort_params": gc.genomon_conf.get("star_align", "samtools_sort_params"),
        "fastq1": input_files[0],
        "fastq2": input_files[1],
        "out_prefix": dir_name + '/' + sample_name + '.'
    }
    bind = [
        rc.run_conf.project_root,
        gc.genomon_conf.get("REFERENCE", "star_genome")
    ]
    if sample_name in sc.sample_conf.fastq_src:
        bind.extend(sc.sample_conf.fastq_src[sample_name])
        
    singularity_params = {
        "image": gc.genomon_conf.get("star_align", "image"),
        "option": gc.genomon_conf.get("star_align", "singularity_option"),
        "bind": bind,
    }
    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    star_align.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name , rc.run_conf.project_root + '/script/' + sample_name, singularity_params)
    os.unlink(input_files[0])
    os.unlink(input_files[1])


@ruffus.follows( link_import_bam )
@ruffus.follows( task_star_align )
@ruffus.transform(chimeric_ctrl_sams, ruffus.formatter(), "{subpath[0][2]}/fusion/{subdir[0][0]}/{subdir[0][0]}.Chimeric.count")
def task_fusion_count(input_file, output_file):

    input_dir_name = os.path.dirname(input_file)
    sample_name = os.path.basename(input_dir_name)
    output_dir_name = os.path.dirname(output_file) 

    arguments = {
        "chimeric_sam": input_file,
        "output": output_file,
        "additional_params": gc.genomon_conf.get("fusion_count_control", "params")
    }
    bind = [rc.run_conf.project_root]
    if sample_name in sc.sample_conf.bam_import_src:
        bind.extend(sc.sample_conf.bam_import_src[sample_name])
    
    singularity_params = {
        "image": gc.genomon_conf.get("fusion_count_control", "image"),
        "option": gc.genomon_conf.get("fusion_count_control", "singularity_option"),
        "bind": bind,
    }
    if not os.path.isdir(output_dir_name): os.mkdir(output_dir_name)
    fusion_count.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name, rc.run_conf.project_root + '/script/' + sample_name, singularity_params)


@ruffus.follows( task_fusion_count )
@ruffus.transform(fusion_control_panel, ruffus.formatter(".+/(?P<NAME>.+).Chimeric.count.list"), "{subpath[0][2]}/fusion/control_panel/{NAME[0]}.merged.Chimeric.count")
def task_fusion_merge(input_file, output_file):

    arguments = {
        "count_list": input_file,
        "output": output_file,
        "additional_params": gc.genomon_conf.get("fusion_merge_control", "params")
    }
    singularity_params = {
        "image": gc.genomon_conf.get("fusion_merge_control", "image"),
        "option": gc.genomon_conf.get("fusion_merge_control", "singularity_option"),
        "bind": [
            rc.run_conf.project_root,
        ],
    }
    fusion_merge.task_exec(arguments, rc.run_conf.project_root + '/log/fusion_merge', rc.run_conf.project_root + '/script/fusion_merge', singularity_params)


@ruffus.follows( task_fusion_merge )
@ruffus.transform(fusionfusion_bams, ruffus.formatter(), "{subpath[0][2]}/fusion/{subdir[0][0]}/{subdir[0][0]}.genomonFusion.result.filt.txt")
def task_fusionfusion(input_file, output_file):

    input_dir_name = os.path.dirname(input_file[0])
    sample_name = os.path.basename(input_dir_name)
    input_chimeric_sam = input_dir_name + '/' + sample_name + ".Chimeric.out.sam"
    output_dir_name = os.path.dirname(output_file) 

    params = ""
    if input_file[1] != None:
        params = "--pooled_control_file " + input_file[1] + " "

    arguments = {
        "ref_fa":gc.genomon_conf.get("REFERENCE", "ref_fasta"),
        "chimeric_sam": input_chimeric_sam,
        "output_prefix": output_dir_name,
        "additional_params": params + gc.genomon_conf.get("fusionfusion", "params"),
        "filt_params": gc.genomon_conf.get("fusionfusion", "filt_params"),
        "sample": sample_name
    }
    singularity_params = {
        "image": gc.genomon_conf.get("fusionfusion", "image"),
        "option": gc.genomon_conf.get("fusionfusion", "singularity_option"),
        "bind": [
            rc.run_conf.project_root,
            gc.genomon_conf.get("REFERENCE", "ref_fasta")
        ],
    }
    if not os.path.isdir(output_dir_name): os.mkdir(output_dir_name)
    fusionfusion.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name, rc.run_conf.project_root + '/script/' + sample_name, singularity_params)


@ruffus.follows( link_import_bam )
@ruffus.follows( task_star_align )
@ruffus.transform(expression_bams, ruffus.formatter(), "{subpath[0][2]}/expression/{subdir[0][0]}/{subdir[0][0]}.genomonExpression.result.txt")
def task_genomon_expression(input_file, output_file):

    input_dir_name = os.path.dirname(input_file)
    sample_name = os.path.basename(input_dir_name)
    output_dir_name = os.path.dirname(output_file)  

    arguments = {
        "input_bam": input_file,
        "output_prefix": output_dir_name + '/' + sample_name,
        "additional_params": gc.genomon_conf.get("genomon_expression", "params")
    }
    singularity_params = {
        "image": gc.genomon_conf.get("genomon_expression", "image"),
        "option": gc.genomon_conf.get("genomon_expression", "singularity_option"),
        "bind": [
            rc.run_conf.project_root,
        ],
    }
    if not os.path.isdir(output_dir_name): os.mkdir(output_dir_name)
    genomon_expression.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name, rc.run_conf.project_root + '/script/' + sample_name, singularity_params)


@ruffus.follows( link_import_bam )
@ruffus.follows( task_star_align )
@ruffus.transform(intron_retention_bams, ruffus.formatter(), "{subpath[0][2]}/intron_retention/{subdir[0][0]}/{subdir[0][0]}.genomonIR.result.txt")
def task_intron_retention(input_file, output_file):

    input_dir_name = os.path.dirname(input_file)
    sample_name = os.path.basename(input_dir_name)
    output_dir_name = os.path.dirname(output_file)  

    arguments = {
        "input_bam": input_file,
        "output_prefix": output_dir_name + "/" + sample_name,
        "additional_params": gc.genomon_conf.get("intron_retention", "params"),
    }
    singularity_params = {
        "image": gc.genomon_conf.get("intron_retention", "image"),
        "option": gc.genomon_conf.get("intron_retention", "singularity_option"),
        "bind": [
            rc.run_conf.project_root,
        ],
    }
    if not os.path.isdir(output_dir_name): os.mkdir(output_dir_name)
    intron_retention.task_exec(arguments, rc.run_conf.project_root + '/log/' + sample_name, rc.run_conf.project_root + '/script/' + sample_name, singularity_params)


@ruffus.active_if(gc.genomon_conf.getboolean("post_analysis", "enable"))
@ruffus.active_if(len(pa_inputs_fusion) > 0)
@ruffus.follows(task_fusionfusion)
@ruffus.follows(task_genomon_expression)
@ruffus.collate(pa_inputs_fusion, ruffus.formatter(), pa_outputs_fusion["outputs"])
def post_analysis_fusion(input_files, output_file):

    arguments = {
        "mode": "fusion",
        "genomon_root": rc.run_conf.project_root,
        "output_dir": rc.run_conf.project_root + "/post_analysis/" + sample_conf_name,
        "sample_sheet": os.path.abspath(rc.run_conf.sample_conf_file),
        "config_file": gc.genomon_conf.get("post_analysis", "config_file"),
        "input_file_case1": ",".join(pa_outputs_fusion["case1"]["samples"]),
        "input_file_case2": ",".join(pa_outputs_fusion["case2"]["samples"]),
        "input_file_case3": "",
        "input_file_case4": ""
    }
    singularity_params = {
        "image": gc.genomon_conf.get("post_analysis", "image"),
        "option": gc.genomon_conf.get("post_analysis", "singularity_option"),
        "bind": [
            rc.run_conf.project_root,
            gc.genomon_conf.get("post_analysis", "config_file"),
            os.path.abspath(rc.run_conf.sample_conf_file)
        ],
    }
    r_post_analysis.task_exec(arguments, rc.run_conf.project_root + '/log/post_analysis', rc.run_conf.project_root + '/script/post_analysis', singularity_params)

@ruffus.active_if(gc.genomon_conf.getboolean("post_analysis", "enable"))
@ruffus.active_if(len(pa_inputs_starqc) > 0)
@ruffus.follows(task_fusionfusion)
@ruffus.follows(task_genomon_expression)
@ruffus.collate(pa_inputs_starqc, ruffus.formatter(), pa_outputs_starqc["outputs"])
def post_analysis_starqc(input_files, output_file):

    arguments = {
        "mode": "starqc",
        "genomon_root": rc.run_conf.project_root,
        "output_dir": rc.run_conf.project_root + "/post_analysis/" + sample_conf_name,
        "sample_sheet": os.path.abspath(rc.run_conf.sample_conf_file),
        "config_file": gc.genomon_conf.get("post_analysis", "config_file"),
        "input_file_case1": ",".join(sc.sample_conf.qc),
        "input_file_case2": "",
        "input_file_case3": "",
        "input_file_case4": ""
    }
                
    singularity_params = {
        "image": gc.genomon_conf.get("post_analysis", "image"),
        "option": gc.genomon_conf.get("post_analysis", "singularity_option"),
        "bind": [
            rc.run_conf.project_root,
            gc.genomon_conf.get("post_analysis", "config_file"),
            os.path.abspath(rc.run_conf.sample_conf_file)
        ],
    }
    r_post_analysis.task_exec(arguments, rc.run_conf.project_root + '/log/post_analysis', rc.run_conf.project_root + '/script/post_analysis', singularity_params)
    
@ruffus.active_if(gc.genomon_conf.getboolean("post_analysis", "enable"))
@ruffus.active_if(gc.genomon_conf.getboolean("paplot", "enable"))
@ruffus.active_if(len(paplot_inputs) > 0)
@ruffus.follows(post_analysis_fusion)
@ruffus.follows(post_analysis_starqc)
@ruffus.collate(paplot_inputs, ruffus.formatter(), rc.run_conf.project_root + '/paplot/' + sample_conf_name + '/index.html')
def paplot(input_file, output_file):
    
    # software version in index.html
    remark = gc.genomon_conf.get("paplot", "remarks")
    remark += "<ul>"
    
    for item in gc.genomon_conf.get("paplot", "software").split(","):
        key = item.split(":")[0].strip(" ").rstrip(" ")
        name = item.split(":")[1].strip(" ").rstrip(" ")
        try:
            version = gc.get_version(key)
        except Exception:
            print ("[WARNING] paplot: %s is not defined." % (key))
            continue
        
        remark += "<li>" + name + " " + version + "</li>"

    remark += "</ul>"

    arguments = {
        "inputs_qc": ",".join(paplot_inputs_starqc),
         "inputs_sv": ",".join(paplot_inputs_fusion),
         "output_dir": rc.run_conf.project_root + "/paplot/" + sample_conf_name,
         "title": gc.genomon_conf.get("paplot", "title"),
         "remarks": remark,
         "config_file": gc.genomon_conf.get("paplot", "config_file"),
    }
    singularity_params = {
        "image": gc.genomon_conf.get("paplot", "image"),
        "option": gc.genomon_conf.get("paplot", "singularity_option"),
        "bind": [
            rc.run_conf.project_root,
            gc.genomon_conf.get("paplot", "config_file")
        ]
    }
    r_paplot.task_exec(arguments, rc.run_conf.project_root + '/log/paplot', rc.run_conf.project_root + '/script/paplot', singularity_params)


