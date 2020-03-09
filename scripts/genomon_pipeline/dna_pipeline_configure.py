import os
import glob

## link the input fastq to project directory
#@ruffus.originate(linked_fastq_list)
def link_input_fastq(sample, output_file, genomon_conf, run_conf, sample_conf):
    fastq_dir = run_conf.project_root + '/fastq/' + sample
    fastq_prefix, ext = os.path.splitext(sample_conf.fastq[sample][0][0])
    # Todo
    # 1. should compare the timestamps between input and linked file
    # 2. check md5sum ?
    for (count, fastq_files) in enumerate(sample_conf.fastq[sample][0]):
        fastq_prefix, ext = os.path.splitext(fastq_files)
        if not os.path.exists(fastq_dir + '/'+str(count+1)+'_1'+ ext):
            os.symlink(sample_conf.fastq[sample][0][count], fastq_dir + '/'+str(count+1)+'_1'+ ext)
        if not os.path.exists(fastq_dir + '/'+str(count+1)+'_2'+ ext):
            os.symlink(sample_conf.fastq[sample][1][count], fastq_dir + '/'+str(count+1)+'_2'+ ext)

# split fastq
#@ruffus.subdivide([bam2fastq,  ], ruffus.formatter(), "{path[0]}/*_*.fastq_split", "{path[0]}")
def split_files(sample_name, input_files, output_files, target_dir, genomon_conf, run_conf, sample_conf):

    import genomon_pipeline.dna_resource.fastq_splitter  as dr_fastq_splitter 
    fastq_splitter = dr_fastq_splitter.Fastq_splitter(genomon_conf.get("split_fastq", "qsub_option"), run_conf.drmaa)

    sample_name = os.path.basename(target_dir)

    for oo in output_files:
        os.unlink(oo)

    split_lines = genomon_conf.get("split_fastq", "split_fastq_line_number")

    input_prefix, ext = os.path.splitext(input_files[0][0])
    arguments = {"lines": split_lines,
                 "fastq_filter": genomon_conf.get("split_fastq", "fastq_filter"),
                 "target_dir": target_dir,
                 "ext": ext}

    bind = [run_conf.project_root]
    if sample_name in sample_conf.fastq_src:
        bind.extend(sample_conf.fastq_src[sample_name])

    singularity_params = {
        "image": genomon_conf.get("split_fastq", "image"),
        "option": genomon_conf.get("split_fastq", "singularity_option"),
        "bind": bind,
    }
    fastq_splitter.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/'+ sample_name, singularity_params, 2)
   
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
#@ruffus.subdivide(split_files, ruffus.formatter(".+/(.+)/1_0000.fastq_split"), ruffus.add_inputs("{subpath[0][2]}/fastq/{subdir[0][0]}/2_0000.fastq_split"), "{subpath[0][2]}/bam/{subdir[0][0]}/{subdir[0][0]}_*.sorted.bam", "{subpath[0][2]}/fastq/{subdir[0][0]}", "{subpath[0][2]}/bam/{subdir[0][0]}")
def map_dna_sequence(input_files, output_files, input_dir, output_dir, genomon_conf, run_conf, sample_conf):
    import genomon_pipeline.dna_resource.bwa_align       as dr_bwa_align
    bwa_align = dr_bwa_align.Bwa_align(genomon_conf.get("bwa_mem", "qsub_option"), run_conf.drmaa)

    sample_name = os.path.basename(output_dir)

    all_line_num = 0
    with open(input_dir + "/fastq_line_num.txt") as in_handle:
        tmp_num = in_handle.read()
        all_line_num = int(tmp_num)
    split_lines = genomon_conf.get("split_fastq", "split_fastq_line_number")

    import math
    ans_quotient = int(math.floor(all_line_num / int(split_lines)))
    ans_remainder = all_line_num % int(split_lines)
    max_task_id = ans_quotient if ans_remainder == 0 else ans_quotient + 1

    arguments = {
        "input_dir": input_dir,
        "output_dir": output_dir,
        "sample_name": sample_name,
        "bwa_params": genomon_conf.get("bwa_mem", "bwa_params"),
        "ref_fa": genomon_conf.get("REFERENCE", "ref_fasta"),
    }
    singularity_params = {
        "image": genomon_conf.get("bwa_mem", "image"),
        "option": genomon_conf.get("bwa_mem", "singularity_option"),
        "bind": [run_conf.project_root, os.path.dirname(genomon_conf.get("REFERENCE", "ref_fasta"))],
    }
    
    bwa_align.task_exec(arguments, run_conf.project_root + '/log/' + sample_name , run_conf.project_root + '/script/' + sample_name, singularity_params, max_task_id) 

    for task_id in range(max_task_id):
        num = str(task_id).zfill(4)
        os.unlink(input_dir +'/1_'+str(num)+'.fastq_split')
        os.unlink(input_dir +'/2_'+str(num)+'.fastq_split')
        os.unlink(output_dir+'/'+sample_name+'_'+str(num)+'.bwa.sam')


# merge sorted bams into one and mark duplicate reads with biobambam
#@ruffus.collate(map_dna_sequence, ruffus.formatter(), "{subpath[0][2]}/bam/{subdir[0][0]}/{subdir[0][0]}.markdup.bam", "{subpath[0][2]}/bam/{subdir[0][0]}")
def markdup(input_files, output_file, output_dir, genomon_conf, run_conf, sample_conf):
    import genomon_pipeline.dna_resource.markduplicates  as dr_markduplicates 
    markduplicates = dr_markduplicates.Markduplicates(genomon_conf.get("markduplicates", "qsub_option"), run_conf.drmaa)

    sample_name = os.path.basename(output_dir)

    output_prefix, ext = os.path.splitext(output_file)

    input_bam_files = ""
    for input_file in input_files:
        input_bam_files = input_bam_files + " I=" + input_file

    arguments = {"out_prefix": output_prefix,
                 "input_bam_files": input_bam_files,
                 "out_bam": output_file}
    
    singularity_params = {
        "image": genomon_conf.get("markduplicates", "image"),
        "option": genomon_conf.get("markduplicates", "singularity_option"),
        "bind": [run_conf.project_root],
    }
        
    markduplicates.task_exec(
            arguments, run_conf.project_root + '/log/' + sample_name , 
            run_conf.project_root + '/script/'+ sample_name, singularity_params)

    for input_file in input_files:
        os.unlink(input_file)
        os.unlink(input_file + ".bai")

def bam_path_to_sample_id(bam_path, genomon_conf, run_conf, sample_conf):
    return os.path.basename(bam_path).replace(".markdup.bam", "")
def get_param(conf, section, param, default="", flag=True):
    if flag == False:
        return default
    if conf.has_option(section, param) == False:
        return default
    return conf.get(section, param)

# identify mutations
#@ruffus.follows( markdup )
#@ruffus.follows( link_import_bam )
#@ruffus.subdivide(markdup_bam_list, ruffus.formatter(), "{subpath[0][2]}/mutation/{subdir[0][0]}/{subdir[0][0]}.genomon_mutation.result.filt.txt", "{subpath[0][2]}/mutation/{subdir[0][0]}")
def identify_mutations(input_file, output_file, output_dir, genomon_conf, run_conf, sample_conf):

    import genomon_pipeline.dna_resource.mutation_call   as dr_mutation_call
    mutation_call = dr_mutation_call.Mutation_call(genomon_conf.get("mutation_call", "qsub_option"), run_conf.drmaa)
    import genomon_pipeline.dna_resource.mutation_merge  as dr_mutation_merge
    mutation_merge = dr_mutation_merge.Mutation_merge(genomon_conf.get("mutation_merge", "qsub_option"), run_conf.drmaa)

    sample_name = os.path.basename(output_dir)

    bind = [
        run_conf.project_root,
        os.path.dirname(genomon_conf.get("REFERENCE", "ref_fasta")),
        os.path.dirname(genomon_conf.get("REFERENCE", "interval_list")),
        os.path.dirname(genomon_conf.get("REFERENCE", "simple_repeat_tabix_db")),
    ]
    databases = []

    active_inhouse_normal_flag = get_param(genomon_conf, "annotation", "active_inhouse_normal_flag", default=False)
    inhouse_normal_tabix_db = get_param(genomon_conf, "REFERENCE", "inhouse_normal_tabix_db", flag=active_inhouse_normal_flag)
    databases.append(inhouse_normal_tabix_db)

    active_inhouse_tumor_flag = get_param(genomon_conf, "annotation", "active_inhouse_tumor_flag", default=False)
    inhouse_tumor_tabix_db = get_param(genomon_conf, "REFERENCE", "inhouse_tumor_tabix_db", flag=active_inhouse_tumor_flag)
    databases.append(inhouse_tumor_tabix_db)

    active_HGMD_flag = get_param(genomon_conf, "annotation", "active_HGMD_flag", default=False)
    HGMD_tabix_db = get_param(genomon_conf, "REFERENCE", "HGMD_tabix_db", flag=active_HGMD_flag)
    databases.append(HGMD_tabix_db)

    active_annovar_flag = get_param(genomon_conf, "annotation", "active_annovar_flag", default=False)
    databases.append(get_param(genomon_conf, "SOFTWARE", "annovar", flag=active_annovar_flag))
    databases.append(get_param(genomon_conf, "annotation", "annovar_database", flag=active_annovar_flag))

    active_hotspot_flag = get_param(genomon_conf, "hotspot","active_hotspot_flag", default=False)
    databases.append(get_param(genomon_conf, "REFERENCE","hotspot_db", flag=active_hotspot_flag))

    active_ExAC_flag = get_param(genomon_conf, "annotation", "active_ExAC_flag", default=False)
    databases.append(get_param(genomon_conf, "REFERENCE", "ExAC_tabix_db", flag=active_ExAC_flag))
        
    active_HGVD_2013_flag = get_param(genomon_conf, "annotation", "active_HGVD_2013_flag", default=False)
    databases.append(get_param(genomon_conf, "REFERENCE", "HGVD_2013_tabix_db", flag=active_HGVD_2013_flag))

    active_HGVD_2016_flag = get_param(genomon_conf, "annotation", "active_HGVD_2016_flag", default=False)
    databases.append(get_param(genomon_conf, "REFERENCE", "HGVD_2016_tabix_db", flag=active_HGVD_2016_flag))
    
    active_hotspot_flag = get_param(genomon_conf, "hotspot", "active_hotspot_flag", default=False)
    databases.append(get_param(genomon_conf, "REFERENCE", "hotspot_db", flag=active_hotspot_flag))

    for item in databases:
        if item == "":
            continue
        if os.path.isdir(item):
            bind.append(item)
        elif os.path.isfile(item):
            bind.append(os.path.dirname(item))
  
    arguments = {
        # fisher mutation
        "fisher_pair_params": genomon_conf.get("fisher_mutation_call", "pair_params"),
        "fisher_single_params": genomon_conf.get("fisher_mutation_call", "single_params"),
        # realignment filter
        "realignment_params": genomon_conf.get("realignment_filter","params"),
        # indel filter
        "indel_params": genomon_conf.get("indel_filter", "params"),
        # breakpoint filter
        "breakpoint_params": genomon_conf.get("breakpoint_filter","params"),
        # simplerepeat filter
        "simple_repeat_db": genomon_conf.get("REFERENCE", "simple_repeat_tabix_db"),
        # EB filter
        "eb_map_quality": genomon_conf.get("eb_filter","map_quality"),
        "eb_base_quality": genomon_conf.get("eb_filter","base_quality"),
        "filter_flags": genomon_conf.get("eb_filter","filter_flags"),
        "control_bam_list": input_file[2],
        # hotspot mutation caller
        "hotspot_database": genomon_conf.get("REFERENCE","hotspot_db"),
        "active_hotspot_flag": genomon_conf.get("hotspot","active_hotspot_flag"),
        "hotspot_params": genomon_conf.get("hotspot","params"),
        # original_annotations
        "active_inhouse_normal_flag": active_inhouse_normal_flag,
        "inhouse_normal_database":inhouse_normal_tabix_db,
        "active_inhouse_tumor_flag": active_inhouse_tumor_flag,
        "inhouse_tumor_database":inhouse_tumor_tabix_db,
        "active_HGVD_2013_flag": genomon_conf.get("annotation", "active_HGVD_2013_flag"),
        "HGVD_2013_database": genomon_conf.get("REFERENCE", "HGVD_2013_tabix_db"),
        "active_HGVD_2016_flag": genomon_conf.get("annotation", "active_HGVD_2016_flag"),
        "HGVD_2016_database": genomon_conf.get("REFERENCE", "HGVD_2016_tabix_db"),
        "active_ExAC_flag": genomon_conf.get("annotation", "active_ExAC_flag"),
        "ExAC_database": genomon_conf.get("REFERENCE", "ExAC_tabix_db"),
        "active_HGMD_flag": active_HGMD_flag,
        "HGMD_database": HGMD_tabix_db,
        # annovar
        "active_annovar_flag": genomon_conf.get("annotation", "active_annovar_flag"),
        "annovar": genomon_conf.get("SOFTWARE", "annovar"),
        "annovar_database": genomon_conf.get("annotation", "annovar_database"),
        "table_annovar_params": genomon_conf.get("annotation", "table_annovar_params"),
        "annovar_buildver": genomon_conf.get("annotation", "annovar_buildver"),
        # commmon
        "ref_fa": genomon_conf.get("REFERENCE", "ref_fasta"),
        "interval_list": genomon_conf.get("REFERENCE", "interval_list"),
        "disease_bam": input_file[0],
        "control_bam": input_file[1],
        "out_prefix": output_dir + '/' + sample_name,
    }
    interval_list = genomon_conf.get("REFERENCE", "interval_list")
    max_task_id = sum(1 for line in open(interval_list))
    
    if sample_name in sample_conf.bam_import_src:
        bind.extend(sample_conf.bam_import_src[sample_name])
    if input_file[1] != None:
        normal_sample = bam_path_to_sample_id(input_file[1])
        if normal_sample in sample_conf.bam_import_src:
             bind.extend(sample_conf.bam_import_src[normal_sample])
    if input_file[2] != None:
        for row in open(input_file[2]).readlines():
            panel_sample = bam_path_to_sample_id(row.rstrip())
            if panel_sample in sample_conf.bam_import_src:
                 bind.extend(sample_conf.bam_import_src[panel_sample])

    bind = list(set(bind))    
    singularity_params = {
        "image": genomon_conf.get("mutation_call", "image"),
        "option": genomon_conf.get("mutation_call", "singularity_option"),
        "bind": bind,
    }
    mutation_call.task_exec(
        arguments, 
        run_conf.project_root + '/log/' + sample_name, 
        run_conf.project_root + '/script/' + sample_name, 
        singularity_params, 
        max_task_id
    )
    
    arguments = {
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
        "pair_params": genomon_conf.get("mutation_util","pair_params"),
        "single_params": genomon_conf.get("mutation_util","single_params"),
        "active_hotspot_flag": genomon_conf.get("hotspot","active_hotspot_flag"),
        "hotspot_database": genomon_conf.get("REFERENCE","hotspot_db"),
        "meta_info_em": genomon_conf.get_meta_info(["mutation_call"]),
        "meta_info_m": genomon_conf.get_meta_info(["mutation_call"]),
        "meta_info_ema": genomon_conf.get_meta_info(["mutation_call"]),
        "meta_info_ma": genomon_conf.get_meta_info(["mutation_call"]),
        "out_prefix": output_dir + '/' + sample_name}

    singularity_params = {
        "image": genomon_conf.get("mutation_call", "image"),
        "option": genomon_conf.get("mutation_call", "singularity_option"),
        "bind": bind,
    }
    mutation_merge.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name, singularity_params)

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
#@ruffus.follows( link_import_bam )
#@ruffus.follows( markdup )
#@ruffus.transform(parse_sv_bam_list, ruffus.formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.junction.clustered.bedpe.gz")
def parse_sv(input_file, output_file, genomon_conf, run_conf, sample_conf):

    import genomon_pipeline.dna_resource.sv_parse        as dr_sv_parse
    sv_parse = dr_sv_parse.SV_parse(genomon_conf.get("sv_parse", "qsub_option"), run_conf.drmaa)

    dir_name = os.path.dirname(output_file)
    if not os.path.isdir(dir_name): os.mkdir(dir_name)
    sample_name = os.path.basename(dir_name)

    arguments = {
                 "input_bam": input_file,
                 "output_prefix": output_file.replace(".junction.clustered.bedpe.gz", ""),
                 "param": genomon_conf.get("sv_parse", "params"),
                }
    bind = [run_conf.project_root]
    if sample_name in sample_conf.bam_import_src:
        bind.extend(sample_conf.bam_import_src[sample_name])
    
    singularity_params = {
        "image": genomon_conf.get("sv_merge", "image"),
        "option": genomon_conf.get("sv_merge", "singularity_option"),
        "bind": bind,
    }
    sv_parse.task_exec(arguments, run_conf.project_root + '/log/' + sample_name , run_conf.project_root + '/script/' + sample_name, singularity_params)


# merge SV
#@ruffus.follows( parse_sv )
#@ruffus.transform(merge_bedpe_list, ruffus.formatter(".+/(?P<NAME>.+).control_info.txt"), "{subpath[0][2]}/sv/non_matched_control_panel/{NAME[0]}.merged.junction.control.bedpe.gz")
def merge_sv(input_files,  output_file, genomon_conf, run_conf, sample_conf):

    import genomon_pipeline.dna_resource.sv_merge        as dr_sv_merge
    sv_merge = dr_sv_merge.SV_merge(genomon_conf.get("sv_merge", "qsub_option"), run_conf.drmaa)

    arguments = {
                 "control_info": input_files[0],
                 "merge_output_file": output_file,
                 "param": genomon_conf.get("sv_merge", "params"),
                }
    bind = [
        run_conf.project_root,
    ]
    if input_files[0] != None:
        for row in open(input_files[0]).readlines():
            panel_sample = row.split("\t")[0]
            if panel_sample in sample_conf.bam_import_src:
                 bind.extend(sample_conf.bam_import_src[panel_sample])
    
    singularity_params = {
        "image": genomon_conf.get("sv_merge", "image"),
        "option": genomon_conf.get("sv_merge", "singularity_option"),
        "bind": bind,
    }
    sv_merge.task_exec(arguments, run_conf.project_root + '/log/sv_merge', run_conf.project_root + '/script/sv_merge', singularity_params)


# filt SV
#@ruffus.follows( merge_sv )
#@ruffus.transform(filt_bedpe_list, ruffus.formatter(), "{subpath[0][2]}/sv/{subdir[0][0]}/{subdir[0][0]}.genomonSV.result.filt.txt")
def filt_sv(input_files,  output_file, genomon_conf, run_conf, sample_conf):

    import genomon_pipeline.dna_resource.sv_filt         as dr_sv_filt
    sv_filt = dr_sv_filt.SV_filt(genomon_conf.get("sv_filt", "qsub_option"), run_conf.drmaa)

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

    arguments = {
                 "input_bam": run_conf.project_root + "/bam/" + sample_name + '/' + sample_name + ".markdup.bam",
                 "output_prefix": run_conf.project_root + "/sv/" + sample_name + '/' + sample_name,
                 "reference_genome": genomon_conf.get("REFERENCE", "ref_fasta"),
                 "param": filt_param,
                 "meta_info": genomon_conf.get_meta_info(["sv_parse", "sv_merge", "sv_filt"]),
                 "sv_utils_param": genomon_conf.get("sv_filt", "sv_utils_params"),
                }
    bind = [
        run_conf.project_root,
        os.path.dirname(genomon_conf.get("REFERENCE", "ref_fasta"))
    ]
    if sample_name in sample_conf.bam_import_src:
        bind.extend(sample_conf.bam_import_src[sample_name])
    if complist[1] != None and complist[1] in sample_conf.bam_import_src: 
        bind.extend(sample_conf.bam_import_src[complist[1]])

    singularity_params = {
        "image": genomon_conf.get("sv_filt", "image"),
        "option": genomon_conf.get("sv_filt", "singularity_option"),
        "bind": bind,
    }
    sv_filt.task_exec(arguments, run_conf.project_root + '/log/' + sample_name, run_conf.project_root + '/script/' + sample_name, singularity_params)

def configure(genomon_conf, run_conf, sample_conf):

    import shutil
    
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
        #linked_fastq_list.append([link_fastq_arr1,link_fastq_arr2])
        link_input_fastq(sample, [link_fastq_arr1,link_fastq_arr2], genomon_conf, run_conf, sample_conf)
        
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
    
    genomon_conf_name, genomon_conf_ext = os.path.splitext(os.path.basename(run_conf.genomon_conf_file))
    sample_conf_name, sample_conf_ext = os.path.splitext(os.path.basename(run_conf.sample_conf_file))
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
    