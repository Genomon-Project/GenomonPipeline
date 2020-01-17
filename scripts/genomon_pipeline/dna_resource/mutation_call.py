#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Mutation_call(st.Stage_task):

    task_name = "mutation_call"

    script_template = """
#!/bin/bash
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv
set -o pipefail

SAMTOOLS=/usr/local/bin/samtools
BLAT=/tools/userApps/bin/blat

REGION=`head -n $SGE_TASK_ID {interval_list} | tail -n 1`

if [ _{control_bam} = "_None" ]; then 
    fisher single -R ${{REGION}} -o {out_prefix}.fisher_mutations.${{SGE_TASK_ID}}.txt --ref_fa {ref_fa} -1 {disease_bam} --samtools_path ${{SAMTOOLS}} {fisher_single_params} || exit $?

    mutfilter realignment --target_mutation_file {out_prefix}.fisher_mutations.${{SGE_TASK_ID}}.txt -1 {disease_bam} --output {out_prefix}.realignment_mutations.${{SGE_TASK_ID}}.txt --ref_genome {ref_fa} --blat_path ${{BLAT}} {realignment_params} || exit $?

    mutfilter simplerepeat --target_mutation_file {out_prefix}.realignment_mutations.${{SGE_TASK_ID}}.txt --output {out_prefix}.simplerepeat_mutations.${{SGE_TASK_ID}}.txt --simple_repeat_db {simple_repeat_db} || exit $?

else
    fisher comparison -R ${{REGION}} -o {out_prefix}.fisher_mutations.${{SGE_TASK_ID}}.txt --ref_fa {ref_fa} -2 {control_bam} -1 {disease_bam} --samtools_path ${{SAMTOOLS}} {fisher_pair_params} || exit $?

    realignment_input={out_prefix}.fisher_mutations.${{SGE_TASK_ID}}.txt

    if [ _{active_hotspot_flag} = "_True" ]; then 

        hotspotCall {hotspot_params} {disease_bam} {control_bam} {out_prefix}.hotspot_mutations.${{SGE_TASK_ID}}.txt {hotspot_database}.${{SGE_TASK_ID}} || exit $?

        mutil merge_hotspot -i {out_prefix}.hotspot_mutations.${{SGE_TASK_ID}}.txt -f {out_prefix}.fisher_mutations.${{SGE_TASK_ID}}.txt -o {out_prefix}.fisher_hotspot_mutations.${{SGE_TASK_ID}}.txt --hotspot_header || exit $?

        realignment_input={out_prefix}.fisher_hotspot_mutations.${{SGE_TASK_ID}}.txt

    fi

    mutfilter realignment --target_mutation_file ${{realignment_input}} -1 {disease_bam} -2 {control_bam} --output {out_prefix}.realignment_mutations.${{SGE_TASK_ID}}.txt --ref_genome {ref_fa} --blat_path ${{BLAT}} {realignment_params} || exit $?

    mutfilter indel --target_mutation_file {out_prefix}.realignment_mutations.${{SGE_TASK_ID}}.txt -2 {control_bam} --output {out_prefix}.indel_mutations.${{SGE_TASK_ID}}.txt --samtools_path ${{SAMTOOLS}} {indel_params} || exit $?

    mutfilter breakpoint --target_mutation_file {out_prefix}.indel_mutations.${{SGE_TASK_ID}}.txt -2 {control_bam} --output {out_prefix}.breakpoint_mutations.${{SGE_TASK_ID}}.txt {breakpoint_params} || exit $?

    mutfilter simplerepeat --target_mutation_file {out_prefix}.breakpoint_mutations.${{SGE_TASK_ID}}.txt --output {out_prefix}.simplerepeat_mutations.${{SGE_TASK_ID}}.txt --simple_repeat_db {simple_repeat_db} || exit $?
fi

if [ _{control_bam_list} != "_None" ]; then 
    EBFilter --loption --region ${{REGION}} -f anno -q {eb_map_quality} -Q {eb_base_quality} --ff {filter_flags} {out_prefix}.simplerepeat_mutations.${{SGE_TASK_ID}}.txt {disease_bam} {control_bam_list} {out_prefix}.ebfilter_mutations.${{SGE_TASK_ID}}.txt || exit $?
else
    cp {out_prefix}.simplerepeat_mutations.${{SGE_TASK_ID}}.txt {out_prefix}.ebfilter_mutations.${{SGE_TASK_ID}}.txt
fi

if [ _{active_inhouse_normal_flag} = "_True" ]; then 
    mutanno mutation -t {out_prefix}.ebfilter_mutations.${{SGE_TASK_ID}}.txt -o {out_prefix}.inhouse_normal.${{SGE_TASK_ID}}.txt -d {inhouse_normal_database} || exit $?
else
    cp {out_prefix}.ebfilter_mutations.${{SGE_TASK_ID}}.txt {out_prefix}.inhouse_normal.${{SGE_TASK_ID}}.txt
fi

if [ _{active_inhouse_tumor_flag} = "_True" ]; then 
    mutanno mutation -t {out_prefix}.inhouse_normal.${{SGE_TASK_ID}}.txt -o {out_prefix}.inhouse_tumor.${{SGE_TASK_ID}}.txt -d {inhouse_tumor_database} || exit $?
else
    cp {out_prefix}.inhouse_normal.${{SGE_TASK_ID}}.txt {out_prefix}.inhouse_tumor.${{SGE_TASK_ID}}.txt
fi

if [ _{active_HGVD_2013_flag} = "_True" ]; then 
    mutanno mutation -t {out_prefix}.inhouse_tumor.${{SGE_TASK_ID}}.txt -o {out_prefix}.HGVD_2013.${{SGE_TASK_ID}}.txt -d {HGVD_2013_database} -c 5 || exit $?
else
    cp {out_prefix}.inhouse_tumor.${{SGE_TASK_ID}}.txt {out_prefix}.HGVD_2013.${{SGE_TASK_ID}}.txt
fi

if [ _{active_HGVD_2016_flag} = "_True" ]; then 
    mutanno mutation -t {out_prefix}.HGVD_2013.${{SGE_TASK_ID}}.txt -o {out_prefix}.HGVD_2016.${{SGE_TASK_ID}}.txt -d {HGVD_2016_database} -c 5 || exit $?
else
    cp {out_prefix}.HGVD_2013.${{SGE_TASK_ID}}.txt {out_prefix}.HGVD_2016.${{SGE_TASK_ID}}.txt
fi

if [ _{active_ExAC_flag} = "_True" ]; then 
    mutanno mutation -t {out_prefix}.HGVD_2016.${{SGE_TASK_ID}}.txt -o {out_prefix}.ExAC.${{SGE_TASK_ID}}.txt -d {ExAC_database} -c 8 || exit $?
else
    cp {out_prefix}.HGVD_2016.${{SGE_TASK_ID}}.txt {out_prefix}.ExAC.${{SGE_TASK_ID}}.txt
fi

if [ _{active_HGMD_flag} = "_True" ]; then 
    mutanno mutation -t {out_prefix}.ExAC.${{SGE_TASK_ID}}.txt -o {out_prefix}.HGMD.${{SGE_TASK_ID}}.txt -d {HGMD_database} -c 7 || exit $?
else
    cp {out_prefix}.ExAC.${{SGE_TASK_ID}}.txt {out_prefix}.HGMD.${{SGE_TASK_ID}}.txt
fi

if [ _{active_annovar_flag} = "_True" ];then
    {annovar}/table_annovar.pl --outfile {out_prefix}_mutations_candidate.${{SGE_TASK_ID}} {table_annovar_params} {out_prefix}.HGMD.${{SGE_TASK_ID}}.txt {annovar_database} || exit $?
else
    cp {out_prefix}.HGMD.${{SGE_TASK_ID}}.txt {out_prefix}_mutations_candidate.${{SGE_TASK_ID}}.{annovar_buildver}_multianno.txt
fi

"""
    def __init__(self, qsub_option, script_dir):
        super(Mutation_call, self).__init__(qsub_option, script_dir)


