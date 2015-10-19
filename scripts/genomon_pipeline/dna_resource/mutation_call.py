#! /usr/bin/env python

from genomon_pipeline.stage_task import *

class Mutation_call(Stage_task):

    task_name = "mutation_call"

    script_template = """
#!/bin/bash
#
#$ -S /bin/bash         # set shell in UGE
#$ -cwd                 # execute at the submitted dir
#$ -e {log}             # log file directory
#$ -o {log}             # log file directory
pwd                     # print current working directory
hostname                # print hostname
date                    # print date
set -xv

# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export LD_LIBRARY_PATH={ld_library_path}
export PYTHONPATH={pythonpath}

TASK_ID={task_id}
REGION=`head -n $TASK_ID {interval_list} | tail -n 1`
OUT_FISHER={out_prefix}.fisher_mutations.${{TASK_ID}}.txt
OUT_REALIGNMENT={out_prefix}.realignment_mutations.${{TASK_ID}}.txt
OUT_INDEL={out_prefix}.indel_mutations.${{TASK_ID}}.txt
OUT_BREAKPOINT={out_prefix}.breakpoint_mutations.${{TASK_ID}}.txt
OUT_SIMPLE_REPEAT={out_prefix}.simplerepeat_mutations.${{TASK_ID}}.txt
OUT_EB={out_prefix}.ebfilter_mutations.${{TASK_ID}}.txt
OUT_MUTATIONS={out_prefix}_mutations_candidate.${{TASK_ID}}

if [ _{control_bam} = "_None" ]; then 
    {fisher} single -R $REGION -o $OUT_FISHER  --ref_fa {ref_fa} --mapping_quality {map_quality} --base_quality {base_quality}  --min_allele_freq {min_allele_freq} --post_10_q {post_10_q} --min_depth {min_depth} -1 {disease_bam} --samtools_path {samtools} || exit $?

    {mutfilter} realignment --tumor_min_mismatch {realign_min_mismatch} --normal_max_mismatch {realign_max_mismatch} --score_difference {realign_score_diff} --window_size {realign_window_size} --max_depth {realign_max_depth} --target_mutation_file $OUT_FISHER -1 {disease_bam} --output $OUT_REALIGNMENT --ref_genome {ref_fa} --blat_path {blat} || exit $?

    {mutfilter} simplerepeat --target_mutation_file $OUT_REALIGNMENT --output $OUT_SIMPLE_REPEAT --simple_repeat_db {simple_repeat_db} || exit $?

else
    {fisher} comparison -R $REGION -o $OUT_FISHER  --ref_fa {ref_fa} --mapping_quality {map_quality} --base_quality {base_quality}  --min_allele_freq {min_allele_freq} --max_allele_freq {max_allele_freq} --min_depth {min_depth} --fisher_value {fisher_thres} -2 {control_bam} -1 {disease_bam} --samtools_path {samtools} || exit $?

    {mutfilter} realignment --tumor_min_mismatch {realign_min_mismatch} --normal_max_mismatch {realign_max_mismatch} --score_difference {realign_score_diff} --window_size {realign_window_size} --max_depth {realign_max_depth} --target_mutation_file $OUT_FISHER -1 {disease_bam} -2 {control_bam} --output $OUT_REALIGNMENT --ref_genome {ref_fa} --blat_path {blat} || exit $?

    {mutfilter} indel --search_length {indel_search_length} --neighbor {indel_neighbor} --base_qual {indel_base_quality} --min_depth {indel_min_depth} --min_mismatch {indel_min_mismatch} --af_thres {indel_min_allele_freq} --target_mutation_file $OUT_REALIGNMENT -2 {control_bam} --output $OUT_INDEL || exit $?

    {mutfilter} breakpoint --max_depth {bp_max_depth} --min_clip_size {bp_min_clip_size} --junc_num_thres {bp_junc_num_thres} --mapq_thres {bp_map_quality} --target_mutation_file $OUT_INDEL -2 {control_bam} --output $OUT_BREAKPOINT || exit $?

    {mutfilter} simplerepeat --target_mutation_file $OUT_BREAKPOINT --output $OUT_SIMPLE_REPEAT --simple_repeat_db {simple_repeat_db} || exit $?

fi

if [ _{control_bam_list} != "_None" ]; then 
    {EBFilter} -f anno -q {eb_map_quality} -Q {eb_base_quality} $OUT_SIMPLE_REPEAT {disease_bam} {control_bam_list} $OUT_EB || exit $?
else
    cp $OUT_SIMPLE_REPEAT $OUT_EB
fi

if [ _{active_annovar_flag} = "_True" ];then
    {annovar}/table_annovar.pl --outfile $OUT_MUTATIONS {table_annovar_params} $OUT_EB {annovar}/humandb || exit $?
else
    cp $OUT_EB {out_prefix}_mutations_candidate.${{TASK_ID}}.hg19_multation.txt
fi

"""
    def __init__(self, qsub_option, script_dir):
        super(Mutation_call, self).__init__(qsub_option, script_dir)

