#! /usr/bin/env python

import genomon_pipeline.stage_task as st

class Mutation_merge(st.Stage_task):

    task_name = "mutation_merge"

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

# set python environment
export PYTHONHOME={pythonhome}
export PATH=$PYTHONHOME/bin:$PATH
export LD_LIBRARY_PATH={ld_library_path}
export PYTHONPATH={pythonpath}

mut_header=""

if [ _{control_bam} = "_None" ]
then
    mut_header="Chr,Start,End,Ref,Alt,depth,variantNum,bases,A_C_G_T,misRate,strandRatio,10%_posterior_quantile,posterior_mean,90%_posterior_quantile,readPairNum,variantPairNum,otherPairNum,10%_posterior_quantile(realignment),posterior_mean(realignment),90%_posterior_quantile(realignment),simple_repeat_pos,simple_repeat_seq,P-value(EBCall)"
else
    mut_header="Chr,Start,End,Ref,Alt,depth_tumor,variantNum_tumor,depth_normal,variantNum_normal,bases_tumor,bases_normal,A_C_G_T_tumor,A_C_G_T_normal,misRate_tumor,strandRatio_tumor,misRate_normal,strandRatio_normal,P-value(fisher)"

    if [ _{active_hotspot_flag} = "_True" ]
    then 
        mut_header="${{mut_header}},score(hotspot)"
    fi

    mut_header="${{mut_header}},readPairNum_tumor,variantPairNum_tumor,otherPairNum_tumor,readPairNum_normal,variantPairNum_normal,otherPairNum_normal,P-value(fisher_realignment),indel_mismatch_count,indel_mismatch_rate,bp_mismatch_count,distance_from_breakpoint,simple_repeat_pos,simple_repeat_seq,P-value(EBCall)"
fi

if [ _{control_bam_list} = "_None" ]
then
    mut_header=`echo ${{mut_header}} | awk -F"," -v OFS="," '{{$NF=""; sub(/.$/,""); print $0}}'` || exit $?
fi

print_header=""

if [ _{active_annovar_flag} = _True ]
then
    tmp_header=`head -n 1 {out_prefix}_mutations_candidate.1.{annovar_buildver}_multianno.txt | awk -F"\t" -v OFS="\t" '{{$NF=""; sub(/.$/,""); print $0}}'` || exit $?
    print_header=${{tmp_header}}

    tmp_header=`echo $mut_header | cut -d "," -f 6- | tr "," "\t"` || exit $?
    print_header=${{print_header}}'\t'${{tmp_header}}
else
    tmp_header=`echo $mut_header | tr "," "\t"` || exit $?
    print_header=${{tmp_header}}
fi

if [ _{active_inhouse_normal_flag} = "_True" ]
then
    inhouse_normal_header="inhouseSNV"
    print_header=${{print_header}}'\t'${{inhouse_normal_header}} || exit $?
fi

if [ _{active_inhouse_tumor_flag} = "_True" ]
then
    inhouse_tumor_header="inhouseCandidateSNV"
    print_header=${{print_header}}'\t'${{inhouse_tumor_header}} || exit $?
fi

if [ _{active_HGVD_2013_flag} = "_True" ]
then
    HGVD_header="HGVD_20131010:#Sample,HGVD_20131010:Filter,HGVD_20131010:NR,HGVD_20131010:NA,HGVD_20131010:Frequency(NA/(NA+NR))"
    tmp_header=`echo $HGVD_header | tr "," "\t"` || exit $?
    print_header=${{print_header}}'\t'${{tmp_header}} || exit $?
fi

if [ _{active_HGVD_2016_flag} = "_True" ]
then
    HGVD_header="HGVD_20160412:#Sample,HGVD_20160412:Filter,HGVD_20160412:NR,HGVD_20160412:NA,HGVD_20160412:Frequency(NA/(NA+NR))"
    tmp_header=`echo $HGVD_header | tr "," "\t"` || exit $?
    print_header=${{print_header}}'\t'${{tmp_header}} || exit $?
fi

if [ _{active_ExAC_flag} = "_True" ]
then
    ExAC_header="ExAC:Filter,ExAC:AC_Adj,ExAC:AN_Adj,ExAC:Frequency(AC_Adj/AN_Adj),ExAC:AC_POPMAX,ExAC:AN_POPMAX,ExAC:Frequency(AC_POPMAX/AN_POPMAX),ExAC:POPMAX"
    tmp_header=`echo $ExAC_header | tr "," "\t"` || exit $?
    print_header=${{print_header}}'\t'${{tmp_header}} || exit $?
fi

if [ _{active_HGMD_flag} = "_True" ]
then
    HGMD_header="HGMD_201504:ID,HGMD_201504:CLASS,HGMD_201504:GENE,HGMD_201504:STRAND,HGMD_201504:DNA,HGMD_201504:PROT,HGMD_201504:PHEN"
    tmp_header=`echo $HGMD_header | tr "," "\t"` || exit $?
    print_header=${{print_header}}'\t'${{tmp_header}} || exit $?
fi

if [ _{control_bam_list} = "_None" ]
then
    if [ _{active_annovar_flag} = "_True" ]
    then
        echo -e "{meta_info_ma}" > {out_prefix}.genomon_mutation.result.txt || exit $?
    else
        echo -e "{meta_info_m}" > {out_prefix}.genomon_mutation.result.txt || exit $?
    fi
else
    if [ _{active_annovar_flag} = "_True" ]
    then
        echo -e "{meta_info_ema}" > {out_prefix}.genomon_mutation.result.txt || exit $?
    else
        echo -e "{meta_info_em}" > {out_prefix}.genomon_mutation.result.txt || exit $?
    fi
fi

echo "$print_header" >> {out_prefix}.genomon_mutation.result.txt || exit $?

for i in `seq 1 1 {filecount}`
do
    if [ _{active_annovar_flag} = "_True" ]
    then
        awk 'NR>1 {{print}}' {out_prefix}_mutations_candidate.${{i}}.{annovar_buildver}_multianno.txt >> {out_prefix}.genomon_mutation.result.txt || exit $?
    else
        cat {out_prefix}_mutations_candidate.${{i}}.{annovar_buildver}_multianno.txt >> {out_prefix}.genomon_mutation.result.txt || exit $?
    fi
done

if [ _{control_bam} = "_None" ]
then 
    if [ _{active_hotspot_flag} = "_True" ]
    then 
        {mutil} filter -i {out_prefix}.genomon_mutation.result.txt -o {out_prefix}.genomon_mutation.result.filt.txt {single_params} --hotspot_db {hotspot_database}
    else
        {mutil} filter -i {out_prefix}.genomon_mutation.result.txt -o {out_prefix}.genomon_mutation.result.filt.txt {single_params}
    fi
else
    if [ _{active_hotspot_flag} = "_True" ]
    then 
        {mutil} filter -i {out_prefix}.genomon_mutation.result.txt -o {out_prefix}.genomon_mutation.result.filt.txt {pair_params} --hotspot_db {hotspot_database}
    else
        {mutil} filter -i {out_prefix}.genomon_mutation.result.txt -o {out_prefix}.genomon_mutation.result.filt.txt {pair_params}
    fi
fi

"""
    def __init__(self, qsub_option, script_dir):
        super(Mutation_merge, self).__init__(qsub_option, script_dir)


