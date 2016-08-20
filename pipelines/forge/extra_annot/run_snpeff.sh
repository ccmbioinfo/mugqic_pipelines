#!/bin/bash

kit_arr=( illumina agilent nimblegen )
updir_arr=( cheo_wbed cheo_indivbed )

base_dir=/hpf/largeprojects/ccmbio/kng/mcgill_jacek

for updir in ${updir_arr[@]}
do
    for kit in ${kit_arr[@]}
    do
	direct=${base_dir}/${updir}/${kit}
	for path in ${direct}/*
	do
	    output_dir=${path}/output
	    if [[ ! -d ${output_dir} ]]
	    then
		echo "Uh oh. Output directory ${output_dir} doesn't exist. Exiting."
		continue;
	    fi
	    sample_name=$(ls ${output_dir}/alignment)
	    snpeff_dir=${output_dir}/job_output/snpeff

	    mkdir -p ${snpeff_dir}

	    qsub -W group_list=ccm \
	    	-l vmem=10g \
	    	-joe -o ${snpeff_dir}/snpeff.o \
	    	-N ${sample_name}_snpeff \
	    	-v sample_name=${sample_name},output_dir=${output_dir} \
	    	./snpeff_call.sh
	done
    done
done