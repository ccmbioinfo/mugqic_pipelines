#!/bin/bash

# Should be .../output
if [[ -z ${output_dir} ]]
then
    echo "Output directory missing. Exiting."
    exit 1
fi

if [[ -z ${sample_name} ]]
then
    echo "Sample name missing. Exiting."
    exit 1
fi

module load java/1.8.0_65
module load snpEff/4.11

java -Xmx4g -Xms2g \
    -jar /hpf/tools/centos6/snpEff/4.11/snpEff.jar \
    -chr chr \
    -v GRCh38.81 \
    -dataDir /hpf/tools/centos6/snpEff/source/data/ \
    -csvStats ${output_dir}/annotation/${sample_name}/${sample_name}.snpeff.stats.csv \
    -stats ${output_dir}/annotation/${sample_name}/${sample_name}.snpeff.stats.html \
    ${output_dir}/variants/${sample_name}.flt.vcf \
    > ${output_dir}/annotation/${sample_name}/${sample_name}.snpeff.vcf

echo "Finished snpeff for ${sample_name}!"