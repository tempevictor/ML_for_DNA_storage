#!/bin/bash


basecall_dir=$1 #directory containing the fasta outputs called out.fasta
ref_file=$2 # reference must be a fasta file
threads=200

# compute the differents performance metrics 
for p in `ls ${basecall_dir}`
do
    result_dir="${basecall_dir}"/"${p}"/"${p}"_result_120
    mkdir -p ${result_dir}
    read_alignment="${result_dir}"/"${p}"_reads.paf
    read_data="${result_dir}"/"${p}"_reads.tsv
    reference="${ref_file}"
    basecall_name="${basecall_dir}"/"${p}"/"out.fasta"  #need change
    echo "reads alignment: minimap2..."
    printf "\n"
    minimap2 -x map-ont -t ${threads} -c ${reference} ${basecall_name} > ${read_alignment}
    
    echo "calculate metrics..."
    printf "\n"
    python read_length_identity.py ${basecall_name} ${read_alignment} > ${read_data}
    python performance.py ${read_data} ${result_dir}/performance.txt ${basecall_dir}/${p}/caller_time.out
    echo "${p} finished!"
done


