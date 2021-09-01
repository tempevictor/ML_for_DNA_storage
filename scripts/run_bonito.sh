#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL # required to send email notifcations
#SBATCH --mail-user=vt520# required to send email notifcations - please replace <your_username> with your college login name or email address

export PATH=/vol/bitbucket/vt520/miniconda3/bin/:$PATH
source activate base
. /vol/cuda/10.2.89-cudnn7.6.4.38/setup.sh
TERM=vt100 # or TERM=xterm

fast5_dir=$1
output_name=$2
chunk_size=$3

bonito download --models
bonito basecaller --fastq --chunksize ${chunk_size} dna_r9.4.1 ${fast5_dir} > ${output_name}
