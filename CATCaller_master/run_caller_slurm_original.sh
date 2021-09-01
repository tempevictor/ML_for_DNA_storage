#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL # required to send email notifcations
#SBATCH --mail-user=vt520# required to send email notifcations - please replace <your_username> with your college login name or email address


# ****TO CHANGE path to conda bin****
export PATH=/vol/bitbucket/vt520/miniconda3/bin/:$PATH 
#*************************************
source activate CATCaller
. /vol/cuda/11.1.0-cudnn8.0.4.30/setup.sh
TERM=vt100 # or TERM=xterm

model_path=$1   
fast5_dir_root=$2  
signal_window_length=$3  # 2048
basecalled_dir=$4 
# ****TO CHANGE, path to store temprary files****
tmp_records_dir="/vol/bitbucket/vt520/tmp_data_dir"
#*************************************
for p in `ls ${fast5_dir_root}`
do
    mkdir -p ${basecalled_dir}/${p}
    fast5_path="${fast5_dir_root}/${p}"
    time_out_path="${basecalled_dir}/${p}/time.txt"

    starttime=`date +'%Y-%m-%d %H:%M:%S'`
    # data preprocessing

    # ****TO CHANGE path to CATCaller_master****
    python /vol/bitbucket/vt520/CATCaller_master/generate_dataset/inference_data_original.py -fast5 ${fast5_path} -records_dir ${tmp_records_dir} -raw_len ${signal_window_length} -threads 40
    #*************************************
    printf "${p} preprocessing done"

    # caller
    # ****TO CHANGE path to CATCaller_master****
    python /vol/bitbucket/vt520/CATCaller_master/caller.py -model ${model_path} -records_dir ${tmp_records_dir} -output ${basecalled_dir}/${p} -half 
    #*************************************
    endtime=`date +'%Y-%m-%d %H:%M:%S'`
    echo "${p} finished!"
    start_seconds=$(date --date="$starttime" +%s);
    end_seconds=$(date --date="$endtime" +%s);
    echo "${p} running time： "$((end_seconds-start_seconds))"s" > ${time_out_path}
    #delete ${tmp_records_dir}
    rm -rf ${tmp_records_dir}
done

