#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL # required to send email notifcations
#SBATCH --mail-user=vt520# required to send email notifcations - please replace <your_username> with your college login name or email address

# ****TO CHANGE path to conda****
export PATH=/vol/bitbucket/vt520/miniconda3/bin/:$PATH
#***********************************

source activate CATCaller
. /vol/cuda/11.1.0-cudnn8.0.4.30/setup.sh
TERM=vt100 # or TERM=xterm

as=$1 #train_signal_path
al=$2 #train_label_path
es=$3 #val_signal_path
el=$4 #val_label_path
model=$5 # base model for training


python /vol/bitbucket/vt520/CATCaller_master/train_litetr.py \
-as ${as} \
-al ${al} \
-es ${es} \
-el ${el}  \
--save_model ./ \
--store_model True > info_train.log 2>&1 \
--from_model ${model} #	/vol/bitbucket/vt520/CATCaller_master/model/model.2048.chkpt
