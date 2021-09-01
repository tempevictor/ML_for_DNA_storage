#!/bin/bash


as=$1 #train_signal_path folder
al=$2 #train_label_path folder
es=$3 #val_signal_path fodler
el=$4 #val_label_path fodler
model=$5 # base model for training

# ****TO CHANGE path to CATCaller_master****
python path/to/CATCaller_master/train_litetr.py \
# ******************************************

-as ${as} \
-al ${al} \
-es ${es} \
-el ${el}  \
--save_model ./ \
--store_model True > info_train.log 2>&1 \
--from_model ${model}
