#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL # required to send email notifcations
#SBATCH --mail-user=vt520# required to send email notifcations - please replace <your_username> with your college login name or email address


# note: replace /vol/bitbucket/vt520/ with apropriate path

# **** CHANGE path to conda ****
export PATH=/vol/bitbucket/vt520/miniconda3/bin/:$PATH
#********************************
source activate CATCaller
source /vol/cuda/11.1.0-cudnn8.0.4.30/setup.sh
TERM=vt100 # or TERM=xterm


# **** CHANGE path to CATCaller_master ****
python /path/to/CATCaller_master/dynamicconv_layer/cuda_function_gen.py
python /path/to/CATCaller_master/dynamicconv_layer/setup.py build
python /path/to/CATCaller_master/dynamicconv_layer/setup.py install --user
#********************************
