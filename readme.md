
# ML for DNA storage
#### Victor Tempe

For the MSc project: Machine Learning for DNA storage. 

Code modified from: https://github.com/lvxuan96/CATCaller
The code can be installed or a read-to-use vesion is directly availble on Imperial doc servers.



Repository oragnisation:
* CATCaller_master: contains the basecaller
* scripts: various python files and notebooks for tests or data manipulation
* performance: compute the performance of the basecaller
* commands_preprocess.txt: example commands to create the dataset from fitlered fast5 fast5 reads (need to change paths to be used)



#### Installation

We had the network work with cuda 11.1.0-cudnn8.0.4.30
The commands create a new conda env and install the necessary packages


```angular2
conda create --name CATCaller
source activate CATCaller
pip install ninja 
pip install fast-ctc-decode
pip install --no-cache-dir tensorflow-gpu==1.13.2
cd CATCaller_master
git clone --recursive https://github.com/parlance/ctcdecode.git
cd ctcdecode && pip install .
cd ..
conda env update --file ./requirements/environment.yaml --prune
```



If the Imeprial servers with slurm is inteded to be used. run install_dynamicconv.sh with sbatch from the dynamicconv_layer folder (it should use cuda 11.1.0-cudnn8.0.4.30), Careful: change the apropriate paths in the file
```angular2
cd dynamicconv_layer
sbatch ./../requirements/install_dynamicconv.sh
cd ..
```

Otherwise, on a computer with gpu and right cuda version

```angular2
cd dynamicconv_layer
python cuda_function_gen.py
python setup.py build
python setup.py install --user
cd ..
```




#### Launch in inference mode
 From CATCaller_master:

```angular2
bash run_caller_trim.sh <model file> <fast5 folder> 2048 <basecalled_dir>
```

`model file`: are stored unde CATCaller_master/model/ Use fine-tuned.chkpt for the fine-tuned version and model.2048.chkpt for the original one.   
`fast5 folder`: directory of fast5 files. Inside the directory, there cannot be fast5 files directly but folders containing fast5 files. The files must be multi-fast5 files. If that is not the case, use https://github.com/nanoporetech/ont_fast5_api
`basecalled_dir`: the output directory that will host the results.  




From Imperial cluster in the CATCaller_master fodler:

```angular2
ssh gpucluster.doc.ic.ac.uk
sbatch run_caller_slurm.sh <model file> <fast5 folder> 2048 <basecalled_dir>
```

To try directly on doc servers
The basecalled result will appear in /vol/bitbucket/vt520/demo/results
```angular2
ssh gpucluster.doc.ic.ac.uk
cd /vol/bitbucket/vt520/demo/
sbatch ./../CATCaller_master/run_caller_slurm.sh ./../CATCaller_master/model/fine-tuned.chkpt ./../data/3XR6/datasets/dataset_resquiggle_256/fast5s_multi/test0/ 2048 ./results/
```


Note:
* To launch the original CATCaller, use the same commands but replacing the bash scripts with run_caller_slurm_original.sh and run_caller_trim_original.sh. The input must be single-fast5 files.
* When launching CATCaller with slurm, make sure beforehand that no conda environment is activated



### Launch fine-tuning
Uses the npy files produced after pre-processing
```angular2
bash train_litetr.sh <as> <al> <es> <el> <model>
```
`as`: train signal path folder
`al`: train label path folder
`es`: val signal path folder
`el`: val label path folder
`model`: base model path


Example on doc servers
```angular2
ssh gpucluster.doc.ic.ac.uk
cd /vol/bitbucket/vt520/demo/fine_tune/
sbatch ./../../CATCaller_master/train_litetr_slurm.sh ./../../data/3XR6/datasets/dataset_resquiggle_256/signal_output_train \
./../../data/3XR6/datasets/dataset_resquiggle_256/label_output_train \
./../../data/3XR6/datasets/dataset_resquiggle_256/signal_output_val \
./../../data/3XR6/datasets/dataset_resquiggle_256/label_output_val \
./../../CATCaller_master/model/model.2048.chkpt
```

During fine-tuning, info_train.log allows to see the progress



#### Performance
Minimap2 must be installed first to run the performance analysis. 
The results directory must have the same architecture as the output directory of the basecaller in inference mode.
To compute the performance of the result of a basecaller use:
```angular2
bash metrics.sh <basecall_dir> <ref_file>
```

`basecall_dir`: the output directory containing the results of the basecaller in inference mode. Inside it, fasta outputs are called out.fasta
`ref_file`: reference file, must be a fasta file


The results of interest are in a file called performance.txt


#### Data

Original data including 3XR6: https://github.com/helixworks-technologies/dos

Toy dataset containing files ready to use for training and fitlered fast5s: https://drive.google.com/drive/folders/10KBVDsusPEa-cnV_SHONaf9jDzupnoc5?usp=sharing




