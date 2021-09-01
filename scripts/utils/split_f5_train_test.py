import numpy as np
from shutil import copyfile
import shutil
from six import b
from statsmodels import robust
import h5py
import os
import sys
from tqdm import tqdm
import glob
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()



''' Dic of reads that align to several references (dic_several)
    and dic reads_id to reference (dic_idToRef), read_id with several
    references are dropped'''

def get_reads_ref(blastnFileList):
    dic_idToRef = {}
    several_id_to_remove = set()
    for blastnFile in blastnFileList:
        with open(blastnFile, 'r') as f:
            current_querry = ''
            for line in f:
                if line[0:6] == "Query=":
                    current_querry = line.split()[1]
                elif line[0] == '>':
                    read_id = line[1:].split()[0]
                    if read_id not in dic_idToRef:
                        dic_idToRef[read_id] = current_querry
                    else:
                        several_id_to_remove.add(read_id)
    # remove read_ids in the dic that have several references linked to it
    for read_id in several_id_to_remove:
        del dic_idToRef[read_id]
    return(dic_idToRef, several_id_to_remove)






''' get the references is a list'''

def get_ref_list(blastnFileList):
    ref_set = set()
    for blastnFile in blastnFileList:
        with open(blastnFile, 'r') as f:
            for line in f:
                if line[0:6] == "Query=":
                    reference = line.split()[1]
                    ref_set.add(reference)
    ref_list = np.array(list(ref_set))
    return(ref_list)





''' get the sequence that corresponds to a reference identity'''
def get_seq_ref_dic(reference_fasta):
    dic_seq_ref = {}
    with open(reference_fasta, 'r') as f:
        lines = f.readlines()
        assert len(lines)%2 == 0
        for i in range(0, len(lines), 2):
            curr_ref = lines[i].strip('>').strip('\n').strip()
            if curr_ref not in dic_seq_ref:
                seq = lines[i+1].strip('\n').strip()
                dic_seq_ref[curr_ref] = seq
    return(dic_seq_ref)





''' get the name of the file without path and witohut extension '''
def get_file_name(file_path):
  file_path = file_path.strip()
  splited_fname = file_path.split('/')
  simple_file_name = splited_fname[-1]
  return simple_file_name
















''' ** creates training and testing sets with fast5 files**
    ( use filtered fast5 files )
    - Shuffle and split the references
    - Create a dic that links references for each read_id 
    - Use it to split the read_ids in lists
    - read fast5 files and write them in appropriate folders (4k reads max per file according to analysis)
   
    
    (*) fast5_global_filtered_folder : contains the fast5 flowcell folders
    (*) blastnFileList      : contains the paths to the blast files
    (*) outputFolder        : folder path where the new fast5s will be saved '''

# NB should be apllied a second time to split train in train + validation
def split(fast5_global_filtered_folder, blastnFileList, outputFolder):
    
    ##### shuffle and split the references
    np.random.seed(42)
    # get the list of references from blastn file
    ref_list = get_ref_list(blastnFileList)
    # shuffle
    np.random.shuffle(ref_list)
    # split 75/25
    ref_len = len(ref_list)
    ref_train_list = ref_list[:int(0.75*ref_len)]
    ref_test_list = ref_list[int(0.75*ref_len):]


    ##### Create dic that links references for each read_id, unique ref to each read_id
    dic_idToRef, several_id_to_remove = get_reads_ref(blastnFileList)


    ##### Split the read_ids in lists
    read_train_list = []
    read_test_list = []
    for read_id in dic_idToRef:
        ref = dic_idToRef[read_id]
        if(ref in ref_train_list):
            read_train_list.append(read_id)
        elif(ref in ref_test_list):
            read_test_list.append(read_id)
        else:
            print("Bizarre (SPLIT), ref does not appear in the train or test lists")


    # paths (TO CHANGE)
    train_path = outputFolder + "/train/"
    test_path = outputFolder + "/test/"
    if not os.path.exists(train_path):
        os.makedirs(train_path)
    else:
        print("train folder already exist at this location")
    if not os.path.exists(test_path):
        os.makedirs(test_path)
    else:
        print("test folder already exist at this location")

    
    folders = glob.glob(fast5_global_filtered_folder + "/*")
    for folder in tqdm(folders):
        if( not os.path.isdir(folder)):
            continue
        filenames = glob.glob(folder + "/*.fast5")
        for fast5File in filenames:
            # try to open fast5 file
            if not fast5File.endswith('fast5'):
                print("Not a fast5 file")
                return
            try:
                fast5_data = h5py.File(fast5File, 'r')
            except IOError:
                print('Error opening file')
                return
            
            #### write
            read = fast5_data['Raw/Reads']
            read_id = read[list(read.keys())[0]].attrs['read_id'].decode('UTF-8')

            if(read_id in read_train_list):
                new_path = train_path + get_file_name(fast5File)
                copyfile(fast5File, new_path)
                

            elif(read_id in read_test_list):
                new_path = test_path + get_file_name(fast5File)
                copyfile(fast5File, new_path)

            elif read_id in several_id_to_remove:
                pass

            else:
                print('read_id does not appear in the test or train list')
            
            fast5_data.close()





if __name__ == "__main__":

    ##### TO CHANGE #####

    # args
    base_path = "/media/victor/USB/MSc_basecall/Data/3XR6/"
    blastn_base = base_path + "ref_3xr6/blastn/"

    fast5_global_filtered_folder = sys.argv[1]
    outputFolder = sys.argv[2]

    blastnFileList = [blastn_base + "fc1s_42k.out",
                  blastn_base + "fc2s_42k.out",
                  blastn_base + "fc3s_42k.out" ]

    # function
    split(fast5_global_filtered_folder, blastnFileList, outputFolder)



            