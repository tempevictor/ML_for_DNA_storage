import numpy as np
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



############### UTILS ###############

WIN_SIZE = 2048
LABEL_SIZE = 120

# Some constants for the dictionaries
#****************************************************
# @Author  : Neng Huang
# @github   : https://github.com/huangnengCSU/SACall-basecaller

SIG_PAD = -10
BLANK = 0
A = 1
T = 2
C = 3
G = 4
PAD = 5   # padding
BASE_DIC = {'A': A, 'T': T, 'C': C, 'G': G}
INT_TO_CHAR = {0: 'BLANK', 1: 'A', 2: 'T', 3: 'C', 4: 'G', 5: ' '}
#****************************************************


''' get the dictionary of the ids in the blastn file with low e-values '''

def get_list_id_to_filter(blastnFile):
    with open(blastnFile, 'r') as f:
        #dictionary of blast ids
        blast_ids = {}
        for line in f:
            if line[0] == ">":
                content = line.strip("\n").strip().strip(">")
                read_id = content.split()[0]
                if read_id not in blast_ids:
                    blast_ids[read_id] = 1
                else:
                    blast_ids[read_id] += 1
    return blast_ids







''' get the name of the file without path and witohut extension '''

def get_just_file_name(file_path):
  file_path = file_path.strip()
  splited_fname = file_path.split('/')
  simple_file_name = splited_fname[-1]
  simple_file_name = os.path.splitext(simple_file_name)[0]
  return simple_file_name







############### CREATE DATASET ###############

''' Keep the signals with an e_value > 1e-14      
        (*) fast5Folder : folder containing fast5 files to filter for a flowcell
        (*) blastnFile  : blastn result of the flowcell
        (*) outputFolder: folder wherre the filtered fast5s will be saved'''

def filter_quality_signal(fast5Folder, blastnFile, outputFolder):

    # get the dic of filtered ids
    blast_ids = get_list_id_to_filter(blastnFile)

    # for each of the fast5 file   
    filenames = glob.glob(fast5Folder + "/*.fast5")
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

        # new fast5 file to write to
        new_f5_path =  outputFolder + "/" + get_just_file_name(fast5File) + "_FILT.fast5"
        new_fast5 = h5py.File(new_f5_path, 'a')
        # copy
        read_list = list(fast5_data.items())
        for i in range(len(read_list)):
            readstep = (read_list[i])
            read = readstep[0]
            raw_id = fast5_data[read+'/Raw'].attrs['read_id'].decode('UTF-8')
            # copy if in filtered ids
            if raw_id in blast_ids:
                fast5_data.copy(fast5_data[read], new_fast5)
        fast5_data.close()
        new_fast5.close()





########## RUN FILTER ##########
# fast5Folder = filedialog.askdirectory()
# blastnFile = filedialog.askopenfilename()
# outputFolder = filedialog.askdirectory()

# filter_quality_signal(fast5Folder, blastnFile, outputFolder)






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





''' translate a string sequence to the corresponding int sequence according to a dic'''
def get_int_ref_seq(ref_sequence, BASE_DIC):
    ref_seq_int = []
    for base in ref_sequence:
        int_base = BASE_DIC.get(base)
        ref_seq_int.append(int_base)
    ref_seq_int = np.array(ref_seq_int)
    return(ref_seq_int)




""" 
Do the pre processing for the singals and split the label
sequence accordingly into segments. If the length is larger 
than the window size (2048 ie WIN_SIZE), split it.
The signal must be half has as long as WIN_SIZE to be accepted"""

def process_signal(signal, int_ref_seq):
    signal_splits = []
    label_splits = []
    #normalise the entire signal
    signal = (signal - np.median(signal)) / np.float32(robust.mad(signal))
    signal_length = len(signal)
    label_length = len(int_ref_seq)
    # accept small signal if its length is between 2048/2 and 2048
    if(signal_length < WIN_SIZE ):
        if(signal_length > WIN_SIZE//2):
            signal = np.pad( signal, (0, WIN_SIZE-signal_length), mode='constant', constant_values=SIG_PAD )
            label = np.pad( int_ref_seq, (0, LABEL_SIZE-label_length), mode='constant', constant_values=PAD )
            signal_splits.append(signal)
            label_splits.append(label)
            return signal_splits, label_splits
        else:
            print('signal not accepted')
            return None
    # else
    label_span = int( WIN_SIZE / signal_length * LABEL_SIZE )
    for i in range(signal_length // WIN_SIZE ):
        signal_splits.append(signal[WIN_SIZE*i: WIN_SIZE*(i+1)])
        label = int_ref_seq[ label_span*i : label_span*(i+1) ]
        label = np.pad( label, (0, LABEL_SIZE-len(label)),
                        mode='constant', constant_values=PAD)
        label_splits.append(label)
    
    # if the last split is long enough
    if(signal_length%WIN_SIZE > WIN_SIZE//2):
        start_last_split = (signal_length // WIN_SIZE)*WIN_SIZE
        last_split = signal[start_last_split : ]
        last_split = np.pad( last_split, (0, WIN_SIZE-len(last_split)), mode='constant', constant_values=SIG_PAD)
        signal_splits.append(last_split)
        last_label_split = int_ref_seq[label_span*(signal_length // WIN_SIZE):]
        last_label_split = np.pad( last_label_split, (0, LABEL_SIZE-len(last_label_split)),
                                    mode='constant', constant_values=PAD)
        label_splits.append(last_label_split)

    return signal_splits, label_splits





''' ** creates training and testing sets **
    ( use filtered fast5 files )
    - Shuffle and split the references
    - Create a dic that links references for each read_id 
    - Use it to split the read_ids in lists
    - read fast5 files and write them in appropriate folders (4k reads max per file according to analysis)
    + write npy labels and read signals for training
    
    (*) fast5_global_filtered_folder : contains the fast5 flowcell folders
    (*) blastnFileList      : containes the paths to the blast files
    (*) reference_fasta     : fasta file with the reference '''

def create_dataset(fast5_global_filtered_folder, blastnFileList, outputFolder, reference_fasta):
    
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
    train_f5_folder_path = outputFolder + "/f5_train/train0/"
    train_npy_folder_path = outputFolder + "/npy_signal_train/"
    train_label_folder_path = outputFolder + "/npy_label_train/"
    test_f5_folder_path = outputFolder + "/f5_test/test0/"
    test_npy_folder_path = outputFolder + "/npy_signal_test/"
    test_label_folder_path = outputFolder + "/npy_label_test/"
    """  #to create parent directories (such as mkdir -p) use the following:
    paths = [train_f5_folder_path, train_npy_folder_path, train_label_folder_path,
             test_f5_folder_path, test_npy_folder_path, test_label_folder_path]
    for path in paths: 
        os.makedirs(path, exist_ok=True) """


    ##### Read all filtered fast5s adn write to fast5 files the result
    #change files to write to when 4000 threshold is reached
    def nex_step_write_train(current_train_file, signal_train_list, label_train_list, file_write_index_train):    
        current_train_file.close()
        current_train_file = h5py.File(train_f5_folder_path + "file_" +  str(file_write_index_train//4000)+".fast5", 'a')
        # save npy files and initialise new list
        signal_train_list = np.array(signal_train_list)
        assert signal_train_list.shape[1] == WIN_SIZE, "not the expected shape"
        np.save(train_npy_folder_path + '/' +  "file_" + str( (file_write_index_train//4000) - 1),
                signal_train_list)
        label_train_list = np.array(label_train_list)   
        assert label_train_list.shape[1] == LABEL_SIZE, "not the expected shape"
        np.save(train_label_folder_path + '/' +  "file_" + str( (file_write_index_train//4000) - 1),
                label_train_list)
        return(current_train_file)
    
    def nex_step_write_test(current_test_file, signal_test_list, label_test_list, file_write_index_test):
        current_test_file.close()
        current_test_file = h5py.File(test_f5_folder_path + "file_" +  str(file_write_index_test//4000)+".fast5", 'a')
        # save npy files and initialise new list
        signal_test_list = np.array(signal_test_list)
        assert signal_test_list.shape[1] == WIN_SIZE, "not the expected shape"
        np.save(test_npy_folder_path + '/' +  "file_" + str( (file_write_index_test//4000) - 1),
                signal_test_list)
        label_test_list = np.array(label_test_list)
        assert label_test_list.shape[1] == LABEL_SIZE, "not the expected shape"
        np.save(test_label_folder_path + '/' +  "file_" + str( (file_write_index_test//4000) - 1),
                label_test_list)
        return(current_test_file)


    # count to open new file every 4000 writes
    file_write_index_train = 0
    file_write_index_test = 0
    # open fast5 files to write to
    current_train_file = h5py.File(train_f5_folder_path + "file_0.fast5", 'a')
    current_test_file = h5py.File(test_f5_folder_path + "file_0.fast5", 'a')
    # lists for npy
    signal_train_list = []
    signal_test_list = []
    label_train_list = []
    label_test_list = []
    # sequence under a reference identity
    dic_seq_ref = get_seq_ref_dic(reference_fasta)
    
    folders = glob.glob(fast5_global_filtered_folder + "/*")
    for folder in tqdm(folders):
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
            read_list = list(fast5_data.items())
            for i in range(len(read_list)):
                readstep = (read_list[i])
                read = readstep[0]
                raw_id = fast5_data[read+'/Raw'].attrs['read_id'].decode('UTF-8')
                
                ## copy in train if raw_id in the right split
                if raw_id in read_train_list:
                    # write in fast5 file
                    fast5_data.copy(fast5_data[read], current_train_file)
                    # add npy files with label
                    ref = dic_idToRef[raw_id]
                    ref_sequence = dic_seq_ref[ref]
                    int_ref_seq = get_int_ref_seq(ref_sequence, BASE_DIC)
                    assert len(int_ref_seq) == LABEL_SIZE, "train label size is not correct"
                    # add npy files with normalised signal
                    signal = np.array(fast5_data.get(read + '/Raw/Signal'), dtype=np.float32)
                    signal_splits, label_splits = process_signal(signal, int_ref_seq)
                    if(signal_splits is not None):
                        for i in range(len(signal_splits)):
                            signal_train_list.append(signal_splits[i])
                            label_train_list.append(label_splits[i])
                            file_write_index_train += 1
                            #reset the writers
                            if( file_write_index_train%4000 == 0 ):
                                current_train_file = nex_step_write_train( current_train_file,
                                                     signal_train_list, label_train_list,
                                                     file_write_index_train )
                                signal_train_list = []
                                label_train_list = []
                    elif(signal_splits is None):
                        print("signal_splits train is None")

                ## copy in test if raw_id in the right split
                elif raw_id in read_test_list:
                    # write in fast5 file
                    fast5_data.copy(fast5_data[read], current_test_file)
                    # add npy files with label
                    ref = dic_idToRef[raw_id]
                    ref_sequence = dic_seq_ref[ref]
                    int_ref_seq = get_int_ref_seq(ref_sequence, BASE_DIC)
                    assert len(int_ref_seq) == LABEL_SIZE, "test label size is not correct"
                    # add npy files with normalised signal
                    signal = np.array(fast5_data.get(read + '/Raw/Signal'), dtype=np.float32)
                    signal_splits, label_splits = process_signal(signal, int_ref_seq)
                    if(signal_splits is not None):
                        for i in range(len(signal_splits)):
                            signal_test_list.append(signal_splits[i])
                            label_test_list.append(label_splits[i])
                            file_write_index_test += 1
                            #reset the writers
                            if( file_write_index_test%4000 == 0):
                                current_test_file =  nex_step_write_test( current_test_file,
                                                     signal_test_list, label_test_list,
                                                      file_write_index_test )  
                                signal_test_list = []
                                label_test_list = []
                    elif(signal_splits is None):
                        print("signal_splits test is None")
                
                ## do nothing if raw_id needs to be removed
                elif raw_id in several_id_to_remove:
                    pass
                    #print("Raw_id in several_id_to_remove", end="\r")
                
                ## there is omething wrong
                else:
                    print("Bizarre (COPY), ref does not appear in the train or test lists nor in the list of ids that have several references")
            
            fast5_data.close()

    #final save
    current_train_file.close()
    current_test_file.close()
    # save if list not empty
    if(len(signal_train_list) > 0 ):
        np.save(train_npy_folder_path + '/' +  "file_" + str( (file_write_index_train//4000)),
                np.array(signal_train_list))
    if(len(label_train_list) > 0):
        np.save(train_label_folder_path + '/' +  "file_" + str( (file_write_index_train//4000)),
                np.array(label_train_list))
    if(len(signal_test_list) > 0):
        np.save(test_npy_folder_path + '/' +  "file_" + str( (file_write_index_test//4000)),
                np.array(signal_test_list))
    if(len(label_test_list) > 0):
        np.save(test_label_folder_path + '/' +  "file_" + str( (file_write_index_test//4000)),
                np.array(label_test_list))










###### RUN CREATE DATASET ######
base_path = "/media/victor/USB/MSc_basecall/Data/3XR6/"
blastn_base = base_path + "ref_3xr6/blastn/"


## for 120 filtered files
print("Beignning building long dataset with 120 bps...")
LABEL_SIZE = 120
fast5_global_filtered_folder = base_path + "fast5s/fast5s_filtered_long/"
blastnFileList = [blastn_base + "fc1s_42k.out",
                  blastn_base + "fc2s_42k.out",
                  blastn_base + "fc3s_42k.out" ]
outputFolder = base_path + "datasets/dataset_long/"
reference_fasta = base_path + "ref_3xr6/reference_full.fasta"
#
create_dataset(fast5_global_filtered_folder, blastnFileList, outputFolder, reference_fasta)
print("end 120 \n\n")


## for 25 filtered files
print("Beignning building short dataset with 25 bps...")
LABEL_SIZE = 25
fast5_global_filtered_folder = base_path + "fast5s/fast5s_filtered_short/"
blastnFileList = [blastn_base + "fc1s_42k_addonly.out",
                  blastn_base + "fc2s_42k_addonly.out",
                  blastn_base + "fc3s_42k_addonly.out" ]
outputFolder = base_path + "datasets/dataset_short/"
reference_fasta = base_path + "ref_3xr6/reference.fasta"
#
create_dataset(fast5_global_filtered_folder, blastnFileList, outputFolder, reference_fasta)
print("end 25")


