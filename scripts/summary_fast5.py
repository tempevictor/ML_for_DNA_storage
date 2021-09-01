from matplotlib.pyplot import get
import numpy as np
import h5py
import os
import sys
import glob
import tkinter as tk
from tqdm import tqdm
from tkinter import filedialog

root = tk.Tk()
root.withdraw()


''' get the name of the file without path and witohut extension '''
def get_just_file_name(file_path):
  file_path = file_path.strip()
  splited_fname = file_path.split('/')
  simple_file_name = splited_fname[-1]
  simple_file_name = os.path.splitext(simple_file_name)[0]
  return simple_file_name


''' get a dictionarry that link a read id to the pass of fail filter status '''
def get_dic_id_pass_fail(sequencing_summary_text):
    dic_pass_fail = {}
    with open(sequencing_summary_text, 'r') as summary:
        first_line = True
        for line in summary:
            if first_line:
                first_line = False
                continue
            line = line.split()
            read_id = line[1]
            filter_status = line[9]
            if read_id not in dic_pass_fail:
                dic_pass_fail[read_id] = filter_status
            else:
                print("read_id already in the dic")
    return dic_pass_fail




''' separate PASS and FAIL fitlered reads from information given in a text file '''
def seprate_pass_fail(sequencing_summary_text, input_fast5_folder, output_folder):
    #get the dictionary that links read to pass or fail folder
    dic_pass_fail = get_dic_id_pass_fail(sequencing_summary_text)
    # wirte index, change file every 4000 write operations
    pass_write_index = 0
    fail_write_index = 0
    # make sub directories
    os.mkdir(output_folder + "/pass")
    os.mkdir(output_folder + "/fail")
    # initiate fast5 files
    current_pass_file = h5py.File(output_folder + "/pass/pass_0.fast5", 'a')
    current_fail_file = h5py.File(output_folder + "/fail/fail_0.fast5", 'a')
    filenames = glob.glob(input_fast5_folder + "/*.fast5")
    for fast5File in tqdm(filenames):
        # try to open fast5 file
        if not fast5File.endswith('fast5'):
            print("Not a fast5 file")
            return
        try:
            fast5_data = h5py.File(fast5File, 'r')
        except IOError:
            print('Error opening file')
            return
        
        read_list = list(fast5_data.items())
        for i in range(len(read_list)):
            readstep = (read_list[i])
            read = readstep[0]
            raw_id = fast5_data[read+'/Raw'].attrs['read_id'].decode('UTF-8')
            filter_status = dic_pass_fail[raw_id]
            if(filter_status == 'TRUE'):
                fast5_data.copy(fast5_data[read], current_pass_file)
                pass_write_index += 1
                if(pass_write_index % 4000 == 0):
                    current_pass_file.close()
                    current_pass_file = h5py.File(output_folder+"/pass/pass_"+str(pass_write_index//4000)+".fast5", 'a')
            elif(filter_status == 'FALSE'):
                fast5_data.copy(fast5_data[read], current_fail_file)
                fail_write_index += 1
                if(fail_write_index % 4000 == 0):
                    current_fail_file.close()
                    current_fail_file = h5py.File(output_folder+"/fail/fail_"+str(fail_write_index//4000)+".fast5", 'a')
            else:
                print("Unknown status:" + filter_status)
                return
    print("[Separate files into pass/fail done]")



if __name__ == "__main__":
    sequencing_summary_text = sys.argv[1]
    input_fast5_folder = sys.argv[2]
    output_folder = sys.argv[3]
    print("****** Begin separation ******")
    seprate_pass_fail(sequencing_summary_text, input_fast5_folder, output_folder)



            
