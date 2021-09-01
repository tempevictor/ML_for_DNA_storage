import numpy as np
import h5py
import os
import argparse
import glob
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()

#def combine_files(dir_e):

## concatenate fastq files in a single fasta file
import sys
import os

''' get the name of the file without path and without extension '''
def get_just_file_name(file_path):
  file_path = file_path.strip()
  splited_fname = file_path.split('/')
  simple_file_name = splited_fname[-1]
  simple_file_name = os.path.splitext(simple_file_name)[0]
  return simple_file_name



'''  Combine fastq/fasta files and turn into one single fasta file '''
def combine_files_to_fasta(folder_path, extension_name, output_path_name = "reference.fasta" ):
    assert extension_name in ["fasta",  "fastq"]
    filenames = glob.glob(folder_path + "/*." + extension_name)
    #filenames = sorted(filenames, reverse=False)
    try:
        with open(output_path_name, 'w') as outfile:
            for f in filenames:
                # if fastq file keep first to lines with arrangements
                if(f[-5:] == "fastq"):
                    transcript = ''
                    with open(f, 'r') as infile:
                        for line_idx, line in enumerate(infile):
                            if(line_idx%4 == 0):
                                line = line.split()[0].strip('@')   
                                transcript += '>' + line + "\n"
                            elif (line_idx%4 == 1):
                                transcript += line           
                        outfile.write(transcript)
                #elif fasta file keep all
                elif(f[-5:] == "fasta"):
                    transcript = ''
                    with open(f, 'r') as infile:
                        for line_idx, line in enumerate(infile):
                            transcript += line         
                        outfile.write(transcript)
                else:
                    print("Wrong format")
        print("Combine fastq/fasta to fasta done")
    except:
        print("Cannot open the files")
        return folder_path


# fastq_dir_path = "/media/victor/USB/MSc_basecall/Data/3XR6/guppy/guppy_FILTERED_all/fail_fasta"
# output_path_name = "/media/victor/USB/MSc_basecall/Data/3XR6/guppy/guppy_FILTERED_all/Final/fail/out.fasta"
# combine_files_to_fasta(fastq_dir_path, "fasta", output_path_name)



'''  Combine fastq files and turn into one single fastq file '''
def combine_files_to_fastq(folder_path, output_path_name):
    filenames = glob.glob(folder_path + "/*.fastq")
    #filenames = sorted(filenames, reverse=False)
    try:
        with open(output_path_name, 'w') as outfile:
            for f in filenames:
                with open(f, 'r') as infile:
                    for line in infile:
                        outfile.write(line)
        print("Combine fastq files done")
    except:
        print("Cannot open the files")
        return folder_path




                

''' store the order the read ids are stored in a fasta file '''
def get_id_order_fasta(fasta_file_name):
    assert fasta_file_name[-5:] == "fasta"
    try:
        id_list = []
        with open(fasta_file_name, 'r') as infile:
            for line_idx, line in enumerate(infile):
                # read one over two lines
                if(line_idx%2 == 0):
                    name_id = line
                    id_list.append(name_id)
        print("id_list done")
        return id_list
    except:
        print("An error occured trying to read the file")
        return None



''' merge all the content under one id which is the file name, all in the order given by an id_list '''
def sorted_fasta(fasta_file_name_bc, output_name, id_list_ref):
    assert fasta_file_name_bc[-5:] == "fasta"
    assert output_name[-5:] == "fasta"
    try:
        with open(output_name, 'w') as outfile:
            with open(fasta_file_name_bc, "r") as infile:
                lines = infile.readlines()
                assert len(lines)%2 == 0
                transcript_list = id_list_ref.copy()
                for index, id_ref in enumerate(id_list_ref):
                    found_read_id = False
                    N = len(lines)//2
                    for i in range(N):
                        read_id = lines[2*i]
                        found_read_id = False
                        if(read_id == id_ref):
                            content = lines[2*i+1]
                            transcript_list[index] = id_ref + content
                            # save the fact we found a corresponding id in the bc fasta
                            found_read_id = True
                            # remove id and content to speed up
                            del lines[2*i+1]
                            del lines[2*i]
                            # break the loop
                            break
                    # if id not found, add an empty line
                    if(not found_read_id):
                        transcript_list[index] = id_ref + "A\n" #add A not to have an empty array
            transcript = "".join(transcript_list)
            outfile.write(transcript)
        print("sorted fasta done")
    except:
        print("An error occured trying to merge sort the file")



''' merge all the content under one id which is the file name, all in the order
     given by an id_list if given'''
def single_elem_fasta(fasta_file_name_bc, output_name, id_list_ref=None):
    assert fasta_file_name_bc[-5:] == "fasta"
    assert output_name[-5:] == "fasta"
    try:
        counter_empty = 0
        with open(output_name, 'w') as outfile:
            with open(fasta_file_name_bc, "r") as infile:
                lines = infile.readlines()
                assert len(lines)%2 == 0

                # no order, just concatenate in that case
                if(id_list_ref == None):
                    transcript_list = []
                    for i in range(len(lines)//2):
                        content = lines[2*i+1]
                        transcript_list.append(content.strip('\n'))
                # in case an order is wanted
                else:
                    transcript_list = id_list_ref.copy()
                    for index, id_ref in enumerate(id_list_ref):
                        found_read_id = False
                        N = len(lines)//2
                        for i in range(N):
                            read_id = lines[2*i]
                            found_read_id = False
                            if(read_id == id_ref):
                                content = lines[2*i+1]
                                transcript_list[index] = content.strip('\n')
                                # save the fact we found a corresponding id in the bc fasta
                                found_read_id = True
                                break
                        # if id not found, add an empty line
                        if(not found_read_id):
                            counter_empty += 1
                            transcript_list[index] = ""
            # write operation
            transcript = "".join(transcript_list)
            outfile.write(">dna_basecalled\n")
            outfile.write(transcript)
        print("final single element fasta done")
        print("Counter empty :" + str(counter_empty))
    except:
        print("An error occured trying to merge the file")





###### CHANGE FILE ######
# fastq_dir_path = "/media/victor/USB/MSc_basecall/Data/3XR6/guppy/guppy_filtered_long/fastq_full/"
# # input_file = filedialog.askopenfilename()
# # input_file = filedialog.askopenfilename()
# for f in ["1a-42k", "2b-42k", "3d-42k"]:
#     path_in = fastq_dir_path + f
#     path_out = fastq_dir_path + "../fastq/" + f + "/out.fastq"
#     combine_files_to_fastq(path_in, path_out)
# # id_list = get_id_order_fasta(fastq_dir_path + "result/reference_tmp.fasta")
# # sorted_fasta("/media/victor/USB/MSc_basecall/Data/shortExamples/fastas/a/out.fasta",
# #                   fastq_dir_path + "result/reference_bc_ins3.fasta",  id_list)
# #########################



''' sanity check '''
def check_same_ids(bc_fasta_file, reference_file):
    with open(bc_fasta_file, "r") as bc:
        with open(reference_file, "r") as ref:
            lines_bc = bc.readlines()
            lines_ref = ref.readlines()
            if(len(lines_bc) != len(lines_ref)):
                print("Not the same number of lines")
                return False
            if(len(lines_bc)%2 != 0):
                print("Not a pair number of lines")
                return False
            it = 0
            while(it<len(lines_bc)//2):
                if(lines_bc[2*it] != lines_ref[2*it]):
                    print("Line " + str(2*it) + " is different")
                    return False
                it += 1
            print("Files have the same ids")
            return True
                
##############
# check_same_ids("/media/victor/USB/MSc_basecall/Data/shortExamples/fastqs/result/a/a_result/",
#               "/media/victor/USB/MSc_basecall/Data/shortExamples/fastqs/result/a/reference_bc.fasta")
##############

# convert tsv file to fasta file
def tsvToFasta(reference_tsv, output_fasta):
    transcript = ''
    with open(reference_tsv, "r") as infile:
        with open(output_fasta, "w") as outfile:
            for line in infile:
                content = line.strip().strip("\n").split("\t")
                transcript += ">"+content[0]+"\n" + content[1]+"\n"
            outfile.write(transcript)

# input_file = filedialog.askopenfilename()
# tsvToFasta(input_file, '/media/victor/USB/MSc_basecall/Data/3XR6/ref_3xr6/reference_full.fasta')



# convert fastq file to fasta file
def fastq2fasta(fastq_file_name, output_fasta_path):
    assert fastq_file_name[-5:] == 'fastq'
    assert output_fasta_path[-5:] == "fasta"
    transcript = ""
    #read the fastq file content
    with open(fastq_file_name, 'r') as fastq:
        for line in fastq:
            if(line[0] == '@'):
                line = line.split()[0].strip('@').strip('\n')   
                transcript += '>' + line + '\n'
            elif (line[0] in ['A', 'C', 'T', 'G']):
                transcript += line           
    # write in fasta format
    with open(output_fasta_path, 'w+') as fasta:
        fasta.write(transcript)
    print("Success: conversion from fastq to fasta")
    
# fastq_file_name = "/media/victor/USB/MSc_basecall/Data/3XR6/fastas/test256_bonito/bontio_resquigg_256.fastq"
# output_fasta_path = "/media/victor/USB/MSc_basecall/Data/3XR6/fastas/test256_bonito/out.fasta"
# fastq2fasta(fastq_file_name, output_fasta_path)

# ####################
import glob
# base_path = "/media/victor/USB/MSc_basecall/Data/3XR6/guppy/guppy_FILTERED_all/"
# file_names = glob.glob(base_path + "pass_fastq/*")
# for f in file_names:
#     name = f.split("/")[-1][:-1] + "a"
#     output_fasta_path = base_path + "pass_fasta/" + name
#     fastq2fasta(f, output_fasta_path)

# for win_l in ["results_2048/", "results_4000/"]:
#     #for ty in ["long/", "short/"]:
#         #for folder in ["1a-42k/", "2b-42k/", "3d-42k/"]:       
#     fastq_file_name = base_path + win_l+ "*.fastq"
#     fastq_file_name = glob.glob(fastq_file_name)[0]
#     print(fastq_file_name)
#     output_fasta_path = base_path + win_l + "out.fasta"
#     fastq2fasta(fastq_file_name, output_fasta_path)

# ################


# convert fasta file to fastq file
def fasta2fastq(fasta_file_name, output_fastq_path):
    assert fasta_file_name[-5:] == 'fasta'
    assert output_fastq_path[-5:] == "fastq"
    transcript = ""
    current_seq_length = 0
    #read the fasta file content
    with open(fasta_file_name, 'r') as fasta:
        for line_idx, line in enumerate(fasta):
            if(line_idx%2 == 0):
                line = line.split()[0].strip('>').strip('\n')
                transcript += '@' + line + '\n'
            else:
                transcript += line
                current_seq_length = len(line.strip().strip('\n'))
                transcript += '+\n'
                transcript += '?'*current_seq_length + '\n'
     # write in fastq format
    with open(output_fastq_path, 'w+') as fastq:
        fastq.write(transcript)
    print("Success: conversion from fastq to fasta")



