# @Author  : Neng Huang
# @github   : https://github.com/huangnengCSU/SACall-basecaller

import os
from multiprocessing import Pool
from itertools import repeat
import argparse
import copy
import time
import h5py
import shutil
from statsmodels import robust
import numpy as np
from trim_raw import trim_and_segment_raw
import constants


def box(array, threshold):
    ori_l = len(array)
    Percentile = np.percentile(array, [0, 25, 50, 75, 100])
    # print(Percentile)
    IQR = Percentile[3] - Percentile[1]
    UpLimit = Percentile[3]+IQR*2
    DownLimit = Percentile[1]-IQR*2.5
    array = array[np.where((array < UpLimit) & (array > DownLimit))]
    new_l = len(array)
    if (ori_l-new_l) > threshold:
        return []
    return array


def raw_split(raw, raw_win_len):
    truly_raw_win_len = raw_win_len
    raw = (raw - np.median(raw)) / np.float(robust.mad(raw))
    raw_len = len(raw)
    wins = list()
    start = 0
    step = truly_raw_win_len
    while True:
        end = start + truly_raw_win_len
        if end > raw_len:
            break
        wins.append(raw[start:end])
        start += step
    last_w = raw[start:raw_len]
    if len(last_w) > 0:
        wins.append(np.pad(last_w, (0, raw_win_len - len(last_w)),
                           mode='constant', constant_values=constants.SIG_PAD))
    return wins


def call_preprocess(filename, argv):
    # print(filename)#lx
    input_file = argv.fast5 + '/' + filename
    if not input_file.endswith('fast5'):
        print("fast5 file not found")
        return filename
    try:
        fast5_data = h5py.File(input_file, 'r')
    except IOError:
        print('Error opening file. Likely a corrupted file.')
        return filename
    ###### modified
    try:
        read_list = list(fast5_data.items())
        raw_dic = dict()

        for i in range(len(read_list)):
            # get the read group and path to signal
            readstep = (read_list[i])
            read = readstep[0]
            signal_path = read + '/Raw/Signal'
            # get the read_id
            read_id = fast5_data[read+'/Raw'].attrs['read_id'].decode('UTF-8')
            # get the signal
            raw_signal = np.array(fast5_data.get(signal_path), dtype=np.float32)
            #add to dictionnary
            raw_dic[read_id]=raw_signal

    except:
        print('missing some component in fast5 file')
        return filename

    ###### modified
    for read_id in raw_dic:
        raw = raw_dic[read_id]
        trimed_raw = raw
        trimed_raw = trim_and_segment_raw(
            raw, argv.trim_start, argv.trim_end, argv.chunk_size, argv.trim_thresh)
        if len(trimed_raw) / len(raw) < 0.6:
            # fail
            #print(read_id + " fail: ratio < 0.6")
            continue
            #return filename
        # detect outlier value
        trimed_raw = box(trimed_raw, argv.outlier_number_threshold)
        # print(trimed_raw)
        if len(trimed_raw) == 0:    # box() return []
            #print(read_id + " timmed is empty")
            continue
            #return filename
        wins = raw_split(trimed_raw, argv.raw_len)
        # print(np.array(wins).shape)
        signal_batch = np.expand_dims(np.array(wins), 2)
        # print(signal_batch.shape)
        np.save(argv.records_dir + '/' + read_id + '.npy', signal_batch)
    return None


if __name__ == "__main__":
    start_t = time.time()
    parse = argparse.ArgumentParser()
    parse.add_argument('-fast5', help="directory to raw fast5", required=True)
    parse.add_argument(
        '-records_dir', help="directory to temporary records files", required=True)
    parse.add_argument('-chunk_size', type=int, default=256)
    parse.add_argument('-trim_thresh', type=float, default=0.9)
    parse.add_argument('-trim_start', type=int, default=150)
    parse.add_argument('-trim_end', type=int, default=150)
    parse.add_argument('-raw_len', type=int, required=True)
    parse.add_argument('-outlier_number_threshold', type=int, default=500)
    parse.add_argument('-threads', type=int, default=8)
    argv = parse.parse_args()

    if os.path.exists(argv.records_dir):
        shutil.rmtree(argv.records_dir)
        os.mkdir(argv.records_dir)
    else:
        os.mkdir(argv.records_dir)
    print('multiprocessing generate records file...')
    filenames = os.listdir(argv.fast5)

    pool = Pool(argv.threads)
    # print(list(zip(filenames, repeat(argv)))) #lx
    results = pool.starmap(call_preprocess, zip(filenames, repeat(argv)))
    end_t = time.time()
    with open(os.path.basename(argv.fast5)+'_fail.txt', 'w') as ffail:
        for v in results:
            if v:
                ffail.write(v)
    print("finish.\n Time cost: {t:.2f} seconds".format(t=(end_t - start_t)))
