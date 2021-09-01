from statsmodels import robust
#import matplotlib.pyplot as plt
import numpy as np
import copy
import h5py
import sys
import os


def extract_raw(input_file):
    if not input_file.endswith('fast5'):
        return None
    try:
        fast5_data = h5py.File(input_file, 'r')
    except IOError:
        print('Error opening file. Likely a corrupted file.')
        return None

    #raw_attr = fast5_data['Raw/Reads/']
    #read_name = list(raw_attr.keys())[0]
    #raw_signal = raw_attr[read_name + '/Signal'][()]

    read_list = list(fast5_data.items())
    signaldata = []
    for i in range(len(read_list)):
        readstep = (read_list[i])
        read = readstep[0]
        signal_path = read + '/Raw/Signal'
        signaldata.append(fast5_data.get(signal_path))
    raw_signal = np.array(signaldata, dtype=np.float32)


    return np.array(raw_signal, dtype=np.float32)


def quantilef(madarr, nchunk, perc, newp):
    perc_lst = [perc]
    # print(perc_lst)#lx
    for i in range(newp):
        assert perc_lst[i] >= 0.0 and perc_lst[i] <= 1.0
    if madarr == []:
        for i in range(newp):
            perc_lst[i] = np.nan
        return perc_lst[0]
    space = copy.deepcopy(madarr)
    space.sort()
    # print('quantilef: space=',space)
    for i in range(newp):
        idx = int(perc_lst[i]*(nchunk-1))
        remf = perc_lst[i]*(nchunk-1)-idx
        # print(idx,remf,perc_lst[i],(nchunk-1),perc_lst[i]*(nchunk-1))
        if (idx < nchunk-1):
            perc_lst[i] = (1.0-remf)*space[idx]+remf*space[idx+1]
        else:
            perc_lst[i] = space[idx]
        # print((1.0-remf)*space[idx]+remf*space[idx+1])
        # print(space[idx])
    return perc_lst[0]


def trim_raw_by_mad(raw, chunk_size, perc):
    assert chunk_size > 1
    assert perc >= 0.0 and perc <= 1.0
    start = 0
    end = len(raw)
    nsample = end-start
    nchunk = nsample // chunk_size
    # print(nsample, chunk_size, nchunk) #lx
    end = nchunk*chunk_size

    madarr = list()
    for i in range(nchunk):
        madarr.append(robust.mad(
            raw[start+i*chunk_size:start+(i+1)*chunk_size]))
    # print('trim_by_mad:', 'nchunk=', nchunk, 'madarr=',madarr)#nchunk=851
    # print('before quantilef:perc=', perc)  #0.9
    perc = quantilef(madarr, nchunk, perc, 1)
    # print('after quantilet: perc=', perc) #97.
    thresh = perc
    for i in range(nchunk):
        if madarr[i] > thresh:
            break
        start += chunk_size
    for i in range(nchunk, 0, -1):
        if madarr[i-1] > thresh:
            break
        end -= chunk_size
    try:
        assert end > start
    except:
        return -1, -1
    return start, end


def trim_and_segment_raw(raw, trim_start, trim_end, varseg_chunk, varseg_thresh):
    start, end = trim_raw_by_mad(raw, varseg_chunk, varseg_thresh)
    # print('trimed_start_end:', start,end)
    if start == -1 and end == -1:
        return []
    n = len(raw)
    if n-start > trim_start:
        t_start = start+trim_start
    else:
        t_start = n
    # print('n, n-start, trim_start,t_start', n, n-start, trim_start,t_start)
    if end > trim_end:
        t_end = end-trim_end
    else:
        t_end = 0
    # print('trim_end, t_end', t_end)
    if t_start >= t_end:
        return []
    return raw[start:end]


if __name__ == "__main__":
    fast5_dir = sys.argv[1]
    for fname in os.listdir(fast5_dir):
        # print(fname)
        raw = extract_raw(fast5_dir+'/'+fname)
        if raw is None:
            continue
        trimed_raw = trim_and_segment_raw(raw, 200, 200, 200, 0.9)
        # print('raw_len:', len(raw), 'trimed_raw_len:', len(trimed_raw))

        import matplotlib.pyplot as plt 

        #fig = plt.figure()
        #ax1 = fig.add_subplot(2,1,1)
        #ax1.plot(range(len(raw)),raw)

        #ax2 = fig.add_subplot(2,1,2)
        #ax2.plot(range(len(trimed_raw)),trimed_raw)
        #plt.tight_layout()

        # Save the full figure...
        #fig.savefig('raw_and_trim' + str(fname) +'.png')
