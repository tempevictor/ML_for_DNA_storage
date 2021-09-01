import pandas as pd
import sys


#sketch file with the results from alignment 
#REQUIRED
read_alignment_file = sys.argv[1] 

# ouputfile 
#REQUIRED
output_file = sys.argv[2] 

# time taken to basecall the folder
#OPTIONAL
try:
    caller_time_file = sys.argv[3] 
except:
    caller_time_file = ''


# dictionnary that stores the results
metrics = {}
# alignment database
df = pd.read_csv(read_alignment_file, header=0, sep='\t')






###### Alignment results
identities = df['Identity']

# the proportion of reads that were not properly aligned
metrics['no_match_ratio'] = len(identities[identities==0]) / len(identities)
# mean identity rate overall
metrics['global_mean_identity'] = identities.mean()
# mean identity rate of matched or aligned sequences
metrics['match_mean_identity'] = identities[identities>0].mean()







###### Length statistics

# extract lengths from the entire and aligned(match) database
global_lengths = df['Length']
match_lengths = df[df['Identity']>0]['Length']

#global
metrics['global_mean_read_length'] = global_lengths.mean()
metrics['global_median_read_length'] = global_lengths.median()
metrics['global_var_read_length'] = global_lengths.var()
metrics['global_mean_relative_lenght'] = df['Relative length'].mean()

#aligned
metrics['match_mean_read_length'] = match_lengths.mean()
metrics['match_median_read_length'] = match_lengths.median()
metrics['match_var_read_length'] = match_lengths.var()
metrics['match_mean_relative_lenght'] = df[df['Identity']>0]['Relative length'].mean()







###### Speeds

# decoding speed that does not consider pre processing. Bases per second 

if caller_time_file != '': # if such a file exists or is given as argument
    with open(caller_time_file, 'r') as inputfile:
        line = inputfile.readline()
        line_split = line.split(':')
        # from string to float in the list
        line_split = list( map(float, line_split) )
        # time in seconds
        time_bc = line_split[0]*3600 + line_split[1]*60 + line_split[2]
    total_bases_number = df['Length'].sum()
    print("TOTAL NB BASES: ", total_bases_number)
    metrics['network_speed'] = total_bases_number / time_bc

    # total basecaller speed, considering preprocessing
    caller_time_path = '/'.join( caller_time_file.split('/')[: -1] )
    total_time_file = caller_time_path + '/time.txt'
    with open(total_time_file, 'r') as inputfile:
        line = inputfile.readline()
        total_time = line.split()[-1].strip('s')
    metrics['total_speed'] = total_bases_number / float( total_time )





###### write in the output file

with open(output_file, 'w') as f:
    for key in metrics:
        f.write(key + ": " + str(metrics[key]) + '\n\n')
