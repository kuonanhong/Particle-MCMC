import sys
import pandas as pd
import numpy as np

# read in data
msgs = pd.read_csv(sys.argv[1], header=None, dtype='str', error_bad_lines=False)

rejected = 0
accepted = 0
nans = 0
nlls = []
ARs = []
for index, row in msgs.iterrows():
    if 'rejected' in str(row):
        rejected += 1
    elif 'accepted' in str(row):
        accepted += 1
    elif 'newLoglike' in str(row):
        if 'nan' in str(row) or 'NaN' in str(row):
            nans += 1
        else:
            my_str = str(row).split(':')[1]
            my_str = my_str.split(' ')[1]
            my_str = my_str.split('\n')[0]
            nlls.append(float(my_str))
#    elif 'AR' in str(row):
#        print( str(row).split(':')[1])

#    # temporary
#    if 'AR' in str(row):
#        ARs.append(float(str(row).split(":")[1].split("\n")[0]))
    

    

print("the acceptance rate was: ", accepted/(accepted+rejected))
if nans > 0:
    print('there were ', nans, ' nans')
#print("the variance of the log-likelihoods was: ", np.var(nlls))    
