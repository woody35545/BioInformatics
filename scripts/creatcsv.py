import pandas as pd
import numpy as np
import os
# print(os.getcwd())

sequence = []
length = []
label = []
for i in range(619, 741):
    f = open(os.path.join(os.getcwd(), '..', 'exon', 'exon_SRR7785{}.txt'.format(i)), 'r')
    data = f.readline()
    data = data.upper()
    sequence.append(data)
    length.append(len(data))
    label.append(1)
    f.close()
    
f = open(os.path.join(os.getcwd(), '..', 'data', 'patient_noise.txt'), 'r')
data = f.readlines()
for d in data:
    sequence.append(d.rstrip('\n'))
    length.append(len(d))
    label.append(1)

f = open(os.path.join(os.getcwd(), '..', 'data', 'normal_data.txt'), 'r')
data = f.readlines()
for line in data:
    sequence.append(line.replace("\n",""))
    length.append(len(line))
    label.append(0)
    f.close()
    
raw_data = {'sequence': sequence,
            'len': length,
            'label': label,}

data = pd.DataFrame(raw_data)
data.to_csv('../data/data.csv', mode='w', index=False)

