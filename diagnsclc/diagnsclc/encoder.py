from itertools import product
import pandas as pd
import numpy as np
from sklearn.preprocessing import OneHotEncoder

def K_mer(data):
	k_mer = []
	for p in product(['A', 'T', 'G', 'C', 'D'],repeat=3):
		k_mer.append(''.join(p))

	k_mer = pd.DataFrame(k_mer)

	ohencoder = OneHotEncoder()
	ohencoder.fit(k_mer)
	k_mer_oh = ohencoder.transform(k_mer).toarray()
	k_mer_dict = dict()
	for i in range(125):
		k_mer_dict[k_mer[0][i]] = k_mer_oh[i]

	result = []
	for i in range(len(data)-3):
		result.append(np.array(k_mer_dict[str(data[i:i+3])]))
	# result = np.array(result)
	return result

def onehot(data):
    mapping = {"A":[1., 0., 0., 0.], "C": [0., 1., 0., 0.], "G": [0., 0., 1., 0.], "T":[0., 0., 0., 1.]}
    result = []
    for i in data:
        result.append(mapping[i] if i in mapping.keys() else [0., 0., 0., 0.])
    # result = np.array(result)
    return result

