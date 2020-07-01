


import pandas as pd
import numpy as np
from collections import Counter

store = pd.HDFStore('GDCw_baa_samples_mutations_cc.h5')
samples_mutations_load = store['GDCw_baa_samples_mutations_cc']


print(np.sum(samples_mutations_load.values))

print(np.sum(samples_mutations_load.values, axis = 1))

print(np.sum(np.sum(samples_mutations_load.values, axis = 0) == 3))

print(len(list(samples_mutations_load.index)))

print(samples_mutations_load.columns[0:10])
x = list(samples_mutations_load.columns)

print('nm:', len(samples_mutations_load.columns))
print('ns:', len(list(samples_mutations_load.index)))

print((np.sum(samples_mutations_load.values == 1))/(np.shape(samples_mutations_load.values)[0] * np.shape(samples_mutations_load.values)[1]))

print(samples_mutations_load.values[0:5, 0:5])

print(list(set([item[0:item.find(':')] for item in x])))

print('------------------------------------------------')

# store = pd.HDFStore('samples_mutations.h5')
# samples_mutations_load = store['samples_mutations']


# print(np.sum(samples_mutations_load.values))

# print(np.sum(samples_mutations_load.values, axis = 1))

# print(np.sum(np.sum(samples_mutations_load.values, axis = 0) == 3))

# print(len(list(samples_mutations_load.index)))

# print(samples_mutations_load.columns[0:10])

# print('nm:',len(samples_mutations_load.columns))

# print((np.sum(samples_mutations_load.values == 1))/(np.shape(samples_mutations_load.values)[0] * np.shape(samples_mutations_load.values)[1]))

# print(samples_mutations_load.values[0:5, 0:5])

# print(list(set([item[0:item.find(':')] for item in list(samples_mutations_load.columns)])))

# # com = []
# # for item in list(samples_mutations_load.columns):
# # 	if item in x:
# # 		com.append(item)


# print('com:', len(set(x) & set(list(samples_mutations_load.columns))))

# # print(len(com))
# # with open('GDC_samples.txt', 'r') as f:
# # 	lines = f.readlines()


# samples = [s.rstrip().split('	') for s in lines[1:]]

# samples = np.asarray(samples)

# print(Counter(list(samples[:, 3])))

# print(Counter(list(samples[:, 4])))
