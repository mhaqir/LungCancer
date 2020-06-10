


import pandas as pd
import numpy as np


store = pd.HDFStore('GDC_samples_mutations.h5')
samples_mutations_load = store['samples_mutations']


print(np.sum(samples_mutations_load.values))

print(np.sum(samples_mutations_load.values, axis = 1))

print(np.sum(np.sum(samples_mutations_load.values, axis = 0) == 3))

print(len(list(samples_mutations_load.index)))

print(samples_mutations_load.columns[0:10])

print(len(samples_mutations_load.columns))

print((np.sum(samples_mutations_load.values == 1))/(np.shape(samples_mutations_load.values)[0] * np.shape(samples_mutations_load.values)[1]))

print(samples_mutations_load.values[0:5, 0:5])


print('------------------------------------------------')

store = pd.HDFStore('samples_mutations.h5')
samples_mutations_load = store['samples_mutations']


print(np.sum(samples_mutations_load.values))

print(np.sum(samples_mutations_load.values, axis = 1))

print(np.sum(np.sum(samples_mutations_load.values, axis = 0) == 3))

print(len(list(samples_mutations_load.index)))

print(samples_mutations_load.columns[0:10])

print(len(samples_mutations_load.columns))

print((np.sum(samples_mutations_load.values == 1))/(np.shape(samples_mutations_load.values)[0] * np.shape(samples_mutations_load.values)[1]))

print(samples_mutations_load.values[0:5, 0:5])