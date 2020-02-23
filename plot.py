
import seaborn as sns
import numpy as np
import sys
import matplotlib.pyplot as plt


with open(sys.argv[1], 'r') as f:
	lines = f.readlines()
lines = [line.rstrip().split('\t') for line in lines]
records = np.asarray(lines[1:])

ALT_REF = records[:,4:6].astype(int)
counts = np.sum(ALT_REF, axis = 1)

sns.set()
fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (15, 15))
sns.distplot(counts, ax = ax, bins = 100)
plt.xlim(0, 500)
fig.savefig( '{}.png'.format('coverage')) 
plt.close(fig)