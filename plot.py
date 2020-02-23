
import seaborn as sns
import numpy as np
import sys
import matplotlib.pyplot as plt


if sys.argv[2] == 'Coverage':
	with open(sys.argv[1], 'r') as f:
		lines = f.readlines()
	lines = [line.rstrip().split('\t') for line in lines]
	records = np.asarray(lines[1:])

	ALT_REF = records[:,4:6].astype(int)
	counts = np.sum(ALT_REF, axis = 1)

	sns.set()
	fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (15, 15))
	sns.distplot(counts, ax = ax, bins = 1000)
	plt.xlim(0, 500)
	fig.savefig( '{}.png'.format('coverage')) 
	plt.close(fig)

elif sys.argv[2] == 'VAF':
	with open(sys.argv[1], 'r') as f:
		lines = f.readlines()
	lines = [line.rstrip().split('\t') for line in lines]
	records = np.asarray(lines[1:])

	ALT_REF = records[:,4:6].astype(int)
	counts = np.sum(ALT_REF, axis = 1)

	ALT = records[:,4].astype(int)

	VAF = (2 * ALT) / counts


	b = list(range(101)) + [1000]
	b = [item/100 for item in b]

	frq = np.zeros((101))
	for i in range(len(VAF)):
		for j in range(len(b)-1):
			if VAF[i] > b[j] and VAF[i] <= b[j + 1]:
				frq[j] += 1
			if VAF[i] == 0:
				frq[0] += 1


	# frq, edges = np.histogram(VAF, bins = b)
	edges = list(b[0:-1])
	edges.append(1.01)


	sns.set()
	fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (15, 15))
	sns.distplot(VAF, ax = ax, bins = edges)
	plt.xlim(0, 1.02)
	fig.savefig( 'plots/{}.png'.format('VAF_dist')) 
	plt.close(fig)


	# fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (15, 15))
	# weights =np.ones_like(VAF)/len(VAF)
	# n, bins, patches = ax.hist(VAF, bins = edges)
	# fig.savefig( '{}.png'.format('VAF_hist'))   
	# plt.close(fig)


	fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (15, 15))
	ax.bar(edges[:-1], frq, width=np.diff(edges), ec="k", align="edge")
	fig.savefig( 'plots/{}.png'.format('VAF_count'))   
	plt.close(fig)

else:
	x = [1, 2, 3, 5, 3, 6]
	frq, edges = np.histogram(x, bins = [1,2,3,5])


	
