
import seaborn as sns
import numpy as np
import sys
import matplotlib.pyplot as plt
import pydotplus
# from graphviz import Source
import pydot


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

elif sys.argv[2] == 'Tree':
	with open(sys.argv[1], 'r') as f:
		lines = f.readlines()
	lines = [line.rstrip().split(' ') for line in lines]
	tree_data = np.asarray(lines).astype(float)

	filename = sys.argv[1][sys.argv[1].index('/') + 1:]
	first_occ = filename.index('_')
	graph_name = filename[0:first_occ + 1]  + filename[first_occ + 1: first_occ + filename[first_occ + 1:].index('_') + 1]

	with open(graph_name + '.dot','w') as out:
	    for line in ('digraph G {','size="16,16";','splines=true;'):
	        out.write('{}\n'.format(line))  
	    for i in range(np.shape(tree_data)[0]):
	        out.write('{} -> {} [ label="{}" ];\n'.format(tree_data[i, 7], tree_data[i, 8], tree_data[i, 9]))
	    out.write('}\n')

	(graph, )= pydot.graph_from_dot_file(graph_name + '.dot')
	graph.write_png(sys.argv[3] + '/' + graph_name + '_tree.png')

elif sys.argv[2] == 'CONETT':
	with open(sys.argv[1], 'r') as f:
		lines = f.readlines()
	lines = [line.rstrip().split(' ') for line in lines[1:]]
	tree_data = np.asarray(lines)

	graph_name = 'CONETT_graph'

	with open(graph_name + '.dot','w') as out:
	    for line in ('digraph G {','size="16,16";','splines=true;'):
	        out.write('{}\n'.format(line))  
	    for i in range(np.shape(tree_data)[0]):
	    	out.write('{} -> {};\n'.format(tree_data[i, 1], tree_data[i, 3]))
	        # out.write('{} -> {} [ label="{}" ];\n'.format(tree_data[i, 1], tree_data[i, 3], tree_data[i, 0]))
	    out.write('}\n')

	(graph, )= pydot.graph_from_dot_file(graph_name + '.dot')
	graph.write_png(graph_name + '_tree.png')


elif sys.argv[2] == 'fr_removed':
	with open(sys.argv[1], 'r') as f:
		lines = f.readlines()
	fr_removed = [float(line.rstrip().split(' ')[0]) for line in lines]

	# print(fr_removed)
	edges = list(range(1, len(fr_removed) + 1))
	fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (15, 15))
	ax.bar(edges, fr_removed, width=1, ec="k", align="edge")
	plt.xlabel('Patients', fontsize = 14)
	plt.ylabel('Fraction of alterations removed', fontsize = 14)
	fig.savefig( 'plots/{}.png'.format('fr_removed'))   
	plt.close(fig)

else:
	x = [1, 2, 3, 5, 3, 6]
	frq, edges = np.histogram(x, bins = [1,2,3,5])


	
