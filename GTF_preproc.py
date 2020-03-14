#

from time import time

GTFDIR = '/home/haghirebrahimm2/refgenome/hg19refGene.gtf.sorted'

def readGTF(GTFDIR):
	with open(GTFDIR, 'r') as f:
		lines = f.readlines()

	GTF = []
	for line in lines[1:]:
		GTF_record = line.rstrip().split('\t')
		GTF.append([GTF_record[0][3:], GTF_record[3], GTF_record[4], GTF_record[8].split(';')[0]])
	GTF = [record[0:3] + [record[3].split(' ')[1][1:-1]] for record in GTF if record[3].split(' ')[0] == 'gene_id']

	GTF_dict = {}
	for record in GTF:
		if record[0] not in GTF_dict.keys():
			GTF_dict[record[0]] = {}
			GTF_dict[record[0]][record[3]] = [[int(record[1]), int(record[2])]]
		else:
			if record[3] not in GTF_dict[record[0]].keys():

				GTF_dict[record[0]][record[3]] = [[int(record[1]), int(record[2])]]
			else:
				GTF_dict[record[0]][record[3]] = GTF_dict[record[0]][record[3]] + [[int(record[1]), int(record[2])]]
	return GTF_dict


def gene(CHROM, POSITION, GTF_dict):
	gene_id = ''
	for k, v in GTF_dict[CHROM].items():
		for item in v:
			if (POSITION >= item[0]) and (POSITION <= item[1]):
				gene_id = k
				break
	return gene_id





# GTF_dict = readGTF(GTFDIR)



# CHROM = '1'
# POSITION = 1666175

# st = time()
# gene_id = gene(CHROM, POSITION, GTF_dict)
# print(gene_id)
# print(time() - st)