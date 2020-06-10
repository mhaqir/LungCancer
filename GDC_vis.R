#install.packages('treeio')
library('treeio')
library('ggtree')
library('rhdf5')
library('lsa')
library('ape')
library('readxl')

d <- h5read('GDC_samples_mutations.h5', '/GDC_samples_mutations')
# clinical <- as.matrix(read_excel('Clinical.xlsx', sheet = 1))

clinical <- read.table('GDC_samples.txt', header = TRUE, sep = '	')
# head(clinical)

samples_mutations <- d$block0_values


gender <- clinical[, 'gender']
histology <- clinical[, 'project_id']
sample_id <- clinical[, 'case_id']


combine <- function(x, y){paste(x, y, sep ='_')}


gender_histology <- mapply(combine, gender, histology)


label <- c()
st <- c()
g <- c()
i <- 1
for (id in d$axis1)
{
	index <- match(id, sample_id)
	label[i] <- gender_histology[index]
	st[i] <- histology[index]
	g[i] <- gender[index]
	i <- i + 1
}


dist1 <- as.matrix(1 - cosine(samples_mutations))
dist2 <- dist1
dist3 <- dist1

dimnames(dist1) <- list(label, label)
dimnames(dist2) <- list(g, g)
dimnames(dist3) <- list(st, st)


tr1 <- nj(dist1)
tr2 <- nj(dist2)
tr3 <- nj(dist3)


for (lo in c('slanted')) { #, 'fan', 'circular', 'equal_angle', 'daylight'

pdf(file = sprintf('plots/GDC_tree1_%s_bl.pdf', lo), width = 20, height = 25)
p1 <- ggtree(tr1, branch.length = 'none', layout = lo)  + geom_tippoint(aes(color = label), size = 1) # geom_tiplab()#
print(p1) 
dev.off()

pdf(file = sprintf('plots/GDC_tree2_%s_bl.pdf', lo), width = 20, height = 25)
p2 <- ggtree(tr2, branch.length = 'none', layout = lo) +  geom_tippoint(aes(color = label), size = 1) # geom_tiplab()# 
print(p2)  
dev.off()

pdf(file = sprintf('plots/GDC_tree3_%s_bl.pdf', lo), width = 20, height = 25)
p3 <- ggtree(tr3, branch.length = 'none', layout = lo) +  geom_tippoint(aes(color = label), size = 1) # geom_tiplab()# 
print(p3)  
dev.off()

}



