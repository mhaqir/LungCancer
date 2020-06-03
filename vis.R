#install.packages('treeio')
library('treeio')
library('ggtree')
library('rhdf5')
library('lsa')
library('ape')
library('readxl')

d <- h5read('samples_mutations.h5', '/samples_mutations')
clinical <- as.matrix(read_excel('Clinical.xlsx', sheet = 1))

# head(clinical)

samples_mutations <- d$block0_values


gender <- clinical[, 'Gender']
histology <- clinical[, 'Histology']
sample_id <- clinical[, 'Tumor_Sample_Barcode']

combine <- function(x, y){paste(x, y, sep ='_')}
# iden <- function(x){c(x)}

gender_histology <- mapply(combine, gender, histology)

# histology <- mapply(iden, hist)
# dim(gender_histology)
# dim(histology)

# typeof(histology)
# typeof(gender_histology)

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
# factor(st)
# label
# length(st)
# length(label)

dist1 <- as.matrix(1 - cosine(samples_mutations))
dimnames(dist1) <- list(label, label)

dist2 <- as.matrix(1 - cosine(samples_mutations))
dimnames(dist2) <- list(g, g)


dist3 <- as.matrix(1 - cosine(samples_mutations))
dimnames(dist3) <- list(st, st)
# dim(dist2)

tr1 <- nj(dist1)
tr2 <- nj(dist2)
tr3 <- nj(dist3)

# gs = factor(label)

# dist1
p1 <- ggtree(tr1, branch.length = 'none', layout = 'daylight')  + geom_tippoint(aes(color = label), size = 5) # geom_tiplab()# 

pdf(file = 'tree1.pdf', width = 10, height = 10)
p1
dev.off()

# v <- rownames(dist2)

# v <- data.frame(st = rownames(dist2))
# row.names(v) <- NULL
#geom_tiplab() , branch.length = 'none', layout = 'daylight'
p2 <- ggtree(tr2, branch.length = 'none', layout = 'daylight') +  geom_tippoint(aes(color = label), size = 5) # geom_tiplab()# 

pdf(file = 'tree2.pdf', width = 10, height = 10)
p2
dev.off()


p3 <- ggtree(tr3, branch.length = 'none', layout = 'daylight') +  geom_tippoint(aes(color = label), size = 5) # geom_tiplab()# 

pdf(file = 'tree3.pdf', width = 10, height = 10)
p3
dev.off()

