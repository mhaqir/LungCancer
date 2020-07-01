#install.packages('treeio')
library('treeio')
library('ggtree')
library('rhdf5')
library('lsa')
library('ape')
library('readxl')

d <- h5read('GDCw_baa_samples_mutations_cc.h5', '/GDCw_baa_samples_mutations_cc')
# clinical <- as.matrix(read_excel('Clinical.xlsx', sheet = 1))

clinical <- read.table('GDC_samples_white_BAA.txt', header = TRUE, sep = '	')
# head(clinical)

samples_mutations <- d$block0_values
# dim(samples_mutations)
root_node <- matrix(0, nrow = dim(samples_mutations)[1])
samples_mutations <- cbind(samples_mutations, root_node)
# dim(samples_mutations)
# head(samples_mutations)
# typeof(samples_mutations)

gender <- clinical[, 'gender']
histology <- clinical[, 'project_id']
sample_id <- clinical[, 'case_id']
# race <- clinical[, 'race']

j <- 1
race <- c()
for (item in clinical[, 'race']){
	if (item == 'black or african american'){
		race[j] <- 'AA'
	}
	else{
		race[j] <- 'white'
	}
	j <- j + 1
}


combine <- function(x, y){paste(x, y, sep ='_')}


gender_histology <- mapply(combine, gender, histology)
race_histology <- mapply(combine, race, histology)
gender_race <- mapply(combine, gender, race)

gh <- c()
h <- c()
g <- c()
r <- c()
rh <- c()
gr <- c()
i <- 1
for (id in d$axis1)
{
	index <- match(id, sample_id)
	gh[i] <- gender_histology[index]
	h[i] <- histology[index]
	g[i] <- gender[index]
	r[i] <- race[index]
	rh[i] <- race_histology[index]
	gr[i] <- gender_race[index]
	i <- i + 1
}

gh[dim(samples_mutations)[2]] <- 'root'
h[dim(samples_mutations)[2]] <- 'root'
g[dim(samples_mutations)[2]] <- 'root'
r[dim(samples_mutations)[2]] <- 'root'
rh[dim(samples_mutations)[2]] <- 'root'
gr[dim(samples_mutations)[2]] <- 'root'

# length(gr)

# a <- matrix(0, nrow = 4)
# a
# b <- matrix(1:10, nrow = 4)
# b
# c <- cbind(b, a)
# c
# dim(samples_mutations)
# length(gh)
# length(h)
# length(g)
# length(r)
# length(rh)
# length(gr)



dist1 <- as.matrix(1 - cosine(samples_mutations))
dist2 <- dist1
dist3 <- dist1
dist4 <- dist1
dist5 <- dist1
dist6 <- dist1

dimnames(dist1) <- list(gh, gh)
dimnames(dist2) <- list(g, g)
dimnames(dist3) <- list(h, h)
dimnames(dist4) <- list(r, r)
dimnames(dist5) <- list(rh, rh)
dimnames(dist6) <- list(gr, gr)
# dist1


tr1 <- njs(dist1)

# is.rooted(tr1)
tr1 <- root(tr1, outgroup ='root', resolve.root = TRUE)
# tr1
# tr1$root.edge <- 0

# is.rooted(tr1)
# tr1
# tr1$root
tr2 <- nj(dist2)
tr2 <- root(tr2, outgroup ='root', resolve.root = TRUE)
tr3 <- nj(dist3)
tr3 <- root(tr3, outgroup ='root', resolve.root = TRUE)
tr4 <- nj(dist4)
tr4 <- root(tr4, outgroup ='root', resolve.root = TRUE)
tr5 <- nj(dist5)
tr5 <- root(tr5, outgroup ='root', resolve.root = TRUE)
tr6 <- nj(dist6)
tr6 <- root(tr6, outgroup ='root', resolve.root = TRUE)

for (lo in c('slanted', 'circular', 'equal_angle', 'daylight')) { #, 'fan', 'circular', 'equal_angle', 'daylight'

pdf(file = sprintf('plots/GDCw_baa_tree1_%s_cc_rooted.pdf', lo), width = 20, height = 25)
p1 <- ggtree(tr1, branch.length = 'none', layout = lo)  + geom_tippoint(aes(color = label), size = 1) # geom_tiplab()#
print(p1) 
dev.off()

pdf(file = sprintf('plots/GDCw_baa_tree2_%s_cc_rooted.pdf', lo), width = 20, height = 25)
p2 <- ggtree(tr2, branch.length = 'none', layout = lo) +  geom_tippoint(aes(color = label), size = 1) # geom_tiplab()# 
print(p2)  
dev.off()

pdf(file = sprintf('plots/GDCw_baa_tree3_%s_cc_rooted.pdf', lo), width = 20, height = 25)
p3 <- ggtree(tr3, branch.length = 'none', layout = lo) +  geom_tippoint(aes(color = label), size = 1) # geom_tiplab()# 
print(p3)  
dev.off()

pdf(file = sprintf('plots/GDCw_baa_tree4_%s_cc_rooted.pdf', lo), width = 20, height = 25)
p4 <- ggtree(tr4, branch.length = 'none', layout = lo) +  geom_tippoint(aes(color = label), size = 1) # geom_tiplab()# 
print(p4)  
dev.off()

pdf(file = sprintf('plots/GDCw_baa_tree5_%s_cc_rooted.pdf', lo), width = 20, height = 25)
p5 <- ggtree(tr5, branch.length = 'none', layout = lo) +  geom_tippoint(aes(color = label), size = 1) # geom_tiplab()# 
print(p5)  
dev.off()

pdf(file = sprintf('plots/GDCw_baa_tree6_%s_cc_rooted.pdf', lo), width = 20, height = 25)
p6 <- ggtree(tr6, branch.length = 'none', layout = lo) +  geom_tippoint(aes(color = label), size = 1) # geom_tiplab()# 
print(p6)  
dev.off()
}



