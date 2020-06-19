#install.packages('treeio')
library('treeio')
library('ggtree')
library('rhdf5')
library('lsa')
library('ape')
library('readxl')

d <- h5read('GDCbaa_BMR_samples_mutations.h5', '/GDCbaa_BMR_samples_mutations')
# d <- h5read('GDC_samples_mutations.h5', '/GDC_samples_mutations')

clinical_GDC <- read.table('GDC_samples_BAA.txt', header = TRUE, sep = '	')
clinical_BMR <- as.matrix(read_excel('Clinical.xlsx', sheet = 1))


samples_mutations <- d$block0_values



gender_BMR <- clinical_BMR[, 'Gender']
gender_GDC <- clinical_GDC[, 'gender']

# print(gender_BMR)
# print(gender_GDC)
# histology <- clinical[, 'project_id']
sample_id_BMR <- clinical_BMR[, 'Tumor_Sample_Barcode']
sample_id_GDC <- clinical_GDC[, 'case_id']
combine <- function(x, y){paste(x, y, sep ='_')}
# gender_histology <- mapply(combine, gender, histology)



race <- c()
gender <- c()
ids <- c()
i <- 1
for (id in d$axis1)
{

	if (substr(id, start =1, stop = 8) == 'combined'){
		race[i] <- 'AA'
		index <- match(id, sample_id_BMR)
		gender[i] <- gender_BMR[index]
		ids[i] <- id

	}
	else{
		race[i] <- 'AA_TCGA'
		index <- match(id, sample_id_GDC)
		gender[i] <- gender_GDC[index]
		ids[i] <- id
	}
	i <- i + 1
}


# gender_race_id <- mapply(combine, race, gender, ids)

gender_race <- mapply(combine, gender, race)

# print(race)
# print(gender)
# print(ids)
# print(gender_race_id)

# length(race)
# length(gender)
# length(gender_race_id)
# length(ids)
# dim(gender_race_id)

dist1 <- as.matrix(1 - cosine(samples_mutations))
dist2 <- dist1

dimnames(dist1) <- list(gender_race, gender_race)
dimnames(dist2) <- list(race, race)

tr1 <- nj(dist1)
tr2 <- nj(dist2)


for (lo in c('slanted', 'fan', 'circular', 'equal_angle', 'daylight')) { #, 'fan', 'circular', 'equal_angle', 'daylight'

pdf(file = sprintf('plots/GDCbaa_BMR_tree1_%s_g_r.pdf', lo), width = 15, height = 14)
p1 <- ggtree(tr1, branch.length = 'none', layout = lo)  + geom_tippoint(aes(color = label), size = 1) # geom_tiplab()#
print(p1) 
dev.off()


pdf(file = sprintf('plots/GDCbaa_BMR_tree1_%s_g_r.pdf', lo), width = 15, height = 14)
p2 <- ggtree(tr2, branch.length = 'none', layout = lo)  + geom_tippoint(aes(color = label), size = 1) # geom_tiplab()#
print(p2) 
dev.off()

}



