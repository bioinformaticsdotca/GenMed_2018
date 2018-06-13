################################################################################
# Canadian Bioinformatics Workshop 2018
# Bioinformatics for Genomic Medicine: Module 6
#
# Created:  12-Jun-2018
# Author: Andrei Turinsky
###############################################################################

library(gplots)
library(lumi)
library(minfi)
library(matrixStats)
library(sva)

##########################################################
# Load data, which is extracted from the GEO dataset GSE52588

dataPath = '.'

pheno = read.csv('samplesheet.csv')
row.names(pheno) = pheno$Sample_Name

betas = read.csv('betas.csv', row.names = 1)
beta_matrix = data.matrix(betas)

anno = read.csv('annotations.csv', row.names = 1)

samplesDisease = rownames(pheno)[ pheno$Sample_Group == 'DS']
samplesControl = rownames(pheno)[ pheno$Sample_Group == 'control']

##########################################################
# Find disease signature

lumi::m2beta
lumi::beta2m

m_matrix = beta2m(beta_matrix)

# Find differentially methylated positions

dmp = dmpFinder(m_matrix, pheno$Sample_Group, type = "categorical")
head(dmp)

adjustMethod = 'bonf'

dmpCpgs = rownames(dmp)
dmp$betaDisease = rowMeans( beta_matrix[dmpCpgs, samplesDisease, drop=F]) 
dmp$betaControl = rowMeans( beta_matrix[dmpCpgs, samplesControl, drop=F]) 
dmp$deltaBeta = dmp$betaDisease - dmp$betaControl
dmp$adjPval = p.adjust(dmp$pval, method = adjustMethod)
dmp$genes = anno[dmpCpgs, "UCSC_RefGene_Name"]
dmp$chr = anno[dmpCpgs, "CHR"]

head(dmp)

# extract signature

pvalThresh = 0.01
deltaBetaThresh = 0.10

sigCpgs = dmpCpgs[ (dmp$adjPval < pvalThresh) & (abs(dmp$deltaBeta) > deltaBetaThresh) ]

cat("Found", length(sigCpgs), "signature CpGs:", "\n")
cat("\t", sum((dmp$pval < pvalThresh), na.rm=T), "CpGs with pval <", pvalThresh, "\n")
cat("\t", sum((dmp$adjPval < pvalThresh), na.rm=T), "CpGs with pval <", pvalThresh, "after", adjustMethod, "adjustment\n")
cat("\t", sum(abs(dmp$deltaBeta) > deltaBetaThresh, na.rm=T), "CpGs with abs(deltaBeta) >", deltaBetaThresh, "\n")
cat("\t", sum((dmp$adjPval < pvalThresh) & (abs(dmp$deltaBeta) > deltaBetaThresh), na.rm=T), "CpGs satisfy both\n")

##########################################################
# hierarchical clustering

beta_matrix_sig = beta_matrix[sigCpgs,]


dist_global = dist( t(beta_matrix) )
hc_global = hclust(dist_global)
plot(hc_global)

dist_signature = dist( t(beta_matrix_sig))
hc_signature = hclust(dist_signature)
plot(hc_signature)

hc_signature = hclust(dist_signature, method="ward.D2")
plot(hc_signature)

dist_signature = dist( t(beta_matrix_sig), method = 'manhattan')
hc_signature = hclust(dist_signature, method="ward.D2")
plot(hc_signature)

# extract cluster membership

clusters_signature = cutree(hc_signature, k = 2)	
clusters_signature

##########################################################
# heat maps

# basic heatmap

heatmap(beta_matrix_sig)

# better heatmaps

require(gplots)

heatmap.2(beta_matrix_sig)
heatmap.2(beta_matrix_sig, trace="none")

# make a pretty heatmap

dendrogram = "column" # c("both","row","column","none")
method = 'manhattan' 
clustMethod = "ward.D2" 

colorRampPalette = colorRampPalette(c("blue", "white",  "orange"))(n = 100)
sampleColors = ifelse(pheno$Sample_Group == 'DS', 'red', 'turquoise')

symmetric = FALSE
reorder = TRUE


heatmap.2( 
		beta_matrix_sig, 
		dendrogram = dendrogram,
		Rowv = reorder,
		Colv=if(symmetric) "Rowv" else reorder,
		symm=symmetric, 
		symkey=F, 
		symbreaks=F,
		trace="none", 
		density.info="density", #"density" "none"
		distfun = function(x) { dist(x, method = method) },
		hclustfun = function(x) hclust(x, method = clustMethod),
		col =  colorRampPalette, 
		key.title = "DNAm scale", 
		ColSideColors = sampleColors,
		labRow = ''
)

##########################################################
# Principal component analysis plots

pca_global = prcomp(t(beta_matrix))
pca_signature = prcomp(t(beta_matrix_sig))

samples = rownames(pheno)
sampleColors = ifelse(pheno$Sample_Group == 'DS', 'red', 'turquoise')
		
par(mfrow = c(1,2))

plot(pca_global$x[,1:2], type='n', main="PCA - global")
text(pca_global$x[,1:2], samples, col = sampleColors)

plot(pca_signature$x[,1:2], type='n', main="PCA - signature")
text(pca_signature$x[,1:2], samples, col = sampleColors)

par(mfrow = c(1,1))

##########################################################
# Classification plot - similar to (Choufani et al., 2015)

# split data into training and test sets, find a training-set signature

set.seed(1234)
trainSamples = sample( colnames(beta_matrix), size = 20)
testSamples = setdiff(colnames(beta_matrix), trainSamples)

dmp_train = dmpFinder(m_matrix[,trainSamples], 
                      pheno[trainSamples,]$Sample_Group, 
                      type = "categorical")
head(dmp_train)

train_sigCpgs = rownames(dmp_train)[ dmp_train$qval < 0.01]

# build DNAm profiles and score the training set

trainingDisease = intersect(trainSamples, samplesDisease)
trainingControl = intersect(trainSamples, samplesControl)

sigProfileDisease = rowMedians(beta_matrix[train_sigCpgs, trainingDisease])
sigProfileControl = rowMedians(beta_matrix[train_sigCpgs, trainingControl])

coords = sapply(testSamples, function(newSample) {
    sigProfile = beta_matrix[train_sigCpgs, newSample]
    x = cor(sigProfile, sigProfileControl)
    y = cor(sigProfile, sigProfileDisease)
    c(x,y)
} )
coords = t(coords)

# plot the classification score

minCorr = .7
plot( NULL, xlim=c(minCorr,1), ylim=c(minCorr,1), 
      xlab = "Similarity to control", ylab = "Similarity to disease")
abline(a=0, b=1, col = 'pink')

colors = ifelse(pheno[testSamples, "Sample_Group"] == 'DS', 'red', 'turquoise')
text(coords, testSamples, col = colors)


##########################################################
# DMRs Bumps

designMatrix = model.matrix(~ Sample_Group + Sex + Age, data = pheno)

# Explore the number of DMRs with different cutoffs, without permutations 

for( cutoff in seq(.1, .5, by = .1)) {
	cat('Cutoff: ', cutoff, '\n')
	dmrs = bumphunter(beta_matrix,
	                  chr = anno$CHR,
	                  pos = anno$MAPINFO,
	                  design = designMatrix, 
	                  cutoff = cutoff, 
	                  B=0, 
	                  type="Beta")
}

# run the bump hunting

date()
dmrs = bumphunter(beta_matrix,
                  chr = anno$CHR,
                  pos = anno$MAPINFO,
                  design = designMatrix, 
                  cutoff = .1, 
                  B=100, 
                  nullMethod = 'bootstrap',
                  type="Beta")
date()

head(dmrs$table)

selected_dmr_table = subset(dmrs$table, L >= 3) 

# plot DMR -  similar to Bioconductor package DMRcate

plotDMR = function(selectedDMR, beta_table, annotationTable, margins = 0) {
	
	dataTable = beta_table
	
	# retrieve array annotations
	
	arrayCpgs = rownames(annotationTable)
	chromosomes = annotationTable$CHR
	positions = annotationTable$MAPINFO
	names(positions) = arrayCpgs
	
	# find matching array probes
	
	chr = selectedDMR[,'chr'] 
	start = selectedDMR[,'start'] - margins
	end = selectedDMR[,'end'] + margins
	
	isMatchingCpg =  (chromosomes == chr) & (positions >= start) & (positions <= end)
	
	# match to the CpGs in the data and order by chromosomal position
	
	matchingCpgs = arrayCpgs [isMatchingCpg]
	matchingCpgs = intersect(matchingCpgs, rownames(dataTable))
	matchingCpgs = matchingCpgs[ order(positions[matchingCpgs]) ]
	
	matchingPositions = positions[matchingCpgs]
	
	# get a list of matching genes
	
	genes = as.character(annotationTable[matchingCpgs,'UCSC_RefGene_Name'])
	geneName = paste(sort(unique(unlist(strsplit(genes, ';')))), collapse = ', ')

	cat("Found ", length(matchingCpgs), "CpGs in the DMR\n")
	cat("REGION:", chr, start, end, '\n', "\tCPGS:", matchingCpgs, '\n')
	cat("\tGenes: ", geneName)
	
	# plot the DMR
	
	colorDisease = rgb(1,0,0, .5)
	colorControl = rgb(0,0.8,0, .5)
		
	plot(NULL, xlim=c(start, end), ylim=c(0,1), xlab  = paste0("DMR on ", chr, ": ", geneName) , ylab = "DNAm")
	
	# plot each CpG
	
	set.seed(1234)
	for(i in 1:length(matchingCpgs)) {
		cpg = matchingCpgs[i]
		
		jitterAmount = 0.01 * (end-start)
		
		pos_jittered = jitter(rep(matchingPositions[i], length(samplesDisease)), amount = jitterAmount)
		points( cbind( pos_jittered, dataTable[cpg, samplesDisease]), pch=19, col=colorDisease)
		
		pos_jittered = jitter(rep(matchingPositions[i], length(samplesControl)), amount = jitterAmount)
		points( cbind( pos_jittered, dataTable[cpg, samplesControl]), pch=4,  col=colorControl)
	}
	
	# plot averages as lines
	
	lines(matchingPositions, rowMeans(dataTable[matchingCpgs, samplesDisease]), col=colorDisease)
	lines(matchingPositions, rowMeans(dataTable[matchingCpgs, samplesControl]), col=colorControl)
}

plotDMR( selected_dmr_table[1,], beta_matrix, anno)
plotDMR( selected_dmr_table[1,], beta_matrix, anno, margins = 2000)

plotDMR( selected_dmr_table[2,], beta_matrix, anno)
plotDMR( selected_dmr_table[2,], beta_matrix, anno, margins = 2000)

##########################################################
# Batch correction

# simulate a contaminant

contaminant = 1

# add 15-20% contamination to randomly chosen samples

set.seed(1234)
compromisedSamples = sample( colnames(beta_matrix), size = 30)

batch_beta = beta_matrix
for(j in compromisedSamples) {
    q = runif(n=1, min = 0.02, max = .05)
    batch_beta[, j] = q * contaminant + (1-q) * batch_beta[,j]
}

# PCA plot

sampleNames = colnames(batch_beta)
batchColors = ifelse(sampleNames %in% compromisedSamples, 'magenta', 'blue')

par(mfrow=c(1,2))

pca_batch_compr = prcomp(t(batch_beta))
plot(pca_batch_compr$x[,1:2], type='n', main="PCA: Batch effect")
text(pca_batch_compr$x[,1:2], sampleNames, col = batchColors)

# batch correction

require(sva)

batch = as.numeric( colnames(beta_matrix) %in% compromisedSamples)
modcombat = model.matrix(~ Sample_Group, data=pheno)
combat_beta = ComBat(dat=batch_beta, batch=batch, mod=modcombat)

pca_combat = prcomp(t(combat_beta))
plot(pca_combat$x[,1:2], type='n', main="PCA: after ComBat")
text(pca_combat$x[,1:2], sampleNames, col = batchColors)
par(mfrow=c(1,1))


# deviations of batch-corrected data from the original data

ylims = c(-1,1) * .20

par(mfrow=c(1,2))
boxplot( batch_beta  - beta_matrix, ylim = ylims, col = batchColors, las=2, 
         main = "Deviations: Batch effect")
boxplot( combat_beta - beta_matrix, ylim = ylims, col = batchColors, las=2, 
         main = "Deviations: after ComBat")
par(mfrow=c(1,1))


##########################################################
# Influence of cell type composition

library(FlowSorted.Blood.450k)

data(FlowSorted.Blood.450k.compTable)
data(FlowSorted.Blood.450k.JaffeModelPars)

head(FlowSorted.Blood.450k.JaffeModelPars)

isSuspicious = p.adjust(FlowSorted.Blood.450k.compTable$p.value, 'bonf') < 0.05
suspiciousBloodCpgs_all = rownames(FlowSorted.Blood.450k.compTable)[ isSuspicious ]

mean( rownames(beta_matrix) %in% suspiciousBloodCpgs_all)
mean( sigCpgs %in% suspiciousBloodCpgs_all)
sum( sigCpgs %in% suspiciousBloodCpgs_all)

sigCpgs_clean = setdiff(sigCpgs, suspiciousBloodCpgs_all)
sigCpgs_clean

# Top 600 CpGs variable used in Jaffe-Irizarry cell subtype model

suspiciousBloodCpgs_top600 = rownames(FlowSorted.Blood.450k.JaffeModelPars)

mean(rownames(beta_matrix) %in% suspiciousBloodCpgs_top600)
mean(sigCpgs %in% suspiciousBloodCpgs_top600)
sum(sigCpgs %in% suspiciousBloodCpgs_top600)

##########################################################
# End of script

