
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(readr)


# sample info dataframe with array_id and chemical columns
samples <- read.csv('group_1_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)


#------Run Limma on Meth
rma.subset <- rma[paste0('X',samples$array_id[c(1:9,48:74)])]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control','3-METHYLCHOLANTHRENE')
  )
)
colnames(design) <- c('Intercept','3-METHYLCHOLANTHRENE')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
meth <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(meth,'meth_limma_results.csv')



#------Run Limma on Clot
rma.subset <- rma[paste0('X',samples$array_id[c(10:38, 48:74)])]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control', 'CLOTRIMAZOLE')
  )
)
colnames(design) <- c('Control','CLOTRIMAZOLE')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
clot <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(clot,'clot_limma_results.csv')



#------Run Limma on Chlor
rma.subset <- rma[paste0('X',samples$array_id[c(39:74)])]

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix(
  ~factor(
    samples$chemical,
    levels=c('Control','CHLOROFORM')
  )
)
colnames(design) <- c('Control','CHLOROFORM')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
chlor <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(chlor,'chlor_limma_results.csv')


#----Filter for differentially expressed 
meth1 <- meth %>%
  filter(adj.P.Val < 0.05) %>%
  count() %>%
  rename(Methylcholanthrene = n)

clot1 <- clot %>%
  filter(adj.P.Val < 0.05) %>%
  count()%>%
  rename(Clotrimazole = n)

chlor1 <- chlor %>%
  filter(adj.P.Val < 0.05) %>%
  count()%>%
  rename(Chloroform = n)

#create table for counts
tab1 <- meth1 %>% 
  merge(clot1) %>%
  merge(chlor1)


#plot table
png("significant.png", height = 50*nrow(tab1), width = 200*ncol(tab1))
grid.table(tab1)

#-----Set row names as sample names
meth$Methylcholanthrene <- row.names(meth)
row.names(meth) <- NULL

clot$Clotrimazole <- row.names(clot)
row.names(clot) <- NULL

chlor$Chloroform <- row.names(chlor)
row.names(chlor) <- NULL

#Select top 10 DE genes, rename to PROBEID
meth2 <- meth %>%
  select(Methylcholanthrene)%>%
  head(10)%>%
  rename(PROBEID = Methylcholanthrene)

clot2 <- clot %>%
  select(Clotrimazole)%>%
  head(10)%>%
  rename(PROBEID = Clotrimazole)

chlor2 <- chlor %>%
  select(Chloroform)%>%
  head(10)%>%
  rename(PROBEID = Chloroform)

#---- Join on Symbols to get the Symbol name to the sample name
meth2_aff <- left_join(meth2, affymap, on=PROBEID)
meth2_aff %>%
  select(GENENAME)%>%
  head(10) 

clot2_aff <- left_join(clot2, affymap, on=PROBEID)
clot2_aff %>%
  select(GENENAME)%>%
  head(10)

chlor2_aff <- left_join(chlor2, affymap, on=PROBEID)
chlor2_aff %>%
  select(GENENAME)%>%
  head(10)

meth2_aff <- left_join(meth2, affymap, on=PROBEID)
methaff <- meth2_aff %>%
  distinct(SYMBOL) %>%
  select(SYMBOL)%>%
  head(10) %>%
  rename(Methylcholanthrene = SYMBOL)

clot2_aff <- left_join(clot2, affymap, on=PROBEID)
clotaff <- clot2_aff %>%
  distinct(SYMBOL) %>%
  select(SYMBOL)%>%
  head(10)%>%
  rename(Clotrimazole = SYMBOL)

chlor2_aff <- left_join(chlor2, affymap, on=PROBEID)
chloraff <- chlor2_aff %>%
  distinct(SYMBOL) %>%
  select(SYMBOL)%>%
  head(10)%>%
  rename(Chloroform = SYMBOL)


# Print the tables
png("10meth.png", height = 50*nrow(methaff), width = 300*ncol(methaff))
grid.table(methaff)

png("10clot.png", height = 50*nrow(clotaff), width = 300*ncol(clotaff))
grid.table(clotaff)

png("10chlor.png", height = 50*nrow(chloraff), width = 300*ncol(chloraff))
grid.table(chloraff)


#----Print histogram of the DE genes
png('histmeth.png',width = 800,height = 300,res = 100)

meth <- meth %>%
  filter(adj.P.Val < 0.05)
p <- ggplot(meth, aes(logFC))+
  geom_histogram()+
  ggtitle("3-Methylcholanthrene significant DE genes")

print(p)

png('histclot.png',width = 800,height = 300,res = 100)

clot <- clot %>%
  filter(adj.P.Val < 0.05)
p <- ggplot(clot, aes(logFC))+
  geom_histogram()+
  ggtitle("Clotrimazole significant DE genes")

print(p)

png('histchlor.png',width = 800,height = 300,res = 100)

chlor <- chlor %>%
  filter(adj.P.Val < 0.05)
p <- ggplot(chlor, aes(logFC))+
  geom_histogram()+
  ggtitle("Chloroform significant DE genes")

print(p)


#-----Print scatterplots of the DE genes
png('scattermeth.png',width = 600,height = 300,res = 100)

p <- ggplot(meth, aes(x = logFC, y= P.Value))+
  geom_point()+
  ggtitle("3-Methylcholanthrene significant DE genes")

print(p)

png('scatterclot.png',width = 600,height = 300,res = 100)

p <- ggplot(clot, aes(x = logFC, y= P.Value))+
  geom_point()+
  ggtitle("Clotrimazole significant DE genes")

print(p)

png('scatterchlor.png',width = 600,height = 300,res = 100)

p <- ggplot(chlor, aes(x = logFC, y= P.Value))+
  geom_point()+
  ggtitle("Chloroform significant DE genes")

print(p)


#Load in affymap data
affymap <- read.csv('refseq_affy_map.csv')
affymap


#Load in DEseq results from programmer
deseq_meth <- read_csv("/projectnb/bf528/users/wheeler_2022/project_3/Analyst/deseq_1.csv")
deseq_clot <- read_csv("/projectnb/bf528/users/wheeler_2022/project_3/Analyst/ 2_deseq_results.csv")
deseq_chlor <- read_csv("/projectnb/bf528/users/wheeler_2022/project_3/Analyst/ 3_deseq_results.csv")
deseq_meth


#-----Meth Map the Gene names using the affymap to both analysis data
meth_limma <- meth %>%
  filter(adj.P.Val < 0.05)
meth_deseq <- deseq_meth %>%
  filter(padj < 0.05)

meth_limma$PROBEID <- row.names(meth_limma)
row.names(meth_limma) <- NULL
meth_deseq <- meth_deseq %>% rename(REFSEQ = ...1)
meth_deseq

meth_limma_aff <- inner_join(meth_limma, affymap, on=PROBEID)

inner_join(meth_limma_aff, meth_deseq, on=REFSEQ) %>% count()




#-----Clot Map the Gene names using the affymap to both analysis data
clot_limma <- clot %>%
  filter(adj.P.Val < 0.05)
clot_deseq <- deseq_clot %>%
  filter(padj < 0.05)

clot_limma$PROBEID <- row.names(clot_limma)
row.names(clot_limma) <- NULL
clot_deseq <- clot_deseq %>% rename(REFSEQ = ...1)

clot_limma_aff <- inner_join(clot_limma, affymap, on=PROBEID)

inner_join(clot_limma_aff, clot_deseq, on=REFSEQ) %>% count()



#-----Chlor Map the Gene names using the affymap to both analysis data
chlor_limma <- chlor %>%
  filter(adj.P.Val < 0.05)
chlor_deseq <- deseq_chlor %>%
  filter(padj < 0.05)

chlor_limma$PROBEID <- row.names(chlor_limma)
row.names(chlor_limma) <- NULL
chlor_deseq <- chlor_deseq %>% rename(REFSEQ = ...1)

chlor_limma_aff <- inner_join(chlor_limma, affymap, on=PROBEID)

inner_join(chlor_limma_aff, chlor_deseq, on=REFSEQ) %>% count()


#------Concordance function
f.x <- function(n0,n1,n2,N) {
  x <- (n0*N-n1*n2)/(n0+N-n1-n2)
  if(x <= n0) {
    x
  } else {
    NA
  }
}


#Calculating and printing concordance 
concordance <- c(f.x(21,58,313,20000), f.x(768,5803,913,20000),f.x(1787,11407,1810,20000))

condf <- do.call(rbind, Map(data.frame, toxin=c("meth", "clot", "chlor"), Concordance=concordance, deseq=c(313, 931, 1810), limma=c(58, 5803, 11407)))
condf

png('deseq_concordance_plot.png',width = 600,height = 250,res = 100)


#---Plotting concordance and # DEGs
p <- ggplot(condf, aes(x = deseq, y=Concordance, label=toxin))+
  geom_point(alpha=.25)+
  geom_text(hjust=0.5, vjust=0)+
  xlab("Treatment effect (number of DEGs)")+
  ylab("Concordance of DEGs")+
  ggtitle("Concordance of DEGs for DESeq analysis")

print(p)

png('limma_concordance_plot.png',width = 600,height = 250,res = 100)

p <- ggplot(condf, aes(x = limma, y=Concordance, label=toxin))+
  geom_point(alpha=.25)+
  geom_text(hjust=0.5, vjust=0)+
  xlab("Treatment effect (number of DEGs)")+
  ylab("Concordance of DEGs")+
  ggtitle("Concordance of DEGs for limma analysis")

print(p)


#----Subsetting data to above and below mean
above_meth <- meth_deseq %>%
  filter(baseMean > median(baseMean))
below_meth <- meth_deseq %>%
  filter(baseMean < median(baseMean))
above_clot <- clot_deseq %>%
  filter(baseMean > median(baseMean))
below_clot <- clot_deseq %>%
  filter(baseMean < median(baseMean))
above_chlor <- chlor_deseq %>%
  filter(baseMean > median(baseMean))
below_chlor <- chlor_deseq %>%
  filter(baseMean < median(baseMean))

#getting the values joined
print(inner_join(meth_limma_aff, above_meth, on=REFSEQ) %>% count())
print(inner_join(meth_limma_aff, below_meth, on=REFSEQ) %>% count())
print(inner_join(clot_limma_aff, above_clot, on=REFSEQ) %>% count())
print(inner_join(clot_limma_aff, below_clot, on=REFSEQ) %>% count())
print(inner_join(chlor_limma_aff, above_chlor, on=REFSEQ) %>% count())
print(inner_join(chlor_limma_aff, below_chlor, on=REFSEQ) %>% count())

c(11,437,1037)
c(10, 329, 750)

#Calculating the different concordances
above_concordance <- c(f.x(11,58,156,20000), f.x(437,5803,465,20000),f.x(1037,11407,905,20000))
below_concordance <- c(f.x(10,58,156,20000), f.x(329,5803,465,20000),f.x(750,11407,905,20000))


allcondf <- do.call(rbind, Map(data.frame, toxin=c("3-Methylcholanthrene", "Clotrimazole", "Chloroform"), median=concordance, above = above_concordance, below = below_concordance))
allcondf

allcondf_1 <- pivot_longer(allcondf, cols=c('median','above', 'below'))


#----Plotting the concordance barplot
png('allconcordance.png',width = 600,height = 400,res = 100)
p <- ggplot(allcondf_1, aes(fill=toxin, y=value, x=name)) +
  geom_bar(position='dodge', stat='identity')+
  xlab("Concordance of DEGs")+
  ylab("Subgroups based on median")+
  ggtitle("Concordance values for subgroups above, below and at baseMean median")
print(p)



meth_d <- deseq_meth %>%
  filter(padj < 0.05) %>%
  count() %>%
  rename(Methylcholanthrene = n)

clot_d <- deseq_clot %>%
  filter(padj < 0.05) %>%
  count()%>%
  rename(Clotrimazole = n)

chlor_d <- deseq_chlor %>%
  filter(padj < 0.05) %>%
  count()%>%
  rename(Chloroform = n)

tab_d <- meth_d %>% 
  merge(clot_d) %>%
  merge(chlor_d)
tab_d

png("significant.png", height = 50*nrow(tab1), width = 200*ncol(tab1))
grid.table(tab1)

meth2 <- deseq_1 %>%
  filter(padj < 0.05) %>%
  rename(Methylcholanthrene =...1)
meth1 <- meth %>%
  filter(adj.P.Val < 0.05)

clot2 <- deseq_clot %>%
  filter(padj < 0.05) %>%
  rename(Clotrimazole =...1)
clot1 <- clot %>%
  filter(adj.P.Val < 0.05)








