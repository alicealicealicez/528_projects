library("ggplot2")
library("tidyverse")

#Part 1 

for (i in c("Ad_1", "Ad_2", "P4_1", "P4_2", "P7_1", "P7_2")){
  n = paste("fpkmtracking/",i,".fpkm_tracking",sep = "")
  name = paste(i,"_df", sep = "")
  assign(name, read_tsv(n))
}

p0_Ad_df <- read_tsv('programmer/cuffdiff_out/genes.fpkm_tracking')

#Create lists for different gene ontology groups 

Sarcomere_GO <- c("Pdlim5", "Pygm", "Myoz2", "Des", "Csrp3", "Tcap", "Cryab")
Mitochondria_GO <- c("Mpc1", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh")
CellCycle_GO <- c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "E2f1", "Cdc27", "Bora", "Cdc45", "Rad51", "Aurkb", "Cdc23")

#Get FPKM values for Sarcomere

p0_Ad_df <- p0_Ad_df %>%
  select(gene_short_name, P0_FPKM)
p0_Ad_sarcomere <- p0_Ad_df %>%
  filter(gene_short_name %in% Sarcomere_GO)

Ad_1_sarcomere <- Ad_1_df %>% 
  select(gene_short_name, FPKM)%>%
  filter(gene_short_name %in% Sarcomere_GO)%>%
  rename(Ad_1 = FPKM) 
Ad_2_sarcomere <- Ad_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% Sarcomere_GO)%>%
  rename(Ad_2 = FPKM) 
P4_1_sarcomere <- P4_1_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% Sarcomere_GO)%>%
  rename(P4_1 = FPKM) 
P4_2_sarcomere <- P4_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% Sarcomere_GO)%>%
  rename(P4_2 = FPKM) 
P7_1_sarcomere <- P7_1_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% Sarcomere_GO)%>%
  rename(P7_1 = FPKM) 
P7_2_sarcomere <- P7_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% Sarcomere_GO)%>%
  rename(P7_2 = FPKM) 

#Merge data to create plot

fpkms <- p0_Ad_sarcomere %>% 
  left_join( Ad_1_sarcomere,by="gene_short_name")%>%
  left_join( Ad_2_sarcomere,by="gene_short_name")%>%
  left_join( P4_1_sarcomere,by="gene_short_name")%>%
  left_join( P4_2_sarcomere,by="gene_short_name")%>%
  left_join( P7_1_sarcomere,by="gene_short_name")%>%
  left_join( P7_2_sarcomere,by="gene_short_name")
fpkms_long <- pivot_longer(fpkms, 
                           cols = P0_FPKM:P7_2,
                           names_to = "Samples",
                           values_to = "FPKM")

#Deliverable 1 (2/3)

png('Sarcomere.png',width = 800,height = 600,res = 100)
p <- ggplot(data = fpkms_long, aes(x = Samples, y = FPKM, color = gene_short_name, group = gene_short_name))+
  geom_point()+
  geom_line()+
  ggtitle("Sacromere")

print(p)
dev.off()

#Filtering for Mitochondrial GO genes

p0_Ad_sarcomere <- p0_Ad_df %>%
  filter(gene_short_name %in% Mitochondria_GO)

Ad_1_sarcomere <- Ad_1_df %>% 
  select(gene_short_name, FPKM)%>%
  filter(gene_short_name %in% Mitochondria_GO)%>%
  rename(Ad_1 = FPKM) 
Ad_2_sarcomere <- Ad_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% Mitochondria_GO)%>%
  rename(Ad_2 = FPKM) 
P4_1_sarcomere <- P4_1_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% Mitochondria_GO)%>%
  rename(P4_1 = FPKM) 
P4_2_sarcomere <- P4_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% Mitochondria_GO)%>%
  rename(P4_2 = FPKM) 
P7_1_sarcomere <- P7_1_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% Mitochondria_GO)%>%
  rename(P7_1 = FPKM) 
P7_2_sarcomere <- P7_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% Mitochondria_GO)%>%
  rename(P7_2 = FPKM) 

fpkms <- p0_Ad_sarcomere %>% 
  left_join( Ad_1_sarcomere,by="gene_short_name")%>%
  left_join( Ad_2_sarcomere,by="gene_short_name")%>%
  left_join( P4_1_sarcomere,by="gene_short_name")%>%
  left_join( P4_2_sarcomere,by="gene_short_name")%>%
  left_join( P7_1_sarcomere,by="gene_short_name")%>%
  left_join( P7_2_sarcomere,by="gene_short_name")
fpkms_long <- pivot_longer(fpkms, 
                           cols = P0_FPKM:P7_2,
                           names_to = "Samples",
                           values_to = "FPKM")

#Deliverable 1 (2/3)

png('Mitochondria.png',width = 800,height = 600,res = 100)
p <- ggplot(data = fpkms_long, aes(x = Samples, y = FPKM, color = gene_short_name, group = gene_short_name))+
  geom_point()+
  geom_line()+
  ggtitle("Mitochondria")

print(p)
dev.off()

#Filtering for Cell Cycle Go genes

p0_Ad_sarcomere <- p0_Ad_df %>%
  filter(gene_short_name %in% CellCycle_GO)

Ad_1_sarcomere <- Ad_1_df %>% 
  select(gene_short_name, FPKM)%>%
  filter(gene_short_name %in% CellCycle_GO)%>%
  rename(Ad_1 = FPKM) 
Ad_2_sarcomere <- Ad_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% CellCycle_GO)%>%
  rename(Ad_2 = FPKM) 
P4_1_sarcomere <- P4_1_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% CellCycle_GO)%>%
  rename(P4_1 = FPKM) 
P4_2_sarcomere <- P4_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% CellCycle_GO)%>%
  rename(P4_2 = FPKM) 
P7_1_sarcomere <- P7_1_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% CellCycle_GO)%>%
  rename(P7_1 = FPKM) 
P7_2_sarcomere <- P7_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% CellCycle_GO)%>%
  rename(P7_2 = FPKM) 

fpkms <- p0_Ad_sarcomere %>% 
  left_join( Ad_1_sarcomere,by="gene_short_name")%>%
  left_join( Ad_2_sarcomere,by="gene_short_name")%>%
  left_join( P4_1_sarcomere,by="gene_short_name")%>%
  left_join( P4_2_sarcomere,by="gene_short_name")%>%
  left_join( P7_1_sarcomere,by="gene_short_name")%>%
  left_join( P7_2_sarcomere,by="gene_short_name")
fpkms_long <- pivot_longer(fpkms, 
                           cols = P0_FPKM:P7_2,
                           names_to = "Samples",
                           values_to = "FPKM")

#Deliverable 1 (3/3)

png('Cell Cycle.png',width = 800,height = 600,res = 100)
p <- ggplot(data = fpkms_long, aes(x = Samples, y = FPKM, color = gene_short_name, group = gene_short_name))+
  geom_point()+
  geom_line()+
  ggtitle("Cell Cycle")

print(p)
dev.off()

#Part 3 download differential expression data

expdiff <- read_tsv('programmer/cuffdiff_out/gene_exp.diff')

#subset for top 1000 expressed genes

expdiff <- expdiff %>%
  select(q_value, gene)%>%
  arrange(q_value)%>%
  head(1000)

#Filter for top expressed genes, join FPKM data for all samples

p0_Ad_sarcomere <- p0_Ad_df %>%
  filter(gene_short_name %in% expdiff$gene)

Ad_1_sarcomere <- Ad_1_df %>% 
  select(gene_short_name, FPKM)%>%
  filter(gene_short_name %in% expdiff$gene)%>%
  rename(Ad_1 = FPKM) 
Ad_2_sarcomere <- Ad_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% expdiff$gene)%>%
  rename(Ad_2 = FPKM) 
P4_1_sarcomere <- P4_1_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% expdiff$gene)%>%
  rename(P4_1 = FPKM) 
P4_2_sarcomere <- P4_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% expdiff$gene)%>%
  rename(P4_2 = FPKM) 
P7_1_sarcomere <- P7_1_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% expdiff$gene)%>%
  rename(P7_1 = FPKM) 
P7_2_sarcomere <- P7_2_df %>% 
  select(FPKM, gene_short_name)%>%
  filter(gene_short_name %in% expdiff$gene)%>%
  rename(P7_2 = FPKM) 

fpkms <- p0_Ad_sarcomere %>% 
  left_join( Ad_1_sarcomere,by="gene_short_name")%>%
  left_join( Ad_2_sarcomere,by="gene_short_name")%>%
  left_join( P4_1_sarcomere,by="gene_short_name")%>%
  left_join( P4_2_sarcomere,by="gene_short_name")%>%
  left_join( P7_1_sarcomere,by="gene_short_name")%>%
  left_join( P7_2_sarcomere,by="gene_short_name")

#plot heatmap

png('diffexp.png',width = 800,height = 1000,res = 100)
p <- as.matrix(fpkms[, -1])
rownames(p) <- fpkms$gene_short_name
print(heatmap(p))

dev.off()

#looking at sections with high FPKM values in p0 samples

hpla <- fpkms %>%
  filter(P0_FPKM > mean(P0_FPKM))

# printing a heatmap for visualization

png('diffexphpla.png',width = 800,height = 1000,res = 100)
p <- as.matrix(hpla[, -1])
rownames(p) <- hpla$gene_short_name
print(heatmap(p))





