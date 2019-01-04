# if having issues installing libraries just type install.packages('nameofprogram',repos="http://cran.cnr.berkeley.edu")
# import phyloseq and dependencies
library("phyloseq")
packageVersion("phyloseq")

library("ggplot2")
packageVersion("ggplot2")

library("scales")
packageVersion("scales")

library("grid")
packageVersion("grid")

library("DESeq2")
packageVersion("DESeq2")
library("biom") #allows us to read the biom table
packageVersion("biom")
x=read_biom("/Users/Pedro_Torres/Desktop/PCOS_Analysis/MGRAST2_3/QIIME.Subsystems.data.biom")
otumat = as(biom_data(x), "matrix")
OTU = otu_table(otumat, taxa_are_rows=TRUE)
taxmat = as.matrix(observation_metadata(x), rownames.force=TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)
# all this might take a while, but 

#prepare to load mapping file
map_file=('/Users/Pedro_Torres/Desktop/PCOS_Analysis/MGRAST2_3/QIIME_PCOS_mappingfile.txt')
bmsd=import_qiime_sample_data(map_file)
class(bmsd)
dim(bmsd)
#merge mapping and Biom
PCOS=merge_phyloseq(physeq,bmsd)
PCOS
sample_sums(PCOSP5vsL5)
#Begin log2 change
#will remove PO andLO just compare P5and L5
PCOSP5vsL5= subset_samples(PCOS, Time !="0")

PCOSdeseq=phyloseq_to_deseq2(PCOSP5vsL5,~Treatment)
#calculating geometric means, this is not arithmetic mean, but will be normalized
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(PCOSdeseq), 1, gm_mean)
PCOSdeseq = estimateSizeFactors(PCOSdeseq, geoMeans = geoMeans)
PCOSdeseq = DESeq(PCOSdeseq, fitType="local")

res = results(PCOSdeseq)

res = res[order(res$padj, na.last=NA), ]
head(res)
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(PCOS)[rownames(sigtab), ], "matrix"))
sigtab
head(sigtab, contrasts=list("P","L"))

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "ontology1", "ontology2", "ontology3", "ontology4")]

library("ggplot2")
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(ontology3))
# Phylum order


x = tapply(sigtabgen$log2FoldChange, sigtabgen$ontology1, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$ontology1 = factor(as.character(sigtabgen$ontology1), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$ontology3, function(x) max(x))
x = sort(x, TRUE)
x
sigtabgen$ontology3 = factor(as.character(sigtabgen$ontology3), levels=names(x))
ggplot(sigtabgen, aes(y=ontology3, x=log2FoldChange, color=ontology1)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  
  theme(axis.text.x = element_text(angle = -90, hj



ust = 0, vjust=0.5))
posigtab















pcos.t5=subset_samples(phylo, time == '5')
pcos.t4=subset_samples(phylo, time == '4')











dim(Time4_Deseq_sigtab)
dim(Time5_Deseq_sigtab)
write.table(x = Time4_Deseq_sigtab, file = "Reviewers/Pubertal/DESeq2/Time4_Deseq_sigtab.csv",sep = ",", row.names = FALSE)
write.table(x = Time5_Deseq_sigtab, file = "Reviewers/Pubertal/DESeq2/Time5_Deseq_sigtab.csv",sep = ",", row.names = FALSE)

theme_set(theme_bw())

scale_fill_discrete <- function(palname = "Set1"){
    scale_fill_brewer(palette = palname)
}

Time4_X <- tapply(Time4_Deseq_sigtab$log2FoldChange, Time4_Deseq_sigtab$Phylum, function(Time4_X) max(Time4_X))
Time4_X <- sort(Time4_X, TRUE)
Time4_Deseq_sigtab$Phylum <- factor(as.character(Time4_Deseq_sigtab$Phylum), levels = names(Time4_X))

Time4_Family_X <- tapply(Time4_Deseq_sigtab$log2FoldChange, Time4_Deseq_sigtab$Family, function(Time4_Family_X) max(Time4_Family_X))
Time4_Family_X <- sort(Time4_Family_X, TRUE)
Time4_Deseq_sigtab$Family <- factor(as.character(Time4_Deseq_sigtab$Family), levels = names(Time4_Family_X))

Time4_Genus_X <- tapply(Time4_Deseq_sigtab$log2FoldChange, Time4_Deseq_sigtab$Genus, function(Time4_Genus_X) max(Time4_Genus_X))
Time4_Genus_X <- sort(Time4_Genus_X, TRUE)
Time4_Deseq_sigtab$Genus <- factor(as.character(Time4_Deseq_sigtab$Genus), levels = names(Time4_Genus_X))

Time4_Species_X <- tapply(Time4_Deseq_sigtab$log2FoldChange, Time4_Deseq_sigtab$Species, function(Time4_Species_X) max(Time4_Species_X))
Time4_Species_X <- sort(Time4_Species_X, TRUE)
Time4_Deseq_sigtab$Species <- factor(as.character(Time4_Deseq_sigtab$Species), levels = names(Time4_Species_X))

#pdf('Reviewers/Adult/DESeq2/T4_log2fold.pdf')
ggplot(Time4_Deseq_sigtab, aes(y=Genus, x=log2FoldChange, color=Family)) +
geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
geom_point(size=6) +

theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ggtitle("Time4 - Placebo vs Letrozole")



ggplot(Time4_Deseq_sigtab, aes(y=Order, x=log2FoldChange, color=Genus)) +
geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
geom_point(size=6) +

theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ggtitle("Time4 - Placebo vs Letrozole")

ggplot(Time4_Deseq_sigtab, aes(y=Phylum, x=log2FoldChange, color=Family)) +
geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
geom_point(size=6) +

theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ggtitle("Time4 - Placebo vs Letrozole")


#dev.off()




