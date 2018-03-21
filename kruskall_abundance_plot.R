'Kruskal-Wallis differential abundance  analysis
#'
#'This function finds the features that are significantly
#'differentially abundant  in the provided taxa abundance data under different conditions
#'using Kruskal-Wallis test. The p-values values generated are corrected for multiple testing
#'using  family wise error rate. Significance is based on the corrected pvalue threshold. Significant features
#'then plotted on a boxplot
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#'        information, this should be one of the components of grouping vector.
#' @param  pvalue.threshold. Cut off p-value for significance of differentially abundant taxa, default is 0.05.
#' @param  Color. What variable in your mapping file should determine the color of your box plots (Phylum,class, order ect..)
#' @param  Facet How your box plots will be split up
#' @return Returns a list of three items: \itemize{
#'         \item  SignfeaturesTable: A \code{data.frame} of taxa with coresponding raw p-values, corrected p-values, family wise error rate and expected abundance
#'         computed by using raw p-values.
#'         \item  importance: A \code{data.frame} of taxa mean decrease in accuracy as obtained by random forest classifier.
#'         \item plotdata: A \code{data.frame} of taxa and corresponding corrected p-values, importance rank organised
#'         in form accepted by ggplot.
#'        }
#' @examples
#' #plot the significant features
#' kruskal_abundance_plot(group12367rarefyna1,'prebiotic',pvalue.threshold=0.1,Color='Genus', Facet='Genus')
#'
#' This code was modified from microbiomeSeq by Pedro J. Torres \email{pjtorres88@gmail.com} from the original code from the authors below to only have the 
#' signigficant kw and a list of the significant taxa, transform my data, and follow it up with normalization and boxplots of those significant taxa.
#' @references \url{http://userweb.eng.gla.ac.uk/umer.ijaz/}, Umer Ijaz, 2015
#' @references \url{http://www.bigre.ulb.ac.be/courses/statistics_bioinformatics/practicals/microarrays_berry_2010/berry_feature_selection.html}
#'
#' @author Pedro J. Torres \email{pjtorres88@gmail.com}, Alfred Ssekagiri \email{assekagiri@gmail.com}, Umer Zeeshan Ijaz \email{Umer.Ijaz@glasgow.ac.uk}
#'
#' @export kruskal_abundance
#'
kruskal_abundance_plot <- function(physeq, grouping_column,pvalue.threshold=0.05,Color,Facet)
{
  #get relative abudnace 
  physeq = transform_sample_counts(physeq, function(x) x/sum(x))
  abund_table <- otu_table(physeq)
  tax=data.frame(tax_table(physeq))
  tax$taxonomy=with(tax,paste0(Phylum,Class,Order,Family,Genus,Species))
  rownames(abund_table)=tax$taxonomy
  abund_table <- t(otu_table(physeq))#needed to transform it for this specific data need samples as rows and features as columns
  
  meta_table <-data.frame(sample_data(physeq))
  meta_table$Groups <- meta_table[,grouping_column]

  kruskal.wallis.table <- data.frame()
  data <- as.data.frame(abund_table)
  for (i in 1:dim(data)[2]){
    ks.test <- kruskal.test(data[,i], g=meta_table$Groups)
    # Store the result in the data frame
    kruskal.wallis.table <- rbind(kruskal.wallis.table,data.frame(id=names(data)[i],p.value=ks.test$p.value))
  }
  kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
  kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
  kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value, decreasing=FALSE), ]
  kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
  kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor#gives qvalue
  rownames(kruskal.wallis.table) <- kruskal.wallis.table$id

  #==================significant feature selection =================================#
  last.significant.element <- max(which(kruskal.wallis.table$q.value <= pvalue.threshold))
  selected <- 1:last.significant.element
  sig_res <-kruskal.wallis.table$id[selected]# gives you the actual variables
  subset.data<-data.frame(data[,as.character(kruskal.wallis.table[rownames(kruskal.wallis.table),"id"])])
  kruskal.wallis.table$id <- colnames(subset.data) #enforce that ids and colnames of subset data remain the same for easy indexing later on
  subset.data <- subset.data[,sig_res]


out <- list("SignfeaturesTable"=kruskal.wallis.table)
#return(out)
out=data.frame(out)
#print(colnames(out))
sig=subset(out, SignfeaturesTable.q.value < pvalue.threshold)
    
sig2=data.frame(sig)
svs=sig2$SignfeaturesTable.id
svs=gsub('X', '', svs)
lsvs=c(svs)
sig = list("SignificantTable"=sig)
    
    
#generate boxplot===========================================
      
#convert physeq to dataframe
mphyseq = psmelt(physeq)
#subset variables that passed p values
mphyseq<- mphyseq[mphyseq$OTU %in% lsvs,]
mphyseq$tax <- paste(mphyseq$OTU,mphyseq$Phylum,mphyseq$Class,mphyseq$Order,mphyseq$Family,mphyseq$Genus,mphyseq$Species)
print('Significant taxa')
print (unique(mphyseq$tax))

p = ggplot(data = mphyseq,  mapping = aes_string(x = grouping_column, y = "Abundance",
                                 color = Color, fill = Color)) +
    geom_boxplot(fill = NA) +
    geom_point(size = 1.5, alpha = 1,
                position = position_jitter(width = 0.3)) +
   facet_wrap(facets = Facet) + ylab('Relative Abundance') +
   scale_y_log10()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
print (p)
    
return (sig)

}

