#'Plot randomforest
#'
#' @author Pedro J. Torres \email{pjtorres88@gmail.com}
#'This function finds the top 10 features from a phyloseq object that have the greatest influence of classification accuracy using 
#randomforest and plots those features in a box plot
#'
#' @param physeq (Required). A \code{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param randomForest (Required)- do random forest
#' @param scales (Required)- will help with scaling of boxplot
#' @param grid (Required)- will help with scaling of boxplot
#' @param knitr (Required)- give a nice output of the important feature in a random foret model
#' @param stringr (Required)- allows modificaiton of important OTUS
#' @param repr (Required)- allows modificaiton of boxplot aka bigger plots!

#' @param grouping_column (Required). Character string specifying name of a categorical variable that is preffered for grouping the information.
#'        information, this should be one of the components of grouping vector.
#' @param  Color. What variable in your mapping file should determine the color of your box plots (Phylum,class, order ect..)
#' @param  Facet How your box plots will be split up

#' @examples
#' #plot the significant features
#' plot_randomforest(physeqobject,grouping_column='Treatment', Color='Family', Facet='Genus')
#'
#' This code was modified by Pedro J. Torres \email{pjtorres88@gmail.com} a number of phyloseq tutorials out there
#'
#'

plot_randomforest=function(physeq, grouping_column,Color,Facet){
    
    #filter out any taxa with a zero or in very low abudnace and prep for random forest, will not make training set here
    physeq_0_filtered <- prune_taxa(taxa_sums(physeq) > 0,physeq)
    physeq_10_filtered <- prune_taxa(taxa_sums(physeq_0_filtered)>5,physeq_0_filtered)
    predictors <- t(otu_table(physeq_10_filtered))
    physeq_10_filtered_sample_Data <- data.frame(sample_data(physeq_10_filtered))
    physeq_10_filtered_sample_Data$Groups <- physeq_10_filtered_sample_Data[,grouping_column]
    response <- as.factor(physeq_10_filtered_sample_Data$Groups)
    rf.data <- data.frame(response,predictors)
    
    # run random forest
    set.seed(123)
    physeq_10_filtered_Classify <- randomForest(response~., data = rf.data, ntree = 1000)
    print(physeq_10_filtered_Classify)
    
    imp <- importance(physeq_10_filtered_Classify)
    imp <- data.frame(predictors = rownames(imp), imp)
    attach(imp)

    # sort by mpg
    imp.sort <- imp[order(-MeanDecreaseGini),] 
    #imp.sort <- arrange(imp, desc(MeanDecreaseGini))
    imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

    imp.20 <- imp.sort[1:10, ]
    
    p=ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
      geom_bar(stat = "identity", fill = "indianred") + 
      coord_flip() +
      ggtitle("Most important OTUs for classifying ")
    print (p)
    
    #-----What are those OTUs?
    otunames <- imp.20$predictors
    otunames=gsub('X', '', otunames)
    r <- rownames(tax_table(physeq_10_filtered)) %in% otunames
    table=kable(tax_table(physeq_10_filtered)[r, ])
    print(table)

    #subset 10 most important variables
    lsvs=c(otunames)

    #generate boxplot===========================================
      
    #convert physeq to dataframe
    physeq_10_filtered = transform_sample_counts(physeq_10_filtered, function(x) x/sum(x))
    mphyseq = psmelt(physeq_10_filtered)
    mphyseq<- mphyseq[mphyseq$OTU %in% lsvs,]

    #plot
    options(repr.plot.width=20, repr.plot.height=15)
    prf=ggplot(data = mphyseq,  mapping = aes_string(x = grouping_column, y = "Abundance",
                                 color = Color, fill = Color)) +
        geom_boxplot(fill = NA) +
        geom_point(size = 1.7, alpha = 1,
                    position = position_jitter(width = 0.3)) +
       facet_wrap(facets = Facet) + ylab('Relative Abundance') +
       scale_y_log10()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    print (prf)

}
