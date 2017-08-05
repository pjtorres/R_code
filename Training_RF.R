#How to train a random forest model using the data from the PArtie paper PARTIE: a partition engine to separate metagenomic and amplicon projects in the Sequence Read Archive. PJ Torres, RA Edwards, KA McNair - Bioinformatics, 2017
Partie.train=read.csv("SRA_used_for_training.csv", header=T)
#SRA_Annotation below can be replaced with whatever category you would like classififed in your data. 
Partie.train$SRA_Annotation=factor(Partie.train$SRA_Annotation)

#randomforest
library(randomForest)
set.seed(20)
xlearn=Partie.train[ ,2:5]
#the 132827 is the number of samples you want to use to trin your dataset. We had a large SRA dataset, but you may not have as big of a dataset
#this is a good tutorial on trining datasets http://machinelearningmastery.com/tune-machine-learning-algorithms-in-r/
i=sample(nrow(Partie.train),132827, rep=F)
ylearn=Partie.train[,6]
xtest=xlearn[i,]
xtrain=xlearn[-i,]
ytest=ylearn[i]
ytrain=ylearn[-i]

#RF of training dataset
rf=randomForest(xtrain, ytrain, xtest, ytest,importance=T, keep.forest=TRUE)

# save the model to disk
saveRDS(rf, "./final_model.rds")
