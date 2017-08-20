# Combat Analysis 
# This analysis takes two inputs having data from microaray o rnaseq or any expression dataset with different batches.
# As two datasets ae coming from different batches, PCA analysis is also carried to check any outlier in the dataset.
# Follwod by need to remove the batch effects using SVA package after combining the two datasets(with equal nuber of gene names should match)

rm(list=ls()) 

combat_analysis <- function(file1,file2,outFileName){
require(sva)
############## Creating Folder structure

wdDir <- getwd()
outDir <- paste0(wdDir,"/combat_output/")
dir.create(outDir, showWarnings = FALSE)


### Reading two different batch effect files (e.g from microarray and  or  )

processFile1 = file1
filename1 = unlist(strsplit(processFile1,".txt"))
process1 = read.delim(processFile1,header = TRUE,row.names = 1);dim(process1)

processFile2 = file2
filename2 = unlist(strsplit(processFile2,"txt"))
process2 = read.delim(processFile2,header = TRUE,row.names = 1);dim(process2)

### PCA plot of file1 
pdf(paste0(outDir,Sys.Date(),'_PCA_',filename1,'.pdf'), width = 12, height = 9)
plot(prcomp(t(apply(process1,2,as.numeric)),scale=TRUE)$x[,1:2],type = 'n')
text(prcomp(t(apply(process1,2,as.numeric)),scale=TRUE)$x[,1:2],colnames(process1))
dev.off()

### PCA plot of file2
pdf(paste0(outDir,Sys.Date(),'_PCA_',filename2,'.pdf'), width = 12, height = 9)
plot(prcomp(t(apply(process2,2,as.numeric)),scale=TRUE)$x[,1:2],type = 'n')
text(prcomp(t(apply(process2,2,as.numeric)),scale=TRUE)$x[,1:2],colnames(process2))
dev.off()

### Combine both the dataset base on matching number of gene names and do the combat anlysis
m<-match(toupper(rownames(process1)),toupper(rownames(process2))) # Match gene name from both file and check dimension of the data
w<-which(!is.na(m))
process1m<-process1[w,]; dim(process1m)
process2m<-process2[m[w],]; dim(process2m)
Data <-cbind(process1m,process2m);dim(Data)
sampleNo = dim(Data)[2]
combined_data_ = paste0(outDir,Sys.Date(),"_combined_data_",outFileName,"_",sampleNo,".txt")
write.table(Data,combined_data_,sep="\t",quote=FALSE,row.names = FALSE)


# Combat analysis

combine_data<-cbind(process1m,process2m[,1:dim(process2m)[2]])
dim(combine_data)
bat<-rep(c(1,2),c((dim(process1m)[2]),(dim(process2m)[2])))
length(bat)
combat<-ComBat(dat=combine_data,batch=bat,mod=NULL, par.prior=TRUE, prior.plots=FALSE)

rownames(combat)<-rownames(combine_data)
combat_genes=cbind(rownames(combat),combat)
colnames(combat_genes) <- c("Genes",colnames(combat))
combat_out = paste0(outDir,Sys.Date(),"_combat_combined_data_",outFileName,"_",sampleNo,".txt")
write.table(combat_genes,combat_out, sep="\t",quote=FALSE,row.names = FALSE)

pdf(paste0(outDir,Sys.Date(),'_PCA_after_combat_combine_data_',outFileName,"_",sampleNo,'.pdf'), width = 12, height = 9)
cl = rep(1:ncol(combat))
plot(prcomp(t(apply(combat,2,as.numeric)),scale=TRUE)$x[,1:2],col=cl,type = 'n')
text(prcomp(t(apply(combat,2,as.numeric)),scale=TRUE)$x[,1:2],colnames(combat),col=cl)
dev.off()


}