args <- commandArgs(TRUE)

# parameters:
mix.mat.file <- args[1]
output <- args[2]
expr <- args[3]
arrays <- args[4]
signame <- args[5]
tumor <- args[6]
mRNAscale <- args[7]
method <- args[8]
prefix <- args[9]
btotalcells <- args[10]
rmgenes <- args[11]

library(preprocessCore)
source("/home/deconvolution/DOCKER_codes.R")

message("\nRunning quanTIseq deconvolution module\n")

# List of genes to be discarded
if (rmgenes=="unassigned" && arrays==TRUE) { # For Microarrays
  rmgenes<-"none"
  
} else if (rmgenes=="unassigned" && arrays==FALSE) { # For RNA-seq
  rmgenes<-"default"
  
}


# Files
listsig<-c("TIL10")
if (signame %in% listsig) {
  
  sig.mat.file<-file.path("/home/deconvolution/",
                          paste0(signame, "_signature.txt"))
  
  mRNA.file<-file.path("/home/deconvolution/",
                       paste0(signame, "_mRNA_scaling.txt"))
  
  fileab<-file.path("/home/deconvolution/",
    paste0(signame,"_TCGA_aberrant_immune_genes.txt"))
  
  if (rmgenes=="default") {
    filerm<-file.path("/home/deconvolution/",
                      paste0(signame,"_rmgenes.txt"))
    
  } else if (rmgenes=="path") {
    filerm<-file.path("/home/deconvolution/",
                      paste0(signame,"rmgenes.txt"))
    
  }

  
} else {
  
  sig.mat.file<-paste0(signame, "_signature.txt")
  mRNA.file<-paste0(signame, "_mRNA_scaling.txt")
 
}

# Load mixture matrix
mix.mat<-read.table(mix.mat.file, header=TRUE, sep="\t", row.names=1)

if(is.numeric(mix.mat[1,1])!= TRUE){
  stop("Wrong input format for the mixture matrix! Please follow the instructions of the documentation.")
}

# Load signature 
sig.mat<-read.table(sig.mat.file, header=TRUE, sep="\t", row.names=1)

# Load normalization factors (set all to 1 if mRNAscale==FALSE)
if (mRNAscale) {

  mRNA<-read.table(mRNA.file, 
                   sep="\t", 
                   header=FALSE, 
                   stringsAsFactors=FALSE)
  colnames(mRNA)<-c("celltype", "scaling")
  mRNA<-as.vector(as.matrix(mRNA$scaling[match(colnames(sig.mat), mRNA$celltype)]))
  
} else {
  
  mRNA<-rep(1, ncol(sig.mat))
  
}

# Preprocess mixture matrix
message(paste0("Gene expression normalization and re-annotation (arrays: ", 
  arrays, ")\n"))
mix.mat<-fixMixture(mix.mat, arrays=arrays)

# Remove noisy genes 
if (rmgenes!="none") {
  
  if (signame %in% listsig) {
    
    lrmgenes<-as.vector(read.table(filerm, header=FALSE, sep="\t")[,1])
    n1<-nrow(sig.mat)
    sig.mat<-sig.mat[!rownames(sig.mat) %in% lrmgenes,, drop=FALSE]
    n2<-nrow(sig.mat)
    message(paste0("Removing ", n1-n2, " noisy genes\n"))
    
  }
}

# Fix tumor data
if (tumor) {
  
  if (signame %in% listsig) {

  abgenes<-as.vector(read.table(fileab, header=FALSE, sep="\t")[,1])
  n1<-nrow(sig.mat)
  sig.mat<-sig.mat[!rownames(sig.mat) %in% abgenes,, drop=FALSE]
  n2<-nrow(sig.mat)
  message(paste0("Removing ", n1-n2, " genes with high expression in tumors\n"))
  
  }
}

# Signature genes present in the mixture
ns<-nrow(sig.mat)
us<-length(intersect(rownames(sig.mat), rownames(mix.mat)))
perc<-round(us*100/ns,digits=2)
message(paste0("Signature genes found in data set: ", 
  us, "/", ns, " (", perc, "%)\n"))

# Run deconvolution
message(paste0("Mixture deconvolution (method: ", method, ")\n"))
results1<-quanTIseq(sig.mat,
                   mix.mat,
                   scaling=mRNA,
                   method=method)
if ("Tregs" %in% colnames(sig.mat) && "T.cells.CD4" %in% colnames(sig.mat) && method %in% c("lsei")) {

  minTregs<-0.02
  i<-which(colnames(sig.mat)=="T.cells.CD4")
  results2<-quanTIseq(sig.mat[,-i],
    mix.mat,
    scaling=mRNA[-i],
    method=method)

  ind<-which(results1[,"Tregs"]<minTregs)
  if (length(ind)>0) {

    results1[ind,"Tregs"]<-(results2[ind,"Tregs"]+results1[ind,"Tregs"])/2
    results1[ind,"T.cells.CD4"]<-pmax(0,results1[ind,"T.cells.CD4"]-(results2[ind,"Tregs"]+results1[ind,"Tregs"])/2)

  }

}
results<-results1
results<-results/apply(results,1,sum)

# Save results using user's output ID
fileout<-paste0(output, prefix, "_cell_fractions.txt")
DCres <- results
results<-cbind(rownames(results), results)
colnames(results)[1]<-"Sample"
write.table(results, 
            sep="\t",
            row.names=FALSE,
            quote=FALSE,
            file=fileout)

if (btotalcells == TRUE){
  
  celldens <- celldensities(DCres)
  fileout2 <- paste0(output, prefix, "_cell_densities.txt")
  celldens<-cbind(rownames(celldens), celldens)
  colnames(celldens)[1]<-"Sample"
  write.table(celldens,
              sep="\t",
              row.names=FALSE,
              quote=FALSE,
              file=fileout2)
}

message("Deconvolution results saved to file!")

