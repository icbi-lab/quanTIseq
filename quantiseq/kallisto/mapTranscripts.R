mapGenes<-function(data) {
  
  # Map transcripts ids to Gene Symbols
  annotfile<-"/opt/quantiseq/kallisto/FullAnnot.txt"
  
  entrezAnnot<-read.csv(annotfile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  ind<-match(data$Feature, entrezAnnot$Transcript)
  data.gene<-entrezAnnot$Gene[ind]
  noNA<-which(!is.na(data.gene))
  data<-data[noNA,]        
  data.gene<-data.gene[noNA]
  dataout<-as.data.frame(matrix(NA, ncol=ncol(data), nrow=length(unique(data.gene))))
  colnames(dataout)<-c("GENE", colnames(data)[-1])
  dataout$GENE<-unique(data.gene)
  for (i in 1:nrow(dataout)) {
    
    ind<-which(data.gene==dataout$GENE[i])
    cdata<-as.matrix(data[ind,-1],nrow=length(ind))
    if (length(ind)>1) {
      cdata<-apply(cdata,2,sum)
    }
    dataout[i,-1]<-cdata
  }
  return(dataout)
}


createMat<-function(cdata,files) {
  #add the filenames twice once for the tpm mapping and once for the transcript cnt mapping
  numfiles <- length(files)
  colnameFromFile <- vector(mode="character", length=2*numfiles)
  for (i in 1:numfiles) {
    filename <- tools::file_path_sans_ext(files[i])
    colnameFromFile[i]<-filename
    colnameFromFile[i+numfiles]<-filename
  }
  
  data<-as.data.frame(matrix(NA, ncol=2*numfiles+1, nrow=nrow(cdata)))
  colnames(data)<-c("Feature", colnameFromFile)
  data$Feature<-cdata$target_id
  
  return(data)
}

fileNameWithPrefix<-function(filepath, prefix) {
 
  return(file.path(dirname(filepath), paste0(basename(filepath), "_", prefix, ".txt")))

}

writeData<-function(dataout, outputfile, subjectcnt) {
  #remove the unneeded cols from the result before storing it to a file
  write.table(dataout[,-((subjectcnt+2):(2*subjectcnt+1)),drop=FALSE], file = fileNameWithPrefix(outputfile,"gene_tpm"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
  dataout[,((subjectcnt+2):(2*subjectcnt+1))] = round(dataout[,((subjectcnt+2):(2*subjectcnt+1))])
  write.table(dataout[,-(2:(subjectcnt+1)),drop=FALSE], file = fileNameWithPrefix(outputfile,"gene_count"), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE) 
}

mapTPM<-function(inputfile, outputfile) {
   cdata<-read.table(inputfile, header=TRUE, stringsAsFactors=FALSE)
   data<-createMat(cdata,c(basename(inputfile)))
   data[,2]<-cdata$tpm[match(data$Feature, cdata$target_id)]
   data[,3]<-cdata$est_counts[match(data$Feature, cdata$target_id)]
   dataout<-mapGenes(data)
   writeData(dataout, outputfile, 1)
}

makeTPMtable<-function(inputdir, outputfile) {
  
  files<-list.files(inputdir)
  init<-FALSE
  for (i in 1:length(files)) {
    
    file<-file.path(inputdir, files[i])
    cdata<-read.table(file, header=TRUE, stringsAsFactors=FALSE)
    
    #just consider files with correct headers
    if(sum(colnames(cdata) == c("target_id","length","eff_length","est_counts","tpm")) != 5) next

    if (!init) {        
      data<-createMat(cdata,files)
      init<-TRUE
    }
    data[,i+1]<-cdata$tpm[match(data$Feature, cdata$target_id)]
    data[,i+length(files)+1]<-cdata$est_counts[match(data$Feature, cdata$target_id)]
    
  }
  dataout<-mapGenes(data)
  writeData(dataout, outputfile, length(files))
}

args <- commandArgs(TRUE)

input <- args[1]
output <- args[2]

if(dir.exists(input)) {
  makeTPMtable(input, output)
} else {
  mapTPM(input, output)
}