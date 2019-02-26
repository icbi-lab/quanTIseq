transcriptexpression <- function(inputfile, outputpath, threads, preproc, preprocpath, avgFragLen, sdFragLen){
  
  input <- as.matrix(read.table(inputfile))
  
   
  for (i in 1:nrow(input)){
    
    # check if single or paired end:
    if (input[i,3] == "None"){
      libtype = paste0("--single -l ", avgFragLen," -s ", sdFragLen)
      
      # check if preprocessing with trimmomatic was performed:
      if (preproc){
        files = paste0(preprocpath, basename(input[i,2]))
      }
      else files = paste0("/opt/quantiseq/Input/", basename(input[i,2]))
    }
    
    else{
      libtype = ""
      
      #check if preprocessing with trimmomatic was performed:
      if (preproc){
        files = paste0(preprocpath, basename(input[i,2]), " ", preprocpath, basename(input[i,3]))
      }
      else files = paste0("/opt/quantiseq/Input/", basename(input[i,2]), " /opt/quantiseq/Input/", basename(input[i,3]))
    }
    
    opath = paste0(outputpath, input[i,1])
    # Run kallisto:
    system(paste("/opt/quantiseq/kallisto/kallisto quant -i /opt/quantiseq/kallisto/hg19_M_rCRS_kallisto.idx -o ", opath," ", libtype, "-t ", threads," ", files))
    
    # move tsv-file to output directory and delete the other directories:
    system(paste0("mv ", opath, "/abundance.tsv ", outputpath, "/", input[i,1], ".tsv"))
    system(paste("rm -r", opath))
  }
  
}

args <- commandArgs(TRUE)

inputfile <- args[1]
outputpath <- args[2]
threads <- args[3]
preproc <- args[4]
preprocpath <- args[5]
avgFragLen <- args[6]
sdFragLen <- args[7]


transcriptexpression(inputfile, outputpath, threads, preproc, preprocpath, avgFragLen, sdFragLen)
