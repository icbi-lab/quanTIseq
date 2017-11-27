library(tools)

preprocessing <- function(inputfile, outputpath, threads, phred, adapters, adapterSeed, palindromeClip, simpleClip, trimLead, trimTrail, minlen, crop)
  {
  
  input <- as.matrix(read.table(inputfile))
  
  # Run Trimmomatic for each line of the input file:
    for (i in 1:nrow(input)){
    
    # check if single or paired end:
    if (input[i,3] == "None"){
      libtype = "SE"
      files = paste0("/home/Input/", basename(input[i,2]), " ", outputpath, basename(input[i,2]))
    }
    
    else{
      libtype = "PE"
      out1 = paste0(outputpath, basename(input[i,2]))
      out2 = paste0(outputpath, basename(input[i,3]))
      files = paste0("/home/Input/", basename(input[i,2])," /home/Input/", basename(input[i,3])," ", out1," ", out1, "_nopair ", out2," ", out2, "_nopair")
    }
    
    system(paste0("java -jar /home/trimmomatic/trimmomatic-0.36.jar ",libtype," -threads ",threads," -phred", phred," ",files," ILLUMINACLIP:",adapters,":",adapterSeed,
                    ":",palindromeClip,":",simpleClip," LEADING:",trimLead," TRAILING:",trimTrail," MINLEN:", minlen," CROP:", crop))
  }
}

args <- commandArgs(TRUE)

inputfile <- args[1]
outputpath <- args[2]
threads <- args[3]
phred <- args[4]
adapters <- args[5]
adapterSeed <- args[6]
palindromeClip <- args[7]
simpleClip <- args[8]
trimLead <- args[9]
trimTrail <- args[10]
minlen <- args[11]
crop <- args[12]

preprocessing(inputfile, outputpath, threads, phred, adapters, adapterSeed, palindromeClip, simpleClip, trimLead, trimTrail, minlen, crop)