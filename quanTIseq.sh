#$ -S /bin/sh

#-------------------------
# quanTIseq pipeline 
#-------------------------

# set initial pipelinesteps to FALSE:
preproc="FALSE"
expr="FALSE"
decon="FALSE"

if [ $help == "TRUE" ]; then
    cat /home/quanTIseq_usage.txt >&2
    exit
fi

# check if inputfile exists:
if [ -f /home/Input/inputfile.txt ]; then
  inputfile=/home/Input/inputfile.txt
else
  echo "ERROR: please specify input files as described in quanTIseq documentation. To inspect this image call: docker run -it --entrypoint=/bin/bash icbi/quantiseq" >&2
  exit
fi

# Get number of available threads from the system and take the minimum of user input and available threads:
maxthreads=$(grep -c ^processor /proc/cpuinfo)
nthreads=$(($maxthreads<$threads?$maxthreads:$threads))

#-----------------------------------
# Preprocessing / Trimmomatic
#-----------------------------------

preprocpath="/home/Output/out_preproc/"

if [ $pipelinestart == "preproc" ]; then

  Rscript /home/trimmomatic/quanTIseq_preproc.R $inputfile $preprocpath $nthreads $phred $adapters $adapterSeed $palindromeClip $simpleClip $trimLead $trimTrail $minlen $crop
  
  rm -f ${preprocpath}*nopair
     
  pipelinestart="expr"
  preproc="TRUE"
  arrays="FALSE"
fi

#--------------------------------------
# Transcript expression / Kallisto
#--------------------------------------

exprpath="/home/Output/out_expr/"

if [ $pipelinestart == "expr" ]; then

  # Run kallisto:
  Rscript /home/kallisto/quanTIseq_expr.R $inputfile $exprpath $nthreads $preproc $preprocpath
  # Map transcripts to human gene symbols:
  Rscript /home/kallisto/mapTranscripts.R $exprpath "${exprpath}${prefix}"
  rm ${exprpath}*.tsv
  
  # Remove rawcount file if not needed:
  if [ $rawcounts == "FALSE" ]; then
    rm ${exprpath}*gene_count.txt
  fi
  
  pipelinestart="decon"
  expr="TRUE"
  arrays="FALSE"
  inputfile="${exprpath}${prefix}_gene_tpm.txt"
  
  cp -r ${exprpath}. /home/user_output/
  chmod a+rwx -R /home/user_output
  
fi

#-------------------------------------
# Deconvolution / R
#-------------------------------------

deconpath="/home/Output/out_decon/"

if [ $pipelinestart == "decon" ]; then
  # Run deconvolution:
  Rscript /home/deconvolution/quanTIseq_decon.R $inputfile $deconpath $expr $arrays $signame $tumor $mRNAscale $method $prefix $btotalcells $rmgenes
  decon="TRUE"
  cp "${deconpath}${prefix}_cell_fractions.txt" /home/user_output/
  
  if [ $btotalcells == "TRUE" ]; then
  cp "${deconpath}${prefix}_cell_densities.txt" /home/user_output/
  
  fi
  
fi
