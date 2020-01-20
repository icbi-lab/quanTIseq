#!/bin/sh
set -eo pipefail
#-------------------------
# quanTIseq pipeline 
#-------------------------

# read the default settings
set -o allexport
source /opt/quantiseq/defaults.conf
set +o allexport

# overwrite with user specified options
for option in "$@"
do
    if [[ $option != --*=* ]]; then
      echo "ERROR: $option is not a recognized option!" >&2
      exit 1
    fi
    param=`echo $option | sed 's/=.*//' | sed 's/-//g'`
    value=`echo $option | sed 's/.*=//'`

    case "$param" in
        inputfile|outputdir|totalcells|rmgenes|prefix|threads|pipelinestart|phred|adapterSeed|palindromeClip|simpleClip|trimLead|trimTrail|minLen|crop|rawcounts|arrays|tumor|mRNAscale|method|avgFragLen|sdFragLen|help|btotalcells)
          export $param=$value
        ;;
        *)
          echo "ERROR: $param is not a valid option!" >&2
          exit 1
        ;;
    esac
done

#-------------------------------#
### check optional parameters:
#-------------------------------#


# starting point of the pipeline (default:preproc)
if [ "$pipelinestart" != "preproc" ] && [ "$pipelinestart" != "expr" ] && [ "$pipelinestart" != "decon" ] && [ "$pipelinestart" != "" ]; then
    echo "ERROR: parameter --pipelinestart ($pipelinestart) is not correct. Choose one of the options: preproc, expr, decon!" >&2
    exit 1
fi

# check if number of threads is an interger
if [[ ! "$threads" =~ ^[+]?[0-9]+$ ]] && [ "$threads" != "" ]; then
    echo "ERROR: --threads ($threads) must be an integer!" >&2
    exit 1
fi

# check if number of min length is an interger
if [[ ! "$minLen" =~ ^[+]?[0-9]+$ ]] && [ "$minLen" != "" ]; then
    echo "ERROR: --minLen ($minLen) must be an integer!" >&2
    exit 1
fi

# check if number of fragment length is an interger
if [[ ! "$avgFragLen" =~ ^[+]?[0-9]+$ ]] && [ "$avgFragLen" != "" ]; then
    echo "ERROR: --avgFragLen ($avgFragLen) must be an integer!" >&2
    exit 1
fi


# check if number of the standard deviation of the fragment length is an interger
if [[ ! "$sdFragLen" =~ ^[+]?[0-9]+$ ]] && [ "$sdFragLen" != "" ]; then
    echo "ERROR: --sdFragLen ($sdFragLen) must be an integer!" >&2
    exit 1
fi

# deconvolution method (default:lsei)
if [ "$method" != "lsei" ] && [ "$method" != "hampel" ] && [ "$method" != "huber" ] && [ "$method" != "bisquare" ] && [ "$method" != "" ]; then
    echo "ERROR: parameter --method ($method) is not correct. Choose one of the options: lsei, hampel, huber, bisquare!" >&2
    exit 1
fi

# check if parameters rawcounts, arrays, tumor and mRNAscale are either TRUE or FALSE
if [ "$rawcounts" != "TRUE" ] && [ "$rawcounts" != "FALSE" ] && [ "$rawcounts" != "" ]; then
    echo "ERROR: --rawcounts ($rawcounts) must be TRUE or FALSE!" >&2
    exit 1
fi

if [ "$arrays" != "TRUE" ] && [ "$arrays" != "FALSE" ] && [ "$arrays" != "" ]; then
    echo "ERROR: --arrays ($arrays) must be TRUE or FALSE!" >&2
    exit 1
fi

if [ "$tumor" != "TRUE" ] && [ "$tumor" != "FALSE" ] && [ "$tumor" != "" ]; then
    echo "ERROR: --tumor ($tumor) must be TRUE or FALSE!" >&2
    exit 1
fi

if [ "$mRNAscale" != "TRUE" ] && [ "$mRNAscale" != "FALSE" ] && [ "$mRNAscale" != "" ]; then
    echo "ERROR: --mRNAscale ($mRNAscale) must be TRUE or FALSE!" >&2
    exit 1
fi


#-------------------------------#
### run pipeline:
#-------------------------------#

# set initial pipelinesteps to FALSE:
preproc="FALSE"
expr="FALSE"
decon="FALSE"

if [ $help == "TRUE" ]; then
    cat /opt/quantiseq/quanTIseq_usage.txt >&2
    exit
fi

# check if inputfile exists:
if [ -f /opt/quantiseq/Input/inputfile.txt ]; then
  inputfile=/opt/quantiseq/Input/inputfile.txt
else
  echo "ERROR: please specify input files as described in quanTIseq documentation. To inspect this image call: docker run -it --entrypoint=/bin/bash icbi/quantiseq" >&2
  exit
fi

# Get number of available threads from the system and take the minimum of user input and available threads:
maxthreads=$(grep -c ^processor /proc/cpuinfo)
nthreads=$(($maxthreads<$threads?$maxthreads:$threads))
 
if [ "$nthreads" != "1" ]; then
    echo
    echo "WARNING: running kallisto in multi-thread mode can produce slightly different results in terms of gene counts and TPMs."
    echo
fi

#-----------------------------------
# Preprocessing / Trimmomatic
#-----------------------------------

preprocpath="/opt/quantiseq/user_output/out_preproc/"
mkdir -p $preprocpath

if [ $pipelinestart == "preproc" ]; then

  Rscript /opt/quantiseq/trimmomatic/quanTIseq_preproc.R $inputfile $preprocpath $nthreads $phred $adapters $adapterSeed $palindromeClip $simpleClip $trimLead $trimTrail $minLen $crop

  rm -f ${preprocpath}*nopair

  pipelinestart="expr"
  preproc="TRUE"
  arrays="FALSE"
fi

#--------------------------------------
# Transcript expression / Kallisto
#--------------------------------------

exprpath="/opt/quantiseq/user_output/out_expr/"
mkdir -p $exprpath

if [ $pipelinestart == "expr" ]; then

  # Run kallisto:
  Rscript /opt/quantiseq/kallisto/quanTIseq_expr.R $inputfile $exprpath $nthreads $preproc $preprocpath $avgFragLen $sdFragLen
  # Map transcripts to human gene symbols:
  Rscript /opt/quantiseq/kallisto/mapTranscripts.R $exprpath "${exprpath}${prefix}"
  rm ${exprpath}*.tsv

  # Remove rawcount file if not needed:
  if [ $rawcounts == "FALSE" ]; then
    rm ${exprpath}*gene_count.txt
  fi

  pipelinestart="decon"
  expr="TRUE"
  arrays="FALSE"
  inputfile="${exprpath}${prefix}_gene_tpm.txt"

  cp -r ${exprpath}. /opt/quantiseq/user_output/

fi

#-------------------------------------
# Deconvolution / R
#-------------------------------------

deconpath="/opt/quantiseq/user_output/out_decon/"
mkdir -p $deconpath

if [ $pipelinestart == "decon" ]; then
  # Run deconvolution:
  Rscript /opt/quantiseq/deconvolution/quanTIseq_decon.R $inputfile $deconpath $expr $arrays $signame $tumor $mRNAscale $method $prefix $btotalcells $rmgenes
  decon="TRUE"
  cp "${deconpath}${prefix}_cell_fractions.txt" /opt/quantiseq/user_output/

  if [ $btotalcells == "TRUE" ]; then
    cp "${deconpath}${prefix}_cell_densities.txt" /opt/quantiseq/user_output/
  fi
fi

rm -fr $preprocpath
rm -fr $exprpath
rm -fr $deconpath
