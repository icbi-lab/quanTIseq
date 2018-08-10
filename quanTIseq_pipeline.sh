#$ -S /bin/sh

## Copyright (c) 2017, Division of Bioinformatics, Innsbruck Medical University
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
## 
## 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
## 
## 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
## 
## 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##

realpath2 ()
{

        relpath=$1
        cpath=`pwd`
        reldir=`dirname $relpath`
        relfile=`basename $relpath`
        cd $reldir
        abspath=`pwd`
        abspath=$abspath"/"$relfile
        cd $cpath
        echo $abspath

}

docker_vars=""
docker_run="docker run --rm"


#---------------------------#
### read input parameters:
#---------------------------#

for option in "$@"
do

        if [[ $option == "--help" ]]; then
            docker run --rm -e help="TRUE" icbi/quantiseq
            exit
        fi
        
	if [[ $option != --*=* ]]; then
		echo "ERROR: $option is not a recognized option!" >&2
		exit 1
	fi

	param=`echo $option | sed 's/=.*//' | sed 's/-//g'`
	value=`echo $option | sed 's/.*=//'`
	
	case "$param" in
	  inputfile|outputdir|totalcells|rmgenes)
	    eval $param=$value
	    ;;
	  prefix|threads|pipelinestart|phred|adapterSeed|palindromeClip|simpleClip|trimLead|trimTrail|minLen|crop|rawcounts|arrays|tumor|mRNAscale|method|avgFragLen|sdFragLen)
	    eval $param=$value
	    docker_vars=""$docker_vars" -e "$param"="$value
	    ;;
	  *)
	    echo "ERROR: $param is not a valid option!" >&2
	    exit 1
	    ;;
	esac
 
done


#---------------------------#
### check parameters:
#---------------------------#

# Optional parameters:

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


# output directory (default: current directory)
if [ "$outputdir" == "" ]; then
  outputdir="$(pwd)/"
elif [ $outputdir == "." ]; then
  outputdir="$(pwd)/"
else
    
    if [ ! -d $outputdir ]; then  
        echo "ERROR: Outputdir ($outputdir) does not exist!" >&2
	exit 1 
    fi
    
    if [ ! -w $outputdir ]; then  
        echo "ERROR: Outputdir ($outputdir) is not writeable!" >&2
	exit 1 
    fi
    
    outputdir=$(realpath2 $outputdir)

fi


# check if file for cell densities exists and is readable (totalcells)
if [ "${totalcells}" != "" ]; then
  
    if [ ! -f "${totalcells}" ]; then
    echo "ERROR: Input file ($totalcells) does not exist!" >&2
    exit 1
    fi
    
    if [ ! -r "${totalcells}" ]; then
    echo "ERROR: Input file ($totalcells) is not readable!" >&2
    exit 1
    fi 
    
    totalcells=" -v "$(realpath2 $totalcells)":/home/deconvolution/totalcells.txt -e btotalcells=TRUE"
fi

# remove genes options
if [ "${rmgenes}" != "" ] && [ "${rmgenes}" != "none" ] && [ "${rmgenes}" != "default" ]; then
  
  if [ ! -f "${rmgenes}" ]; then
    echo "ERROR: Input file ($rmgenes) does not exist!" >&2
    exit 1
    fi
    
    if [ ! -r "${rmgenes}" ]; then
    echo "ERROR: Input file ($rmgenes) is not readable!" >&2
    exit 1
    fi 
    
    rmgenesfile=" -v "$(realpath2 $rmgenes)":/home/deconvolution/rmgenes.txt"
    rmgenes="path"
else
    rmgenesfile=""
fi

if [ "${rmgenes}" == "" ]; then
  rmgenes="unassigned"
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

# Mandatory parameter:

if [ "$inputfile" == "" ]; then
  echo "ERROR: parameter --inputfile is mandatory!" >&2
  exit 1
fi
# check if inputfile exists and is readable
if [ ! -f "${inputfile}" ]; then
    echo "ERROR: Input file ($inputfile) does not exist!" >&2
    exit 1
fi
if [ ! -r "${inputfile}" ]; then
    echo "ERROR: Input file ($inputfile) is not readable!" >&2
    exit 1
fi

if [ "$pipelinestart" != "decon" ]; then

  docker_run=" "$docker_run" -v "$(realpath2 $inputfile)":/home/Input/inputfile.txt"
  
  # check if rnaSeq-files exist:
  inputfiles=(`cut -f 2,3 $inputfile | tr "\t" "\n" | grep -v None`)
    
    for file in ${inputfiles[@]}
    do
        if [ ! -f "${file}" ]; then
             echo "$file not found! Check input files!" >&2
             exit 1
        
        else
	    filename=$(basename $file)
	    docker_run=" "$docker_run" -v "$(realpath2 $file)":/home/Input/"$filename
        fi
    done
else
    docker_run=" "$docker_run" -v "$(realpath2 $inputfile)":/home/Input/inputfile.txt"
fi

#---------------------------#
### run Docker:
#---------------------------#

# define volumes:
docker_run=""$docker_run" -v "$outputdir":/home/user_output/ "$rmgenesfile" "$totalcells

# define parameters as environment variables:
docker_run=$docker_run" "$docker_vars" -e rmgenes="$rmgenes" icbi/quantiseq"

$docker_run
