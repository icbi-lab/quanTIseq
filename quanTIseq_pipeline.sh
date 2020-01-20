#!/bin/bash
set -eo pipefail
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

RUNUID=`id -u`
RUNGID=`id -g`
# get operating system:
UNAME=`uname` 

# check OS to define if the pipeline will be run in Docker or Singularity:
if [ $UNAME == "Linux" ] && [ -z $FORCE_DOCKER ]; then
    container="s"
    echo "Starting quanTIseq pipeline with Singularity"
else
    container="d"
    echo "Starting quanTIseq pipeline with Docker"
fi

# get realpath function, which works on Linux and Mac OS
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

case $container in
  d)
    # check if docker is installed:
    docker --version > /dev/null 2>&1
    if [ "$?" != "0" ]; then
      echo "Error: Docker is not installed!"
      exit 1
    fi
    
    docker_opts=""
    docker_run="docker run --rm --user ${RUNUID}:${RUNGID}"
    ;;
  s)
    # check if singularity is installed:
    singularity --version > /dev/null 2>&1
    if [ "$?" != "0" ]; then
      echo "Error: Singularity is not installed!"
      exit 1
    fi
    # get docker image and create new singularity image:
    if [ ! -e ./quantiseq2.img ]; then
    echo "Building quantiseq singularity image"
    singularity build quantiseq2.img docker://icbi/quantiseq2
    fi
    singularity_vars=""
    singularity_run="singularity run"
    bindMount="-B"
    ;;
esac


#---------------------------#
### read input parameters:
#---------------------------#

for option in "$@"
do

    if [ $option == "--help" ] || [ $option == "--h" ]; then
        case "$container" in
        d)
	  docker run --rm icbi/quantiseq2 --help="TRUE"
	  exit
	  ;;
	s)
	  singularity run ./quantiseq2.img --help="TRUE"
	  exit
	  ;;
	esac
    fi

    param=`echo $option | sed 's/=.*//g' | sed 's/-//g'`
    value=`echo $option | sed 's/.*=//g'`

    if [ $container == "d" ]; then
      case "$param" in
	  inputfile|outputdir|totalcells|rmgenes)
	      eval $param=$value
	      ;;
	  prefix|threads|pipelinestart|phred|adapterSeed|palindromeClip|simpleClip|trimLead|trimTrail|minLen|crop|rawcounts|arrays|tumor|mRNAscale|method|avgFragLen|sdFragLen)
	      eval $param=$value
	      docker_opts=""$docker_opts" --"$param"="$value
	      ;;
	  *)
	      echo "ERROR: \"$param\" is not a valid option! See --help for further details." >&2
	      exit 1
	      ;;
      esac
    else
      case "$param" in
	  inputfile|outputdir|totalcells|rmgenes|prefix|threads|pipelinestart|phred|adapterSeed|palindromeClip|simpleClip|trimLead|trimTrail|minLen|crop|rawcounts|arrays|tumor|mRNAscale|method|avgFragLen|sdFragLen)
	      singularity_vars="$singularity_vars $option"
	      eval $param=$value
	      ;;
	  *)
	      echo "ERROR: \"$param\" is not a valid option! See --help for further details." >&2
	      exit 1
	      ;;
      esac
    fi
    
      
done


#---------------------------#
### check input parameters:
#---------------------------#


# output directory
if [ "$outputdir" == "" ]; then
   echo "ERROR: parameter --outputdir is mandatory!" >&2
   exit 1
    
elif [ $outputdir == "." ]; then
    outputdir="$(pwd)/quantiseqResults_$$"
    mkdir -p $outputdir || {
        echo "ERROR: can not create outputdir ($outputdir)!" >&2
        exit 1
    }
else

    if [ ! -d $outputdir ]; then
        echo "WARN: Outputdir ($outputdir) does not exist, trying to create it!" >&2
        mkdir -p $outputdir || {
            echo "ERROR: can not create outputdir ($outputdir)!" >&2
            exit 1
        }
    fi

    if [ ! -w $outputdir ]; then
        echo "ERROR: Outputdir ($outputdir) is not writeable!" >&2
        exit 1
    fi

    outputdir=$(realpath2 $outputdir)

fi


# check if cell densities file exists and is readable (totalcells)
btotalcells="FALSE"
if [ "${totalcells}" != "" ]; then

    if [ ! -f "${totalcells}" ]; then
        echo "ERROR: Input file ($totalcells) does not exist!" >&2
        exit 1
    fi

    if [ ! -r "${totalcells}" ]; then
        echo "ERROR: Input file ($totalcells) is not readable!" >&2
        exit 1
    fi
    
    case "$container" in
      d)
	totalcells=" -v "$(realpath2 $totalcells)":/opt/quantiseq/deconvolution/totalcells.txt"
	btotalcells="TRUE"
	;;
      s)
	totalcells=$(realpath2 $totalcells)":/opt/quantiseq/deconvolution/totalcells.txt"
	singularity_vars="$singularity_vars  --btotalcells=TRUE"
	;;
    esac
	
fi

# check if "remove genes file" exists and is readable (rmgenes)
if [ "${rmgenes}" != "" ] && [ "${rmgenes}" != "none" ] && [ "${rmgenes}" != "default" ]; then

    if [ ! -f "${rmgenes}" ]; then
        echo "ERROR: Input file ($rmgenes) does not exist!" >&2
        exit 1
    fi

    if [ ! -r "${rmgenes}" ]; then
        echo "ERROR: Input file ($rmgenes) is not readable!" >&2
        exit 1
    fi
    
    case "$container" in
      d)
	rmgenesfile=" -v "$(realpath2 $rmgenes)":/opt/quantiseq/deconvolution/rmgenes.txt"
	;;
      s)
	rmgenesfile=$(realpath2 $rmgenes)":/opt/quantiseq/deconvolution/rmgenes.txt"
	;;
    esac
	
    rmgenes="path"
else
    rmgenesfile=""
fi

if [ "${rmgenes}" == "" ]; then
    rmgenes="unassigned"
fi

# check if inputfile exists and is readable
if [ "$inputfile" == "" ]; then
    echo "ERROR: parameter --inputfile is mandatory!" >&2
    exit 1
fi

if [ ! -f "${inputfile}" ]; then
    echo "ERROR: Input file ($inputfile) does not exist!" >&2
    exit 1
fi
if [ ! -r "${inputfile}" ]; then
    echo "ERROR: Input file ($inputfile) is not readable!" >&2
    exit 1
fi

if [ "$pipelinestart" != "decon" ]; then
    
    case "$container" in
      d)
	docker_run=" "$docker_run" -v "$(realpath2 $inputfile)":/opt/quantiseq/Input/inputfile.txt"
	;;
      s)
	bindMount="$bindMount "$(realpath2 $inputfile)":/opt/quantiseq/Input/inputfile.txt"
	;;
    esac
	
    # check if rnaSeq-files exist:
    inputfiles=(`cut -f 2,3 $inputfile | tr "\t" "\n" | grep -v None`)

    for file in ${inputfiles[@]}
    do
        if [ ! -f "${file}" ]; then
            echo "$file not found! Check input files!" >&2
            exit 1

        else
            filename=$(basename $file)
            
            case "$container" in
	      d)
		docker_run=" "$docker_run" -v "$(realpath2 $file)":/opt/quantiseq/Input/"$filename
		;;
	      s)
		bindMount="$bindMount,"$(realpath2 $file)":/opt/quantiseq/Input/"$filename
		;;
	    esac
        fi
    done
else
    case "$container" in
    d)
      docker_run=" "$docker_run" -v "$(realpath2 $inputfile)":/opt/quantiseq/Input/inputfile.txt"
      ;;
    s)
      bindMount="$bindMount "$(realpath2 $inputfile)":/opt/quantiseq/Input/inputfile.txt"
      ;;
    esac
fi


#----------------------------------#
### run Docker or Singularity:
#----------------------------------#


case "$container" in
    d)
      # define volumes and add parameters:
      docker_run=""$docker_run" -v "$outputdir":/opt/quantiseq/user_output/ "$rmgenesfile" "$totalcells" icbi/quantiseq2 --rmgenes="$rmgenes" --btotalcells="$btotalcells" "$docker_opts
      echo $docker_run
      $docker_run
      ;;
    s)
      bindMount=$bindMount","$outputdir":/opt/quantiseq/user_output/"
      if [ "$rmgenesfile" != "" ]; then
	  bindMount=$bindMount","$rmgenesfile
      fi
      if [ "$totalcells" != "" ]; then
	  bindMount=$bindMount","$totalcells
      fi
      # define parameters as environment variables:
      singularity_run="$singularity_run $bindMount ./quantiseq2.img $singularity_vars --rmgenes=$rmgenes"
      # echo $singularity_run
      $singularity_run
      ;;
esac

