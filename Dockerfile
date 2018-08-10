FROM centos
 
LABEL maintainer="christina.plattner@i-med.ac.at,francesca.finotello@i-med.ac.at"

RUN yum clean all && yum -y update

################# BEGIN INSTALLATION ######################


################## TRIMMOMATIC #Requires Java
RUN yum install -y java

ADD trimmomatic/trimmomatic-0.36.jar /usr/local/bin/trimmomatic
ADD trimmomatic/ /home/trimmomatic/

##################  R 
RUN yum install -y epel-release
RUN yum install -y R-core R-devel
 
ADD dependencies.R /tmp/
RUN Rscript /tmp/dependencies.R

ADD deconvolution/ /home/deconvolution/

################## KALLISTO #Requires R

ADD kallisto/ /home/kallisto/

################## main
ADD quanTIseq.sh /home/
ADD quanTIseq_usage.txt /home/
ADD LICENSE.txt /home/
ADD Output/ /home/Output/
 
##################### INSTALLATION END #####################



####### DEFAULT PARAMETERS as ENVIRONMENT VARIABLES ########

# TRIMMOMATIC:
ENV phred=33 adapters="home/trimmomatic/TruSeqAdapt.fa" adapterSeed=2 palindromeClip=30 simpleClip=10 trimLead=20 trimTrail=20 minlen=36 crop=10000

# KALLISTO:
ENV rawcounts="FALSE" avgFragLen=50 sdFragLen=20

# DECONVOLUTION:
ENV arrays="FALSE" signame="TIL10" tumor="FALSE" mRNAscale="TRUE" method="lsei" btotalcells="FALSE" rmgenes="unassigned"

# GENERAL:
ENV threads=1 prefix="quanTIseq" pipelinestart="preproc" help=FALSE

# Run main shell script when the container launches
CMD sh /home/quanTIseq.sh --threads=$threads --prefix=$prefix --pipelinestart=$pipelinestart --phred=$phred --adapters=$adapters --adapterSeed=$adapterSeed --palindromeClip=$palindromeClip --simpleClip=$simpleClip --trimLead=$trimLead --trimTrail=$trimTrail --minlen=$minlen --crop=$crop --rawcounts=$rawcounts --avgFragLen=$avgFragLen --sdFragLen=$sdFragLen --arrays=$arrays --signame=$signame --tumor=$tumor --mRNAscale=$mRNAscale --method=$method --btotalcells=$btotalcells --rmgenes=$rmgenes

#Clean up
RUN rm -rf /tmp/* /var/tmp/* ~/.cache/*
