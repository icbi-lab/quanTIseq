FROM centos:7

LABEL maintainer="christina.plattner@i-med.ac.at,francesca.finotello@i-med.ac.at,dietmar.rieder@i-med.ac.at"

RUN yum clean all && yum -y update

################# BEGIN INSTALLATION ######################


################## Java
RUN yum install -y java

##################  R (> version 3.4.3)
RUN yum install -y epel-release
RUN yum install -y R-core R-devel

#### utils

RUN yum install -y dos2unix
RUN yum install -y mc


ADD dependencies.R /tmp/
RUN Rscript /tmp/dependencies.R

################## main

ADD quantiseq /opt/quantiseq

##################### INSTALLATION END #####################

#Clean up
RUN rm -rf /tmp/* /var/tmp/* ~/.cache/*

# entrypoint
ENTRYPOINT ["/opt/quantiseq/quanTIseq.sh"]
