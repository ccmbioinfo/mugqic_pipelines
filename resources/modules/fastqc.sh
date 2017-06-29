#!/bin/bash

###################
################### FASTQC
###################
VERSION="0.11.2"
INSTALL_PATH=/hpf/projects/brudno/daniel/mugqic_pipelines/software/fastqc/fastqc_v"$VERSION" # where to install..
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v"$VERSION".zip
unzip fastqc_v"$VERSION".zip
rm fastqc_v"$VERSION".zip
chmod +x FastQC/fastqc
chmod -R g+w $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - FASTQC \"
}
module-whatis \"MUGQIC -FASTQC \"
                      
set             root                \$::env /hpf/projects/brudno/daniel/mugqic_pipelines/software/fastqc/fastqc_v"$VERSION"/FastQC
prepend-path    PATH               \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p /hpf/projects/brudno/daniel/mugqic_pipelines/modulefiles/mugqic/fastqc
mv .version $VERSION /hpf/projects/brudno/daniel/mugqic_pipelines/modulefiles/mugqic/fastqc/


