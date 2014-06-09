
###################
################### STAR
###################
VERSION="2.3.0e"
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/star
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download and extract
wget https://rna-star.googlecode.com/files/STAR_$VERSION.tgz
tar zxvf STAR_$VERSION.tgz
chmod -R 775 STAR_$VERSION

##INSTALL
cd STAR_$VERSION
make
cd ..
chmod -R 775 STAR_$VERSION

# Module file
echo "#%Module1.0
proc ModulesHelp { } {STAR RNA-aligner \"
}
module-whatis \"MUGQIC - star \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME_DEV)/software/star/STAR_$VERSION
prepend-path    PATH               \$root
" > $VERSION

# Version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/star
mv .version $VERSION $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/star/
rm STAR_$VERSION.tgz
