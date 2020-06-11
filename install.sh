#!/bin/bash
#########################
### run ./install.sh ####
#########################

### get the working directory ######
BASEDIR=$(pwd)
#echo $BASEDIR
####################################


###### create additional directories ######
mkdir lib;
mkdir include;
mkdir bin;
echo -e "\ngenerating LocHap binary ...\n\n";
###########################################

###### create libz ######
tar -xzf zLib.tar.gz
cd zLib;
./configure > ../out.log 2>&1;
make >> ../out.log 2>&1;
cp libz.a $BASEDIR/lib;
cp zlib.h $BASEDIR/include;
cp zconf.h $BASEDIR/include;
cd ..;
echo -e "libz installed and copied to ./lib\n";
#########################

###### create libbam ######
tar -xzf bamLib.tar.gz
cd bamLib;
make libbam.a >> ../out.log 2>&1;
cp libbam.a $BASEDIR/lib;
cp bam.h $BASEDIR/include;
cp bgzf.h $BASEDIR/include;
cd ..;
echo -e "libbam installed and copied to ./lib\n";
###########################

###### create libgmp ##############################
tar -xzf gmpLib.tar.gz;
cd gmpLib;
mkdir GMP_INSTALL;
./configure --prefix=$BASEDIR/gmpLib/GMP_INSTALL >> ../out.log 2>&1;
make >> ../out.log 2>&1;
make check >> ../out.log 2>&1;
make install >> ../out.log 2>&1;
cp GMP_INSTALL/lib/libgmp.a $BASEDIR/lib;
cp GMP_INSTALL/include/gmp.h $BASEDIR/include;
cd ..;
echo -e "libgmp installed and copied to ./lib\n";
####################################################

#### generate LocHap ####
cd src;
#echo -e "\ngenerating LocHap binary ...\n\n";
#echo -e "-------------------------------\n";
make >> ../out.log 2>&1;
if [ -e "LocHap" ] && [ -s "LocHap" ] && [ -x "LocHap" ]
then 
	echo -e "\nLocHap binary generated !\n";
	mv LocHap ../bin/;
	echo -e "\n-------------------------------------\n";
	echo -e "Lochap installed and copied to bin ";
	echo -e "\n-------------------------------------\n";
else
	echo -e "\nError in LocHap binary generation !\n";
fi 
cd ..;
#########################

#### generate Plugin binaries ####
cd plugin;
#echo -e "\ngenerating additional binaries ...\n\n";
#echo -e "-------------------------------\n";
make >> ../out.log 2>&1;
mv addHomozygous ../bin/;
mv multiSampleRun ../bin/;
mv HCF2BED ../bin/;
cd ..;
#########################


################################
### run ./filter_install.sh ####
################################
cd filter;
./filter_install.sh
cd ..;
#####################################################

