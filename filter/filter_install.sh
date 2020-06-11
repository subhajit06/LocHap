#!/bin/bash
################################
### run ./filter_install.sh ####
################################

##### generate binary executable filter in bin folder #####

#####################################################
### generating libbamtools - bamtools C++ library ###
echo -e "\ngenerating filter binary ...\n";
tar -xzf bamtools_CPP.tar.gz;
cd ./bamtools_CPP;
mkdir build;
cd ./build;
cmake .. > ../out.log 2>&1;
make >> ../out.log 2>&1;
cd ..;
cd ..;
echo -e "\nlibbamtools installed\n";
#####################################################
### compile filter src file ###
make >> ../out.log 2>&1;

if [ -e "filter" ] && [ -s "filter" ] && [ -x "filter" ]
then
	echo -e "\nfilter binary generated !\n";
	mv filter ../bin;
else
	echo -e "\nError in filter binary generation !\n";
fi
#####################################################


