## Compile with this:
## run make in this folder
all: multiSampleRun addHomozygous HCF2BED
multiSampleRun: multiSampleRun.cpp multiSampleRun.h	
	g++ -g -I. -Wall -o multiSampleRun multiSampleRun.cpp
addHomozygous:	addHomozygous.cpp
	g++ -g -Wall -o addHomozygous addHomozygous.cpp
HCF2BED: HCF2BED.cpp
	g++ -g -Wall -o HCF2BED HCF2BED.cpp	