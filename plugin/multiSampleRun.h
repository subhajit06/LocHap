/*#############################################################################
  
  # LocHap : Local-Haplotype Variant Calling Software.
  # Copyright (C) <2016>  <Subhajit Sengupta and Kamalakar Gulukota>
  # This file is part of LocHap.

  #  LocHap is free software: you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation, either version 3 of the License, or
  #  (at your option) any later version.

  #  LocHap is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #  GNU General Public License for more details.

  #  You should have received a copy of the GNU General Public License
  #  along with LocHap. If not, see <http://www.gnu.org/licenses/>.
  ############################################################################  


##############################################################################*/

/**************************************************************
***************************************************************

	#### multiSampleRun.h
 	#### It is Plug-in that is part of LocHap-ver2.0
	#### This wrapper is for running a multisample VCF file automatically
	#### Find all the samples in the VCF file and generate a shell script
	#### that contains the commands to run all the samples with appropriate BAM files
 	#### BAM file names should be named as <sample_name>.bam
 	#### LocHap binary, BAM files and VCF files need to be in the same directory
 	#### in order to run the shell script   	
	
	
 	#### Created on: Feb 6, 2016 as an add-on to LocHap
 	####  It need an exsting hcf and vcf file to generate the output
 	####      Author: subhajit sengupta
 	#### Assumes an sorted input VCF file !!
 
	#### Usage:: 
	#### ./muttiSampleRun --vcf VCF_FILE_NAME --out OUTPUT_FILE_NAME [--sig --size BLOCK_SIZE --qual Min_QUAL]
 	
	#### In this version of LocHap:	
	#### This is a bug fix for last line of VCF not parsing
	#### This code is for JOINT LIKELIHOOD 
	#### This code is fixed for dynamic memory allocation
	#### update for existing code
	
	LocHap: Calling LHVs for a sample where   
	input: bam + bai and associated vcf
	output: hcf or igv 
	Developed by: Subhajit Sengupta and Kamalakar Gulukota
	July 2014
	
	Details of the commands, output format are described
	in the Quick Manual. 

	Methods, results can be found in the Main paper and 
	Online methods titled "Ultra-Fast Local-Haplotype
	Variant Calling Using Paired-end DNA-Sequencing Data 
	Reveals Somatic Mosaicism in Tumor and Normal Blood Samples" 
	by Sengupta et al.

	citation: 
	Sengupta, S., Gulukota, K., Zhu, Y., Ober, C., Naughton, 
	K., Wentworth-Sheilds, W., Ji, Y. (2015), 
	Ultra-Fast Local-Haplotype Variant Calling Using 
	Paired-end DNA-Sequencing Data Reveals Somatic Mosaicism 
	in Tumor and Normal Blood Samples 
	(Nucleic Acids Research; http://dx.doi.org/10.1093/nar/gkv953)

		
***************************************************************
***************************************************************/

#ifndef MULTISAMPLERUN_H_
#define MULTISAMPLERUN_H_

#include "getopt.h"
#include "iostream"
#include "sstream"
#include "cstdio"
#include "cstring"
#include "cstdlib"
#include "new"
#include "vector"
#include "set"
#include "utility"

// constants
#define MAX_READ_LEN    50000
#define MAX_LINE_LEN    1000
#define EXIT_CMD_ERR    1
#define EXIT_FILE_ERR   2
#define EXIT_VCF_ERR    3
#define EXIT_BAM_ERR    4

#define MAX_STR_LEN    200

using namespace std;

/* Functions for command line parsing */
// variable for function pointer added @subhajit
int parseCommandLine (int ac, char **av, char *vn, char* on, int *bs, int *q, int* sf, int *ig);
void printSyntax (char *command_name);

/* Utility functions */
FILE *myopen (const char name[], const char md[]);
int my_fieldNumber (char *ln, const char *lfor, const char *delim);
char *getField (char *str, int n, const char *delim);

/* Functions for parsing VCF */
int my_readNextLine(char** pBuff,FILE* fp,int* TMP_READ_LEN);/// longer line read for VCF (adaptive). 17thNov, 2014
int my_skipComments (FILE *fp, char **dl, char **hl,int* TMP_READ_LEN);///modified 17th Nov,2014
int my_errorCheck_getSamples(int rf, int af, int fkf, char* hl, int* nSamples, vector<string> &sampleNames);


#endif /* MULTISAMPLERUN_H_ */

