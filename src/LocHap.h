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
 
 
############################################################################*/
/**************************************************************
***************************************************************

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

/************************
 LocHap.h: header file 
required for LocHap.cpp
*************************/

#ifndef LOCHAP_H_
#define LOCHAP_H_

#include "getopt.h"
#include "iostream"
#include "sstream"
#include "cmath"
#include "cstdio"
#include "cstring"
#include "cstdlib"
#include "new"
#include "vector"
#include "set"
#include "algorithm"
#include "map"
#include "utility"

#include "gmp.h"  /// we need multiple precision library

#include "bam.h"  /// we need bam library

// constants
#define MAX_READ_LEN    50000
#define MAX_LINE_LEN    1000
#define EXIT_CMD_ERR    1
#define EXIT_FILE_ERR   2
#define EXIT_VCF_ERR    3
#define EXIT_BAM_ERR    4

#define MAX_STR_LEN    200
#define MAX_SNPS	3
#define MAX_HAPLO_CALL	8
#define FDR_THRESHOLD	0.01
#define NUM_BASES	4
#define BIN_BASES	2

#define HYPER_PARAM_ALPHA 0.05
#define HYPER_PARAM_BETA	1.00	

#define MIN_LEN 3 ///required for JL

using namespace std;

/// data structures
typedef struct _baseTable {
	int cnt;
	char id[NUM_BASES];
  } baseTable;

typedef struct _HCF_struct	{
	int MinQ;
	double error_rate;
	char sampleName[MAX_STR_LEN];
	FILE *OUT;

	/* Defining the SNPs (number, chrom, posn's and ref_seq */
	int nSNP;
	char chr[MAX_STR_LEN];
	long pos[MAX_SNPS]; // should be long @subhajit
	char ref[MAX_SNPS+1];

	// Frequencies and Read-based haplotypes
	int numBlankReads;	/* # of reads that are all ____ */
	int numDiscrepantReads;	/* paired end reads that disagree about one or more bases */
	map <string, int> Freq;
	
	// required for bayesModule
	int numInformativeReads;// number of reads excluding all "___" @subhajit
	vector<pair<double,string> > Haplo_call;
	vector<double> Haplo_fdr;
	int numSigHaploCall;// how many are significant calls for the given fdr // @subhajit

	//optional fields
	int missing_arr[MAX_SNPS+1];
	int distinct_noMiss_cnt;
 
	////JL@subhajit 
	vector<pair<double,string> > JHL_call;
	int nJL;///should be same as JHL_call.size()

	//// fix the bug for analyzing and printing the last block
	int analyzed_printed;	


	} HCF_struct;

/* Functions for command line parsing */
// variable for function pointer added @subhajit
int ParseCommandLine (int ac, char **av, char *vn, char *bn, char *sn, char *on, int *bs, int *q, int* sf, int *ig);
void printSyntax (char *command_name);

/* Utility functions */
FILE *myopen (const char name[], const char md[]);
int fieldNumber (const char *ln, const char *lfor, const char *delim);
int my_fieldNumber (char *ln, const char *lfor, const char *delim);
char *getField (char *str, int n, const char *delim);
string Combine (char *nw, string ld);
int countD_ (string s, int *D);

/* Functions for parsing VCF */
int readNextLine (char *l, FILE *f);
void skipComments (FILE *fp, char *dl, char *hl);
int my_skipComments (FILE *fp, char **dl, char **hl,int* TMP_READ_LEN);///modified 17th Nov,2014
void createParseString (char *str, int f);
////@subhajit
int errorCheck (int rf, int af, int fkf, int svf,char* hl);
int my_errorCheck (int rf, int af, int fkf, int svf,char* hl);

int altLen (char *a, int w);

/* Utilities for dealing with HCF Structure */
void addSNP (HCF_struct *b, int p, char *rf);
int lastPosition (HCF_struct *b);
int printBlock (HCF_struct *b);
void firstInitialize (HCF_struct *b, const char *s, const char *o, int q);
void initialize (HCF_struct *tb, char *c, int p, char *r);
int incrementFreqs (HCF_struct *tb);

/* bam related functions and macros */
/* 1 for A, 2 for C, 4 for G, 8 for T and 15 for N.  */
#define base(i)	("0AC3G567T901234N"[(i)])
#define cigar_op(i)	("MIDNSHP=X"[(i)])
void CIGAR_Operate (uint32_t FullCig, int *readC, int *genC);
static int HaploOfAlign (const bam1_t *a, void *d);
int getBamIDforChrom (const char *c, bam_header_t *b);

//functions for printing HCF
void printHCFHeader(FILE* in_fpOut);
void printHeaderLine (FILE *f, char *sample_name);
void printCommand (FILE *f, int ac, char **av);
void printBamContigs (FILE *f, bam_header_t *b);
void printSig2Igv (HCF_struct* in_hcfVar);
void print2Igv (HCF_struct* in_hcfVar);
void printTrackLine (FILE *f, char *sample_name);

//// bayes module functions
int runBayesModel(HCF_struct* in_hcfVar);
int dataPopulate(HCF_struct* in_hcfVar,char** out_data,int* out_freq,char** out_sym,char** out_binStr);
//now commented //int fillHCFStruct(HCF_struct* out_HCFVar);
int fillOptStruct(HCF_struct* in_hcfVar);
void printHCFStruct(HCF_struct* in_hcfVar);
void printToHCFFile(HCF_struct* in_hcfVar);
void printDataToHCFFile(HCF_struct* in_hcfVar);
void printSigToHCFFile(HCF_struct* in_hcfVar);
int preScreen(int in_RL,int in_N,char** in_data,char** in_sym,int* out_hK,int *out_bSet);
void calculatePrDataPriorTable(int in_RL,char** in_data,int in_N,char** in_hSet,int in_hK,
		int** in_mat_a,int** in_mat_d,int* in_arr_m,double** in_mat_e,char** in_sym,
		double** out_P_data,double* out_prior_table,int** out_lambdaStr);
void genSymTab(int L,char** str_arr,char** binStr);
void convertToBase(int in_v, int in_L,int in_B,char* out_str);
int buildBinaryStrings(int M,int **lambdaStr);
void buildErrMatrix(int in_N,int in_RL,double in_err, double** out_mat_e);
void buildMatrices(char** in_data,int in_N,int in_RL,char** in_hSet,int in_hK,
		int** out_mat_a,int** out_mat_d, int* out_arr_m);
void indices_minus_j(int in_hK,int in_j,int* out_ids);
void fillLambdaStr(int in_hK,int* in_ids,int* in_otherLambdaStr,int* out_lambdaStr);

//// utils functions
int bin2dec(int* in_bStr,int in_hK);
int findID_Str(int in_L, char* in_arr_a, char in_elem, int* out_id);
int find_missing_positions(char* in_s,int in_RL,int* out_ids); ////@subhajit: 07.16:14
int fdr_cumsum_cal(double *in_d,int in_K,double* out_d);
bool pairCompare(const pair<double, string>& A, const pair<double, string>& B);
double calcGammaFactor(double a,double b);
int my_readNextLine(char** pBuff,FILE* fp,int* TMP_READ_LEN);/// longer line read for VCF (adaptive). 17thNov, 2014

////JL functions
void my_dec2bin(int i,int L,int *bin_str);
int isSubset(int *a,int L1, int *b,int L2);
double find_all_summable_string_sum(int *bin_str,int hK,int zero_cnt, mpf_t *MM_gmp,int L3,mpf_t r_sumMM);
void printJLToHCFFile(HCF_struct* in_hcfVar);
void printINTstring(int hK,int *S);
#endif /* LOCHAP_H_ */

