/*#############################################################################
  
  # LocHap : Local-Haplotype Variant Calling Software.
  # Copyright (C) <2014>  <Subhajit Sengupta and Kamalakar Gulukota>
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

/********************************
	filter.h: header file used by
	filter.cpp
*********************************/

#include "iostream"
#include "fstream"
#include "cstdlib"
#include "cstdio"
#include "cstring"
#include "cmath"

#include "stdint.h"
#include "inttypes.h"
// Fisher 2-by-2 exact test 

#define SMALLISH_EPSILON 0.00000000003
#define SMALL_EPSILON 0.0000000000001
// This helps us avoid premature floating point overflow.
#define EXACT_TEST_BIAS 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125

/**
	to import the header of the bamtools library
*/
#include "api/BamReader.h"
#include "api/BamIndex.h"
#include "api/BamConstants.h"
#include "api/BamAlignment.h"
#include "api/BamAlgorithms.h"
#include "api/BamAux.h"
#include "api/SamHeader.h"
#include "api/SamSequence.h"
#include "api/SamSequenceDictionary.h"

// common constant for string manipulation
#define MAX_STR_LEN    	200
#define MAX_LINE_LEN    10000

/// mapping and base quality
#define MIN_MAP_QUAL    30
#define MIN_BASE_QUAL	30

/// nHaplo and nClus for HCF filter
#define N_HAPLO			3
#define	N_CLUSTER				N_HAPLO

#define MAX_SNP     			3

/// constants for BAM level filter
#define END_DIST_RATIO	0.10	/// 10% of the read length
#define	ALLELE_MIN_CNT	1
/// used in (1-p_value) so more the value more strand bias (for Fisher test) 
#define FISHER_P_VAL_THRESHOLD	0.95

/// not used right now
#define P_VAL_SIG_LEVEL	0.05  	/// used for strand-bias test


/// constants for VCF file part

#define PHRED_QUAL_OFFSET		33 

#define EXIT_VCF_ERR 			3
#define EXIT_FILE_ERR  			2

#define HOMOZYGOUS				0
#define HETEROZYGOUS			1
#define INDEL					2

/// Constants for filter types
#define TAIL_FILTER				0
#define LEG_FILTER				1
#define ELEPHANT_FILTER			2


using namespace std;
using namespace BamTools;

// definition of structures

typedef struct _SNP_readInfo
{
	int posSNP;
	char refSeq;
	char altSeq; // the top one other than the ref seq character 
	int n_total_mapped_read;
	int n_total_hQ_mapped_read;///passed mapping quality
	int n_total_no_indel_read;///
	int n_total_hQ_base_read;///
	int n_total_not_end_read;///
	int n_total_good_read;
	//int n_low_baseQual_read;
	int n_pos_direction_read;
	int n_neg_direction_read;


	// for strand direction-test	
	int n_tot_pos_direction_read; // n1
	int n_mut_pos_direction_read; // x1
	int n_tot_neg_direction_read; // n2
	int n_mut_neg_direction_read; // x2

	int n_nearby_SNP;
	int n_nearby_indel;
	//int n_end_readBase;
	bool goodSNP;

}SNP_readInfo;

typedef struct _LH_block
{
	string blockStr;
	string chrName;
	int nSNP;
	SNP_readInfo readInfoSNP[MAX_SNP];
	bool goodBlock;

}LH_block;

typedef struct _seqFreq
{
	char seq;
	int freq;
}seqFreq;

typedef struct _idPos
{
	int chrID;/// 0 to 21 for chromosome 1-22 and 23 and 24 for X and Y, respectively
	int pos;// with relative to that chromosome
}idPos;

typedef struct _vcfInfo
{
	string chrom;
	idPos iP;
	int infoType;/// 0 if homozygous, 1 if heterozygous, 2 if indel	
}vcfInfo;


/// function prototype

void initBlockStat(LH_block& L);
void initSNPReadInfo(SNP_readInfo& S);
  
void print_SNP_readInfo(SNP_readInfo& S);
void print_LH_block(LH_block& L); 

void print_SNP_readInfo(SNP_readInfo& S,ofstream& outFs);
void print_LH_block(LH_block& L,ofstream& outFs); 

bool indelExistBAM(BamAlignment& B);
int returnRelativePosOnRead(BamAlignment &B,int snpPos,int short_RL);

/* Utility functions */ 
FILE* vcfOpen (const char name[], const char md[]);
int fieldNumber (const char *ln, const char *lfor, const char *delim);
char *getField (char *str, int n, const char *delim);
/* Functions for parsing VCF */
int readNextLine (char *l, FILE *f);
void skipComments (FILE *fp, char *dl, char *hl);
void createParseString (char *str, int f);
void errorCheck (int rf, int af, int fkf, int svf);
int altLen (char *a, int w);

int getBamIDforChrom (string str1, SamHeader sH);
void getVCFInfo(const char vcfFileName[],const char sampleName[],SamHeader sH, vector<vcfInfo>& vcfVector);

bool passed_two_proportions_z_test(int x1,int n1,int x2,int n2);
bool passed_Fisher_exact_test(int x1,int n1,int x2,int n2);
double compute_quantile(int* x,int n,double f);

int determine_short_RL_from_BAM(BamReader& bR);

//// main Filtering functions - Level-1,Level-2 and Level-3
bool Level_1_Filtering_HCF(char* tStr,char* chrName,int* nSNP,long *posSNP,char* refSeq,int* tmpType,int short_RL,int filterType);
bool Level_2_Filtering_VCF(idPos tmpIdPos,vector<idPos>& F1,vector<idPos>& F2,int* tmpType,int short_RL,int filterType);
bool Level_3_Filtering_BAM(BamReader& bR, BamAlignment& al, LH_block& LH_blk,int i,long posSNP[],char refSeq[],int short_RL,
										int& level_3_subfilter_rem_cnt_1,int& level_3_subfilter_rem_cnt_2,
										int& level_3_subfilter_rem_cnt_3,int& level_3_subfilter_rem_cnt_4);


// cdf normal functions
double calculate_z_value(double p_val_sig);
double NormalCDFInverse(double p);
double RationalApproximation(double t);
int n_choose_k(int n,int k);


//// for fisherExact.cpp
double fisher2by2(uint32_t m11, uint32_t m12, uint32_t m21, uint32_t m22, uint32_t midp);


