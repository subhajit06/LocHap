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

	#### addHomozygous.cpp
 	#### It is Plug-in that is part of LocHap-ver2.0
	#### It generates an HCF file from an existing HCF file
	#### with adding Homozygous SNVs that exist between or within a 
	#### distance(specified by user) from the Heterozygous SNVs. 

 	#### Created on: Feb 6, 2016 as an add-on to LocHap
 	####  It need an exsting hcf and vcf file to generate the output
 	####      Author: subhajit sengupta
 	#### Assumes an sorted input VCF file !!
 
	#### Usage:: 
	#### ./addHomozygous <HCF_file_name> <VCF_file_name> <Max_len_either_side> <sample> <output_file_prefix>
 	
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



#include "stdio.h"
#include "assert.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "cstring"
#include "vector"

using namespace std;

#define MAX_LINE_LEN	50000
#define MAX_FLD_LEN		500
#define SAMPLE_OFFSET	9
#define FORMAT_INDX		8
#define VCF_CHR_NAME_INDX	0
#define VCF_POS_INDX	1
#define	VCF_REF_INDX	3
#define	VCF_ALT_INDX	4

#define HCF_CHR_INDX	1
#define HCF_POS_INDX	2
#define HCF_REF_INDX	3
#define HCF_N_SIG_INDX	4
#define HCF_SIG_HAP_INDX	5
#define HCF_ALL_HAP_INDX	6

typedef struct _SNV_VCF
{
	char chr[50];
	int pos;
	char base;
	char refBase;
	int hom; // 1 if homozygous 0 otherwise

}SNV_VCF;

int readNextLine (char *l, FILE *f)
{
	if (fgets (l, MAX_LINE_LEN, f) == NULL)
	{
		/* NULL.  EOF was read.  reset l, return */
		l[0] = '\0';
	   	return 1;
	}
	else
	{
		if (strlen (l) == MAX_LINE_LEN -1)
		{
			/* Line too long - do some malloc stuff.
	 		 * For now, just terminate and return. */
 	 		//error("Encountered line longer than %d:\n%s\n\nEnding ...\n", MAX_LINE_LEN, l);
 	 		fprintf(stdout,"Encountered line longer than %d:\n%s\n\nEnding ...\n", MAX_LINE_LEN, l);
	     	l[0] = '\0';
	     	return 2;
		}
		else
		{
			/* Proper line read. Return 0 */
	    	return 0;
	   	}
	}
}

void VCF_parse(FILE* fp_vcf, char* sample_name, vector<SNV_VCF>& vcf_rcd, int* VCF_n_rcd)
{

	char vcf_line[MAX_LINE_LEN], header_line[MAX_LINE_LEN];
	char parse_str[MAX_LINE_LEN];
	char* pch,*pch1;
	//char *pch2;
	int cnt1;
	int GT_fld;
	char fmt_fld[MAX_FLD_LEN];
	char sample_GT_str[MAX_FLD_LEN];
	int e_flag;

	int header_line_cnt = 0;

	char chrName[50];
	char refStr[50];
	char altStr[50];

	int posSNV;

	int vcf_line_cnt = 0;
	int sample_idx = 0;

	SNV_VCF info;

	while (readNextLine(vcf_line, fp_vcf) == 0)
	{
		/* If it is a comment, keep track as header line */
		if (vcf_line[0] == '#')
		{
			header_line_cnt++;
			strcpy (header_line, vcf_line);
		}
		/* If not a comment, break out of this loop */
		else
		{
			vcf_line_cnt++;
			if(vcf_line_cnt == 1)
			{
				// we can now get index of the sample
				cnt1=0;
				pch = strtok (header_line,"\t\n");
				while (pch != NULL)
				{
					pch = strtok (NULL, "\t\n ");
					if(strcmp(sample_name,pch) == 0)
						break;
					cnt1++;
				}

				sample_idx = cnt1+1;
			}
			/*  Chr Pos ID Ref Alt Qua Fil Inf Fmt */
			/*   0   1   2  3   4   5   6   7   8  */
			strcpy (parse_str, "%s %d %*s %s %s %*s %*s %*s %s");
			for (int i = SAMPLE_OFFSET; i < sample_idx; ++i)
			{
				strcat (parse_str, " %*s");
			}
			 /* read sample field as string */
			strcat (parse_str, " %s");

			sscanf (vcf_line, parse_str, chrName, &posSNV, refStr, altStr, fmt_fld, sample_GT_str);
			/*
			cnt1=0;
			pch = strtok(vcf_line,"\t\n");
			while (pch != NULL)
			{
				if(cnt1==VCF_CHR_NAME_INDX)
					strcpy(chrName,pch);
				if(cnt1==VCF_POS_INDX)
					sscanf(pch,"%d",&posSNV);
				if(cnt1==VCF_REF_INDX)
					strcpy(refStr,pch);
				if(cnt1==VCF_ALT_INDX)
					strcpy(altStr,pch);
				if(cnt1==FORMAT_INDX)
					strcpy(fmt_fld,pch);
				if(cnt1 == sample_idx)
				{
					strcpy(sample_GT_str,pch);
					break;
				}
				pch = strtok (NULL, "\t\n");
				cnt1++;
			}
			*/

			cnt1=0;
			pch = strtok (fmt_fld,":");
			while (pch != NULL)
			{
				if(strcmp("GT",pch)==0)
				{
					GT_fld=cnt1;
					break;
				}
				pch = strtok (NULL, ":");
				cnt1++;
			}

			e_flag=0;

			if( (strlen(refStr) != 1) || (strlen(altStr) != 1) )
				e_flag=1;
			if(strcmp(sample_GT_str,"./.")==0 || strcmp(sample_GT_str,".|.")==0 )
				e_flag=1;

			if(e_flag == 0)
			{
				int j=0;
				pch1 = strtok(sample_GT_str,":");
				while (pch1 != NULL)
				{
					if(j==GT_fld)
					{
						if ((strcmp(pch1,"./.")==0) || (strcmp(pch1,".|.")==0))
						{
							e_flag = 1;
							break;
						}
						if ((strcmp(pch1,"0/1")==0) || (strcmp(pch1,"0|1")==0) || (strcmp(pch1,"1/0")==0) || (strcmp(pch1,"1|0")==0))
						{
							strcpy(info.chr,chrName);
							info.pos = posSNV;
							info.hom = 0;
							info.refBase = refStr[0];
							info.base = refStr[0];
							vcf_rcd.push_back(info);
							(*VCF_n_rcd)++;
							break;
						}
						if((strcmp(pch1,"0/0")==0) || (strcmp(pch1,"0|0")==0))
						{
							//homozygous wildtype
							strcpy(info.chr,chrName);
							info.pos = posSNV;
							info.hom = 1;
							info.refBase = refStr[0];
							info.base = refStr[0];
							vcf_rcd.push_back(info);
							(*VCF_n_rcd)++;
							break;
						}

						if((strcmp(pch1,"1/1")==0) || (strcmp(pch1,"1|1")==0))
						{
							//homozygous variant
							strcpy(info.chr,chrName);
							info.pos = posSNV;
							info.hom = 1;
							info.refBase = refStr[0];
							info.base = altStr[0];
							vcf_rcd.push_back(info);
							(*VCF_n_rcd)++;
							break;

						}
					}
					pch1 = strtok (NULL, ":");
					j++;
				}//while

			} // if e_flag==0


		} //else
	}//while file

}

void padHomozygousSNPs(char* origStr,char* maskedStr,char* returnedStr)
{
	int hJ = 0;
	// find heterozygous locations on input str
	unsigned int i=0;
	for(i=0;i<strlen(maskedStr);i++)
	{
		if(maskedStr[i] == 'H')
		{
			returnedStr[i] = origStr[hJ];
			hJ++;
		}
		else
		{
			returnedStr[i] = maskedStr[i];
		}
	}

	returnedStr[i] = '\0';

	//fprintf(stderr,"\ninput1 = %s; input2 = %s; output1 = %s\n",origStr,maskedStr,returnedStr);

}


void form_masked_and_ref_str(int nSNV,char* chrName,int* pos, vector<SNV_VCF> vcf_rcd,int VCF_n_rcd,int boundary_len,
		char* maskedStr,char *oldRefStr,char* newRefStr,char* newPosStr,int* oldVCFIndx)
{
    int start_pos,end_pos;
    start_pos = pos[0]-boundary_len-1;
    end_pos = pos[nSNV-1]+boundary_len+1;
    int base_cnt = 0;
    vector<char> ref_base_arr;
    vector<char> hap_base_arr;
    strcpy(newPosStr,"");
    char tempPosStr[50];
    char newPosStr1[MAX_FLD_LEN];
    strcpy(newPosStr1,"");
    int newVCFIndx;

    ref_base_arr.clear();
    hap_base_arr.clear();

    for(int i=(*oldVCFIndx);i<VCF_n_rcd;i++)
    {
        if((strcmp(vcf_rcd[i].chr,chrName)==0) && (vcf_rcd[i].pos > start_pos ) && (vcf_rcd[i].pos < end_pos ) )
        {

            if((nSNV==2 && (vcf_rcd[i].pos != pos[0]) && (vcf_rcd[i].pos != pos[1]) )
                || (nSNV==3 && (vcf_rcd[i].pos != pos[0]) && (vcf_rcd[i].pos != pos[1]) && (vcf_rcd[i].pos != pos[2]) ))
            {
            	if(vcf_rcd[i].hom == 1)
            	{
            		ref_base_arr.push_back(vcf_rcd[i].refBase);
            		hap_base_arr.push_back(vcf_rcd[i].base);
            		sprintf(tempPosStr,"%d,",vcf_rcd[i].pos);
            		strcat(newPosStr1,tempPosStr);
            	}
            }

            if(vcf_rcd[i].pos == pos[0])
            {
                ref_base_arr.push_back(oldRefStr[0]);
                hap_base_arr.push_back('H');
                sprintf(tempPosStr,"%d,",vcf_rcd[i].pos);
                strcat(newPosStr1,tempPosStr);

            }
            if(vcf_rcd[i].pos == pos[1])
            {
            	ref_base_arr.push_back(oldRefStr[1]);
            	hap_base_arr.push_back('H');
            	newVCFIndx = i;
            	sprintf(tempPosStr,"%d,",vcf_rcd[i].pos);
            	strcat(newPosStr1,tempPosStr);
            }
            if(nSNV == 3 && vcf_rcd[i].pos == pos[2])
            {
                ref_base_arr.push_back(oldRefStr[2]);
                hap_base_arr.push_back('H');
                newVCFIndx = i;
                sprintf(tempPosStr,"%d,",vcf_rcd[i].pos);
                strcat(newPosStr1,tempPosStr);
            }

            base_cnt++;
        } // end if
        if(vcf_rcd[i].pos > end_pos )
        	break;
    } // end for

    int L = strlen(newPosStr1);
    strncpy(newPosStr,newPosStr1,L-1); // remove "," at the end
    newPosStr[L-1] = '\0';
    if((int)ref_base_arr.size() != base_cnt)
        fprintf(stderr,"problem!!\n");
    else
    {	int j;
        for(j=0;j<base_cnt;j++)
        {
            maskedStr[j] = hap_base_arr[j];
            newRefStr[j] = ref_base_arr[j];
        }
        maskedStr[j] = '\0';
        newRefStr[j] = '\0';
    }

    *oldVCFIndx = newVCFIndx;
    ref_base_arr.clear();
    hap_base_arr.clear();
}

int main(int argc, char* argv[])
{
	if(argc < 6)
	{
		fprintf(stdout,"Usage:: ./addHomozygous <HCF_file_name> <VCF_file_name> <Max_len_either_side> <sample> <output_file_prefix> \n\n");
		return -1;
	}

	fprintf(stdout,"\nGenerating HCF file with Homozygous loci ... \n");
	FILE* fp_hcf,*fp_vcf;
	FILE* fp_out;
	fp_hcf = fopen(argv[1],"r");
	fp_vcf = fopen(argv[2],"r");
	char outFileName[MAX_FLD_LEN];
	sprintf(outFileName,"%s_homozygous.hcf",argv[5]);
	fp_out = fopen(outFileName,"w");

	if(fp_hcf == NULL || fp_vcf == NULL || fp_out == NULL)
	{
		fprintf(stderr,"File open error\n");
		return -1;
	}
	char sample_name[MAX_FLD_LEN];
	sprintf(sample_name,"%s",argv[4]);
	int boundary_len = atoi(argv[3]);

	//VCF part
	vector<SNV_VCF> vcf_rcd;
	int VCF_n_rcd = 0;
	VCF_parse(fp_vcf, sample_name, vcf_rcd, &VCF_n_rcd);
    fprintf(stderr,"\nVCF parsing done with %d records !!\n",VCF_n_rcd);
    /// print VCF records
	/*
    for(int i=0;i<VCF_n_rcd;i++)
    {
        fprintf(stdout,"%s\t%d\t%d\t%c\t%c\n",vcf_rcd[i].chr, vcf_rcd[i].pos,
            vcf_rcd[i].hom, vcf_rcd[i].refBase, vcf_rcd[i].base);
    }
	*/

	//HCF part
	char tmpLine[MAX_LINE_LEN];
	char tmpLine1[MAX_LINE_LEN];

	char *pch1,*pch2,*pch3, *pch4;
	int cnt1=0,cnt2=0,cnt3=0,cnt4=0;
	char chrName[MAX_FLD_LEN];
	char posStr[MAX_FLD_LEN];
	char posStr1[MAX_FLD_LEN];
	char hcfRefStr[10];
	char hapStrArr[MAX_FLD_LEN];
	char sigHapStrArr[MAX_FLD_LEN];
	char hapStr[MAX_FLD_LEN];
	char newHapStr[MAX_FLD_LEN];
	char newHapStrArr[MAX_FLD_LEN];
	char newSigHapStrArr[MAX_FLD_LEN];

	int nLine = 0;
	int pos[3];
	int nSNV;
	char nSigHapStr[10];
	char *newHCFStr;
	char* hcfTailStr;

	int oldVCFIndx = 0;

	while(1)
	{
		readNextLine(tmpLine, fp_hcf);
		strcpy(tmpLine1,tmpLine);
		if(feof(fp_hcf))
			break;
		if((tmpLine1[0] != '#') && (strlen(tmpLine1) > 1))
		{
			nLine++;
			cnt1 = 0;
			cnt2 = 0;
			cnt3 = 0;
			cnt4 = 0;
			pch1 = NULL;
			pch2 = NULL;
			pch3 = NULL;
			pch4 = NULL;

			int newLen = strlen(tmpLine1)+MAX_FLD_LEN;
			newHCFStr = (char*) (calloc(newLen,sizeof(char)));
			hcfTailStr = (char*) (calloc(newLen,sizeof(char)));
			strcpy(hcfTailStr,"");

			//int L = strlen(tmpLine)+1;
			//tmpLine[L] = '\0';
			//fprintf(stdout,"LINE:: %s\n",tmpLine);
			strcpy(chrName,"");
			strcpy(posStr,"");
			strcpy(hapStrArr,"");
			pch1 = strtok(tmpLine1,"\t\n");
			while (pch1 != NULL)
			{
				//printf("%s\n",pch1);
				cnt1++;
				if(cnt1 == HCF_CHR_INDX)
					strcpy(chrName,pch1);
				if(cnt1 == HCF_POS_INDX)
					strcpy(posStr,pch1);
				if(cnt1 == HCF_REF_INDX)
					strcpy(hcfRefStr,pch1);
				if(cnt1 == HCF_N_SIG_INDX)
					strcpy(nSigHapStr,pch1);
				if(cnt1 == HCF_SIG_HAP_INDX)
					strcpy(sigHapStrArr,pch1);
				if(cnt1 == HCF_ALL_HAP_INDX)
					strcpy(hapStrArr,pch1);
				if(cnt1 > HCF_ALL_HAP_INDX)
				{
					strcat(hcfTailStr,pch1);	// last part of new hcf line
					strcat(hcfTailStr,"\t");
				}
				pch1 = strtok(NULL, "\t\n");
			}

			int L1 = strlen(hcfTailStr);
			hcfTailStr[L1-1] = '\0';
			// tail part
			//fprintf(stdout,"%s\n",hcfTailStr);

			strcpy(posStr1,posStr);
			pch2 = strtok(posStr1,",");

			while (pch2 != NULL)
			{
				//printf("%s\n",pch2);
				pos[cnt2] = atoi(pch2);
				cnt2++;
				pch2 = strtok(NULL,",");
			}

			nSNV = cnt2;
			/// form the masked str
			/// which has H is the Heterozygous position and correct base in Homozygous position
			char maskedStr[MAX_FLD_LEN];
			char newRefStr[MAX_FLD_LEN];
			char newPosStr[MAX_FLD_LEN];

			strcpy(maskedStr,"");
			strcpy(newRefStr,"");
			strcpy(newPosStr,"");

			//oldVCFIndx = 0;
			form_masked_and_ref_str(nSNV,chrName,pos,vcf_rcd,VCF_n_rcd,boundary_len,maskedStr,hcfRefStr,newRefStr,newPosStr,&oldVCFIndx);

			sprintf(newHCFStr,"%s\t%s\t%s\t%s\t",chrName,newPosStr,newRefStr,nSigHapStr); /// first part of new hcf line
			//fprintf(stdout,"%s\n",newHCFStr);

			int nSigHap = atoi(nSigHapStr);
			/// sig haplotype str part
			strcpy(newSigHapStrArr,"");
			pch3 = strtok(sigHapStrArr,",");
         	char tempStr1[MAX_FLD_LEN];
            char tempStr2[MAX_FLD_LEN];
            char tempStr3[MAX_FLD_LEN];
            char tempStr4[MAX_FLD_LEN];

			while (pch3 != NULL)
			{
				//printf("%s\n",pch3);
               strcpy(tempStr1,"");
               strcpy(tempStr2,"");
               strcpy(tempStr3,pch3);
			   strncpy(hapStr,pch3,nSNV);
               int L = strlen(tempStr3);
               strncpy(tempStr2,pch3+nSNV,L-nSNV);
               tempStr2[L-nSNV] = '\0';
               padHomozygousSNPs(hapStr,maskedStr,newHapStr);
               strcpy(tempStr1,newHapStr);
               strcat(tempStr1,tempStr2);
               if(cnt3 < (nSigHap-1))
                   strcat(tempStr1,",");
               cnt3++;
			   pch3 = strtok(NULL,",");
			   strcat(newSigHapStrArr,tempStr1);
			}

			//fprintf(stdout,"%s\n",newSigHapStrArr);

			//// all haplotype str
			int nAllHap = 4;
			if(nSNV == 3)
				nAllHap = 8;
			strcpy(newHapStrArr,"");
			pch4 = strtok(hapStrArr,",");
			while (pch4 != NULL)
			{
				//printf("%s\n",pch3);
               strcpy(tempStr1,"");
               strcpy(tempStr2,"");
               strcpy(tempStr4,pch4);
	       strncpy(hapStr,pch4,nSNV);
	       int L = strlen(tempStr4);
               strncpy(tempStr2,pch4+nSNV,L-nSNV);
               tempStr2[L-nSNV] = '\0';
               padHomozygousSNPs(hapStr,maskedStr,newHapStr);
               strcpy(tempStr1,newHapStr);
               strcat(tempStr1,tempStr2);
               if(cnt4 < (nAllHap-1))
                   strcat(tempStr1,",");

			   cnt4++;
			   pch4 = strtok(NULL,",");

			   strcat(newHapStrArr,tempStr1);
			}

			//fprintf(stdout,"%s\n\n",newHapStrArr);

			/// form the new line for hcf file

			strcat(newHCFStr,newSigHapStrArr);
			strcat(newHCFStr,"\t");
			strcat(newHCFStr,newHapStrArr);
			strcat(newHCFStr,"\t");
			strcat(newHCFStr,hcfTailStr);

			///write back to file
			//fprintf(stdout,"%s\n",newHCFStr); //add \n ?
			fprintf(fp_out,"%s\n",newHCFStr); //add \n ?

			free(newHCFStr);
			free(hcfTailStr);
		}
		else
		{
			fprintf(fp_out,"%s",tmpLine); //add \n ?
		}
	}// while file

    fclose(fp_hcf);
    fclose(fp_vcf);
    fclose(fp_out);
	
    fprintf(stdout,"\nnew HCF file %s is generated !!\n\n",outFileName); 
}

