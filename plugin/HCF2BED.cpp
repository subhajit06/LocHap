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

	#### HCF2BED.cpp
 	#### It is Plug-in that is part of LocHap-ver2.0
	#### It generates an BED file from an existing HCF file
	
 	#### Created on: Feb 6, 2016 as an add-on to LocHap
 	####  It need an exsting hcf and vcf file to generate the output
 	####      Author: subhajit sengupta
 	#### Assumes an sorted input VCF file !!
 
	#### Usage:: 
	#### ./HCF2BED <HCF_file_name> <output_BED_file_prefix> 
 	
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
#include "string.h"
#include "stdlib.h"

#define MAX_LINE_LEN 10000
int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		fprintf(stdout,"Usage::  ./HCF2BED <HCF_or_VCF_file_name> <BED_file_prefix>\n\n");
		return -1;
	}
	
	FILE* fp,*ft;
	fp = fopen(argv[1],"r");
	char outFileName[500];
	sprintf(outFileName,"%s.bed",argv[2]);
	ft = fopen(outFileName,"w");
	if(fp == NULL || ft == NULL)
	{
		fprintf(stderr,"File open error\n");
		return -1;
	}

	char tmpLine[MAX_LINE_LEN];
	char *pch1,*pch2;
	int cnt1=0,cnt2=0;
	char chrName[500];
	char posStr[500];
	char refStr[500];

	int nLine = 0;
	int pos[3];
	int refLEN=0;

	while(1)
	{
		fgets(tmpLine,MAX_LINE_LEN,fp);
		if(feof(fp))
			break;
		if((tmpLine[0] != '#') && (strlen(tmpLine) > 1))
		{
			nLine++;
			cnt1 = 0;
			cnt2 = 0;
			pch1 = NULL;
			pch2 = NULL;
			//int L = strlen(tmpLine)+1;
			//tmpLine[L] = '\0';
			//fprintf(stdout,"LINE:: %s\n",tmpLine);
			strcpy(chrName,"");		
			strcpy(posStr,"");
			strcpy(refStr,"");
			pch1 = strtok(tmpLine,"\t");
  			while (pch1 != NULL)
  			{
    			//printf("%s\n",pch1);
    			cnt1++;
				if(cnt1 == 1)
					strcpy(chrName,pch1);
				if(cnt1 == 2)
				{
					strcpy(posStr,pch1);
					//break;
				}
				if(cnt1==4)
				{
					strcpy(refStr,pch1);
					refLEN = strlen(refStr);
					break;
				}
				pch1 = strtok(NULL, "\t");
			}
			pch2 = strtok(posStr,",");
  			while (pch2 != NULL)
  			{
    			//printf("%s\n",pch2);
				pos[cnt2] = atoi(pch2);
    			cnt2++;
				pch2 = strtok(NULL,",");
			}
			
			//fprintf(stdout,"%d SNPs\n",cnt2);
	
			if(cnt2 == 1) //this will be used for VCF file
			{	
				fprintf(ft,"%s\t%d\t%d\tComments: Line = %d\n",chrName,pos[0],pos[0]+refLEN,nLine);
			}
			if(cnt2 == 2)
			{	
				fprintf(ft,"%s\t%d\t%d\tComments: Line = %d\n",chrName,pos[0],pos[0]+1,nLine);
				fprintf(ft,"%s\t%d\t%d\tComments: Line = %d\n",chrName,pos[1],pos[1]+1,nLine);
			}
			if(cnt2 == 3)
			{	
				fprintf(ft,"%s\t%d\t%d\tComments: Line = %d\n",chrName,pos[0],pos[0]+1,nLine);
				fprintf(ft,"%s\t%d\t%d\tComments: Line = %d\n",chrName,pos[1],pos[1]+1,nLine);
				fprintf(ft,"%s\t%d\t%d\tComments: Line = %d\n",chrName,pos[2],pos[2]+1,nLine);
			}
			
		}
	}

	fclose(fp);
	fclose(ft);

	fprintf(stdout,"\nBED file conversion is done !\n");

}
