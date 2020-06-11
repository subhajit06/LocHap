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

	#### multiSampleRun.cpp
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


#include "multiSampleRun.h"

int main (int argc, char **argv)
{
	/* Parameters to be set by command line options */
	char vcf_name[MAX_STR_LEN]; 
	//char bam_name[MAX_STR_LEN], sample_name[MAX_STR_LEN], out_name[MAX_STR_LEN];
	char shell_file_name[MAX_STR_LEN]; 
	int block_size, MinQ;
	// sigFlag = 1 means only significant called haplos will be in the HCF file
	// it is an optional flag in command line
	// by default it is 0; so NONE_SIGNIFICANT will be in the HCF file 
	// if it is 1; NONE_SIGNIFICANT will NOT be in the HCF file 
	int sigFlag = 0, igvFlag = 0;

	/* Parse the command line params, catch errors and echo appropriately */
	int iRet;
	iRet = parseCommandLine (argc, argv, vcf_name, shell_file_name, &block_size, &MinQ, &sigFlag, &igvFlag);
	if(iRet != 0)
	{
		fprintf(stderr,"\n");
		fprintf(stderr,">>>>> Error: User command error !\n");
		exit (EXIT_CMD_ERR);	
	}

	/* Print out a happy to analyze message */
	fprintf (stderr, "\nAnalyzing multisample VCF file (%s).\nBlock size = (%d).\nMin Quality = (%d).\nSig Flag = (%d).\nIGV Flag = (%d).\n",
									vcf_name,  block_size,	MinQ,	sigFlag, igvFlag);

	FILE *VCFP;
	//int i; 
	int ref_field, alt_field, fmt_key_field, hlen; 
	//int sample_val_field;
	//char vcf_line[MAX_LINE_LEN], header_line[MAX_LINE_LEN];

	int TMP_READ_LEN = MAX_READ_LEN; /// this TMP_READ_LEN is dynamic
	char* vcf_line = (char*)(calloc(TMP_READ_LEN,sizeof(char))); 
	char* header_line = (char*)(calloc(TMP_READ_LEN,sizeof(char)));

	if((vcf_line == NULL) || (header_line == NULL))
	{
		fprintf(stderr,">>>>> Error: Memory allocation !\n");
		return -1;
	}
	//int vcf_read_flag;

	VCFP = myopen (vcf_name, "r");

	/* Get the "FORMAT" and the needed sample field numbers from header line.  Check for errors in
 	* the header line */
	//skipComments (VCFP, vcf_line, header_line);
	iRet = my_skipComments (VCFP, &vcf_line, &header_line, &TMP_READ_LEN);
	if(iRet != 0)
	{
		fprintf(stderr,"\n");
		fprintf(stderr,">>>>> Error: VCF file memory allocation error !\n");
		exit(EXIT_VCF_ERR);
	}

	hlen = strlen (header_line);
	if (header_line[hlen-1] == '\n')
	{	
		header_line[hlen-1] = '\0';	
	}
	
	ref_field = my_fieldNumber (header_line, "REF", "\t");
	alt_field = my_fieldNumber (header_line, "ALT", "\t");
	fmt_key_field = my_fieldNumber (header_line, "FORMAT", "\t");
	//sample_val_field = my_fieldNumber (header_line, this_block.sampleName, "\t");

	int nSamples;
	vector<string> sampleNames;

	sampleNames.clear();
	iRet = my_errorCheck_getSamples (ref_field, alt_field, fmt_key_field, header_line, &nSamples, sampleNames);

	if(iRet != 0)
	{
		fprintf(stderr,"\n");
		fprintf(stderr,">>>>> Error: VCF file error !\n");
		exit(EXIT_VCF_ERR);
	}

	fprintf(stderr,"\nNumber of samples in the VCF = %d.\n\n",nSamples);
	int sCnt = 0;
	vector<string>::iterator it;
	for(it = sampleNames.begin(); it != sampleNames.end(); ++it)
    	fprintf(stderr,"Sample #%d:\t%s\n",++sCnt,(*it).c_str());

	fprintf(stderr,"\n\n");
		

	// form one command line for LocHap run and put it in the shell script for all the samples

	FILE* multiSample_fp;	
	
	multiSample_fp = myopen(shell_file_name, "w");

	fprintf(stderr,"Writing LocHap run commands to \"%s\" ... \n\n",shell_file_name);

	char* locahap_command_line = (char*)(calloc(TMP_READ_LEN,sizeof(char))); 


	fprintf(multiSample_fp,"#This shell script will execute LocHap on all the samples found in the VCF file.\n");
	fprintf(multiSample_fp,"#This file needs to be in the same directory as LocHap binary executable.\n\n");
	fprintf(multiSample_fp,"#Number of samples in the VCF = %d.\n",nSamples);
	for(int i = 0; i < nSamples; i++)
    	fprintf(multiSample_fp,"#Sample %d:\t%s\n",i,sampleNames[i].c_str());

	fprintf(multiSample_fp,"#########################################\n\n");
	

	for(int i=0;i<nSamples;i++)
	{
		fprintf(multiSample_fp,"echo \"\"\n");
		fprintf(multiSample_fp,"echo \"Running LocHap for sample %s.....\"\n",sampleNames[i].c_str());
			
		sprintf(locahap_command_line,"./LocHap --vcf %s --bam %s.bam --sample %s --out %s --size %d --qual %d",
								vcf_name, sampleNames[i].c_str(), sampleNames[i].c_str(), sampleNames[i].c_str(), block_size,MinQ);
		if(sigFlag == 1)
			sprintf(locahap_command_line,"%s --sig",locahap_command_line);
		
		if(igvFlag == 1)
			sprintf(locahap_command_line,"%s --igv",locahap_command_line);

		//fprintf(stderr,"%s\n",locahap_command_line);

		fprintf(multiSample_fp,"%s\n\n",locahap_command_line);

	}

	free(locahap_command_line);
	free(header_line);
	free(vcf_line);
	fclose(multiSample_fp);

	/* make the shell script executable */
	char cmd1[MAX_STR_LEN];
	sprintf(cmd1,"chmod +x %s",shell_file_name);
	iRet = system(cmd1);
	if(iRet != 0)
		fprintf(stderr,"Couldn't make the script executable; Please run 'chmod +x %s' on your console",shell_file_name);

	return (0);
}

int my_fieldNumber (char *ln, const char *lfor, const char *delim)
{/* Tokenizes ln on delim and looks for lfor.  Returns the field number for lfor
	return -1 if not found */
int tk = 0;
char *token; 
//char dt[MAX_LINE_LEN];
int ln_len = strlen(ln)+1;
char* dt = (char*)(calloc(ln_len,sizeof(char))); 

if(dt == NULL)
{
	fprintf(stderr,">>>>> Error: Memory allocation!!\n");
	return -2;
}
strcpy (dt, ln);

//fprintf(stderr,"HEADER_LINE(%s)::%s\n(FLD):%s\n(DELIM):%s\n",__FUNCTION__,dt,lfor,delim);
token = strtok (dt, delim);
while (token != NULL)	
{
	if (strcmp (token, lfor) == 0)
	{
		//fprintf(stderr,"tk = %d\n",tk);	
		return (tk);	
	}
	else	
	{
		token = strtok (NULL, delim);
		++tk;
	}
  }
/* All tokens exhausted.  lfor not found.  Return -1 */
return (-1);
}

char *getField (char *str, int n, const char *delim)
{
/* split str by delim.  Return a pointer to the n'th field */
char *t;
int i = 0;
t = strtok (str, delim);
while (t != NULL)	{
	if (i == n)
	  {	return (t);	}
	else	{
		t = strtok (NULL, delim);
		++i;
	  }
  }
/* all done before n't field could be reached
 * t right now is NULL.  Return it.  */
return (t);
}

FILE *myopen (const char name[], const char md[])
{
/* opens file named in name with mode md.  If open fails, 
 * prints an error msg and exits.  Else, returns the
 * file pointer
 */
FILE *FP;
if ((FP = fopen (name, md)) == NULL)	{
	fprintf (stderr, "File\n<%s>\ncannot be opened for mode <%s>\n", name, md);
	exit (EXIT_FILE_ERR);
  }
return FP;
}

int my_readNextLine(char** pBuff,FILE* fp,int* TMP_READ_LEN)
{
    if (fgets (*pBuff, *TMP_READ_LEN, fp) == NULL) 
    {   
        /* NULL.  EOF was read.  reset l, return */
        //fprintf(stderr,"EOF reached!!\n");
        (*pBuff)[0] = '\0';
        return 1; 
    }   
    else 
    {  
		//fprintf(stderr,"~~%s~~\n",(*pBuff)); 
        if ((int)strlen((*pBuff)) == (*TMP_READ_LEN) -1) 
        {   
            //fprintf(stderr,"limit reached !!\n");
            fseek( fp, -(*TMP_READ_LEN-1), SEEK_CUR);  
            (*TMP_READ_LEN) = 2*(*TMP_READ_LEN);  
            (*pBuff)[0] = '\0';
            (*pBuff) = (char*)(realloc((*pBuff),(*TMP_READ_LEN)*sizeof(char)));    
            if((*pBuff) == NULL)
			{
				fprintf(stderr,">>>>> Error: Reallocation memory !\n");
				return -1;
			}
			//fprintf(stderr,"success: reallocation !!\n");
            
			return 2;  
        }   
        else
        {   
            //(*lCnt)++;
            //int L = strlen((*pBuff));   
            //(*pBuff)[L-1] = '\0';
            //fprintf(stderr,"(%d):%s\n",(*TMP_READ_LEN),(*pBuff));
            //if(L>1)
            //	fprintf(ft,"%s\n",(*pBuff));
            return 0;
        }   
    }   
}

int my_skipComments (FILE *fp, char **dl, char **hl,int* TMP_READ_LEN)
{ /* Reads file pointed by fp, line by line.  Each comment line in it
	starts with a '#'.  This sub reads (and skips all comment lines.
	the last comment line is the header line - puts it in hl.
	first non-comment line is the data line - puts it in dl.
	Returns error code  */
	//fprintf(stderr,"inside %s !!\n",__FUNCTION__);

	int vcf_read_flag;
	vcf_read_flag = my_readNextLine (dl, fp, TMP_READ_LEN);	
	while ( (vcf_read_flag == 0) || (vcf_read_flag == 2) )	
	{
		//fprintf(stderr,"inside while\n");
		if(vcf_read_flag == 0)
		{
			/* If it is a comment, keep track as header line */
			if ((*dl)[0] == '#')
			{	
				int dl_len = strlen((*dl))+1;
    			(*hl) = (char*)(realloc((*hl),(dl_len)*sizeof(char)));    
        		if((*hl) == NULL)
				{
					fprintf(stderr,">>>>> Error: Memory reallocation(for VCF file comments) !\n");
					return -1;	
				}
				strcpy ((*hl), (*dl));	
			}
			/* If not a comment, break out of this loop */
			else
	  			break;
  		}
		//else // my_readNextLine returned 2
		//{
			//fprintf(stderr,"my_readNextLine returned 2\n");	
			//fprintf(stderr,"STR: %s; TMP_READ_LEN: %d\n",*dl,*TMP_READ_LEN);
		//}
		vcf_read_flag = my_readNextLine (dl, fp, TMP_READ_LEN);	
	}
  	return 0;
}

int my_errorCheck_getSamples (int rf, int af, int fkf, char* hl, int *nSamples, vector<string> &sampleNames)
{
/* rf is reference field, af is alt_field, fkf is "FORMAT" field
 * Here we do error checking.  If
 * format fails, exit.  Else return 
 * number of samples and their names 
 * iRet - return value to indicate error */

///@subhajit
int iRet = 0;
int tk = 0;
char *token;

int TMP_LEN = strlen(hl)+1;
char* hLine = (char*)(calloc(TMP_LEN,sizeof(char)));
if(hLine==NULL)
{	
	fprintf(stderr,"memory allocation error!!\n");
	return -1;
}
strcpy (hLine, hl);
(*nSamples) = 0;
/////////////////////////

	if ((rf != 3) || (af != 4) || (fkf != 8))	
	{
		fprintf (stderr, "VCF file appears mal-formatted:\n\
			ref is in %d th field.  ALT is in %d the field.  FORMAT is in %d the field.\n", rf+1, af+1, fkf+1);
	
		/* Other fields: ref, alt and format were not in the right place */
		if (rf != 3)
	  	{	fprintf (stderr, "Header line of VCF does not have REF as 4th field.\n");	}
		if (af != 4)
	  	{	fprintf (stderr, "Header line of VCF does not have ALT as 5th field.\n");	}
		if (fkf != 8)
	  	{	fprintf (stderr, "Header line of VCF does not have FORMAT as 9th field.\n");	}

		iRet = -1;
		//exit (EXIT_VCF_ERR);
	}
	else
	{
		token = strtok (hLine, "\t");
		//fprintf(stderr,"\nThe following samples were found in the VCF:\n");
		while (token != NULL)
		{
			token = strtok (NULL, "\t");
 			if(tk >= 8 && token != NULL)
			{
				//fprintf(stderr,"%s; ",token);
				sampleNames.push_back(string(token));
				(*nSamples)++;
			}
			++tk;
		}

		//fprintf(stderr,"\nNumber of samples in the VCF = %d\n",*nSamples);
			
		fprintf(stderr,"\n");
	}


  free(hLine);		
  return iRet;
}

int parseCommandLine (int ac, char **av, char *vn, char *on, int *bs, int *q, int *sf, int *igv)
{
/* Parses ac command line arguments in av.  Fills up the values:
 * vn - VCF file name (required)
 * *bs - block size (default = 500; must be > 50 and < 1000)
 * *q - min quality score (default 30; must be >=15)
 * *sf - variable to select print routine: Set to 1 if only
 * 	blocks with significant haplos are to be printed
 * *igv - variable to select print routine: Set to 1 if igv-compatible
 *	output is desired in the file (default : 0 print in hcf format)

 * iRet - return value denotes if the given command is successful

 */
int iRet = 0;

vn[0] = '\0';
on[0] = '\0';
*bs = 500;
*q = 30;
*sf = 0;
static struct option long_options[] = {
        /* These options don't set a flag. */
        {"sig", no_argument,       sf, 1},
        {"igv", no_argument,       igv, 1},
	/* These options set a flag.
	 * We distinguish them by their indices. */
	/* longopt	type		    ptr_t  return_value */
        {"vcf",		required_argument,	0, 'v'},
        {"out",		required_argument,	0, 'o'},
	/* Optional arguments:  out defaults to "sample".HCF; block_size to 500
	 * and MinQ to 30 */
        {"size",	required_argument,	0, 'z'},
        {"qual",	required_argument,	0, 'q'},
        {0, 0, 0, 0}
    };
int option_index = 0, c;

/* getopt_long stores the option index here. */
while ((c = getopt_long_only (ac, av, "v:o:z:q:", long_options, &option_index)) != -1)	{
	switch (c) {
	case 'v': strcpy (vn, optarg); break;
	case 'o': strcpy (on, optarg); break;
	case 'z': sscanf (optarg, "%d", bs); break;
	case 'q': sscanf (optarg, "%d", q); break;

	case 0:		case '?':	default:
		break;
   }
  }

/* If any command line arguments are not known or if any req params not set  */
if ((optind < ac) || (vn[0] == '\0') || (on[0] == '\0') 
		|| (*bs < 50) || (*bs > 1000)
		|| (*q < 15))	{
	/* Unknown options */
	if (optind < ac)	{
		fprintf (stderr, ">>>>>Error: Unknown options: ");
		while (optind < ac)
		  {	fprintf (stderr, "%s ", av[optind++]);	}
		fprintf (stderr, "\n");
		iRet = -1;
	  }
	if (vn[0] == '\0')
	  {	fprintf (stderr, ">>>>> Error:  Please specify name of VCF File\n");}
	if (on[0] == '\0')
	  {	fprintf (stderr, ">>>>> Error:  Please specify name of output script File\n");}
	if ((*bs < 50) || (*bs > 1000))
	  {	fprintf (stderr, ">>>>> Error:  Please set the block size (--size) to between 50 and 1000 (default = 500)\n");}
	if (*q < 15)
	  {	fprintf (stderr, ">>>>> Error:  Please set the quality threshold (--qual) to greater than 15 (default = 30)\n");}
	printSyntax (av[0]);
	iRet = -1;
	//exit (EXIT_CMD_ERR);
  }

/* If outname is not set, set it now based on sample name */
//if (on[0] == '\0')
// {	sprintf (on, "%s.%s", sn, (*igv) ? "igv" : "hcf");	}
//else
//  {	strcat (on, (*igv) ? ".igv" : ".hcf");	}

return iRet;
}

void printSyntax (char *command_name)
{
/* Print the syntax.  command_name is name of this command */
fprintf (stderr, "\nMulti-Sample LocHap Script Generation:\n%s --vcf VCF_FILE_NAME --out OUTPUT_FILE_NAME [--sig --size BLOCK_SIZE --qual Min_QUAL]\n", command_name);

fprintf (stderr, "\nMulti-Sample LocHap Script Generation:\n%s --vcf VCF_FILE_NAME --out OUTPUT_FILE_NAME [--sig --size BLOCK_SIZE --qual Min_QUAL --igv]\n", command_name);
fprintf (stderr, "\n\
Option	Type	Meaning\n\
------  ----    -----------------------------------------------------\n\
vcf	string	(required) Name of VCF file to analyze for identifying blocks\n\n\
out	string	(required) Name of output script file to run multi-sample analysis\n\n\
sig	Flag	(optional) Sets a flag for printing only blocks with at least\n\
		one significant haplotyope called\n\
		defaults to ==> 0: all blocks will be printed\n\n\
size	int	(optional) size of blocks while analyzing VCF\n\
		(default = 500 NTs; must be between 50 and 1000)\n\n\
qual	int	(optional) ignore all reads with Phred-scaled Mapping Quality less than this.\n\
		Ignore all bases with Base Quality less than this\n\
		(default = 30; must be greater than 15)\n\n\
igv	Flag	(optional) Sets a flag for printing in format suitable for loading into igv\n\
		defaults to ==> 0: prints in native hcf format\n\
		Output file name is <sample>.hcf or <sample>.igv, based on this choice\n\n\
");
return;
}
