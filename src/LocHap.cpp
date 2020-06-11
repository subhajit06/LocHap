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

#include "LocHap.h"

map <string, string> ReadHaplo;
double PRIOR_PROB;

int main (int argc, char **argv)
{
/* Parameters to be set by command line options */
char vcf_name[MAX_STR_LEN], bam_name[MAX_STR_LEN], sample_name[MAX_STR_LEN], out_name[MAX_STR_LEN];
int block_size, MinQ;
// sigFlag = 1 means only significant called haplos will be in the HCF file
// it is an optional flag in command line
// by default it is 0; so NONE_SIGNIFICANT will be in the HCF file 
// if it is 1; NONE_SIGNIFICANT will NOT be in the HCF file 
int sigFlag = 0, igvFlag = 0;

/* Parse the command line params, catch errors and echo appropriately */
int iRet;
iRet = ParseCommandLine (argc, argv, vcf_name, bam_name, sample_name, out_name, &block_size, &MinQ, &sigFlag, &igvFlag);
if(iRet != 0)
{
	fprintf(stderr,"\n");
	fprintf(stderr,">>>>> Error: User command error !\n");
	exit (EXIT_CMD_ERR);	
}

/* Print out a happy to analyze message */
fprintf (stderr, "\nAnalyzing VCF (%s) for sample (%s) using bam (%s).\nOutput to => (%s)\nBlock size = (%d).\nMin Quality = (%d).\nSig Flag = (%d)\n",
				vcf_name,  sample_name,		bam_name,	out_name,	block_size,		MinQ,	sigFlag);


// decide the print routine by using a function pointer
// by default it prints NONE_SIGNiFICANT line
// function pointer for print routine
void (*print_function)(HCF_struct* in_hcfVar);
if ((sigFlag == 0) && (igvFlag == 0))
  {	print_function = &printToHCFFile;	}
else if ((sigFlag == 1) && (igvFlag == 0))
  {	print_function = &printSigToHCFFile;	}
else if ((sigFlag == 0) && (igvFlag == 1))
  {	print_function = &print2Igv;	}
else if ((sigFlag == 1) && (igvFlag == 1))
  {	print_function = &printSig2Igv;	}


FILE *VCFP;
int i, ref_field, alt_field, fmt_key_field, sample_val_field, hlen;
//char vcf_line[MAX_LINE_LEN], header_line[MAX_LINE_LEN];

int TMP_READ_LEN = MAX_READ_LEN; /// this TMP_READ_LEN is dynamic
char* vcf_line = (char*)(calloc(TMP_READ_LEN,sizeof(char))); 
char* header_line = (char*)(calloc(TMP_READ_LEN,sizeof(char)));

if((vcf_line == NULL) || (header_line == NULL))
{
	fprintf(stderr,">>>>> Error: Memory allocation !\n");
	return -1;
}
int vcf_read_flag;

/* Open the VCF file and skip all comments in it.  (Retain the last comment line) */
HCF_struct this_block;
firstInitialize (&this_block, sample_name, out_name, MinQ);
ReadHaplo.clear();

VCFP = myopen (vcf_name, "r");

/* Make sure the bam file is openable */
bam_header_t *BamHeader;
bamFile bamFP;
bam_index_t *bidx;

if ((bamFP = bam_open (bam_name, "r")) == NULL)		{
	fprintf (stderr, "Cannot open bam file: %s\n", bam_name);
	exit (EXIT_FILE_ERR);
  }
BamHeader = bam_header_read (bamFP);
bidx = bam_index_load (bam_name);
if(bidx == NULL)
{
	fprintf(stderr,"Cannot open the index file of the bam file. \n");
	exit (EXIT_FILE_ERR);
}

/* Get the "FORMAT" and the needed sample field numbers from header line.  Check for errors in
 * the header line */
//skipComments (VCFP, vcf_line, header_line);
iRet = my_skipComments (VCFP, &vcf_line, &header_line,&TMP_READ_LEN);
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
sample_val_field = my_fieldNumber (header_line, this_block.sampleName, "\t");
iRet = my_errorCheck (ref_field, alt_field, fmt_key_field, sample_val_field, header_line);
//@subhajit
if(iRet != 0)
{
	fprintf(stderr,"\n");
	fprintf(stderr,">>>>> Error: VCF file error !\n");
	exit(EXIT_VCF_ERR);
}
//// @subhajit: all inputs are ok, open the HCF file for writing
this_block.OUT = myopen (out_name, "w");
printHCFHeader (this_block.OUT);
printCommand (this_block.OUT, argc, argv);
printBamContigs (this_block.OUT, BamHeader);
if (igvFlag)
  {	printTrackLine (this_block.OUT, sample_name);	}
else
  {	printHeaderLine (this_block.OUT, sample_name);	}

//////////////////////////////////////////////////////////

/* Now vcf_line is the "next" (first) vcf_line.  Process them one by one till end*/
int posn, line_num = 1, num_printed = 0, target_id = -1, same;
int Freq_Sig_Haplo[MAX_HAPLO_CALL+1], not_analyzed = 0;
char parse_string[MAX_LINE_LEN], chrom[MAX_LINE_LEN], ref_seq[MAX_LINE_LEN], alt_seq[MAX_LINE_LEN], fmt_key[MAX_LINE_LEN], fmt_val[MAX_LINE_LEN];
createParseString (parse_string, sample_val_field);

for (i = 0; i <= MAX_HAPLO_CALL; ++i)
  {	Freq_Sig_Haplo[i] = 0;	}

fflush (this_block.OUT);
/// based on beta-bernouli prior
PRIOR_PROB = calcGammaFactor(HYPER_PARAM_ALPHA,HYPER_PARAM_BETA);

/*fprintf(stdout,"gamma-factor = %0.12f\n",PRIOR_PROB);*/
fprintf(stdout,"\n\nPerforming Local Haplotype Analysis and generating %s file ..... \n\n",out_name);

//TMP_READ_LEN = 100; // to test the functionlity in the main loop

/* vcf_line gets blank when EOF is reached */
while (vcf_line[0] != '\0')	{
	/* Get the chromosome, position, format_key and format_vals for our sample */
	sscanf (vcf_line, parse_string, chrom, &posn, ref_seq, alt_seq, fmt_key, fmt_val);

	/* Get the format field number from the key and retrieve that field from the values */
	int f;
	char *v;
	/* Ignore lines where there is no GT field in the format */
	if ((f = my_fieldNumber (fmt_key, "GT", ":")) >= 0)	{
		/* Get the GT value */
		v = getField (fmt_val, f, ":");
		int lp = lastPosition (&this_block);

		/* Can't read the genotype string, exit */
		if ((strlen (v) != 3) || (v[1] != '/'))		{
			fprintf (stderr, "Genotype %s not parseable on line %d of file (%s)\n", v, line_num, vcf_name);
			exit (EXIT_VCF_ERR);
		  }

		/* if its homozygous, or if indel (ref or alt > 1 NT long):  do nothing (skip it)*/
		else if ((v[0] == v[2])
			|| (strlen (ref_seq) != 1)
			|| (altLen (alt_seq, v[0] - '0') != 1) || (altLen (alt_seq, v[2] - '0') != 1))
		  {	;	}
		/* If in old Chromosome, within distance:  We are within the old block.  Add a SNP */
		else if (((same = strcmp (chrom, this_block.chr)) == 0)
					&& ((posn - lp) < block_size))
		  {	addSNP (&this_block, posn, ref_seq);	  }
		else	{
			/* Starting, a new block: First Print the old block if it has 2 - MAX_SNPS snps*/
			if ((this_block.nSNP > 1) && (this_block.nSNP <= MAX_SNPS))	{
				/* Fetch all alignments that overlay this block */
				bam_fetch (bamFP, bidx, target_id, this_block.pos[0], lp, &this_block, HaploOfAlign);

				incrementFreqs (&this_block);
				/* Call Bayes Module here: fill rest of the fields
				 * in this_block */

				/* printBlock (&this_block); */
				/* 	runBayesModel gets called and the results are printed in the HCF file */
				if (runBayesModel(&this_block) >= 0)	{
					Freq_Sig_Haplo[this_block.numSigHaploCall] +=1;
					(*print_function)(&this_block);
					this_block.analyzed_printed = 1;
					++num_printed;
				  }
				//else
				//	fprintf(stdout,"1.Bayes module throws an error !!\n");	
			  }

			if (this_block.nSNP > MAX_SNPS)
			  {	++not_analyzed;	}

			/* set target_id if going to next chromosome */
			if (same != 0)	{
				target_id = getBamIDforChrom (chrom, BamHeader);
				if (target_id < 0)	{
					fprintf (stderr, "Chromosome %s (Line %d of %s) not found in BAM %s\n",
						chrom, line_num, vcf_name, bam_name);
					exit (EXIT_BAM_ERR);
				  }
			  }

			//fprintf (stderr, "Started new block.  (%d printed) from %d lines\n", num_printed, line_num);
			/* initialize the block for next block */
			initialize (&this_block, chrom, posn, ref_seq);
			ReadHaplo.clear ();
		  }
	  }
	//readNextLine (vcf_line, VCFP);
	
	vcf_read_flag = my_readNextLine (&vcf_line, VCFP, &TMP_READ_LEN);
	//fprintf(stderr,"out vcf_read_flag: %d ###\n",vcf_read_flag);
	while(vcf_read_flag == 2)
	{
		vcf_read_flag = my_readNextLine (&vcf_line, VCFP, &TMP_READ_LEN);
		//fprintf(stderr,"in vcf_read_flag: %d ###\n",vcf_read_flag);
		if(vcf_read_flag == 0)
			break;
	}
	++line_num;
  }

/* 
	This is a bug fix !! after parsing the entire vcf file, check the final segment 
	if that is analyzed before based on the flag analyzed_printed in the structure.
	in the following code, the final block is getting analyzed and printed.
*/

if ((this_block.analyzed_printed == 0) && (this_block.nSNP > 1) && (this_block.nSNP <= MAX_SNPS)) 
{
	/// we need to get this last position which is saved inside the  while loop before
	int lp = lastPosition (&this_block);
	//fprintf(stdout,"here here lp = %d!!\n",lp);
	
	/* Fetch all alignments that overlay this block */
    bam_fetch (bamFP, bidx, target_id, this_block.pos[0], lp, &this_block, HaploOfAlign);
  
    incrementFreqs (&this_block);
    /* Call Bayes Module here: fill rest of the fields
    * in this_block */

	//printBlock(&this_block);
    if (runBayesModel(&this_block) >= 0)	
    {
		Freq_Sig_Haplo[this_block.numSigHaploCall] +=1;
		(*print_function)(&this_block);
		this_block.analyzed_printed = 1;
		++num_printed;
    }
	//else
	//	fprintf(stdout,"2.Bayes module throws an error !!\n");
	
}

if(igvFlag == 0)
{
	fprintf (this_block.OUT, "\n# Analyzed %d Variants in %s\n", line_num-1, vcf_name);
	for (i = 0; i <= MAX_HAPLO_CALL; ++i)
  	{	fprintf (this_block.OUT, "# Number of Blocks with %d significant haplotypes = %d\n", i, Freq_Sig_Haplo[i]);	}
	fprintf (this_block.OUT, "#\n# Number of Blocks with more than %d Variants = %d (Not analyzed)\n", MAX_SNPS, not_analyzed);
}
fclose (this_block.OUT);
bam_close(bamFP);

free(header_line);
free(vcf_line);

return (0);
}

void initialize (HCF_struct *tb, char *c, int p, char *r)
{/* initializes all the values in strucnt pointed by tb.
 * sets chromosome to c, posn[0] to p, initiates refseq with r */

/* initialize the block */
tb->nSNP = 1;
tb->pos[0] = p;
strcpy (tb->chr, c);
strcpy (tb->ref, r);

/* ALl other ints set to 0 */
tb->numBlankReads = tb->numInformativeReads
	= tb->numSigHaploCall
	= tb->distinct_noMiss_cnt
	= tb->numDiscrepantReads
	= 0;
tb->analyzed_printed = 0;

//fprintf(stderr,"call 1\n");

int i;
for (i = 0; i <=MAX_SNPS; ++i)
  {	tb->missing_arr[i] = 0;		}

/* Clear all maps */
tb->Freq.clear ();
tb->Haplo_call.clear ();
tb->Haplo_fdr.clear ();
////JL@subhajit
tb->nJL = 0;
tb->JHL_call.clear ();

return;
}

void firstInitialize (HCF_struct *b, const char *s, const char *o, int q)
{
strcpy (b->sampleName, s);
b->chr[0] = b->ref[0] = '\0';
//@subhajit
//b->OUT = myopen (o, "w");

b->MinQ = q;
b->error_rate = pow (10.0, ((double) -q / 10.0));
b->nSNP = 0;
b->analyzed_printed = 0;

return;
}

int fieldNumber (const char *ln, const char *lfor, const char *delim)
{/* Tokenizes ln on delim and looks for lfor.  Returns the field number for lfor
	return -1 if not found */
int tk = 0;
char *token, dt[MAX_LINE_LEN];
strcpy (dt, ln);

token = strtok (dt, delim);
while (token != NULL)	{
	if (strcmp (token, lfor) == 0)
	  {	return (tk);	}
	else	{
		token = strtok (NULL, delim);
		++tk;
	  }
  }
/* All tokens exhausted.  lfor not found.  Return -1 */
return (-1);
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

int readNextLine (char *l, FILE *f)
{
/* Reads teh next line from f and puts it into l.  Returns chars read.
 * Need a util on top of fgets to deal with the case where the line might
 * be too long.  We'll leave it alone for now, but a malloc call may be needed
 */
if (fgets (l, MAX_LINE_LEN, f) == NULL)	{
	/* NULL.  EOF was read.  reset l, return */
	l[0] = '\0';
	return 1;
  }
else if (strlen (l) == MAX_LINE_LEN -1)	{
	/* Line too long - do some malloc stuff.
	 * For now, just terminate and return. */
	fprintf (stderr, "Encountered line longer than %d:\n%s\n\nEnding ...\n", MAX_LINE_LEN, l);
	l[0] = '\0';
	return 2;
  }
else	{
	/* Proper line read. Return 0 */
	return 0;
  }
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

void skipComments (FILE *fp, char *dl, char *hl)
{ /* Reads file pointed by fp, line by line.  Each comment line in it
	starts with a '#'.  This sub reads (and skips all comment lines.
	the last comment line is the header line - putsit in hl.
	first non-comment line is the data line - puts it in dl.
	Returns nothing */
while (readNextLine (dl, fp) == 0)	{
	/* If it is a comment, keep track as header line */
	if (dl[0] == '#')
	  {	strcpy (hl, dl);	}
	/* If not a comment, break out of this loop */
	else
	  {	break;	}
  }
return;
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

void createParseString (char *str, int f)
{
/* Creates a format string for use in sscanf.  We want to read
 * 0'th field as string - chromosome
 * 1'th field as int - position
 * 3'th field as string - reference seq
 * 4'th field as string - alternate sequence
 * 8'th field as string - format_key_string
 * f'th field as string - val_string for sample
 */
	/*  Chr Pos ID Ref Alt Qua Fil Inf Fmt */
    /*   0   1   2  3   4   5   6   7   8  */
strcpy (str, "%s %d %*s %s %s %*s %*s %*s %s");

/* f is guaranteed to be > 8 (prev error checking) */
/* Skip all fields except f */
int i;
/* ignore all fields till f'th */
for (i = 9; i < f; ++i)
  {	strcat (str, " %*s");	}

/* read f'th field as string */
strcat (str, " %s");

return;
}

int errorCheck (int rf, int af, int fkf, int svf,char* hl)
{
/* rf is reference field, af is alt_field, fkf, if "FORMAT" field and
 * svf is sample's format values field.  Here we do error checking.  If
 * format fails, exit.  Else return 
 * iRet - return value to indicate error */

///@subhajit
int iRet = 0;
int tk = 0;
char *token;
char hLine[MAX_LINE_LEN];
strcpy (hLine, hl);

/////////////////////////

if ((svf < 8) || (rf != 3) || (af != 4) || (fkf != 8))	{
	fprintf (stderr, "VCF file appears mal-formatted:\n\
		ref is in %d th field.  ALT is in %d the field.  FORMAT is in %d the field.  Sample is in %d the field\n", rf+1, af+1, fkf+1, svf+1);
	/* Sample was not found or is among the first 9 fields */
	if (svf < 0)
	{	
		fprintf (stderr, "Our Sample was not found in the VCF.\n");	
		/////@subhajit
		token = strtok (hLine, "\t");
		fprintf(stderr,"\nThe following samples were found in the VCF:\n");
		while (token != NULL)
		{
			token = strtok (NULL, "\t");
 			if(tk >= 8 && token != NULL)
				fprintf(stderr,"%s; ",token);
			++tk;
		}
		fprintf(stderr,"\n");
		//////////////////////
	}
	else if (svf < 8)
	  {	fprintf (stderr, "Sample field is among the mandatory first 9 fields of the VCF.\n");	}

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
	
  return iRet;
}

int my_errorCheck (int rf, int af, int fkf, int svf,char* hl)
{
/* rf is reference field, af is alt_field, fkf, if "FORMAT" field and
 * svf is sample's format values field.  Here we do error checking.  If
 * format fails, exit.  Else return 
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

/////////////////////////

if ((svf < 8) || (rf != 3) || (af != 4) || (fkf != 8))	{
	fprintf (stderr, "VCF file appears mal-formatted:\n\
		ref is in %d th field.  ALT is in %d the field.  FORMAT is in %d the field.  Sample is in %d the field\n", rf+1, af+1, fkf+1, svf+1);
	/* Sample was not found or is among the first 9 fields */
	if (svf < 0)
	{	
		fprintf (stderr, "Our Sample was not found in the VCF.\n");	
		/////@subhajit
		token = strtok (hLine, "\t");
		fprintf(stderr,"\nThe following samples were found in the VCF:\n");
		while (token != NULL)
		{
			token = strtok (NULL, "\t");
 			if(tk >= 8 && token != NULL)
				fprintf(stderr,"%s; ",token);
			++tk;
		}
		fprintf(stderr,"\n");
		//////////////////////
	}
	else if (svf < 8)
	  {	fprintf (stderr, "Sample field is among the mandatory first 9 fields of the VCF.\n");	}

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

  free(hLine);		
  return iRet;
}

int altLen (char *a, int w)
{/* a is alternate sequence (potentially has multiple comma sep values). 
  * w is like 1, 2, 3, ... such that we have to take the 1'th, 2'th etc of the 
  * multiple values.   Here we return the length of the w'th element.  If w is 0,
  * return 1 */

if (w == 0)
  {	return 1;	}

/* Which element do we want?  w'th */
int len = 0, ele = 0;
while (*a != '\0')	{
	if (*a == ',')	{
		++ele;
		if (ele == w)
		  {	return (len);	}
		len = 0;
	  }
	else
	  {	++len;	}
	++a;
  }
++ele;
if (ele == w)
  {	return (len);	}
else
  {	return -1;	}
}

void addSNP (HCF_struct *b, int p, char *rf)
{
/* Adds a SNP to block pointed by b.  New position is p and ref seq at p is rf */

/* Have to add the ref seq to the ref_seq string - but only if it is within
 * Max allowed length */
if (b->nSNP >= MAX_SNPS)
  {	b->pos[MAX_SNPS-1] = p;	}
else	{
	b->ref[b->nSNP] = *rf;
	b->ref[b->nSNP+1] = '\0';
	b->pos[b->nSNP] = p;
  }

/* increment the number of SNPs in the block */
++b->nSNP;

return;
}

int lastPosition (HCF_struct *b)
{
/* Returns the position of the last SNP in block b */
if (b->nSNP > MAX_SNPS)
  {	return b->pos[MAX_SNPS - 1];	}
else if (b->nSNP >0)
  {	return (b->pos[b->nSNP-1]);	}
else
  {	return -50000;	}
}

int printBlock (HCF_struct *b)
{
/* Prints block b to file pointer b->OUT */
/* If only 1 SNP, nothing to do */
fprintf (b->OUT, "NumSNP = %d\n\
Chrom = %s\n\
Ref = %s\n\
numBlankReads = %d\n\
numDiscrepantReads = %d\n\
numInformativeReads = %d\n", b->nSNP, b->chr, b->ref, b->numBlankReads, b->numDiscrepantReads, b->numInformativeReads);

int i;
for (i = 0; i < b->nSNP; ++i)	{
	fprintf (b->OUT, "Pos[%d] = %ld\tmissing_arr[%d] = %d\n", i, b->pos[i], i, b->missing_arr[i]);
  }
fprintf (b->OUT, "\t\tmissing_arr[%d] = %d\n", i, b->missing_arr[i]);

if (b->nSNP < 2)
  {	return 0;	}
else	{
	/*if This block has too many SNPs.  Don't print */
	if (b->nSNP > MAX_SNPS)
	  {	return 1;	}

	/* OK to print */
	else	{
		if (2 == b->nSNP)
		  {	fprintf (b->OUT, "%s\t%ld,%ld\t%s\t2\n", b->chr, b->pos[0], b->pos[1], b->ref);	}
		if (3 == b->nSNP)
		  {	fprintf (b->OUT, "%s\t%ld,%ld,%ld\t%s\t3\n", b->chr, b->pos[0], b->pos[1], b->pos[2], b->ref);	}

		/* Now to print the Freq table */
		std::map <string, int>::iterator it;
		for (it = b->Freq.begin(); it != b->Freq.end(); ++it)
		  {	fprintf (b->OUT, "\t\t\t%s\t%d\n", it->first.c_str(), it->second);	}
		return 2;
	  }
  }
fprintf (b->OUT, "Done with printBlock\n");
}

int ParseCommandLine (int ac, char **av, char *vn, char *bn, char *sn, char *on, int *bs, int *q, int *sf, int *igv)
{
/* Parses ac command line arguments in av.  Fills up the values:
 * vn - VCF file name (required)
 * bn - BAM file name (required)
 * sn - sample name (required)
 * on - output file name (default = "sn".HCF)
 * *bs - block size (default = 500; must be > 50 and < 1000)
 * *q - min quality score (default 30; must be >=15)
 * *sf - variable to select print routine: Set to 1 if only
 * 	blocks with significant haplos are to be printed
 * *igv - variable to select print routine: Set to 1 if igv-compatible
 *	output is desired in the file (default : 0 print in hcf format)

 * iRet - return value denotes if the given command is successful

 */
int iRet = 0;

vn[0] = bn[0] = sn[0] = on[0] = '\0';
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
        {"bam",		required_argument,	0, 'b'},
        {"sample",	required_argument,	0, 's'},

	/* Optional arguments:  out defaults to "sample".HCF; block_size to 500
	 * and MinQ to 30 */
        {"out",		required_argument,	0, 'o'},
        {"size",	required_argument,	0, 'z'},
        {"qual",	required_argument,	0, 'q'},
        {0, 0, 0, 0}
    };
int option_index = 0, c;

/* getopt_long stores the option index here. */
while ((c = getopt_long_only (ac, av, "v:b:s:o:z:q:", long_options, &option_index)) != -1)	{
	switch (c) {
	case 'v': strcpy (vn, optarg); break;
	case 'b': strcpy (bn, optarg); break;
	case 's': strcpy (sn, optarg); break;
	case 'o': strcpy (on, optarg); break;
	case 'z': sscanf (optarg, "%d", bs); break;
	case 'q': sscanf (optarg, "%d", q); break;

	case 0:		case '?':	default:
		break;
   }
  }

/* If any command line arguments are not known or if any req params not set  */
if ((optind < ac) || (vn[0] == '\0') || (bn[0] == '\0') || (sn[0] == '\0')
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
	if (bn[0] == '\0')
	  {	fprintf (stderr, ">>>>> Error:  Please specify name of BAM File\n");}
	if (sn[0] == '\0')
	  {	fprintf (stderr, ">>>>> Error:  Please specify name of Sample to analyze\n");}
	if ((*bs < 50) || (*bs > 1000))
	  {	fprintf (stderr, ">>>>> Error:  Please set the block size (--size) to between 50 and 1000 (default = 500)\n");}
	if (*q < 15)
	  {	fprintf (stderr, ">>>>> Error:  Please set the quality threshold (--qual) to greater than 15 (default = 30)\n");}
	printSyntax (av[0]);
	iRet = -1;
	//exit (EXIT_CMD_ERR);
  }

/* If outname is not set, set it now based on sample name */
if (on[0] == '\0')
  {	sprintf (on, "%s.%s", sn, (*igv) ? "igv" : "hcf");	}
else
  {	strcat (on, (*igv) ? ".igv" : ".hcf");	}

return iRet;
}

void printSyntax (char *command_name)
{
/* Print the syntax.  command_name is name of this command */
fprintf (stderr, "\nSYNOPSIS (HCF format):\n%s --vcf VCF_FILE_NAME --bam BAM_FILE_NAME --sample SAMPLE_NAME [--sig --out OUTFILE --size BLOCK_SIZE --qual Min_QUAL]\n", command_name);

fprintf (stderr, "\nSYNOPSIS (IGV format):\n%s --vcf VCF_FILE_NAME --bam BAM_FILE_NAME --sample SAMPLE_NAME [--sig --out OUTFILE --size BLOCK_SIZE --qual Min_QUAL --igv]\n", command_name);
fprintf (stderr, "\n\
Option	Type	Meaning\n\
------  ----    -----------------------------------------------------\n\
vcf	string	(required) Name of VCF file to analyze for identifying blocks\n\n\
bam	string	(required) Name of BAM file for getting read based haplotypes\n\n\
sample	string	(required) Name of sample to examine in the VCF file\n\
		The BAM file must correspond to this sample\n\
		This sample must be included in the VCF file\n\n\
sig	Flag	(optional) Sets a flag for printing only blocks with at least\n\
		one significant haplotyope called\n\
		defaults to ==> 0: all blocks will be printed\n\n\
out	string	(optional) Prefix of output file.  (Defaults to <sample>)\n\n\
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

void CIGAR_Operate (uint32_t FullCig, int *readC, int *genC)
{/* Takes the CIGAR operatoe in FullCig.  Increments read coordinate (readC) and/or
  * Genome coordinate (genC) as needed */
int op = FullCig & BAM_CIGAR_MASK;
int len = FullCig >> BAM_CIGAR_SHIFT;
/* Check what the operation is.  If M, increment both genome and read coord.
 * 'I', increment only read.  'D', increment only genome */
/* #define cigar_op(i)	("MIDNSHP=X"[(i)]) */

     /* 'M'	     '='	  'X' : THose are the cases */
if ((0 == op) || (7 == op) || (8 == op))	{
	/* In all these cases, read is aligned to the genome.
	 * Increment both gc and rc */
	*readC += len;
	*genC += len;
  }
	/*  'I' (ins)   'S' (softclip).  Increment only rc */
else if ((1 == op) || (4 == op))
  {	*readC += len;		}
	/*  'H' (hardclip)   'P' (padding).  do nothing */
else if ((5 == op) || (6 == op))
  {	;	/* do nonthing */	}
	/*  'D' (Del)   'N' (intron).  Increment only gc */
else if ((2 == op) || (3 == op))
  {	*genC += len;		}
else	{
	fprintf (stderr, "Unknown CIGAR Operator : %d of length %d\n", op, len);
	exit (10);
  }

return;
}

static int HaploOfAlign (const bam1_t *a, void *tb)
{/* Fills Haplo with bases from alignment a which matches to genome coordinates 
  * in HCF struct of this block in tb. (Remember these coords are 1-based while
  * bam file coords are 0-based).  If this read was already seen (checks the
  * ReadHaplo map, combines the old with the new and updates the Freqs map.
  * If this read has not yet been seen, stores it in the ReadHaplo map.
  * Returns 0 if alignment is added.  1 other wise*/

int iw;
unsigned char *q = bam1_qual (a);
char Haplo[MAX_SNPS+1];
HCF_struct *blo;
blo = (HCF_struct *) tb;

/* If mapping score is less than required Min, do nothing */
if (a->core.qual < blo->MinQ)
  {	return 0;	}

for (iw = 0; iw < blo->nSNP; ++iw)	{
	/* Set Haplo to default '_'.  We'll set to other value, if needed, in the rest*/
	Haplo[iw] = '_';
	int i, rc = 0, gc = a->core.pos, wantedC;

	/* The '- 1' takes care of the fact that VCF coords are 1-based while BAM coordinates are 0-based */
	wantedC = blo->pos[iw] - 1;
	/* If this position is to left of gc, it remains default _
	 * If this coord is to the right of gc, possible that we get a base from this */
	if (wantedC >= gc)	{
		uint32_t *cig = bam1_cigar (a);

		for (i = 0; i < a->core.n_cigar; ++i)	{
			uint32_t l = cig[i];
			/* CIGAR Operate on gc and rc for this */
			CIGAR_Operate (l, &rc, &gc);

			/* If we crossed wanted coordinate genome, get the base */
			if (gc > wantedC)	{
				/* Crossed wanted point.  go back to wanted point */
				rc = rc - (gc - wantedC);
				/* Take the base only if Quality is acceptable*/
				if (q[rc] > blo->MinQ)
				  {	Haplo[iw] = base(bam1_seqi (bam1_seq(a), rc));	}
				break;
			  }
		  } /* All Cigar operations done*/
	  } 
   }
/* Close out Haplo */
Haplo[iw] = '\0';

string rn = string (bam1_qname(a));

/* Put into ReadHaplo a combination of this Haplo and what was already there for this read's pair */
ReadHaplo[rn] = Combine (Haplo, ReadHaplo[rn]);
return 1;
}

int incrementFreqs (HCF_struct *tb)
{
/*tb is a pointer to struct of this block.  THis function
 * takes all singletons in it and increments the appropriate 
 * frequencies.  Returns the number of singletons incremented */

int dne = 0, b, d;
std::map <string, string>::iterator it;
for (it = ReadHaplo.begin(); it != ReadHaplo.end(); ++it)	{
	++dne;
	b = countD_ (it->second, &d);
	tb->missing_arr[b] += 1;

	if (b == tb->nSNP)	{
		tb->numBlankReads += 1;
		continue;
	  }
	
	/* Check if it is Discrepant */
	if (d > 0)
	  {	tb->numDiscrepantReads += 1;	}
	/* Else It is a good read */
	else	{
		tb->Freq[it->second] += 1;
		tb->numInformativeReads += 1;
	  }
  }
return dne;
}

int countD_ (string s, int *D)
{ /*Fills into *D and *U the number of D's and Underscores in s */
int u;
*D = u = 0;

for (std::string::iterator c = s.begin (); c != s.end (); ++c)	{
	if (*c == '_')
	  {	++u;	}
	else if (*c == 'D')
	  {	*D = *D + 1;	}
  }

return u;
}

int getBamIDforChrom (const char *c, bam_header_t *b)
{
/* looks for chromosome c in the bam header b.  Returns its id (for use in 
 * bam parsing.  If not found, returns -1 */
int i;
for (i = 0; i < b->n_targets; ++i)	{
	if (strcmp (c, b->target_name[i]) == 0)
	  {	return i;	}
  }
/* Not found.  Return -1 */
return -1;
}

string Combine (char *nw, string ld)
{
/* n and old are strings (haplos of the new and old read of a red ID.
 * Here we combine them.  into n.  return a pointer to n */

char *n;
n = nw;
string::iterator old = ld.begin ();
/* If old string is empty, just return teh new string */
if (*old == '\0')
  {	return string (nw);	}

/* Till we reach the end of n */
while (*n != '\0')	{
	/* if old reaches end before n: throw error and quit */
	if (*old == '\0')	{
		cerr << "Haplos of different sizes.  OLD (" << ld << ") ended while NEW (" << nw << ") remains\n";
		exit (5);
	  }
	/* if new is blank (underscore), retain the old */
	else if (*n == '_') {	*n = *old;	}

	/* If old is blank (underscore), nothing to do: n already has the base */
	else if (*old == '_')	{	;	}

	/* Neither n nor old is blank.  Check if they are discrepant */
	else if (*n != *old)	{	*n = 'D';	}

	/* if they are not discrepant, both have same base: Again nothing to do*/

	/* increment both n and old for the next char */
	++n;
	++old;
  }

/* Now n is \0.  old also has to be.  If it is - all is well, return nw (the edited new ). */
if (*old == '\0')
  {	return string (nw);	}

/* If not, old still has string remaining.  Not right.  throw error and exit */
else	{
	cerr << "Haplos of diff size.  NEW (" << nw << ") ended while OLD (" << ld << ") remains\n";
	exit (5);
  }
}

void printHCFHeader(FILE* in_fpOut)
{
fprintf(in_fpOut,
"##fileformat=HCFv1.0\n\
##COLUMNS\n\
##CHROM=Chromosome name\n\
##POS=List of coordinated on chromosome (comma separated)\n\
##REF=Sequence at POS coordinates according to the reference genome\n\
##NumSig=Number of significant haplotypes (called based on FDR)\n\
##HAP_Call=List of NumSig haplotypes (posterior probability in paratheses), comma separated\n\
##All_HAP=List of all possible haplotypes (\"posterior probability;FDR\" in paratheses), comma separated(the first NumSig of these were called based on FDR threshold)\n\
##DATA=<ID=nSNP,TYPE=int,Description=\"Number of SNPS in the block\">\n\
###DATA=<ID=nTot,TYPE=int,Description=\"Total Number of Reads overlapping the block\">\n\
###DATA=<ID=nACGT,TYPE=int,Description=\"Number of Reads overlapping the block with an informative base (A, C G or T) at at least one position\">\n\
###DATA=<ID=nBlank,TYPE=int,Description=\"Number of Reads overlapping the block with no information at any position\">\n\
###DATA=<ID=nDisc,TYPE=int,Description=\"Number of Reads overlapping the block where the two reads of a pair disagree on the base at a position\">\n\
###DATA=<ID=nM0,TYPE=int,Description=\"Number of Reads which cover all SNPS in the block\">\n\
###DATA=<ID=nM1,TYPE=int,Description=\"Number of Reads which cover all but 1 SNPS in the block\">\n\
###DATA=<ID=nMk,TYPE=int,Description=\"Number of Reads which cover all but k SNPS in the block\">\n\
###DATA=<ID=nClus,TYPE=int,Description=\"Number of distinct haplotypes which have direct read-based evidence covering all SNPs\">\n");
return;
}

void printHeaderLine (FILE *f, char *sample_name)
{
fprintf (f,
"##CHROM\tPOS\t\
REF\t\
NumSig\t\
HAP_Call\t\
All_HAP\t\
DataForSample=%s\n", sample_name);
return;
}

void printTrackLine (FILE *f, char *sample_name)
{
fprintf (f, "#track name=\"%s\" description=\"Local Haplotype Analysis for %s\" viewLimits=-8:2 color=0,0,255 altColor=255,0,0\n", sample_name, sample_name);
return;
}

void printCommand (FILE *f, int ac, char **av)
{
fprintf (f, "##User_Command=");
int i;
for (i = 0; i < ac; ++i)
  {	fprintf (f, " %s", av[i]);	}
fprintf (f, "\n");
return;
}

void printBamContigs (FILE *f, bam_header_t *b)
{
int i;
fprintf (f, "##LIST_CHROM:List of Chromosome Contigs in the bam file:\n");
for (i = 0; i < b->n_targets; ++i)
  {	fprintf (f, "###%s (Len=%d)\n", b->target_name[i], b->target_len[i]);	}

return;
}

int runBayesModel(HCF_struct* in_hcfVar)
{
	/*
		in_hcfVar: input HCF_struct that contains data as a map
		iRet: output that returns the success (0) or falure (-1).

		This function is the main function to implement the statistical model
		It finds out the probabilities of all the haplotypes and hence the
		final haplotypes after FDR thesholding and fill in the appropriate fields
		of the HCF structure that is returned.
	*/
	int iRet = 0;
		
	if (in_hcfVar->nSNP == 0)
	  	return -1;	

	iRet = fillOptStruct(in_hcfVar);
	if(iRet < 0)
		return iRet;
	

	int nFreqMapSize = in_hcfVar->Freq.size();
	char **data = new char*[nFreqMapSize];
	int *freq = new int[nFreqMapSize];

	for(int i=0;i<nFreqMapSize;i++)
	{
		data[i] = new char[MAX_SNPS+1];
		data[i][MAX_SNPS] = '\0';
	}

	/////////
	int RL = in_hcfVar->nSNP;

	int L1 = pow((double)NUM_BASES,RL);
	int L2 = pow((double)BIN_BASES,RL);

	char** sym = new char*[L1];
	for(int i=0;i<L1;i++)
	{
		sym[i] = new char[MAX_SNPS+1];
		sym[i][MAX_SNPS] = '\0';
	}
	char** binStr = new char*[L2];
	for(int i=0;i<L2;i++)
	{
		binStr[i] = new char[MAX_SNPS+1];
		binStr[i][MAX_SNPS] = '\0';
	}
	int N = dataPopulate(in_hcfVar,data,freq,sym,binStr);
	if(N <= 0)
	{
		//fprintf(stdout,"error in data generation !!!\n\n");
		
		for(int i=0;i<nFreqMapSize;i++)
			delete[] data[i];
		delete[] data;
		delete[] freq;

		for(int i=0;i<L1;i++)
			delete[] sym[i];
		delete[] sym;

		for(int i=0;i<L2;i++)
			delete[] binStr[i];
		delete[] binStr;
		
		iRet = -1;
		return iRet;	
	}

	int hK;
	int* bSet = new int[L1];

	iRet = preScreen(RL,N,data,sym,&hK,bSet);

	if(iRet < 0)
	{	
		//fprintf(stderr,"error in screening the data!!\n");
		
		for(int i=0;i<nFreqMapSize;i++)
			delete[] data[i];
		delete[] data;
		delete[] freq;

		for(int i=0;i<L1;i++)
			delete[] sym[i];
		delete[] sym;

		for(int i=0;i<L2;i++)
			delete[] binStr[i];
		delete[] binStr;

		delete[] bSet;
		
		iRet = -1;	
		return iRet;
	}

	char** hSet = new char*[hK];
	for(int i=0;i<hK;i++)
	{
		hSet[i] = new char[MAX_SNPS+1];
		hSet[i][in_hcfVar->nSNP] = '\0';
	}
	int hCnt = 0;
	for(int i=0;i<L1;i++)
	{
		if(bSet[i] == RL)
		{
			strcpy(hSet[hCnt++],sym[i]);
		}
	}

	double** mat_e = new double*[N];
	int** mat_a = new int*[N];
	int** mat_d = new int*[N];
	int* arr_m = new int[N];

	for(int i=0;i<N;i++)
	{
		mat_e[i] = new double[RL];
		mat_a[i] = new int[hK];
		mat_d[i] = new int[hK];
	}

	buildErrMatrix(N,RL,in_hcfVar->error_rate, mat_e);
	buildMatrices(data,N,RL,hSet,hK,mat_a,mat_d,arr_m);

	double** P_data = new double*[N];
	for(int i=0;i<N;i++)
		P_data[i] = new double[hK];

	int L3 = pow((double)BIN_BASES,hK);
	double* prior_table = new double[L3];

	int** lambdaStr = new int*[L3];
	for(int i=0;i<L3;i++)
		lambdaStr[i] = new int[hK];

	calculatePrDataPriorTable(RL,data,N,hSet,hK,mat_a,mat_d,arr_m,mat_e,sym,P_data,prior_table,lambdaStr);
	//// running the model

	double* pr_z = new double[hK];

	double sum_lambda_str;

	int id;

	//////////////// main model computation /////////////////
	/// GMP part
	// Intialize the mpf variables
	mpf_t r1,r2,r3,r4,r_sumMM;
	mpf_init (r1);
	mpf_init (r2); 
	mpf_init (r3);
	mpf_init (r4);
	mpf_init (r_sumMM);

	mpf_t *MM_gmp = new mpf_t[L3];
	for(int i=0;i<L3;i++)
		mpf_init(MM_gmp[i]);
	//////////////////
	
	mpf_set_d(r_sumMM,0.0);
	for (int i=0;i<L3;i++)
	{
		//lambdaStr_rev = fliplr(lambdaStr(i,:));
		sum_lambda_str = 0.0;
		for(int j=0;j<hK;j++)
			sum_lambda_str += lambdaStr[i][j];

	    if(sum_lambda_str == 0.0)
		{
			for(int j=0;j<hK;j++)
	           	pr_z[j] = 0.0;
		}
	    else
		{
			for(int j=0;j<hK;j++)
	        	pr_z[j] = lambdaStr[i][j]/sum_lambda_str;
	    }

		mpf_set_d (r4, 1.0); 
		for(int ii=0;ii<N;ii++)
	    {
			mpf_set_d(r3, 0.0);
			for(int jj=0;jj<hK;jj++)
			{	
				mpf_set_d (r1, P_data[ii][jj]); 
				mpf_set_d (r2, pr_z[jj]); 
				mpf_mul (r1, r1,r2);
				mpf_add (r3,r3,r1); 
			}
			/// model computation update for using map
			/// multiply freq number of data for one haplotype
			mpf_pow_ui (r1, r3, freq[ii]);
			mpf_mul (r4, r4, r1);
		
		}

		id = bin2dec(lambdaStr[i],hK);
		mpf_set_d (r1, prior_table[id]); 
		mpf_mul (r2, r4, r1); 
		mpf_set(MM_gmp[id],r2);
		mpf_add(r_sumMM,r_sumMM,MM_gmp[id]);
	}

	double *lP = new double[hK];

	int L4 = L3;
	int L5 = L4/2;
	int blockSize;
	int n,startId,tmp_id;
	int *indx0 =  new int[L5];
	//calculation of lP variables
	for(int j=0;j<hK;j++)
	{
		blockSize = L4/2;
	    n = 0;
	    startId = 0;
		mpf_set_d(r1,0.0);
	    while(n < (L3-1) )
		{
	    	for(int ii=0;ii<blockSize;ii++)
	    	{
				indx0[startId+ii] = n+ii;
				tmp_id = indx0[startId+ii];
				mpf_add(r1,r1,MM_gmp[tmp_id]);
			}
	        n = n + 2*blockSize;
	        startId = startId + blockSize;
	    }
		mpf_sub(r2,r_sumMM,r1);
 		mpf_div(r3,r2,r_sumMM);		

		lP[j] = mpf_get_d (r3);
		if(isnan(lP[j]) == 1)
			fprintf(stdout,"NAN error !!!!\n");

	    L4 = blockSize;
	}
	


	///// computation done !!! ///////
	string seqStr;
	for(int i=0;i<hK;i++)
	{
		//hSet[i][in_hcfVar->nSNP] = '\0';
		seqStr = string(hSet[i]);
		in_hcfVar->Haplo_call.push_back(make_pair(lP[i],seqStr));
	}

	////@subhajit_JL ////////////////////
	vector<pair<double,string> > Haplo_call_unsorted(in_hcfVar->Haplo_call);	
	
	//fprintf(stdout,"before sort:\n");
	//for(unsigned int i=0;i<Haplo_call_unsorted.size();i++)
	//	fprintf(stdout,"%s(%5.3f),",Haplo_call_unsorted[i].second.c_str(),Haplo_call_unsorted[i].first);
	//fprintf(stdout,"\n");
	
	sort(in_hcfVar->Haplo_call.begin(), in_hcfVar->Haplo_call.end(),pairCompare); //stable_sort ? @subhajit_JL

	//fprintf(stdout,"after sort:\n");
	//for(unsigned int i=0;i<Haplo_call_unsorted.size();i++)
	//	fprintf(stdout,"%s(%5.3f),",in_hcfVar->Haplo_call[i].second.c_str(),in_hcfVar->Haplo_call[i].first);
	//fprintf(stdout,"\n");
	double *sort_v = new double[hK];

	for(int j=0;j<hK;j++)
		sort_v[j] = in_hcfVar->Haplo_call[j].first;

	int fdr_indx;
	double *d_out = new double[hK];
	fdr_indx = fdr_cumsum_cal(sort_v,hK,d_out);
	in_hcfVar->numSigHaploCall = fdr_indx+1;
	for(int j=0;j<hK;j++)
		in_hcfVar->Haplo_fdr.push_back(d_out[j]);

	//////@subhajit_JL
	int *bin_str = new int[hK];
	int *sig_hap_original_order_bin =  new int[hK];
	for(int i=0;i<hK;i++)
	{
		sig_hap_original_order_bin[i] = 0;
		for(int j=0;j<in_hcfVar->numSigHaploCall;j++)
			if(strcmp(in_hcfVar->Haplo_call[j].second.c_str(),Haplo_call_unsorted[i].second.c_str()) == 0) 
				sig_hap_original_order_bin[i] = 1;

	}
	//fprintf(stdout,"bit string: ");
	//for(int i=0;i<hK;i++)
	//	fprintf(stdout,"%d",sig_hap_original_order_bin[i]);
	//fprintf(stdout,"\tn_sig = %d\n",in_hcfVar->numSigHaploCall);

	int cnt1 = in_hcfVar->numSigHaploCall;
	if(cnt1 >= MIN_LEN)
	{
		int zerocnt = 0;
		for(int i=0;i<hK;i++)
		{	
			if(sig_hap_original_order_bin[i] == 0)
				zerocnt++;
		}		
		int *pos1 = new int[zerocnt];
		int jj=0;
		for(int i=0;i<hK;i++)
		{
			if(sig_hap_original_order_bin[i] == 0)
				pos1[jj++] = i;
		}

		int arr_size = cnt1-MIN_LEN+1;
		int *cnt_arr = new int[arr_size];
		for(int i=0;i<arr_size;i++)
			cnt_arr[i] = 0;

		int *min_str = new int[hK];
		for(int i=0;i<hK;i++)
			min_str[i] = 0;	
		int cnt0 = 0;
		for(int i=hK-1;i>=0;i--)
		{
			if(sig_hap_original_order_bin[i] == 1)
			{
				min_str[i] = sig_hap_original_order_bin[i];
				cnt0++;
			}
			if(cnt0 == MIN_LEN)
				break;
		}
		int MIN_VAL = bin2dec(min_str,hK);
		int MAX_VAL = bin2dec(sig_hap_original_order_bin,hK); 

		//fprintf(stdout,"min = %d; max = %d\n",MIN_VAL,MAX_VAL);
		int cnt3 = 0;
		int *pos2 = new int[hK];
	
 
		for(int i=MIN_VAL;i<=MAX_VAL;i++)
		{
			my_dec2bin(i,hK,bin_str);
			//printINTstring(hK,bin_str);
			int cnt4=0;
			for(int j=0;j<hK;j++)
			{
				pos2[j] = -1;
				if(bin_str[j] == 0)
				{	
					pos2[cnt4] = j;
					cnt4++;
				}
			}	
			int cnt2 = hK-cnt4;
			//printINTstring(zerocnt,pos1);
			//printINTstring(hK,pos2);
			int isSubset_flag = isSubset(pos1,zerocnt,pos2,hK);
			if(isSubset_flag == 1)
			{
				if(cnt2 >= MIN_LEN && cnt2 <= cnt1)
				{
					//for(int j=0;j<hK;j++)
					//	fprintf(stdout,"%d",bin_str[j]);
					//fprintf(stdout," = %d\n",i);
					int tmp_idx = cnt2-MIN_LEN;
					cnt_arr[tmp_idx]++;

					string seqJLStr;
					int cnt6 = 0;				
					for(int k=0;k<hK;k++)
					{
						if(bin_str[k] == 1)
						{	
							cnt6++;
							if(cnt6 > 1)	
								seqJLStr.append(",");
							seqJLStr.append(Haplo_call_unsorted[k].second);
						}
					}

					// this bin string is a candidate string	
					double dd1 = find_all_summable_string_sum(bin_str,hK,cnt4,MM_gmp,L3,r_sumMM);
					//fprintf(stdout,"sum is = %5.6f\n",dd1);
					in_hcfVar->JHL_call.push_back(make_pair(dd1,seqJLStr));
	
					//fprintf(stdout,"--------\n");
					cnt3++;
				}
			}		
		}	

		int KK = in_hcfVar->numSigHaploCall; 
		in_hcfVar->nJL = pow((double)(BIN_BASES),(int)KK) - (int)(0.5*(KK*KK+KK+2));
		if(cnt3 != in_hcfVar->nJL)
		{
			fprintf(stdout,"\tcnt3 = %d nJL = %d\n",cnt3,in_hcfVar->nJL);
			fprintf(stdout,"Error!!!!!\n\n");
		}
		else
		{
			//fprintf(stdout,"\t%d\t",in_hcfVar->nJL);
			sort(in_hcfVar->JHL_call.begin(), in_hcfVar->JHL_call.end(),pairCompare); //stable_sort ? @subhajit_JL
			//for(unsigned int i=0;i<in_hcfVar->JHL_call.size();i++)
			//	fprintf(stdout,"%s(%5.6f);",in_hcfVar->JHL_call[i].second.c_str(),in_hcfVar->JHL_call[i].first);
			//fprintf(stdout,"\n");
		}

		delete[] cnt_arr;
		delete[] min_str;
		delete[] pos1; 
		delete[] pos2; 
	}
	else
	{
		;//fprintf(stdout,"\t0\tNA\n");
	}
	/////////////////////////////////////////@subhajit_JL

	// clearing mpf variables
	mpf_clear(r1);
	mpf_clear(r2);
	mpf_clear(r3);
	mpf_clear(r4);
	for(int i=0;i<L3;i++)
		mpf_clear(MM_gmp[i]);


	/////////// delete allocated memory

	for(int i=0;i<nFreqMapSize;i++)
		delete[] data[i];
	delete[] data;
	delete[] freq;

	for(int i=0;i<L1;i++)
		delete[] sym[i];
	delete[] sym;

	for(int i=0;i<L2;i++)
		delete[] binStr[i];
	delete[] binStr;

	delete[] bSet;
	
	for(int i=0;i<hK;i++)
		delete[] hSet[i];
	delete[] hSet;
	
	delete[] arr_m;
	for(int i=0;i<N;i++)
	{
		delete[] mat_e[i];
		delete[] mat_a[i];
		delete[] mat_d[i];
	}
	delete[] mat_e;
	delete[] mat_a;
	delete[] mat_d;

	for(int i=0;i<N;i++)
		delete[] P_data[i];
	delete[] P_data;

	delete[] prior_table;

	for(int i=0;i<L3;i++)
		delete[] lambdaStr[i];
	delete[] lambdaStr;
	delete[] pr_z;
	
	delete[] MM_gmp;
	
	delete[] lP;
	delete[] indx0;
	delete[] sort_v;
	delete[] d_out;

	//////@subhajit_JL	
	delete[] bin_str;
	delete[] sig_hap_original_order_bin;

	////////////// deallocation done !!!
	return iRet;
}


int dataPopulate(HCF_struct* in_hcfVar,char** out_data,int* out_freq,char** out_sym,char** out_binStr)
{
	/// now data has no_all_missing strings
	int nCnt = 0;
	string tmpStr;
	//int nTotData = 0;
	map<string,int>::iterator freq_it;
	for(freq_it = in_hcfVar->Freq.begin(); freq_it != in_hcfVar->Freq.end();freq_it++ )
	{
		tmpStr.assign(freq_it->first);
		strcpy(out_data[nCnt],tmpStr.c_str());
		out_freq[nCnt] = freq_it->second;
		nCnt++;
	}
	/// data is ready as before
	genSymTab(in_hcfVar->nSNP,out_sym,out_binStr);

	return nCnt;
}

int fillOptStruct(HCF_struct* in_hcfVar)
{
	int iRet = 0;
	if(in_hcfVar->numBlankReads != in_hcfVar->missing_arr[in_hcfVar->nSNP])
	{
		fprintf(stdout,"error in getting missing read sequence!!\n");
		iRet = -1;
		return iRet;
	}
	
	in_hcfVar->distinct_noMiss_cnt = 0;
	string tmpStr;
	map<string,int>::iterator freq_it;
	unsigned int bFlag = 0;
	for(freq_it = in_hcfVar->Freq.begin(); freq_it != in_hcfVar->Freq.end();freq_it++ )
	{
		bFlag = 0;
		tmpStr.assign(freq_it->first);
		for(int j=0;j<in_hcfVar->nSNP;j++)
		{
			if(tmpStr.at(j) == '_' || tmpStr.at(j) == 'N')
			{
				bFlag = 1;
				break;
			}
		}
		if(bFlag == 0)
			in_hcfVar->distinct_noMiss_cnt++;
	}

	return iRet;
}

void printHCFStruct(HCF_struct* in_hcfVar)
{
	fprintf(stdout,"=========================\n");
	fprintf(stdout,"Printing HCF Structure\n");
	fprintf(stdout,"=========================\n");
	fprintf(stdout,"Error rate = %f\n",in_hcfVar->error_rate);
	fprintf(stdout,"FDR threshold = %f\n",FDR_THRESHOLD);
	fprintf(stdout,"sample = %s\n",in_hcfVar->sampleName);	
	fprintf(stdout,"output file = %s.HCF",in_hcfVar->sampleName);
	fprintf(stdout,"chromosome = %s\n",in_hcfVar->chr);
	fprintf(stdout,"SNP = %d\n",in_hcfVar->nSNP);
	fprintf(stdout,"Pos:: ");
	for(int i=0;i<in_hcfVar->nSNP;i++)
		fprintf(stdout,"%ld ",in_hcfVar->pos[i]);
	fprintf(stdout,"\n");
	fprintf(stdout,"total reads = %d\n",in_hcfVar->numBlankReads+in_hcfVar->numInformativeReads);
	fprintf(stdout,"total informative reads = %u\n",in_hcfVar->numInformativeReads);
	fprintf(stdout,"ref = %s\n",in_hcfVar->ref);
	fprintf(stdout,"Haplotype results:: ");
	for(unsigned int j=0;j<in_hcfVar->Haplo_call.size();j++)
		fprintf(stdout,"%s:%f:%f ",in_hcfVar->Haplo_call[j].second.c_str(),in_hcfVar->Haplo_call[j].first,in_hcfVar->Haplo_fdr[j]);
	fprintf(stdout,"\n");
	fprintf(stdout,"significant haplotype call number = %d\n",in_hcfVar->numSigHaploCall);
	fprintf(stdout,"missing info:: ");
	for(int i=0;i<in_hcfVar->nSNP+1;i++)
		fprintf(stdout,"%d:%u ",i,in_hcfVar->missing_arr[i]);
	fprintf(stdout,"\n");
	fprintf(stdout,"distinct no miss seq (from data) = %u\t",in_hcfVar->distinct_noMiss_cnt);
	fprintf(stdout,"\n");
	fprintf(stdout,"PE array length = %d\n",(int)(in_hcfVar->Freq.size()) );
	fprintf(stdout,"PE array:: ");
	map<string,int>::iterator freq_it;
	for(freq_it = in_hcfVar->Freq.begin(); freq_it != in_hcfVar->Freq.end();freq_it++ )
	{
		fprintf(stdout,"%s:%d ",freq_it->first.c_str(),freq_it->second);
	}
	fprintf(stdout,"\n");

	fprintf(stdout,"=========================\n");
	fprintf(stdout,"=========================\n");
}

void printSig2Igv (HCF_struct* in_hcfVar)
{
	if(in_hcfVar->numSigHaploCall > 0)
		print2Igv (in_hcfVar);
	else
		return;
}

void printSigToHCFFile(HCF_struct* in_hcfVar)
{
	if(in_hcfVar->numSigHaploCall > 0)
		printToHCFFile(in_hcfVar);
	else
		return;
}

void print2Igv (HCF_struct* in_hcfVar)
{ /* Prints in_hcfVar structure in igv format */
FILE *OUT = in_hcfVar->OUT;
int nS = in_hcfVar->nSNP, i;

if (nS < 2)
  {	return;	}

fprintf (OUT, "%s\t%ld\t%ld\tREF=%s<br>NumHaplos=%d<br>", in_hcfVar->chr, in_hcfVar->pos[0],
		(nS == 2) ? in_hcfVar->pos[1] : in_hcfVar->pos[nS-1],
		in_hcfVar->ref, in_hcfVar->numSigHaploCall);

if (nS == 3)
  {	fprintf (OUT, "midCrd=%ld<br>", in_hcfVar->pos[1]);	}

fprintf (OUT, "SigHaplos=");
if(in_hcfVar->numSigHaploCall > 0) {
	for(i=0;i<in_hcfVar->numSigHaploCall-1;i++)
		fprintf(OUT, "%s(%5.3f),", in_hcfVar->Haplo_call[i].second.c_str(),in_hcfVar->Haplo_call[i].first);
		fprintf(OUT,"%s(%5.3f)<br>", in_hcfVar->Haplo_call[i].second.c_str(),
								in_hcfVar->Haplo_call[i].first);
	}
else
  {	fprintf (OUT, "NONE_SIGNIFICANT<br>All_Haplos=");	}
int possibleHCnt = pow((double)(BIN_BASES),(int)in_hcfVar->nSNP);
for(i=0;i<possibleHCnt-1;i++)	{
	fprintf(OUT,"%s(%5.3f;%5.3f),", in_hcfVar->Haplo_call[i].second.c_str(),
			in_hcfVar->Haplo_call[i].first,
			in_hcfVar->Haplo_fdr[i]);
  }
fprintf(OUT,"%s(%5.3f;%5.3f)<br>", in_hcfVar->Haplo_call[i].second.c_str(),
			in_hcfVar->Haplo_call[i].first,
			in_hcfVar->Haplo_fdr[i]);

fprintf (OUT, "nSNP=%d<br>nTot=%d<br>nACGT=%d<br>nBlank=%d<br>nDisc=%d<br>", in_hcfVar->nSNP,
		in_hcfVar->numInformativeReads + in_hcfVar->numBlankReads + in_hcfVar->numDiscrepantReads,
		in_hcfVar->numInformativeReads, in_hcfVar->numBlankReads, in_hcfVar->numDiscrepantReads);
for(i=0;i<=in_hcfVar->nSNP;i++)
  {	fprintf(OUT,"nM%d=%d<br>", i, in_hcfVar->missing_arr[i]);	}

fprintf(OUT,"nClus=%d\t", in_hcfVar->distinct_noMiss_cnt);
fprintf(OUT,"%d\n", (in_hcfVar->numSigHaploCall <= 2) ? in_hcfVar->numSigHaploCall : (-1 * in_hcfVar->numSigHaploCall));

return;
}

void printToHCFFile(HCF_struct* in_hcfVar)
{
	int i=0;
	FILE *in_fpOut = in_hcfVar->OUT;
	char positions[MAX_LINE_LEN];

	if (in_hcfVar->nSNP < 2)
	  {	return;	}

if (in_hcfVar->nSNP == 2)
  {	sprintf (positions, "%ld,%ld", in_hcfVar->pos[0], in_hcfVar->pos[1]);	}
if (in_hcfVar->nSNP == 3)
  {	sprintf (positions, "%ld,%ld,%ld", in_hcfVar->pos[0], in_hcfVar->pos[1], in_hcfVar->pos[2]);	}


fprintf(in_fpOut, "%s\t%s\t%s\t%d\t", in_hcfVar->chr, positions, in_hcfVar->ref, in_hcfVar->numSigHaploCall);

if(in_hcfVar->numSigHaploCall > 0) {
		for(i=0;i<in_hcfVar->numSigHaploCall-1;i++)
			fprintf(in_fpOut,"%s(%5.3f),",in_hcfVar->Haplo_call[i].second.c_str(),in_hcfVar->Haplo_call[i].first);
		fprintf(in_fpOut,"%s(%5.3f)\t", in_hcfVar->Haplo_call[i].second.c_str(),
								in_hcfVar->Haplo_call[i].first);
	}
else
  {	fprintf (in_fpOut, "NONE_SIGNIFICANT\t");	}
int possibleHCnt = pow((double)(BIN_BASES),(int)in_hcfVar->nSNP);
for(i=0;i<possibleHCnt-1;i++)	{
	fprintf(in_fpOut,"%s(%5.3f;%5.3f),", in_hcfVar->Haplo_call[i].second.c_str(),
			in_hcfVar->Haplo_call[i].first,
			in_hcfVar->Haplo_fdr[i]);
  }
fprintf(in_fpOut,"%s(%5.3f;%5.3f)\t", in_hcfVar->Haplo_call[i].second.c_str(),
			in_hcfVar->Haplo_call[i].first,
			in_hcfVar->Haplo_fdr[i]);

fprintf (in_fpOut, "nSNP=%d;nTot=%d;nACGT=%d;nBlank=%d;nDisc=%d;", in_hcfVar->nSNP,
		in_hcfVar->numInformativeReads + in_hcfVar->numBlankReads + in_hcfVar->numDiscrepantReads,
		in_hcfVar->numInformativeReads, in_hcfVar->numBlankReads, in_hcfVar->numDiscrepantReads);
for(i=0;i<=in_hcfVar->nSNP;i++)
  {	fprintf(in_fpOut,"nM%d=%d;", i, in_hcfVar->missing_arr[i]);	}

fprintf(in_fpOut,"nClus=%d;", in_hcfVar->distinct_noMiss_cnt);
printJLToHCFFile(in_hcfVar);
//printDataToHCFFile(in_hcfVar);
fprintf(in_fpOut,"\n");

}

void printDataToHCFFile(HCF_struct* in_hcfVar)
{
	FILE *in_fpOut = in_hcfVar->OUT;
	///extra: print data set
	fprintf(in_fpOut,"\tPE_length=%d\t",(int)(in_hcfVar->Freq.size()) );
	map<string,int>::iterator freq_it;
	for(freq_it = in_hcfVar->Freq.begin(); freq_it != in_hcfVar->Freq.end();freq_it++ )
	{
		fprintf(in_fpOut,"%s=%d;",freq_it->first.c_str(),freq_it->second);
	}
}

void printJLToHCFFile(HCF_struct* in_hcfVar)
{
	FILE *in_fpOut = in_hcfVar->OUT;
	///extra: print data set
	fprintf(in_fpOut,"\tJL_length=%d\t",(int)(in_hcfVar->nJL) );

	if(in_hcfVar->nJL > 0)
	{
		for(unsigned int i=0;i<in_hcfVar->JHL_call.size();i++)
		{
			fprintf(in_fpOut,"%s(%5.6f);", in_hcfVar->JHL_call[i].second.c_str(),
										   in_hcfVar->JHL_call[i].first);
		}
	}
	else
	{
		fprintf(in_fpOut,"NO_JL_CALL;");
	}
}

int preScreen(int in_RL,int in_N,char** in_data,char** in_sym,int* out_hK,int *out_bSet)
{
	int iRet = 0;
	char s[NUM_BASES]={'A','C','G','T'};

	int** base_freq = new int*[NUM_BASES];
	for(int k=0;k<NUM_BASES;k++)
		base_freq[k] = new int[in_RL];

	char** data_one_position = new char*[in_RL];
	for(int j=0;j<in_RL;j++)
		data_one_position[j] = new char[in_N];

	int* ids = new int[in_N];
	int *cnt_h = new int[in_RL];

	baseTable* baseTbl = new baseTable[in_RL];

	int nHaplo,t_sum,id_cnt;

	nHaplo = 1;

	for(int j=0;j<in_RL;j++)
	{
		for(int i=0;i<in_N;i++)
			data_one_position[j][i] = in_data[i][j];

		t_sum = 0;
		id_cnt = 0;
		baseTbl[j].cnt = 0;
		for(int k=0;k<NUM_BASES;k++)
		{
			base_freq[k][j] = findID_Str(in_N,data_one_position[j],s[k],ids);
			if(base_freq[k][j] == 0)
				t_sum++;
			else
			{
				baseTbl[j].cnt++;
				baseTbl[j].id[id_cnt++] = s[k];
			}
		}

		cnt_h[j] = NUM_BASES - t_sum;
		nHaplo *= cnt_h[j];
	}

	*out_hK = nHaplo;

	if(nHaplo != pow((double)BIN_BASES,in_RL) )
	{
		//fprintf(stdout,"not every read position has 2 variants !!!\n");
		iRet = -1;			
	}

	int L1 = pow((double)NUM_BASES,in_RL);

	int bFlag = 0;

	for(int i=0;i<L1;i++)
	{
		out_bSet[i] = 0;
		for(int j=0;j<in_RL;j++)
		{
			bFlag = 0;
			for(int k=0;k<baseTbl[j].cnt;k++)
			{
				if(in_sym[i][j]==baseTbl[j].id[k])
				{
					bFlag = 1;
					break;
				}
			}
			if(bFlag == 1)
				out_bSet[i] = out_bSet[i]+1;
		}
	}

	//// deallocation
	for(int j=0;j<NUM_BASES;j++)
		delete[] base_freq[j];
	delete[] base_freq;

	for(int j=0;j<in_RL;j++)
		delete[] data_one_position[j];
	delete[] data_one_position;

	delete[] ids;
	delete[] cnt_h;
	delete[] baseTbl;

	return iRet;
}

void calculatePrDataPriorTable(int in_RL,char** in_data,int in_N,char** in_hSet,int in_hK,
		int** in_mat_a,int** in_mat_d,int* in_arr_m,double** in_mat_e,char** in_sym,
		double** out_P_data,double* out_prior_table,int** out_lambdaStr)
{
	/*
		This function does all the necessary calculation 
		for the prior tables, prior error calculation in a matrix format		 
	*/

	char* y = new char[MAX_SNPS+1];
	y[MAX_SNPS] = '\0';
	char* y_tmp = new char[MAX_SNPS+1];
	y_tmp[MAX_SNPS] = '\0';

	int *ids = new int[in_RL];
	int L1,strt_indx;
	double p;
	double p_sum;
	int M;

	for(int i=0;i<in_N;i++)
	{
		strcpy(y,in_data[i]);

		M = find_missing_positions(y,in_RL,ids); ////@subhajit: 07.16.14
		L1 = pow((double)NUM_BASES,in_arr_m[i]);

		for(int l=0;l<in_hK;l++)
		{
			p_sum = 0.0;
			for(int jj=0;jj<L1;jj++)
			{
				strcpy(y_tmp,y);
				strt_indx = in_RL - in_arr_m[i];
				for(int ll=0;ll<in_arr_m[i];ll++)
					y_tmp[ids[ll]] = in_sym[jj][strt_indx+ll];
				p = 1;
				for(int r=0;r<in_RL;r++)
				{
					if(y_tmp[r] == in_hSet[l][r])
						p = p*(1-in_mat_e[i][r]);
					else
						p = p*(in_mat_e[i][r]/3);
				}
				p = p*pow(0.5,M); ///@subhajit
				p_sum = p_sum+p;
			}

			out_P_data[i][l] = p_sum;
		}
	}

	//// creating a table for prior
	int *a_dot_col = new int[in_hK];
	int *d_dot_col = new int[in_hK];
	for(int j=0;j<in_hK;j++)
	{
		a_dot_col[j] = 0;
		for(int i=0;i<in_N;i++)
		{
			a_dot_col[j] += in_mat_a[i][j];
			d_dot_col[j] += in_mat_d[i][j];

		}
	}

	int T = buildBinaryStrings(in_hK,out_lambdaStr);

	double* prior_lambda_1 = new double[in_hK];
	double* prior_lambda_0 = new double[in_hK];

	for(int j=0;j<in_hK;j++)
	{
		prior_lambda_1[j] = PRIOR_PROB; //// using beta-binomial prior beta(0.05,1) bernoulli  with 1
		prior_lambda_0[j] = 1-prior_lambda_1[j];
	}


	for(int i=0;i<T;i++)
	{
		out_prior_table[i] = 1.0;
		for(int j=0;j<in_hK;j++)
		{
			if(out_lambdaStr[i][j] == 1)
				out_prior_table[i] *= prior_lambda_1[j];
			else
				out_prior_table[i] *= prior_lambda_0[j];
		}
	}


	delete[] y;
	delete[] y_tmp;
	delete[] ids;
	delete[] a_dot_col;
	delete[] d_dot_col;
	delete[] prior_lambda_0;
	delete[] prior_lambda_1;

}

void genSymTab(int L,char** str_arr,char** binStr)
{
	/// length of the sequence
    /// this length determines the symTab entries
    int B = NUM_BASES;
    /// base is 4
    /// A = 0; C = 1; G = 2; T = 3; base are sorted

	for(int i=0;i< pow((double)B,L);i++)
	{
		convertToBase(i,L,B,str_arr[i]);
		for(int j=0;j<L;j++)
		{
			switch(str_arr[i][j])
			{
				case '0': str_arr[i][j] = 'A'; break;
				case '1': str_arr[i][j] = 'C'; break;
				case '2': str_arr[i][j] = 'G'; break;
				case '3': str_arr[i][j] = 'T'; break;
			}
		}
		str_arr[i][L] = '\0'; // terminator needed 
	}

	B = BIN_BASES;
	for(int i=0;i< pow((double)B,L);i++)
	{
		convertToBase(i,L,B,binStr[i]);
	}

}

void convertToBase(int in_v, int in_L,int in_B,char* out_str)
{
	int x,y;char c;
	for(int i=0;i<in_L;i++)
	{
		out_str[i] = '0';
	}
	for(int i=0;i<in_L;i++)
	{
		x = floor(in_v/in_B);
		y = in_v%in_B;
		c = char(y+'0');
		out_str[in_L-1-i] = c;
		in_v = x;
	}
}

int buildBinaryStrings(int M,int **lambdaStr)
{
	/* 
		to get all the binary strings from 0 to (2^hK)-1
	*/
	int d,cnt;
	int L = pow(2.0,M);

	for(int i=0;i<L;i++)
	{
		cnt = 0;
   		for (int c=M-1;c >= 0;c--)
   		{
      		d = i >> c;
 			if ( d & 1 )
         		lambdaStr[i][cnt] = 1;
      		else
         		lambdaStr[i][cnt] = 0;

      		cnt++;
   		}
	}
	return L;
}

void buildErrMatrix(int in_N,int in_RL,double in_err, double** out_mat_e)
{
	/*
		build error matrix
	*/
	for(int i=0;i<in_N;i++)
		for(int j=0;j<in_RL;j++)
			out_mat_e[i][j] = in_err;

}

void buildMatrices(char** in_data,int in_N,int in_RL,char** in_hSet,int in_hK,
		int** out_mat_a,int** out_mat_d, int* out_arr_m)
{
	/*
		 construct match and mismatch matrix AND missing array
    */
	for(int i=0;i<in_N;i++)
	{
		out_arr_m[i] = 0;
		for(int j=0;j<in_RL;j++)
		{
			if(in_data[i][j] == '_')
				out_arr_m[i] += 1;
		}
	}

	for(int i=0;i<in_N;i++)
	{
		for(int j=0;j<in_hK;j++)
		{
			out_mat_a[i][j] = 0;
			out_mat_d[i][j] = 0;
			for(int k=0;k<in_RL;k++)
			{
				if(in_data[i][k] == in_hSet[j][k])
					out_mat_a[i][j] += 1;
			}
			out_mat_d[i][j] = in_RL - out_mat_a[i][j] - out_arr_m[i];
		}
	}
}

int bin2dec(int* in_bStr,int in_hK)
{
	/* 
		find decimal equivalent of an input binary string
    */
	int ret = 0;
    int val = 1;

    for(int j=in_hK-1;j>=0;j--)
	{
       if (in_bStr[j] == 1)
			ret = ret + val;
       val = val*2;
    }
    return ret;
}

void indices_minus_j(int in_hK,int in_j,int* out_ids)
{
	int cnt = 0;
	for(int i=0;i<in_hK;i++)
	{
		if(i != in_j)
			out_ids[cnt++] = in_hK-1-i;
	}
}

void fillLambdaStr(int in_hK,int* in_ids,int* in_otherLambdaStr,int* out_lambdaStr)
{
	for(int i=0;i<(in_hK-1);i++)
	{
		out_lambdaStr[in_ids[i]] = in_otherLambdaStr[in_hK-2-i];
	}
}


/***** few utility functions *****/
int findID_Str(int in_L, char* in_arr_a, char in_elem, int* out_id)
{
	int cnt = 0;

	for(int i=0;i<in_L;i++)
	{
		if(in_arr_a[i] == in_elem)
		{
			out_id[cnt++] = i;
		}
	}

	return cnt;
}

///// @subhajit: 07.16:14
int find_missing_positions(char* in_s,int in_RL,int* out_ids)
{
	int cnt = 0;
	int M = 0;
	for(int j=0;j<in_RL;j++)
	{
		if(in_s[j] == '_')
		{
			out_ids[cnt++] = j;
			M++;
		}
	}
	return M;
}

int fdr_cumsum_cal(double *in_d,int in_K,double* out_d)
{
	int indx = -1;

	for(int i=0;i<in_K;i++)
	{
		out_d[i] = (1-in_d[i]);
		for(int j=i-1;j>=0;j--)
			out_d[i] += (1-in_d[j]);

		out_d[i] = out_d[i]/(i+1);
	}

	for(int i=0;i<in_K;i++)
		if(out_d[i] < FDR_THRESHOLD)
			indx = i;

	return indx;
}

bool pairCompare(const pair<double, string>& A, const pair<double, string>& B)
{
  return (A.first > B.first);
}

double calcGammaFactor(double a,double b)
{
	/* calculation of gammaFactor once in the beginning */
	double dRet;
	int n = 1;
	dRet = (tgamma(1.0+a)*tgamma(b)*tgamma(a+b))/(tgamma(n+a+b)*tgamma(a)*tgamma(b));
	return dRet;
}


////// new functions for JL part
void my_dec2bin(int n,int L,int *bin_str)
{
	for(int i=0;i<L;i++)
		bin_str[i] = 0;

    int oneorzero;
    for(int i=L-1;i>=0;i--) 
	{
    	oneorzero = n % 2;
        if (oneorzero == 1) 
        	bin_str[i] = 1;
        else 
        	bin_str[i] = 0;
            
        n /= 2;   
    }
}

int isSubset(int *a,int L1, int *b,int L2)
{
	int iRet = 0;
	int cnt = 0;

	for(int i=0;i<L1;i++)
	{
		for(int j=0;j<L2;j++)
			if(a[i] == b[j])
				cnt++;
	}

	if(cnt == L1)
		iRet = 1;

	return iRet;
}

double find_all_summable_string_sum(int *bin_str,int hK,int zero_cnt,mpf_t *MM_gmp,int L3,mpf_t r_sumMM)
{
	double sum_d = 0.0;
	//mpf_t r_sumMM;
	//mpf_init (r_sumMM);
	//mpf_set_d(r_sumMM,0.0);

	mpf_t sum_d_gmp;
	mpf_init (sum_d_gmp);
	mpf_set_d(sum_d_gmp,0.0);

	mpf_t res_d_gmp;
	mpf_init (res_d_gmp);
	mpf_set_d(res_d_gmp,0.0);


	int *pos3 = new int[zero_cnt];
	int idx_size = pow((double)(BIN_BASES),(int)zero_cnt);
	int j=0;
	for(int i=0;i<hK;i++)
		if(bin_str[i] == 0)	
			pos3[j++] = i;
	
	int L1 = zero_cnt;	
	int *tmp_bin_str = new int[L1];
	int *new_bin_str = new int[hK];
	for(int i=0;i<hK;i++)
		new_bin_str[i] = bin_str[i];

	int cnt5 = 0;
	//printINTstring(hK,bin_str);
	//printINTstring(zero_cnt,pos3);
	//for(int i=0;i<L3;i++)
	//{
	//	mpf_add(r_sumMM,r_sumMM,MM_gmp[i]);
	//}
	for(int i=0;i<idx_size;i++)
	{
		my_dec2bin(i,L1,tmp_bin_str);
		//printINTstring(L1,tmp_bin_str);
		for(int j=0;j<L1;j++)
			new_bin_str[pos3[j]] = tmp_bin_str[j];
	
		//fprintf(stdout,"\t");
		//for(int j=0;j<hK;j++)
		//	fprintf(stdout,"%d",new_bin_str[j]);

		int tmp_idx = bin2dec(new_bin_str,hK);
		//fprintf(stdout," = %d\n",tmp_idx);
		mpf_add(sum_d_gmp,sum_d_gmp,MM_gmp[tmp_idx]);
		cnt5++;
	}

	//fprintf(stdout,"\n\tneed to sum %d strings\n",cnt5);

	delete[] pos3;
	delete[] tmp_bin_str;
	delete[] new_bin_str;

	mpf_div(res_d_gmp,sum_d_gmp,r_sumMM);
	sum_d = mpf_get_d(res_d_gmp);
	return sum_d;

}			

void printINTstring(int hK,int *S)
{
	fprintf(stdout,"\n");
	for(int i=0;i<hK;i++)
		fprintf(stdout,"%d",S[i]);
	fprintf(stdout,"\n");
}

