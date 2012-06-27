/* Small tool for simulating sequence reads from a reference genome.
 * It writes the true read coordinates as well as the number of
 * polymorphisms and sequencing errors in the read name.
 * Compilation
 * =====================================================
 * gcc -o simulator -lz -lm -std=gnu99 simulator_v1.c
 * Usage
 * =====================================================
 * Command line for simulation:
 * ./simulator [<parameters>] <in_seq1.fa> <out_file1.fq> <out_file2.fq>
 */

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include "kseq.h"
#include "genran.h"
#include <math.h>
#include <zlib.h>
#include <stdint.h>
#include <unistd.h>

#define MAX 3
#define MIN 0
#define PACKAGE_VERSION "0.0.1"

KSEQ_INIT(gzFile, gzread);

static double ERR_RATE = 0.02; //base calling error rate - default
static double MUT_RATE = 0.001; //mutation rate - default
static double INDEL_FRAC = 0.15; //rate of indels - default
static int GAP_SIZE = 1;//size of gap - default
static int SEED = -1;//seed for random number generator - default
static double INDEL_EXT = 0.3;//probability that new random nuclotides are inserted into fragment - default

typedef unsigned short mut_t;

typedef struct{
	int l,m;
	mut_t *s;
} mutseq_t;

char *tot_seq,*read_f,*iter;
int n_sub[2],n_err[2],n_ind[2];
gzFile fp,fpc;

 static double_t frequency_A;
 static double_t frequency_T;
 static double_t frequency_G;
 static double_t frequency_C;
 static double_t frequency_N;
 
 static uint64_t nA;
 static uint64_t nT;
 static uint64_t nG;
 static uint64_t nC;
 static uint64_t nN;

 
/* ran_normal() - normal distribution approximation
 * this function is copied from genran.c
 * normally disturbed random number is returned*/
 
 
double ran_normal()
{ 
    static int iset = 0; 
    static double gset; 
    double fac, rsq, v1, v2; 
    if (iset == 0) {
        do { 
            v1 = 2.0 * ran_uniform() - 1.0;
            v2 = 2.0 * ran_uniform() - 1.0; 
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq); 
        gset = v1 * fac; 
        iset = 1;
        return v2 * fac;
    } else {
        iset = 0;
        return gset;
    }
}


/*poisson_random_number() - Poisson distribution approximation
 * input parameter is lambda (number of base calling errors)
 * number of iterations is set to 1000
 * random number is returned*/


const int poisson_random_number(const double lambda) //Poisson distribution
{
  int k=0;                        
  const int max_k = 1000;           
  double p = ran_uniform(); 
  double P = exp(-lambda);          
  double sum=P;                     
  if (sum>=p) return 0;             
  for (k=1; k<max_k; ++k) {         
           P*=lambda/(double)k;            
           sum+=P;                         
           if (sum>=p) break;              
  }
  return k;                         
}


/* swap_base() - changes the value of the base given as an input parameter
 * value of the new random base is calculated with uniform distribution 
 * (function drand48()) the new value of the base is returned or
 * the same base (if value of input parameter was N - unkown nucleotide)
 * */


char swap_base(char base){
	char new_base;
	double value;
	value = drand48()*MAX; //value is a number between 0 and 4
	switch (base){
		case 'A':
		     if(value >= 1) return 'C';
		     if((value > 1) && (value <= 2)) return 'G';
		     else return 'T';
		     break;
	    case 'T':
	         if(value >= 1) return 'C';
	         if((value > 2) && (value <= 3)) return 'G';
		     else return 'A';
		     break;
	    case 'G':
	         if(value >= 1) return 'C';
	         if((value > 2) && (value <= 3)) return 'T';
		     else return 'A';
		     break;
		case 'C':
		     if(value >= 1) return 'C';
	         if((value > 2) && (value <= 3)) return 'G';
		     else return 'A';
		     break;
	    default:
	         return base; 
	   }	
} 


/* get_complement - calculates complement of the given parameter (nucleotide)
 * A - T
 * T - A
 * G - C
 * C - G
 * returns complement or the input base if value of base was N - unkown nucleotide
 * */


char get_complement(char base){ //get complement A-T G-C
	if (base == 'A') return 'T';
	else if (base == 'T') return 'A';
	else if (base == 'G') return 'C';
	else if (base == 'C') return 'G';
	else return base;
}


/*simulate_BCER() - simulates base calling errors, input parameters are:
 * read_lenght - total length of the read
 * read_local - read in which errors would be added
 * n_errors - total number of errors calculated with poisson distribution
 * returns read_local with errors added*/


char *simulate_BCER(int read_length, char *read_local, int *n_errors){
	double value;
	int i,index;
	*n_errors=(int)(poisson_random_number((read_length*ERR_RATE))); //calculates total number of errors
	for (i=0;i<*n_errors;i++){
		value = drand48()*(read_length - 1);
		index = trunc((int)value); //calculates index. the base on this index will be randomly changed
		*(read_local+index) = swap_base(*(read_local+index)); //change base
    }
    return read_local;
}


/*generate mutations() - generates mutations on input string, input parameters are:
 * dist - distance between two reads
 * size_l - size of the first read
 * size_r - size of the second read
 * function makes mutations on orginial string, so nothing is returned
 * */
  

void generate_mutations(int dist, int size_l,int size_r){
	int i,pos,n_mut;
	double r;
	n_mut = dist * MUT_RATE; //calculate the total number of substitutions
	for (i=0;i<n_mut;i++){
		r=drand48(); //uniform distribution r = [0,1>
		pos = (int)(trunc(r * dist));
		if(pos<=size_l){n_sub[0]++;} //calculate index for mutation
		if(pos>=(dist-size_r)){n_sub[1]++;}
		if(*(read_f+pos)!='N') //if value of base is N - skip!
		    *(read_f+pos)=swap_base(*(read_f+pos));
	}
}


/*generate_gaps() - generates gaps on fragment by adding \0 character to it
 * input parameters are:
 * gap_pos - position of the gap
 * gap size - size of gap (number of nucleotides that will be overwritten)
 * */


void generate_gaps(int gap_pos, int gap_size){
	uint64_t i;
	read_f[gap_pos]='\0';
	strcat(read_f,(read_f + gap_pos + gap_size)); //overwrites the gap_size nucleotides beginning from the gap_pos 
}


/*core() - main function, calls every other function, if
 * everything was fine during the simulation process returns 0 otherwise
 * -1 , input parameters are:
 * fout1 - file in which the fastq result of simulation for the first read will be stored
 * fout2 - file in which the fastq result of simulation for the second read will be stored
 * std_dev - standard deviation for fragments generation
 * size_l - size of the first read
 * size_r - size of the second read
 * N - total number of reads
 * dist - distance between two reads
 * */


int core(FILE *fout1,FILE *fout2,char *argv, int std_dev, int size_l, int size_r, uint64_t N,int dist){ //most important function
	mutseq_t *ret[2];
	gzFile fp;
	uint64_t total_len;
	kseq_t *seq;
	int l,n_ref,max_size,Q,n_errors;
	uint64_t i,j,counter_a,counter_b;
	char *q_string,*q2_string;
	fp = gzopen(argv, "r");
	seq = kseq_init(fp);
	total_len = n_ref = 0;
	frequency_A=frequency_T=frequency_G=frequency_C=0.;
	nA=nT=nG=nC=j=0;
	max_size = size_l > size_r? size_l : size_r; //determine maximum size of the read
	Q = (ERR_RATE == 0.0)? 'I' : (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33; //calculates quality score
	q_string = (char *)malloc((size_l+1)*sizeof(char)); //first quality string
	q2_string = (char *)malloc((size_r+1)*sizeof(char)); // second quality string
	for(int k=0;k<size_l;k++){q_string[k]=Q;};//prepare quality strings
	for(int k=0;k<size_r;k++){q2_string[k]=Q;};
	fprintf(stderr, "[%s] calculating the total length of the sequnce...\n",__func__);
	
	while ((l = kseq_read(seq)) >= 0){ //prints out basic info and calculates the total length of sequence
		     printf("[%s] name: %s\n",__func__,seq->name.s);
		     if(seq->comment.l) printf("[%s] comment: %s\n",__func__,seq->comment.s);
		     total_len+=l;
		     ++n_ref;
		     if (seq->qual.l) printf ("qual: %s\n",seq->qual.s);
	}
	
	if (total_len == 0){ //in case of empty file simulation is terminated, -1 is returned
		     printf("[%s] input file is empty!\n",__func__);
		     return -1;
	}
	
	fprintf(stderr, "[%s] %d sequences, total length: %llu\n", __func__, n_ref, (long long)total_len);
	kseq_destroy(seq);
	gzclose(fp);
	fp = gzopen(argv, "r");
	seq = kseq_init(fp);
	
	
	while ((l = kseq_read(seq)) >= 0){
		if (l < dist + 3*std_dev){
			fprintf(stderr,"[%s] ERROR sequence to short for given parametars!\n",__func__); //check if sequence is too short for given parameters
			return -1;
		}
		for (i=0;i<l;i++){ //calculate the total number of A,T,G,C
		if(seq->seq.s[i]== 'A')
		     nA++;
		else if(seq->seq.s[i] == 'T')
		     nT++;
		else if(seq->seq.s[i] == 'G')
		     nG++;
		else if(seq->seq.s[i] == 'C')
		     nC++;
		else 
		     nN++;
    } 
	frequency_A = (double_t)nA/l; //calculates frequency of nucleotides
	frequency_C = (double_t)nC/l;
	frequency_G = (double_t)nG/l;
	frequency_T = (double_t)nT/l;
	frequency_N = (double_t)nN/l;
	j++;	  
    printf("[%s] frequency per sequence [%llu/%llu] A - %f | C - %f | G - %f | T - %f | *N - %f - unknown nucleotides (percentage) \n",__func__,(long long)j,(long long)n_ref,frequency_A*100,frequency_C*100,frequency_G*100,frequency_T*100,frequency_N*100);
    }
    
    
    //tot_seq = (char *)malloc((total_len+1)*sizeof(char)); //allocates enough memory space for sequence
	tot_seq = seq->seq.s;
	printf("[%s] transferring sequence into memory and generating errors...\n",__func__);
	counter_a=counter_b=0; //set counters to zero
	for(i=0;i<N;i++){
		double ran;
		int d,pos;
		uint64_t fp1_b,fp1_e,fp2_b,fp2_e,begin,end;
		char *read1,*read2;
		n_sub[0]=n_sub[1]=n_err[0]=n_err[1]=n_ind[0]=n_ind[1]=0; //set the number of substitutions, errors and indels to zero (this would be printed in header of the output files)
		
		do{
			ran = ran_normal();
			ran = ran * std_dev + dist;
			d = (int)(ran + 0.5); //calculates size of fragment
			d = d > max_size ? d : max_size; //avoid boundary failure
			pos = (int)((total_len-d+1)*drand48()); //calculates the position (first index) of the fragment with normal distribution
		}while(pos < 0 || pos >= total_len || (pos + d - 1) >= total_len);
		
		read_f = (char *)malloc((d+1+(int)((INDEL_EXT+0.5)*d))*sizeof(char));counter_a++; //allocates enough memory for the fragment
		//read_f = (char *)malloc(2000*sizeof(char));
		read_f[d+1]='\0'; //adds \0 to the end of the fragment
		strncpy(read_f,tot_seq+pos,d); //copies part of the sequence into fragment
		begin=pos; end=pos+d; //store the beginning and the ending of the fragment
		generate_mutations(d,size_l,size_r); //generates mutations
		int n_n=0;
		int n_indel = (int)(INDEL_FRAC * d); //total number of indels
		int curr_dist=d; //distance between two reads
		if((n_indel >= d) || ((GAP_SIZE*n_indel)>=d)){
			fprintf(stderr,"[%s] ERROR sequence too short for given parametars!\n",__func__); //checks if sequence is too short for given parameters
			return -1;
		}
		do{
			int type_indel,gap_size;
			char *fragment;
			type_indel = (drand48()>=INDEL_EXT)?1:0; //if 0 generates gaps otherwise adds new random nucleotides
			if (type_indel == 1){
				pos = (int)trunc(drand48()*(curr_dist-1)); //calculates the position of gap
				if (pos<size_l) n_ind[0]++;
				if (pos >= (d-size_r)) n_ind[1]++;
				gap_size = poisson_random_number(GAP_SIZE); //generates gaps
				if (gap_size){
			        generate_gaps(pos,gap_size);
				}
			    curr_dist=curr_dist-gap_size; //reduce size of fragment
		    }
			else{
				char fragment[4]={0x0}; //new nucleotides will be stored here
				char keeper[1000]={0x0}; //temporary string
				int pos;
				double r;
				char base;
				r = drand48();
				base=(r < 0.25)?'A':((r>=0.25 && r<0.5)?'T':(r>=0.25 && r<0.75)?'C':'G'); //generates new random nucleotide
				fragment[0]=base;fragment[1]='\0'; //adds nucleotide to fragment
				pos = (int)trunc(drand48()*(curr_dist-1)); //calculates the position for insertion
				if (pos<size_l) n_ind[0]++;
				if (pos >= (d-size_r)) n_ind[1]++;
				strcat(keeper,read_f+pos); 
				read_f[pos]='\0';
				strcat(read_f,fragment);
				strcat(read_f,keeper);//inserts nucleotides
				curr_dist++;
				
			}
			n_n++;
		}while(n_n<n_indel);
		read1=(char *)malloc((size_l+1)*sizeof(char)); //generates two pair end reads 
		n_errors=(int)(INDEL_FRAC*size_l); //calculates the total number of base calling errors
		read2=(char *)malloc((size_r+1)*sizeof(char));
		strncpy(read1,read_f,size_l);read1[size_l+1]='\0';
		int internal_counter=0;
		for(int k=(curr_dist-1);internal_counter<size_r;k--){
			*(read2+internal_counter)=get_complement(*(read_f + k)); //makes complement of first fread
			internal_counter++;
		}
		read2[internal_counter]='\0';
		int err_1,err_2;
		read1=simulate_BCER(size_l,read1,&err_1); //generates base calling errors
		read2=simulate_BCER(size_r,read2,&err_2);
		fprintf(fout1,"@%s_%llu_%llu_%d:%d:%d_%d:%d:%d_%llx/%d\n",seq->name.s,(long long)begin,(long long)end,err_1,n_sub[0],n_ind[0],(int)(size_r*ERR_RATE),n_sub[1],n_ind[1],(long long)counter_a,1); //prints results into the file
		fprintf(fout1,"%s\n+\n%s\n",read1,q_string);
		fprintf(fout2,"@%s_%llu_%llu_%d:%d:%d_%d:%d:%d_%llx/%d\n",seq->name.s,(long long)begin,(long long)end,err_2,n_sub[0],n_ind[0],(int)(size_r*ERR_RATE),n_sub[1],n_ind[1],(long long)counter_a,2);
		fprintf(fout2,"%s\n+\n%s\n",read2,q2_string);
		free(read1);
		free(read2);
		free(read_f);
	}
    kseq_destroy(seq);
	gzclose(fp);
	free(q_string);
	free(q2_string);
}


/* simu_usage() - displays menu on standard error*/


static int simu_usage(){
	fprintf(stderr,"**********************************************************\n");
	fprintf(stderr,"Program: simulator (short read simulator)\n");
	fprintf(stderr,"Version %s\n",PACKAGE_VERSION);
	fprintf(stderr,"\nUsage: ./simulator [options] <in_seq.fa> <out_read1.fq> <out_read2.fq>\n\n");
	fprintf(stderr,"Options: -r FLOAT rate of mutations [default 0.001]\n");
	fprintf(stderr,"         -e FLOAT base error rate [default 0.02]\n");
	fprintf(stderr,"         -1 INT length of the first read [default 70bp]\n");
	fprintf(stderr,"         -2 INT length of the second read [default 70bp]\n");
	fprintf(stderr,"         -N INT number of read pairs [default 1000000]\n");
	fprintf(stderr,"         -R FLOAT fraction of indels [default 0.15]\n");
	fprintf(stderr,"         -S INT seed for random generator [default -1]\n");
	fprintf(stderr,"         -d INT outer distance between the two ends [default 500pb]\n");
	fprintf(stderr,"         -g INT average gap size [default 1]\n");
	fprintf(stderr,"         -D INT standard deviation [default 50pb]\n");
	fprintf(stderr,"         -X FLOAT probability that indel is extended [default 0.3]\n");
	fprintf(stderr,"\n**********************************************************\n");
	return 1;
}



int main(int argc, char *argv[])
{
	int dist, std_dev, size_l, size_r,c,ind;
	uint64_t N;
	FILE *fout1, *fout2;
	clock_t start = clock();
	char flag[10];
	N = 1000000; dist = 500; std_dev = 50; size_l = size_r = 70;
	flag[0]='O';flag[1]='K';flag[2]='\0'; ind = 0;
	while ((c = getopt(argc, argv, "e:N:1:2:r:R:S:d:g:D:X:")) >= 0) {
		switch (c) {
		case 'N': N = atoi(optarg); break; //number of reads
		case '1': size_l = atoi(optarg); break;//length of the first read
		case '2': size_r = atoi(optarg); break; //length of the second read
		case 'e': ERR_RATE = atof(optarg); break; //bcer
		case 'r': MUT_RATE = atof(optarg); break; //rate of mutations
		case 'R': INDEL_FRAC = atof(optarg); break;//del.
		case 'S': SEED = atoi(optarg);break;//seed for random number gen.
		case 'd': dist=atoi(optarg);break;//distance between two reads
		case 'g': GAP_SIZE=atoi(optarg);break;//average gap size
		case 'D': std_dev=atoi(optarg);break;//standard deviation
		case 'X': INDEL_EXT=atof(optarg);break;//probability that indel is extended
		}
	}
	if(argc - optind < 3) return simu_usage(); //checks the number of arguments
	fout1 = fopen(argv[optind+1], "w");
	fout2 = fopen(argv[optind+2], "w");
	if (!fout1 || !fout2) {
		fprintf(stderr, "[%s] file open error\n",__func__);
		return 1;
	}
	if (SEED <= 0) SEED = time(0)&0x7fffffff; //sets up seed
	srand48(SEED);
	fprintf(stderr,"[%s] simulator seed = %d\n",__func__,SEED);
	if(core(fout1,fout2,argv[optind],std_dev,size_l,size_r,N,dist)==-1){ind = -1; strcpy(flag,"FAIL");}
    printf("[%s] return value: %d %s\n",__func__, ind,flag);
	printf ( "[%s] Total time taken: %f sec\n",__func__, ( (double)clock() - start ) / CLOCKS_PER_SEC );
	fclose(fout1); fclose(fout2);
	return 0;
}


