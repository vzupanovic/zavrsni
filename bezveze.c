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

#define MAX 4
#define MIN 0

KSEQ_INIT(gzFile, gzread);

static double ERR_RATE = 0.02;
static double MUT_RATE = 0.001;

typedef unsigned short mut_t;

typedef struct{
	int l,m;
	mut_t *s;
} mutseq_t;

char *iter;

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
 
 
double ran_normal() //normal distribution copied from generan.c
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


char swap_base(char base){ //uniform distribution
	char new_base;
	double value;
	srand(time(NULL));
	value = rand()/(RAND_MAX+1.0)*(MAX-MIN)+MIN;
	if ((value >=0 ) && (value < 1 )) return 'A';
	else if ((value >=1 ) && (value < 2 )) return 'T';
	else if ((value >=2 ) && (value < 3 )) return 'G';
	else return 'C';
} 

void generate_mutations(char *argv){ //nesto isprobavam
	mutseq_t *ret[2];
	gzFile fp,fpc;
	uint64_t total_len, iterator;
	uint8_t *tmp_seq[2];
	kseq_t *seq, *copy;
	int l,n_ref;
	uint64_t i,j;
	fp = gzopen(argv, "r");
	//fpc = gzopen("proba.fa","w");
	seq = kseq_init(fp);
	//copy = kseq_init(fpc);
	//copy = (char*)calloc(tot_len+1, 1);
	//tmp_seq[0] = (uint8_t*)calloc(tot_len+2, 1);
	//tmp_seq[1] = (uint8_t*)calloc(tot_len+2, 1);
	while ((l = kseq_read(seq)) >= 0){
		printf("bla\n");
		total_len+=l;
		++n_ref;
		//if (seq->seq.s[i] == 'T') printf("T\n");
		
		iter = seq->seq.s;
		//printf("seq: %c %c \n", *iter, *(iter+1));
	}
	printf("tu ide: \n");
    for (iterator=0; iterator < total_len; iterator++)
    printf("%c",*(iter+iterator));
    kseq_destroy(seq);
	//kseq_destroy(copy);
	gzclose(fp);
	//gzclose(fpc);
	
}

	
void core(char *argv){
	mutseq_t *ret[2];
	gzFile fp;
	uint64_t total_len;
	kseq_t *seq;
	int l,n_ref;
	uint64_t i,j;
	fp = gzopen(argv, "r");
	seq = kseq_init(fp);
	total_len = n_ref = 0;
	frequency_A=frequency_T=frequency_G=frequency_C=0.;
	nA=nT=nG=nC=j=0;
	fprintf(stderr, "[%s] calculating the total length of the sequnce...\n",__func__);
	while ((l = kseq_read(seq)) >= 0){
		printf("[%s] name: %s\n",__func__,seq->name.s);
		if(seq->comment.l) printf("[%s] comment: %s\n",__func__,seq->comment.s);
		total_len+=l;
		++n_ref;
		//if (seq->seq.s[i] == "N") printf("N\n");
		iter = seq->seq.s;
		//printf("seq: %c %c \n", *iter, *(iter+1));
		if (seq->qual.l) printf ("qual: %s\n",seq->qual.s);
	}
	fprintf(stderr, "[%s] %d sequences, total length: %llu\n", __func__, n_ref, (long long)total_len);
	kseq_destroy(seq);
	gzclose(fp);
	fp = gzopen(argv, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0){
		for (i=0;i<l;i++){
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
	
	frequency_A = (double_t)nA/l;
	frequency_C = (double_t)nC/l;
	frequency_G = (double_t)nG/l;
	frequency_T = (double_t)nT/l;
	frequency_N = (double_t)nN/l;
	j++;	  
    printf("[%s] frequency per sequence [%llu/%llu] A - %f | C - %f | G - %f | T - %f | *N - %f - unknown nucleotides (percentage) \n",__func__,(long long)j,(long long)n_ref,frequency_A*100,frequency_C*100,frequency_G*100,frequency_T*100,frequency_N*100);
    }
   
    kseq_destroy(seq);
	gzclose(fp);
	
}

int main(int argc, char *argv[])
{
	clock_t start = clock();
	if (argc == 1){
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		printf("[%s] return value: %d FAIL\n",__func__, 1);
		return 1;
	}
	core(argv[1]);
    printf("[%s] return value: %d OK\n",__func__, 0);
	printf ( "[%s] Total time taken: %f sec\n",__func__, ( (double)clock() - start ) / CLOCKS_PER_SEC );
	swap_base('A');
	generate_mutations(argv[1]);
	return 0;
}
