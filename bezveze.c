//Prevoditi s: gcc -lz -lm -std=gnu99 bezveze.c
//Pokretati s npr. : ./a.out Mus_musculus.NCBIM37.61.dna_rm.chromosome.1.fa 30

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

int generate_mut_index(uint64_t tot_len){
	uint64_t mut_index;
	double drand48();
	double r;
	srand48(time(NULL));
	r = drand48();
	mut_index = (long long)(trunc(r * tot_len));
	return mut_index;
}

uint64_t *get_mut_index_array(uint64_t tot_len, double N_rate){
	uint64_t i;
	uint64_t *array_index;
	array_index = (uint64_t *)malloc(tot_len * sizeof(uint64_t));
	for (i=0;i<N_rate;i++){
		*(array_index + i) = generate_mut_index(tot_len);
	}
	return array_index;
}
		
	
void generate_mutations(char *argv, float m_rate,uint64_t total){ //fali parametara
	mutseq_t *ret[2];
	FILE *fp_outm;
	gzFile fp,fpc;
	uint64_t total_len, iterator, N_rate, *array_index;
	uint8_t *tmp_seq[2];
	kseq_t *seq, *copy;
	int l,n_ref;
	uint64_t i,j;
	N_rate = m_rate * total;
	for (iterator=0; iterator<N_rate;iterator++){
		double ran;
		int d, pos;
		ran = ran_normal();
	}
	fp = gzopen(argv, "r");
	fp_outm = fopen("output_mut1.fa","wb");
	if (!fp_outm){
		fprintf(stderr,"[%s] file open error\n",__func__);
	}
	seq = kseq_init(fp);
	array_index = get_mut_index_array(total, N_rate);
	while ((l = kseq_read(seq)) >= 0){
		total_len+=l;
		++n_ref;
		fprintf(fp_outm,"%s",seq->seq.s);
		//fputs(seq->seq.s,fp_outm);
		//fwrite(seq->seq.s,1,sizeof(seq->seq.s),fp_outm);
		iter = seq->seq.s;
	}
	fclose(fp_outm);
	kseq_destroy(seq); 
	gzclose(fp); 
	
	
}

	
void core(char *argv, char *m_ratec){
	mutseq_t *ret[2];
	gzFile fp;
	uint64_t total_len;
	kseq_t *seq;
	int l,n_ref;
	uint64_t i,j;
	float mut_rate;
	fp = gzopen(argv, "r");
	seq = kseq_init(fp);
	total_len = n_ref = 0;
	mut_rate = atof(m_ratec);
	frequency_A=frequency_T=frequency_G=frequency_C=0.;
	nA=nT=nG=nC=j=0;
	fprintf(stderr, "[%s] calculating the total length of the sequnce...\n",__func__);
	while ((l = kseq_read(seq)) >= 0){
		printf("[%s] name: %s\n",__func__,seq->name.s);
		if(seq->comment.l) printf("[%s] comment: %s\n",__func__,seq->comment.s);
		total_len+=l;
		++n_ref;
		iter = seq->seq.s;
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
	generate_mutations(argv,mut_rate,total_len);
}

int main(int argc, char *argv[])
{
	int dist, std_dev, size_l, size_r;
	int64_t N;
	FILE *fout1, *fout2;
	clock_t start = clock();
	N = 1000000; dist = 500; std_dev = 50; size_l = size_r = 70;
	if (argc == 1){
		fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
		printf("[%s] return value: %d FAIL\n",__func__, 1);
		return 1;
	}
	core(argv[1],argv[2]);
    printf("[%s] return value: %d OK\n",__func__, 0);
	printf ( "[%s] Total time taken: %f sec\n",__func__, ( (double)clock() - start ) / CLOCKS_PER_SEC );
	
	
	return 0;
}
