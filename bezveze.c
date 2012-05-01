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
#include <unistd.h>

#define MAX 4
#define MIN 0

KSEQ_INIT(gzFile, gzread);

static double ERR_RATE = 0.02;
static double MUT_RATE = 0.001;
static double INDEL_FRAC = 0.15;
static double SEED = -1;

typedef unsigned short mut_t;

typedef struct{
	int l,m;
	mut_t *s;
} mutseq_t;

char *iter;
char *tot_seq;
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

const int poisson_random_number(const double lambda, uint64_t seed_gen) //poisson distribution - indel length - lambda: average len.
{
    int k=0;                         
    const int max_k = 1000;           
    double P = exp(-lambda);          
    double sum=P;
    double drand48(); 
    double p;
    srand48(seed_gen);
    p = drand48();                    
    if (sum>=p) return 0;             
        for (k=1; k<max_k; ++k) {         
             P*=lambda/(double)k;            
             sum+=P;                         
             if (sum>=p) break;             
    }

    return k;                       
}

char swap_base(char base,uint64_t i){ //uniform distribution
	char new_base;
	double value;
	srand(i);
	value = rand()/(RAND_MAX+1.0)*(MAX-MIN)+MIN;
	if ((value >=0 ) && (value < 1 )) return 'A';
	else if ((value >=1 ) && (value < 2 )) return 'T';
	else if ((value >=2 ) && (value < 3 )) return 'G';
	else return 'C';
} 

int generate_mut_index(uint64_t tot_len,uint64_t i){
	uint64_t mut_index;
	double drand48();
	double r;
	srand48(i);
	r = drand48();
	mut_index = (long long)(trunc(r * tot_len));
	//printf("ovo je mut index %llu\n",(long long)mut_index);
	return mut_index;
}

uint64_t *get_mut_index_array(uint64_t tot_len, uint64_t N_rate){
	uint64_t i;
	uint64_t *array_index;
	array_index = (uint64_t *)malloc(tot_len * sizeof(uint64_t));
	for (i=0;i<N_rate;i++){
		*(array_index + i) = generate_mut_index(tot_len,i);
	}
	return array_index;
}
	
int check_index(uint64_t *array_of_int, uint64_t N_rate, uint64_t current_index){
	uint64_t i;
	for (i=0;i<N_rate;i++){
		if (*(array_of_int) == current_index) return 1;
	}
	return 0;
}	

char simulate_BCER(char base, uint64_t i){
	double drand48();
	double r;
	srand48(i);
	r = drand48();
	if (r <= ERR_RATE) base=swap_base(base,i);
	//printf("%c\n",base);
	return base;
}
	
void generate_mutations(char *argv, float m_rate,uint64_t total){ //fali parametara
	mutseq_t *ret[2];
	FILE *fp_outm,*proba;
	//gzFile fp,fpc;
	uint64_t total_len, iterator, N_rate, *array_index;
	uint8_t *tmp_seq[2];
	kseq_t *seq, *copy;
	int l,n_ref;
	uint64_t i,j;
	//char *tot_seq;
	N_rate = m_rate * total;
	fp = gzopen(argv, "r");
	proba = fopen(argv,"r");
	fp_outm = fopen("output_mut1.fa","w");
	if (!fp_outm){
		fprintf(stderr,"[%s] file open error\n",__func__);
	}
	seq = kseq_init(fp);
	array_index = get_mut_index_array(total, N_rate);
	tot_seq = (char *)malloc((total+1)* sizeof(char)); n_ref = 0;total_len = 0;
	tot_seq[total]='\0';
	while ((l = kseq_read(seq)) >= 0){
		total_len+=l;
		++n_ref;
		//fprintf(fp_outm,"%s",seq->seq.s);
		//iter = seq->seq.s;
		tot_seq = seq->seq.s;
	}
	//printf("normalna %s\n",tot_seq);
	for(i=0;i<N_rate;i++){
		*(tot_seq + *(array_index+i)) = swap_base(*(tot_seq + *(array_index+i)),i);
	}
	for(i=0;i<total;i++){
		*(tot_seq + i) = simulate_BCER(*(tot_seq+i),i);
		
	}
	//printf("normalna %s\n",tot_seq);
	fclose(fp_outm);
	//kseq_destroy(seq); 
	gzclose(fp); 
	printf("[core] done...\n[core] simulating BCER...\n");
}

void generate_gaps(uint64_t gap_pos, int gap_size){
	uint64_t i;
	tot_seq[gap_pos+1]='\0';
	strcat(tot_seq,(tot_seq + gap_pos + gap_size));
	//printf("%d %llu %llu\n",gap_size, (long long)gap_pos,(long long)strlen(tot_seq));
	//printf("%s\n",tot_seq);
}

void get_gaps(char *g_ratec,float a_len, uint64_t total){
	float g_rate;
	int g_size;
	uint64_t tot_gaps, gap_pos,i,j;
	double drand48();
	double r;
	g_rate = atof(g_ratec);
	tot_gaps = g_rate * total;
	//printf("%llu\n",(long long)tot_gaps);
	//generate_gaps(1,11);
	for (i=0;i<tot_gaps;i++){
		srand48(i);
		r = drand48();
		g_size = poisson_random_number(a_len,i);
		gap_pos = (long long)(trunc(r * (total-g_size)));
		//printf("%d %llu\n",g_size, (long long)gap_pos);
		for (j=gap_pos;j<((gap_pos + g_size)<total?(gap_pos+g_size):total);j++){
			*(tot_seq+j)='_';
	}
}	
}
void core(char *argv, double mut_rate){
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
		//iter = seq->seq.s;
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
	printf("[%s] transferring sequence into memory and generating mutations...\n",__func__);
	generate_mutations(argv,mut_rate,total_len);
	//printf("ovo:\n%s\n",tot_seq);
	printf("[%s] done...\n",__func__);
	printf("[%s] generating gaps...\n",__func__);
	get_gaps("0.05",3,total_len);
	printf("[%s] done...\n",__func__);
	//printf("ovo:\n%s\n",tot_seq);
}

static int simu_usage(){
	fprintf(stderr,"**********************************************************\n");
	fprintf(stderr,"\nUsage: ./a.out [options] <in_seq.fa> <out_read1.fq> <out_read2.fq>\n\n");
	fprintf(stderr,"Options: -r FLOAT rate of mutations\n");
	fprintf(stderr,"         -e FLOAT base error rate (default 0.02)\n");
	fprintf(stderr,"         -1 INT length of first read (default 70bp)\n");
	fprintf(stderr,"         -2 INT length of second read (default 70bp)\n");
	fprintf(stderr,"         -N INT number of read pairs (default 100000)\n");
	fprintf(stderr,"         -R FLOAT fraction of indels\n");
	fprintf(stderr,"         -S INT seed for random generator (default -1)\n");
	fprintf(stderr,"\n**********************************************************\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int dist, std_dev, size_l, size_r,c;
	int64_t N;
	FILE *fout1, *fout2;
	clock_t start = clock();
	N = 1000000; dist = 500; std_dev = 50; size_l = size_r = 70;
	while ((c = getopt(argc, argv, "e:N:1:2:r:R:S")) >= 0) {
		switch (c) {
		case 'N': N = atoi(optarg); break; //broj pair end readova
		case '1': size_l = atoi(optarg); break;//length of first read
		case '2': size_r = atoi(optarg); break; //length of second read
		case 'e': ERR_RATE = atof(optarg); break; //bcer
		case 'r': MUT_RATE = atof(optarg); break; //rate of mutations
		case 'R': INDEL_FRAC = atof(optarg); break;//del.
		case 'S': SEED = atoi(optarg);break;//seed for random number gen.
		}
	}
	if(argc - optind < 3) return simu_usage();
	fout1 = fopen(argv[optind+1], "w");
	fout2 = fopen(argv[optind+2], "w");
	if (!fout1 || !fout2) {
		fprintf(stderr, "[%s] file open error\n",__func__);
		return 1;
	}
	if (SEED <= 0) SEED = time(0)&0x7fffffff;
	core(argv[optind],MUT_RATE);
    printf("[%s] return value: %d OK\n",__func__, 0);
	printf ( "[%s] Total time taken: %f sec\n",__func__, ( (double)clock() - start ) / CLOCKS_PER_SEC );
	
	fclose(fout1); fclose(fout2);
	return 0;
}
