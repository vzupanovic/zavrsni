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

#define MAX 3
#define MIN 0

KSEQ_INIT(gzFile, gzread);

static double ERR_RATE = 0.02;
static double MUT_RATE = 0.001;
static double INDEL_FRAC = 0.15;
static int SEED = -1;

typedef unsigned short mut_t;

typedef struct{
	int l,m;
	mut_t *s;
} mutseq_t;

char *tot_seq,*read_f,*iter;

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
    srand48(seed_gen+SEED);
    p = drand48();                    
    if (sum>=p) return 0;             
        for (k=1; k<max_k; ++k) {         
             P*=lambda/(double)k;            
             sum+=P;                         
             if (sum>=p) break;             
    }

    return k;                       
}

char swap_base(char base){ //uniform distribution
	char new_base;
	double value;
	value = rand()/(RAND_MAX+1.0)*(MAX-MIN)+MIN;
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
	         return -1; //ovo se ne bi smjelo desit
	   }
	
} 

int generate_mut_index(uint64_t tot_len,uint64_t i){
	uint64_t mut_index;
	double drand48();
	double r;
	srand48(i+SEED);
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
	srand48(i+SEED);
	r = drand48();
	if (r <= ERR_RATE) base=swap_base(base);
	//printf("%c\n",base);
	return base;
}
	
void generate_mutations(int dist){ //fali parametara
	int i,pos,n_mut;
	double r;
	n_mut = dist * MUT_RATE;
	srand48(SEED);
	for (i=0;i<n_mut;i++){
		r=drand48();
		pos = (int)(trunc(r * dist));
		*(read_f+pos)=swap_base(*(read_f+pos));
	}
}

void generate_gaps(int gap_pos, int gap_size){
	uint64_t i;
	read_f[gap_pos]='\0';
	strcat(read_f,(read_f + gap_pos + gap_size));
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
	tot_gaps = INDEL_FRAC * total;
	//printf("%llu\n",(long long)tot_gaps);
	//generate_gaps(1,11);
	for (i=0;i<tot_gaps;i++){
		srand48(SEED+i);
		r = drand48();
		g_size = poisson_random_number(a_len,i);
		gap_pos = (long long)(trunc(r * (total-g_size)));
		//printf("%d %llu\n",g_size, (long long)gap_pos);
		for (j=gap_pos;j<((gap_pos + g_size)<total?(gap_pos+g_size):total);j++){
			*(tot_seq+j)='_';
	}
}	
}
int core(char *argv, int std_dev, int size_l, int size_r, uint64_t N,int dist){
	mutseq_t *ret[2];
	gzFile fp;
	uint64_t total_len;
	kseq_t *seq;
	int l,n_ref,max_size,Q;
	uint64_t i,j;
	fp = gzopen(argv, "r");
	seq = kseq_init(fp);
	total_len = n_ref = 0;
	frequency_A=frequency_T=frequency_G=frequency_C=0.;
	nA=nT=nG=nC=j=0;
	max_size = size_l > size_r? size_l : size_r;
	Q = (ERR_RATE == 0.0)? 'I' : (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;
	srand48(SEED);
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
		if (l < dist + 3*std_dev){
			fprintf(stderr,"[%s] ERROR sequence to short for given parametars!\n",__func__);
			return -1;
		}
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
    tot_seq = (char *)malloc((total_len+1)*sizeof(char));
	tot_seq = seq->seq.s;
	printf("[%s] transferring sequence into memory and generating mutations...\n",__func__);
	for(i=0;i<N;i++){
		double ran;
		int d,pos;
		do{
			ran = ran_normal();
			ran = ran * std_dev + dist;
			d = (int)(ran + 0.5);
			d = d > max_size ? d : max_size;
			pos = (int)((total_len-d+1)*drand48());
		}while(pos < 0 || pos >= total_len || (pos + d - 1) >= total_len);
		read_f = (char *)malloc((dist+1)*sizeof(char));
		read_f[dist+1]='\0';
		strncpy(read_f,tot_seq+pos,dist);
		//printf("%s\n",read_f);
		generate_mutations(dist);
		//printf("S: %s\n",read_f);
		//printf("pozicija %d\n",pos);
		int n_n=0;
		int n_indel = (int)(INDEL_FRAC * dist);
		if(n_indel >= dist){
			fprintf(stderr,"[%s] ERROR sequence to short for given parametars!\n",__func__);
			return -1;
		}
		do{
			int type_indel;
			char *fragment;
			type_indel = (rand()/(RAND_MAX+1.0)*2>=1)?1:0;
			if (type_indel == 1){pos = (int)trunc(rand()/(RAND_MAX+1.0)*(dist-1));generate_gaps(pos,1);}
			/*else{
				char fragment[4];
				double r;
				char base,*temp;//pazi kad produzujes read koji je vec 500 ne valja sporo
				temp = (char *)malloc(dist*sizeof(char));
				r = rand()/(RAND_MAX+1.0);
				base=(r < 0.25)?'A':((r>=0.25 && r<0.5)?'T':(r>=0.25 && r<0.75)?'C':'G');
				fragment[0]=base;fragment[1]='\0';
				pos = (int)trunc(rand()/(RAND_MAX+1.0)*dist);//insert(pos,fragment);
				strcpy(temp,read_f+pos);
				read_f[pos]='\0';
				strcat(read_f,fragment);
				strcat(read_f,temp);
				}*/
			n_n++;
		}while(n_n<n_indel);
		//printf("N: %s\n",read_f);
	}
    //kseq_destroy(seq);
	gzclose(fp);
	//printf("lalalalala %s\n",tot_seq);
	//generate_mutations(argv,mut_rate,total_len);
	//printf("ovo:\n%s\n",tot_seq);
	//printf("[%s] done...\n",__func__);
	//printf("[%s] generating gaps...\n",__func__);
	//get_gaps("0.05",3,total_len);
	//printf("[%s] done...\n",__func__);
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
	fprintf(stderr,"         -d INT outer distance between the two ends (default 500)\n");
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
	while ((c = getopt(argc, argv, "e:N:1:2:r:R:S:d:")) >= 0) {
		switch (c) {
		case 'N': N = atoi(optarg); break; //broj pair end readova
		case '1': size_l = atoi(optarg); break;//length of first read
		case '2': size_r = atoi(optarg); break; //length of second read
		case 'e': ERR_RATE = atof(optarg); break; //bcer
		case 'r': MUT_RATE = atof(optarg); break; //rate of mutations
		case 'R': INDEL_FRAC = atof(optarg); break;//del.
		case 'S': SEED = atoi(optarg);break;//seed for random number gen.
		case 'd': dist=atoi(optarg);break;//distance between two reads
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
	if(core(argv[optind],std_dev,size_l,size_r,N,dist)==-1){ind = -1; strcpy(flag,"FAIL");}
    printf("[%s] return value: %d %s\n",__func__, ind,flag);
	printf ( "[%s] Total time taken: %f sec\n",__func__, ( (double)clock() - start ) / CLOCKS_PER_SEC );
	
	fclose(fout1); fclose(fout2);
	return 0;
}

