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
static int GAP_SIZE = 1;
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

const int PoissonRandomNumber(const double lambda)
{
  int k=0;                          //Counter
  const int max_k = 1000;           //k upper limit
  double p = ran_uniform(); //uniform random number
  double P = exp(-lambda);          //probability
  double sum=P;                     //cumulant
  if (sum>=p) return 0;             //done allready
  for (k=1; k<max_k; ++k) {         //Loop over all k:s
    P*=lambda/(double)k;            //Calc next prob
    sum+=P;                         //Increase cumulant
    if (sum>=p) break;              //Leave loop
  }

  return k;                         //return random number
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
	         return base; //ovo se ne bi smjelo desit
	   }
	
} 

char get_complement(char base){
	if (base == 'A') return 'T';
	else if (base == 'T') return 'A';
	else if (base == 'G') return 'C';
	else if (base == 'C') return 'G';
	else return base;
}

int generate_mut_index(uint64_t tot_len,uint64_t i){
	uint64_t mut_index;
	double value;
	double r;
	//srand48(i+SEED);
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

char *simulate_BCER(int read_length, char *read_local){
	double value;
	int i,index,n_errors;
	n_errors=(int)(read_length*ERR_RATE);
	srand(SEED);
	for (i=0;i<n_errors;i++){
		value = rand()/(RAND_MAX+1.0)*(read_length-1);
		index = trunc((int)value);
		*(read_local+index) = swap_base(*(read_local+index));
    }
    return read_local;
	
}

void generate_mutations(int dist){ //fali parametara
	int i,pos,n_mut;
	double r;
	n_mut = dist * MUT_RATE;
	//srand48(SEED);
	for (i=0;i<n_mut;i++){
		r=drand48();
		pos = (int)(trunc(r * dist));
		if(*(read_f+pos)!='N')
		    *(read_f+pos)=swap_base(*(read_f+pos));
	}
}

void generate_gaps(int gap_pos, int gap_size){
	uint64_t i;
	read_f[gap_pos]='\0';
	strcat(read_f,(read_f + gap_pos + gap_size));
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
		//srand48(SEED+i);
		r = drand48();
		//g_size = poisson_random_number(a_len,i);
		gap_pos = (long long)(trunc(r * (total-g_size)));
		//printf("%d %llu\n",g_size, (long long)gap_pos);
		for (j=gap_pos;j<((gap_pos + g_size)<total?(gap_pos+g_size):total);j++){
			*(tot_seq+j)='_';
	}
}	
}
int core(FILE *fout1,FILE *fout2,char *argv, int std_dev, int size_l, int size_r, uint64_t N,int dist){
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
	max_size = size_l > size_r? size_l : size_r;
	Q = (ERR_RATE == 0.0)? 'I' : (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;
	q_string = (char *)malloc((size_l+1)*sizeof(char));
	q2_string = (char *)malloc((size_r+1)*sizeof(char));
	for(int k=0;k<size_l;k++){q_string[k]=Q;};
	for(int k=0;k<size_r;k++){q2_string[k]=Q;};
	//srand48(SEED);
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
	printf("[%s] transferring sequence into memory and generating errors...\n",__func__);
	counter_a=counter_b=0;
	//printf("ovo je size_l Size_r d %d %d\n",size_l,size_r);
	for(i=0;i<N;i++){
		double ran;
		int d,pos;
		uint64_t fp1_b,fp1_e,fp2_b,fp2_e,begin,end;
		char *read1,*read2;
		do{
			ran = ran_normal();
			ran = ran * std_dev + dist;
			d = (int)(ran + 0.5);
			d = d > max_size ? d : max_size;
			pos = (int)((total_len-d+1)*drand48());//sporno uvijek read??? to ne smije!
			//pos=(int)(rand()/(RAND_MAX+1.0)*(total_len-d+1));//neka sad bude ovo...
		}while(pos < 0 || pos >= total_len || (pos + d - 1) >= total_len);
		read_f = (char *)malloc((d+1+(int)((INDEL_FRAC+0.5)*d))*sizeof(char));counter_a++;//zajeb, nije valjalo dist -> d
		read_f[d+1]='\0';
		strncpy(read_f,tot_seq+pos,d);
		begin=pos; end=pos+d;
		//printf("------------------------------------\n");
		//printf("pocetak %llu kraj %llu duljina %d %s\n",(long long)begin,(long long)end,d,read_f);
		generate_mutations(d);
		int n_n=0;
		int n_indel = (int)(INDEL_FRAC * d);
		int curr_dist=d;
		//printf("ovo je d %d\n",curr_dist);
		if((n_indel >= d) || ((GAP_SIZE*n_indel)>=d)){
			fprintf(stderr,"[%s] ERROR sequence to short for given parametars!\n",__func__);
			return -1;
		}
		//printf("prije: %s\n",read_f);
		do{
			int type_indel;
			char *fragment;
			type_indel = (rand()/(RAND_MAX+1.0)*2>=1)?1:0;
			//type_indel = 1;
			if (type_indel == 1){pos = (int)trunc(rand()/(RAND_MAX+1.0)*(curr_dist-1));generate_gaps(pos,GAP_SIZE);curr_dist=curr_dist-GAP_SIZE;}
			else{
				//printf("ub\n");
				char fragment[4]={0x0};
				char keeper[1000]={0x0};
				int pos;
				double r;
				char base;//pazi kad produzujes read koji je vec 500 ne valja sporo
				r = rand()/(RAND_MAX+1.0);
				base=(r < 0.25)?'A':((r>=0.25 && r<0.5)?'T':(r>=0.25 && r<0.75)?'C':'G');
				fragment[0]=base;fragment[1]='\0';
				pos = (int)trunc(rand()/(RAND_MAX+1.0)*(curr_dist-1));//insert(pos,fragment);
				strcat(keeper,read_f+pos);
				read_f[pos]='\0';
				strcat(read_f,fragment);
				strcat(read_f,keeper);
				curr_dist++;
				
			}
			n_n++;
		}while(n_n<n_indel);
		//printf("lala read %s\n",read_f);
		//printf("curr DIST %d\n",curr_dist);
		read1=(char *)malloc((size_l+1)*sizeof(char));
		n_errors=(int)(INDEL_FRAC*size_l);
		read2=(char *)malloc((size_r+1)*sizeof(char));
		strncpy(read1,read_f,size_l);read1[size_l+1]='\0';
		int internal_counter=0;
		for(int k=(curr_dist-1);internal_counter<size_r;k--){
			*(read2+internal_counter)=get_complement(*(read_f + k));
			internal_counter++;
		}
		read2[internal_counter]='\0';
		read1=simulate_BCER(size_l,read1);
		read2=simulate_BCER(size_r,read2);
		fprintf(fout1,"@%s_%llu_%llu_0:0:0_0:0:0_%llx/%d\n",seq->name.s,(long long)begin,(long long)end,(long long)counter_a,1);
		fprintf(fout1,"%s\n+\n%s\n",read1,q_string);
		fprintf(fout2,"@%s_%llu_%llu_0:0:0_0:0:0_%llx/%d\n",seq->name.s,(long long)begin,(long long)end,(long long)counter_a,2);
		fprintf(fout2,"%s\n+\n%s\n",read2,q2_string);
		//printf("read %s\n",read_f);
		//printf("read 1 %s\n",read1);
		//printf("read 2 %s\n",read2);
		free(read1);
		free(read2);
		//printf("%s\n",read_f);
		//printf("pozicija %d\n",pos);
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
	fprintf(stderr,"         -e FLOAT base error rate [default 0.02]\n");
	fprintf(stderr,"         -1 INT length of first read [default 70bp]\n");
	fprintf(stderr,"         -2 INT length of second read [default 70bp]\n");
	fprintf(stderr,"         -N INT number of read pairs [default 100000]\n");
	fprintf(stderr,"         -R FLOAT fraction of indels\n");
	fprintf(stderr,"         -S INT seed for random generator [default -1]\n");
	fprintf(stderr,"         -d INT outer distance between the two ends [default 500]\n");
	fprintf(stderr,"         -g INT average gap size [default 1]\n");
	fprintf(stderr,"         -D INT standard deviation [default 50]\n");
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
	while ((c = getopt(argc, argv, "e:N:1:2:r:R:S:d:g:D:")) >= 0) {
		switch (c) {
		case 'N': N = atoi(optarg); break; //broj pair end readova
		case '1': size_l = atoi(optarg); break;//length of first read
		case '2': size_r = atoi(optarg); break; //length of second read
		case 'e': ERR_RATE = atof(optarg); break; //bcer
		case 'r': MUT_RATE = atof(optarg); break; //rate of mutations
		case 'R': INDEL_FRAC = atof(optarg); break;//del.
		case 'S': SEED = atoi(optarg);break;//seed for random number gen.
		case 'd': dist=atoi(optarg);break;//distance between two reads
		case 'g': GAP_SIZE=atoi(optarg);break;//average gap size
		case 'D': std_dev=atoi(optarg);break;//standard deviation
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
	srand48(SEED);
	if(core(fout1,fout2,argv[optind],std_dev,size_l,size_r,N,dist)==-1){ind = -1; strcpy(flag,"FAIL");}
    printf("[%s] return value: %d %s\n",__func__, ind,flag);
	printf ( "[%s] Total time taken: %f sec\n",__func__, ( (double)clock() - start ) / CLOCKS_PER_SEC );
	
	fclose(fout1); fclose(fout2);
	return 0;
}

