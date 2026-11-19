#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>
static double cpu(void){struct timespec t;clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&t);return t.tv_sec+1e-9*t.tv_nsec;}
int main(void){
  int64_t ns[7]={12,13,15,16,17,20,23}; mp_limb_t a[32],b[32],r[64];
  printf("gmp %s\n",gmp_version);
  for(int k=0;k<7;k++){ int64_t n=ns[k], it=4000000;
    for(int i=0;i<n;i++){a[i]=0x9e3779b97f4a7c15ULL*(i+1);b[i]=0xc2b2ae3d27d4eb4fULL*(i+3);}
    double best_m=1e9,best_s=1e9;
    for(int rep=0;rep<3;rep++){
      double t=cpu(); for(int64_t i=0;i<it;i++){mpn_mul_n(r,a,b,n);a[0]^=r[n];} t=cpu()-t; if(t<best_m)best_m=t;
      double u=cpu(); for(int64_t i=0;i<it;i++){mpn_sqr(r,a,n);a[0]^=r[n];} u=cpu()-u; if(u<best_s)best_s=u; }
    printf("limbs %2ld: mul_n %.1f ns  sqr %.1f ns\n",(long)n,1e9*best_m/it,1e9*best_s/it); }
  int64_t bits[3]={768,1030,1263};
  for(int k=0;k<3;k++){ mpz_t x,e,m,y; mpz_inits(x,e,m,y,NULL); mpz_set_ui(m,1); mpz_mul_2exp(m,m,bits[k]); mpz_sub_ui(m,m,187); mpz_set_ui(x,3); mpz_sub_ui(e,m,1);
    double best=1e9; for(int rep=0;rep<3;rep++){ double v=cpu(); for(int i=0;i<400;i++) mpz_powm(y,x,e,m); v=cpu()-v; if(v<best)best=v; }
    printf("powm %ld bits: %.4f ms\n",(long)bits[k],1e3*best/400); }
  return 0;}
