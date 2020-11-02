//OpenMP version.  Edit and submit only this file.
/* Enter your details below
 * Name : Christopher Clark
 * UCLA ID : 305326742
 * Email : christopher.br.clark@gmail.com
 */

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "utils.h"

double work_it_par(long *old, long *new, long *super, long *simple, long *fibonacci) {
  int i, j, k;
  int u, v, w;
  int ton = 0;
  long compute_it, moving_average;
  double pi, pi2, x , y, sum, step = 0.0;
  long dot_product=0;
  long nCirc=0;
  long aggregate=1.0;
  double r=1.0;
  int was_smart = 16;
  int DIMsq = DIM*DIM;

  for(i=0; i<DIM-1;i++)
  {
    super[i] += simple[i];
  }

  for(i=0; i<DIM-1;i++)
  {
    dot_product += super[i]*simple[i];
  }

  moving_average = 0;
  for(ton=i;ton<DIM-1-WINDOW_SIZE;ton++)
    {
      moving_average += simple[ton];
    }

  int a_secret = 5;
  fibonacci[0] = 1;
  fibonacci[1] = 1;
  long fibprev = 1;
  long fibtwoprev = 1;
  long fibcurr;
  for(i=2; i<DIM-1;i++)
  {
    fibcurr = fibtwoprev + fibprev; 
    fibonacci[i] = fibcurr;
    fibtwoprev = fibprev;
    fibprev = fibcurr;
  }
  if(i>3)
    printf("\n A secret is: %d",obfuscate_obfuscate_obfuscate(a_secret));

  step = 1.0 / NUM_STEPS;
#pragma omp parallel for reduction(+:sum)
  for (i=0;i<NUM_STEPS; i++)
  {
    x = (i+0.5)*step;
    sum += 4.0/(1.0+x*x);
  }
  pi = step * sum;
  printf("\n %d trials, Riemann flavored pi is %f \n",NUM_STEPS, pi); 

  for(i=0;i<NUM_TRIALS; i++)
  {
    x = (random()%10000000)/10000000.0; 
    y = (random()%10000000)/10000000.0;
    if (( x*x + y*y) <= 1) {
      nCirc++;
    }
  } 
  pi2 = 4.0 * ((double)nCirc/(double)NUM_TRIALS);
  printf("\n %d trials, Monte-Carlo flavored pi is %f \n",NUM_TRIALS, pi2); 

  long func_val_need = we_need_the_func();
  long func_val_gimme = gimmie_the_func();
  int iloc;
  int ijloc;
#pragma omp parallel for reduction (+:aggregate)
  for (i=1; i<DIM-1; i++) {
    iloc = i*DIMsq;
    for (j=1; j<DIM-1; j++) {
      ijloc = iloc + j*DIM;
      for (k=1; k<DIM-1-8; k+=8) {
	aggregate += old[ijloc+k]*func_val_need/func_val_gimme;
	aggregate += old[ijloc+k+1]*func_val_need/func_val_gimme;
	aggregate += old[ijloc+k+2]*func_val_need/func_val_gimme;
	aggregate += old[ijloc+k+3]*func_val_need/func_val_gimme;
	aggregate += old[ijloc+k+4]*func_val_need/func_val_gimme;
	aggregate += old[ijloc+k+5]*func_val_need/func_val_gimme;
	aggregate += old[ijloc+k+6]*func_val_need/func_val_gimme;
	aggregate += old[ijloc+k+7]*func_val_need/func_val_gimme;	
      }
      for (; k<DIM-1; k++) {
	aggregate+= old[ijloc+k]*func_val_need/func_val_gimme;
	}
      
    }
  }

  printf("AGGR:%ld\n",aggregate);

  long holderval, holderval2, holderval3;
  int loc1i, loc2i, loc3i;
  long bigtemp, bigtemp2, bigtemp3;
  for (i=1; i<DIM-1; i++) {
    loc1i = (i-1)*DIMsq;
    loc2i = i*DIMsq;
    loc3i = (i+1)*DIMsq;
    for (j=1; j<DIM-1; j++) {
      for (k=1; k<DIM-1-4; k+=3) {
	holderval=0;
	holderval2 = 0;
	holderval3 = 0;
	//set 1
	holderval+=old[loc1i+(j-1)*DIM+(k-1)];
	bigtemp=old[loc1i+(j-1)*DIM+(k)];
	bigtemp2=old[loc1i+(j-1)*DIM+(k+1)];
	holderval += bigtemp + bigtemp2;
	bigtemp3 = old[loc1i+(j-1)*DIM+(k+2)];
	holderval2 += bigtemp + bigtemp2 + bigtemp3;
	holderval3 += bigtemp2 + bigtemp3;
	holderval3 += old[loc1i + (j-1)*DIM + (k+3)];
	
	holderval+=old[loc1i+(j+0)*DIM+(k-1)];
	bigtemp=old[loc1i+(j+0)*DIM+(k)];
	bigtemp2=old[loc1i+(j+0)*DIM+(k+1)];
	holderval += bigtemp + bigtemp2;
	bigtemp3 = old[loc1i+j*DIM+(k+2)];
	holderval2 += bigtemp + bigtemp2 + bigtemp3;
	holderval3 += bigtemp2 + bigtemp3;
	holderval3 += old[loc1i + j*DIM + (k+3)];
	
	holderval+=old[loc1i+(j+1)*DIM+(k-1)];
	bigtemp=old[loc1i+(j+1)*DIM+(k)];
	bigtemp2=old[loc1i+(j+1)*DIM+(k+1)];
	holderval += bigtemp + bigtemp2;
	bigtemp3 = old[loc1i + (j+1)*DIM + (k+2)];
	holderval2 += bigtemp + bigtemp2 + bigtemp3;
	holderval3 += bigtemp2 + bigtemp3;
	holderval3 += old[loc1i + (j+1)*DIM + (k+3)];
	
	//set 2
	holderval+=old[loc2i+(j-1)*DIM+(k-1)];
	bigtemp=old[loc2i+(j-1)*DIM+(k)];
	bigtemp2=old[loc2i+(j-1)*DIM+(k+1)];
	holderval += bigtemp + bigtemp2;
	bigtemp3 = old[loc2i + (j-1)*DIM + (k+2)];
	holderval2 += bigtemp + bigtemp2 + bigtemp3;
	holderval3 += bigtemp2 + bigtemp3;
	holderval3 += old[loc2i + (j-1)*DIM + (k+3)];
	
	holderval+=old[loc2i+(j+0)*DIM+(k-1)];
	bigtemp=old[loc2i+(j+0)*DIM+(k)];
	bigtemp2=old[loc2i+(j+0)*DIM+(k+1)];
	holderval += bigtemp + bigtemp2;
	bigtemp3 = old[loc2i + j*DIM + (k+2)];
	holderval2 += bigtemp + bigtemp2 + bigtemp3;
	holderval3 += bigtemp2 + bigtemp3;
	holderval3 += old[loc2i + j*DIM + (k+3)];
	
	holderval+=old[loc2i+(j+1)*DIM+(k-1)];
	bigtemp=old[loc2i+(j+1)*DIM+(k)];
	bigtemp2=old[loc2i+(j+1)*DIM+(k+1)];
	holderval += bigtemp + bigtemp2;
	bigtemp3 = old[loc2i + (j+1)*DIM + (k+2)];
	holderval2 += bigtemp + bigtemp2 + bigtemp3;
	holderval3 += bigtemp2 + bigtemp3;
	holderval3 += old[loc2i + (j+1)*DIM + (k+3)];
	
	//set 3
	holderval+=old[loc3i+(j-1)*DIM+(k-1)];
	bigtemp=old[loc3i+(j-1)*DIM+(k)];
	bigtemp2=old[loc3i+(j-1)*DIM+(k+1)];
	holderval += bigtemp + bigtemp2;
	bigtemp3 = old[loc3i + (j-1)*DIM + (k+2)];
	holderval2 += bigtemp + bigtemp2 + bigtemp3;
	holderval3 += bigtemp2 + bigtemp3;
	holderval3 += old[loc3i + (j-1)*DIM + (k+3)];
	
	holderval+=old[loc3i+(j+0)*DIM+(k-1)];
	bigtemp=old[loc3i+(j+0)*DIM+(k)];
	bigtemp2=old[loc3i+(j+0)*DIM+(k+1)];
	holderval += bigtemp + bigtemp2;
	bigtemp3 = old[loc3i + j*DIM + (k+2)];
	holderval2 += bigtemp + bigtemp2 + bigtemp3;
	holderval3 += bigtemp2 + bigtemp3;
	holderval3 += old[loc3i + j*DIM + (k+3)];
	
	holderval+=old[loc3i+(j+1)*DIM+(k-1)];
	bigtemp=old[loc3i+(j+1)*DIM+(k)];
	bigtemp2=old[loc3i+(j+1)*DIM+(k+1)];
	holderval += bigtemp + bigtemp2;
	bigtemp3 = old[loc3i + (j+1)*DIM + (k+2)];
	holderval2 += bigtemp + bigtemp2 + bigtemp3;
	holderval3 += bigtemp2 + bigtemp3;
	holderval3 += old[loc3i + (j+1)*DIM + (k+3)];
	
	holderval/=27;
	holderval2/=27;
	holderval3/=27;
	
	new[loc2i+j*DIM+k]=holderval;
	new[loc2i+j*DIM+k+1]=holderval2;
	new[loc2i+j*DIM+k+2]=holderval3;
      }
      for(; k < DIM-1; k++) {
	holderval=0;
	holderval+=old[loc1i+(j-1)*DIM+(k-1)];
	holderval+=old[loc1i+(j-1)*DIM+(k)];
	holderval+=old[loc1i+(j-1)*DIM+(k+1)];
	holderval+=old[loc1i+(j)*DIM+(k-1)];
	holderval+=old[loc1i+(j)*DIM+(k)];
	holderval+=old[loc1i+(j)*DIM+(k+1)];
	holderval+=old[loc1i+(j+1)*DIM+(k-1)];
	holderval+=old[loc1i+(j+1)*DIM+(k)];
	holderval+=old[loc1i+(j+1)*DIM+(k+1)];

	holderval+=old[loc2i+(j-1)*DIM+(k-1)];
	holderval+=old[loc2i+(j-1)*DIM+(k)];
	holderval+=old[loc2i+(j-1)*DIM+(k+1)];
	holderval+=old[loc2i+(j)*DIM+(k-1)];
	holderval+=old[loc2i+(j)*DIM+(k)];
	holderval+=old[loc2i+(j)*DIM+(k+1)];
	holderval+=old[loc2i+(j+1)*DIM+(k-1)];
	holderval+=old[loc2i+(j+1)*DIM+(k)];
	holderval+=old[loc2i+(j+1)*DIM+(k+1)];

	holderval+=old[loc3i+(j-1)*DIM+(k-1)];
	holderval+=old[loc3i+(j-1)*DIM+(k)];
	holderval+=old[loc3i+(j-1)*DIM+(k+1)];
	holderval+=old[loc3i+(j)*DIM+(k-1)];
	holderval+=old[loc3i+(j)*DIM+(k)];
	holderval+=old[loc3i+(j)*DIM+(k+1)];
	holderval+=old[loc3i+(j+1)*DIM+(k-1)];
	holderval+=old[loc3i+(j+1)*DIM+(k)];
	holderval+=old[loc3i+(j+1)*DIM+(k+1)];

	holderval/=27;
	new[loc2i+j*DIM+k]=holderval;
      }
    }
  }
  long u0, u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12, u13, u14, u15;
  for (i=1; i<DIM-1; i++) {
    iloc = i*DIMsq;
    for (j=1; j<DIM-1; j++) {
      ijloc = iloc + j*DIM;
      for (k=1; k<DIM-1-16; k+=16) {
        u0=(new[ijloc+k]/100);
	u1=(new[ijloc+k+1]/100);
	u2=(new[ijloc+k+2]/100);
	u3=(new[ijloc+k+3]/100);
	u4=(new[ijloc+k+4]/100);
	u5=(new[ijloc+k+5]/100);
	u6=(new[ijloc+k+6]/100);
	u7=(new[ijloc+k+7]/100);
	u8=(new[ijloc+k+8]/100);
	u9=(new[ijloc+k+9]/100);
	u10=(new[ijloc+k+10]/100);
	u11=(new[ijloc+k+11]/100);
	u12=(new[ijloc+k+12]/100);
	u13=(new[ijloc+k+13]/100);
	u14=(new[ijloc+k+14]/100);
	u15=(new[ijloc+k+15]/100);

	if (u0<=0) {u0=0;}
	if (u0>=9) u0=9;
	
	if (u1<=0) {u1=0;}
	if (u1>=9) u1=9;
	
	if (u2<=0) {u2=0;}
	if (u2>=9) u2=9;

	if (u3<=0) {u3=0;}
	if (u3>=9) u3=9;

	if (u4<=0) {u4=0;}
	if (u4>=9) u4=9;

	if (u5<=0) {u5=0;}
	if (u5>=9) u5=9;

	if (u6<=0) {u6=0;}
	if (u6>=9) u6=9;

	if (u7<=0) {u7=0;}
	if (u7>=9) u7=9;

	if (u8<=0) {u8=0;}
	if (u8>=9) u8=9;

	if (u9<=0) {u9=0;}
	if (u9>=9) u9=9;

	if (u10<=0) {u10=0;}
	if (u10>=9) u10=9;

	if (u11<=0) {u11=0;}
	if (u11>=9) u11=9;

	if (u12<=0) {u12=0;}
	if (u12>=9) u12=9;

	if (u13<=0) {u13=0;}
	if (u13>=9) u13=9;

	if (u14<=0) {u14=0;}
	if (u14>=9) u14=9;

	if (u15<=0) {u15=0;}
	if (u15>=9) u15=9;
	
        histogrammy[u0]++;
	histogrammy[u1]++;
	histogrammy[u2]++;
	histogrammy[u3]++;
	histogrammy[u4]++;
	histogrammy[u5]++;
	histogrammy[u6]++;
	histogrammy[u7]++;
	histogrammy[u8]++;
	histogrammy[u9]++;
	histogrammy[u10]++;
	histogrammy[u11]++;
	histogrammy[u12]++;
	histogrammy[u13]++;
	histogrammy[u14]++;
	histogrammy[u15]++;
      }
      for(; k < DIM-1; k++) {
	u=(new[iloc+j*DIM+k]/100);
	if(u<=0) u = 0;
	if(u>=9) u = 9;
	histogrammy[u]++;
      }
    }
  }
  return (double) (dot_product+moving_average+pi+pi2);


}
