#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void psum1(float a[], float p[], long n)
{
  long i;
  p[0] = a[0];
  for(i = 1; i < n; i++)
    p[i] = p[i-1] + a[i];
}

void psum2(float a[], float p[], long n)
{
  long i;
  p[0] = a[0];
  for(i = 1; i < n-1; i+=2){
    float mid_val = p[i-1] + a[i];
    p[i] = mid_val;
    p[i+1] = mid_val + a[i+1];
  }
  if(i < n)
    p[i] = p[i-1] + a[i];
}

//5.12
void psum1a(float a[], float p[], long n)
{
  long i;
  /*last_val holds p[i-1], val holds p[i]*/
  float last_val, val;
  last_val = p[0] = a[0];
  for(i = 1; i < n; i++) {
    val = last_val + a[i];
    p[i] = val;
    last_val = val;
  }
}

//5.19 now
void psum3(float a[], float p[], long n)
{
  long i;
  float acc, val1, val2;
  acc = p[0] = a[0];
  for(i = 1; i < n-1; i+=2)
    {
      val1 = a[i];
      val2 = a[i+1];
      p[i] = acc + val1;
      p[i+1] = (acc += (val1 + val2));
    }
  for(; i < n; i++)
    {
      p[i] = (acc += a[i]);
    }
}

void psum3_v2(float a[], float p[], long n)
{
  long i;
  float acc, val1, val2, val3;
  acc = p[0] = a[0];
  for(i = 1; i < n-2; i+=3)
    {
      val1 = a[i];
      val2 = a[i+1];
      val3 = a[i+2];
      acc += val1;
      p[i] = acc;
      p[i+1] = (acc += val2);
      p[i+2] = (acc += val3);
    }
  for(; i < n; i++)
    {
      p[i] = (acc += a[i]);
    }
}

void psum3_v3(float a[], float p[], long n)
{
  long i;
  float acc, val1, val2, val3, val4;
  acc = p[0] = a[0];
  for(i = 1; i < n-3; i+=4)
    {
      val1 = a[i];
      val2 = a[i+1];
      val3 = a[i+2];
      val4 = a[i+3];
      acc += val1;
      p[i] = acc;
      p[i+1] = (acc + (val1 + val2));
      p[i+2] = (acc + (val1 + val2 + val3));
      p[i+3] = (acc += (val1 + val2 + val3 + val4));
    }
  for(; i < n; i++)
    {
      p[i] = (acc += a[i]);
    }
}

int main(){
  float destArr[100000];
  float arr1[100000];
  float arr2[100000];
  float arr3[100000];
  float arr4[100000];
  float arr5[100000];
  float arr6[100000];
  for(long i = 0; i < 100000; i++) {
    arr1[i] = 1;
    arr2[i] = 1;
    arr3[i] = 1;
    arr4[i] = 1;
    arr5[i] = 1;
    arr6[i] = 1;
  }


  clock_t start1, end1, start2, end2;

  start1 = clock();
  //  printf("start 1, %ld\n", start1);
  psum1(arr1, destArr, 100000);
  // for(int i = 0; i < (1<<31); i++) {
  // }
  end1 = clock();
  //printf("end 1, %ld\n", end1);
  printf("total 1: %ld\n", end1-start1);
  //printf("random thing: %d\n", destArr[1000]);

  start1 = clock();
  //printf("start 2, %ld\n", start1);
  psum2(arr2, destArr, 100000);
  // for(int i = 0; i < (1<<31); i++) {
  // }
  end1 = clock();
  //printf("end 2, %ld\n", end1);
  printf("total 2: %ld\n", end1-start1);
  //printf("random thing: %d\n", destArr[1000]);

  start1 = clock();
  //printf("start 3, %ld\n", start1);
  psum1a(arr3, destArr, 100000);
  // for(int i = 0; i < (1<<31); i++) {
  // }
  end1 = clock();
  //printf("end 3, %ld\n", end1);
  printf("total 3: %ld\n", end1-start1);
  //printf("random thing: %d\n", destArr[1000]);
  
  start2 = clock();
  //printf("start 4, %ld\n", start2);
  psum3(arr4, destArr, 100000);
  // for(int i = 0; i < (1<<31); i++) {
  // }
  end2 = clock();
  //printf("end 4, %ld\n", end2*10000);
  printf("total 4: %ld\n", end2-start2);
  //printf("random thing: %d\n", destArr[1000]);

  start2 = clock();
  //printf("start 5, %ld\n", start2);
  psum3_v2(arr5, destArr, 100000);
  // for(int i = 0; i < (1<<31); i++) {                                                                            
  // }                                                                                                             
  end2 = clock();
  //printf("end 5, %ld\n", end2*10000);
  printf("total 5: %ld\n", end2-start2);
  //printf("random thing: %d\n", destArr[1000]);

  start2 = clock();
  //printf("start 6, %ld\n", start2);
  psum3_v3(arr6, destArr, 100000);
  // for(int i = 0; i < (1<<31); i++) {                                                                           \
                                                                                                                   
  // }                                                                                                            \
                                                                                                                   
  end2 = clock();
  //printf("end 6, %ld\n", end2*10000);
  printf("total 6: %ld\n", end2-start2);
  //printf("random thing: %d\n", destArr[1000]);
}
