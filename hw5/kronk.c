#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//#define SIZE 40
//#define KSIZE SIZE*SIZE

int main()
{
  int SIZE = 30;
  int KSIZE = SIZE*SIZE;
  int a[SIZE][SIZE];
  int b[SIZE][SIZE];
  //  int *c = malloc(sizeof(*c)*KSIZE*KSIZE);
  int c[KSIZE][KSIZE];
  clock_t start, end;

  int i, j, x, y;
  int c1, c2;
  printf("%d\n", SIZE);
  for(i = 0; i < SIZE; i++) {
    for(j = 0; j < SIZE; j++) {
      a[i][j] = 1;
      b[i][j] = 1;
    }
  }
  
  for(i = 0; i < KSIZE; i++) {
    for(j = 0; j < KSIZE; j++) {
      c[i][j] = 1;
    }
  }
  
  start = clock();
  for (i=0; i<SIZE; i++)
    for (j=0; j<SIZE; j++)
      for (x=0; x<SIZE; x++)
	for (y=0; y<SIZE; y++)
	  {
	    //c1 = i*SIZE + x;
	    //c2 = j*SIZE + y;
	    c[i*SIZE + x][j*SIZE + y]=a[i][j]*b[x][y]; 
	  }
  end = clock();
  printf("ijxy: %ld\n", end-start);
  
  start= clock();
  for (j=0; i<SIZE; i++)
    for (i=0; j<SIZE; j++)
      for (x=0; x<SIZE; x++)
	for (y=0; y<SIZE; y++)
          c[i*SIZE+x][j*SIZE+y]=a[i][j]*b[x][y];

  end =clock();
  printf("jixy: %ld\n",end-start);

  start= clock();
  for (x=0; i<SIZE; i++)
    for (i=0; j<SIZE; j++)
      for (y=0; x<SIZE; x++)
	for (j=0; y<SIZE; y++)
          c[i*SIZE+x][j*SIZE+y]=a[i][j]*b[x][y];

  end =clock();
  printf("xiyj: %ld\n",end-start);
  printf("this shit broken]n");
  start= clock();
  for (i=0; i<SIZE; i++)
    for (x=0; j<SIZE; j++)
      for (j=0; x<SIZE; x++)
	for (y=0; y<SIZE; y++)
	  {
	    printf("hi");
	    printf("%d ", i);
          c[i*SIZE+x][j*SIZE+y]=a[i][j]*b[x][y];
	  }
  end =clock();
  printf("ixjy: %ld\n",end-start);

  start= clock();
  for (j=0; i<SIZE; i++)
    for (y=0; j<SIZE; j++)
      for (i=0; x<SIZE; x++)
	for (x=0; y<SIZE; y++)
          c[i*SIZE+x][j*SIZE+y]=a[i][j]*b[x][y];

  end =clock();
  printf("jyix: %ld\n",end-start);
}
