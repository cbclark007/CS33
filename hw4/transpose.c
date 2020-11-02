//Christopher Clark
//305326742
//Tested this with sizes 2, 4, 8, 16, and 32. 16 was the fastest

void transpose_5(int *dst, int *src, int dim)
{
  int i, j;
  for(i = 0; i < dim-16; i += 16) {
    for(j = 0; j < dim-16; j+= 16) {
      for(int x = i; x < i+16; x++) {
        for(int y = j; y < j+16; y++) {
          dst[y * dim + x] = src[x*dim + y];
        }
      }
    }
  }

  for(; i < dim; i++) {
    for(j = 0; j < dim-16; j+=16) {
      dst[j*dim +i] = src[i*dim + j];
      dst[(j+1)*dim + i] = src[i*dim + j+1];
      dst[(j+2)*dim + i] = src[i*dim + j+2];
      dst[(j+3)*dim + i] = src[i*dim + j+3];
      dst[(j+4)*dim + i] = src[i*dim + j+4];
      dst[(j+5)*dim + i] = src[i*dim + j+5];
      dst[(j+6)*dim + i] = src[i*dim + j+6];
      dst[(j+7)*dim + i] = src[i*dim + j+7];
      dst[(j+8)*dim + i] = src[i*dim + j+8];
      dst[(j+9)*dim + i] = src[i*dim + j+9];
      dst[(j+10)*dim + i] = src[i*dim + j+10];
      dst[(j+11)*dim + i] = src[i*dim + j+11];
      dst[(j+12)*dim + i] = src[i*dim + j+12];
      dst[(j+13)*dim + i] = src[i*dim + j+13];
      dst[(j+14)*dim + i] = src[i*dim + j+14];
      dst[(j+15)*dim + i] = src[i*dim + j+15];
    }
  }

  for(i=0; i < dim; i++) {
    for(j = dim-16; j < dim; j++) {
      dst[j*dim +i] = src[i*dim + j];
    }
  }
}
