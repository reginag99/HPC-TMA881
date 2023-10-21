#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <threads.h>
#include <complex.h>
#include <stdint.h>

#define PI 3.14159 //265358979323846;
typedef int8_t TYPE_ATTR;
typedef uint8_t TYPE_CONV;

TYPE_ATTR ** attractors;
TYPE_CONV ** convergences;

void FunctionDerivate(float complex *z,float complex *derivate, int d);
void GetRoots(float ** roots,  int d);
void GetColors(char**color, int d);
void GetGrayScale(char**grayscale);

typedef struct {
  int val;
  char pad[60];
} int_padded;

typedef struct {
  TYPE_ATTR **attractors;
  TYPE_CONV **convergences;
  float**roots;
  int numThreads;
  int d;
  int lines;
  int tx;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_t;

typedef struct {
  TYPE_ATTR **attractors;
  TYPE_CONV **convergences;
  int lines;
  int d;
  int numThreads;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_check_t;

int main_thrd(void *args)
{
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  TYPE_ATTR **attractors = thrd_info->attractors;
  TYPE_CONV **convergences = thrd_info->convergences;
  float**roots = thrd_info->roots;
  const int numThreads = thrd_info->numThreads;
  int d = thrd_info->d;
  const int lines = thrd_info->lines;
  const int tx = thrd_info->tx;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;

  int i, j, k;
  float ix, jx;
  float realRoot, imagRoot;
  float complex z,derivate,derivateConj,functionValue,derivateAbs;

  float c0 = 1. - 1./d;
  float c1 = 1./d;

  for (i = tx;  i < lines; i += numThreads ) {
    ix = 2.0 - (4.0 * i /(lines - 1));

    attractors[i] = (TYPE_ATTR*) malloc(lines*sizeof(TYPE_ATTR));
    convergences[i] = (TYPE_CONV*) malloc(lines*sizeof(TYPE_CONV));
    TYPE_ATTR*attractor = attractors[i];
    TYPE_CONV*convergence = convergences[i];

    for ( size_t cx = 0; cx < lines; ++cx ) {
      attractor[cx] = 11;
      convergence[cx] = 128;
    }

    for (j = 0; j < lines; ++j){
      jx = -2.0 + (4.0 * j / (lines - 1));
      z = jx + ix * I;
      for (k = 0; k < 128; k++){
	if (crealf(z)*crealf(z) + cimagf(z)*cimagf(z) < 0.000001){
	 convergence[j] = k;
	  break;}
	if (crealf(z) > 10000000000.||crealf(z) < -10000000000.){
	  break;}
	if (cimagf(z) > 10000000000.||cimagf(z) < -10000000000.){
	  break;}

	FunctionDerivate(&z,&derivate, d);
	functionValue = derivate * z - 1;

	if (crealf(functionValue)*crealf(functionValue) + cimagf(functionValue)*cimagf(functionValue) < 0.0000001){
	  convergence[j] = k;
	  for(int ixd = 0; ixd < d; ixd++){
	     if((crealf(z) <=  (roots[ixd][0] + 0.001) && crealf(z) >=  (roots[ixd][0] - 0.001)) &&
                (cimagf(z) <=  (roots[ixd][1] + 0.001) && cimagf(z) >=  (roots[ixd][1] - 0.001))){
                           attractor [j] = ixd;
                           break;
            }
	  }
	  break;
	}
	derivateConj = conj(derivate);
	derivateAbs = derivate * derivateConj;

	z = c0*z + c1 * derivateConj/derivateAbs;
      }
    }

    mtx_lock(mtx);
    status[tx].val = i + numThreads;
    mtx_unlock(mtx);
    cnd_signal(cnd);
  }
  return 0;
}

int
main_thrd_check(void *args)
{
  const thrd_info_check_t *thrd_info = (thrd_info_check_t*) args;
  TYPE_ATTR **attractors = thrd_info->attractors;
  TYPE_CONV **convergences = thrd_info->convergences;
  const int lines = thrd_info->lines;
  const int d = thrd_info->d;
  const int numThreads = thrd_info->numThreads;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;

  int grayIndex, colorIndex;

  char**color = (char**) malloc(sizeof(char*)*12);

  GetColors(color, d);

  char**grayscale = (char**) malloc(sizeof(char*)*129);

  GetGrayScale(grayscale);

  char filenameA[30];
  char filenameC[30];

  snprintf(filenameA, sizeof(filenameA), "newton_attractors_x%d.ppm", d);
  snprintf(filenameC, sizeof(filenameC), "newton_convergence_x%d.ppm", d);

  FILE *fileColor = fopen(filenameA, "w");
  if (fileColor == NULL) {
    perror("Error opening the file");
    exit(EXIT_FAILURE);
  }

  fprintf(fileColor, "P3\n");
  fprintf(fileColor, "%d %d \n", lines, lines);
  fprintf(fileColor, "255\n");

  FILE *fileGray = fopen(filenameC, "w");
  if (fileGray == NULL) {
    perror("Error opening the file");
    exit(EXIT_FAILURE);
  }

  fprintf(fileGray, "P3\n");
  fprintf(fileGray, "%d %d \n", lines, lines);
  fprintf(fileGray, "129\n");

  char * stringGray = (char*)malloc(12*lines*sizeof(char));
  char * stringColor = (char*)malloc(12*lines*sizeof(char));

  for ( int ix = 0, ibnd; ix < lines; ) {
    for ( mtx_lock(mtx); ; ) {
      ibnd = lines;

      for ( int tx = 0; tx < numThreads; ++tx )
	if ( ibnd > status[tx].val )
	  ibnd = status[tx].val;
      if ( ibnd <= ix )
	cnd_wait(cnd,mtx);
      else {
	mtx_unlock(mtx);
	break;
      }
    }

    for ( ; ix < ibnd; ++ix ) {
      TYPE_ATTR * attractor_line = attractors[ix];
      TYPE_CONV * convergence_line = convergences[ix];
      for(int jx = 0; jx < lines; ++jx){
	memcpy(stringColor + jx *12, color[attractor_line[jx]], sizeof(char)*12);
	memcpy(stringGray + jx *12, grayscale[convergence_line[jx]], sizeof(char)*12);
      }
      stringColor[12*lines-1] = '\n';
      stringGray[12*lines-1] = '\n';

      fwrite(stringColor, sizeof(char), lines * 12,fileColor);
      fwrite(stringGray, sizeof(char), lines*12, fileGray);

      free(attractor_line);
      free(convergence_line);
    }
  }

  fclose(fileGray);
  fclose(fileColor);
  free(color);
  free(grayscale);

  return 0;
}

  int
main(int argc, char *argv[])
{
  int lines, numThreads, d;
  if (argc != 3 && argc != 4) {
    printf("Two or three arguments expected");
    exit(EXIT_FAILURE);
  }

  for (int i = 1; i < 4; i++) {
    if (strncmp(argv[i], "-t", 2) == 0) {
      numThreads = atoi(argv[i] + 2);
    } else if (strncmp(argv[i], "-l", 2) == 0) {
      lines = atoi(argv[i] + 2);
    } else if (strncmp(argv[i], "", 0) == 0) {
      d = atoi(argv[i] + 0);
    } else {
      printf("Wrong format\n");
      return 1;
    }
  }

  float *rootsEntries = (float*) malloc(sizeof(float)* d * 2);
  float ** roots = (float**) malloc(sizeof(float*)*d);

  for ( size_t ix = 0, jx = 0; ix < d; ++ix, jx+=2)
    roots[ix] = rootsEntries + jx;

  for ( size_t ix = 0; ix < d; ++ix )
    for ( size_t jx = 0; jx < 2; ++jx )
      roots[ix][jx] = 0;

  GetRoots(roots, d);

  float stepSize = 1/lines;

  TYPE_ATTR**attractors = (TYPE_ATTR **) malloc(lines*sizeof(TYPE_ATTR *));
  TYPE_CONV**convergences = (TYPE_CONV **) malloc(lines*sizeof(TYPE_CONV *));

  thrd_t thrds[numThreads];
  thrd_info_t thrds_info[numThreads];

  thrd_t thrd_check;
  thrd_info_check_t thrd_info_check;

  mtx_t mtx;
  mtx_init(&mtx, mtx_plain);

  cnd_t cnd;
  cnd_init(&cnd);

  int_padded status[numThreads];

  for ( int tx = 0; tx < numThreads; ++tx ) {
    thrds_info[tx].attractors = attractors;
    thrds_info[tx].convergences = convergences;
    thrds_info[tx].roots=roots;
    thrds_info[tx].numThreads = numThreads;
    thrds_info[tx].d = d;
    thrds_info[tx].lines = lines;
    thrds_info[tx].tx = tx;
    thrds_info[tx].mtx = &mtx;
    thrds_info[tx].cnd = &cnd;
    thrds_info[tx].status = status;
    status[tx].val = -1;

    int r = thrd_create(thrds+tx, main_thrd, (void*) (thrds_info+tx));

    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
  }

  {
    thrd_info_check.attractors = attractors;
    thrd_info_check.convergences = convergences;
    thrd_info_check.lines = lines;
    thrd_info_check.d = d;
    thrd_info_check.numThreads = numThreads;
    thrd_info_check.mtx = &mtx;
    thrd_info_check.cnd = &cnd;
    thrd_info_check.status = status;

    int r = thrd_create(&thrd_check, main_thrd_check, (void*) (&thrd_info_check));
    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
  }

  {
    int r;
    thrd_join(thrd_check, &r);
  }

  free(attractors);
  free(convergences);
  mtx_destroy(&mtx);
  cnd_destroy(&cnd);

  return 0;

}

void GetRoots(float ** roots, int d) {
  for( size_t ix = 0; ix < d; ix++) {
    float theta = (ix*2* PI) / d ;
    roots[ix][0] = cos(theta);
    roots[ix][1] = sin(theta);
  }
}

void FunctionDerivate(float complex *z,float complex *derivate, int d){
  switch (d) {
    case 1:
      *derivate = 1;
      break;
    case 2:
      *derivate = (*z);
      break;
    case 3:
      *derivate = (*z)*(*z);
      break;
    case 4:
      *derivate = (*z)*(*z)*(*z);
      break;
    case 5:
      *derivate = ((*z)*(*z))*((*z)*(*z));
      break;
    case 6:
      *derivate = ((*z)*(*z)*(*z))*((*z)*(*z));
      break;
    case 7:
      *derivate = ((*z)*(*z)*(*z))*((*z)*(*z)*(*z));
      break;
    case 8:
      *derivate = ((*z)*(*z)*(*z))*((*z)*(*z)*(*z)*(*z));
      break;
    case 9:
      *derivate = ((*z)*(*z)*(*z)*(*z))*((*z)*(*z)*(*z)*(*z));
      break;
    case 10:
      *derivate = ((*z)*(*z)*(*z)*(*z))*((*z)*(*z)*(*z)*(*z)*(*z));
      break;

    default:
      fprintf(stderr, "unexpected degree\n");
      exit(1);
  }
}

void GetColors(char**color, int d){
  color[0] = "210 116 188 ";
  color[1] = "154 116 210 ";
  color[2] = "166 160 210 ";
  color[3] = "201 146 116 ";
  color[4] = "116  58 155 ";
  color[5] = "245 179   0 ";
  color[6] = "0   204 245 ";
  color[7] = "245  82   0 ";
  color[8] = "255 172 221 ";
  color[9] = "172 255 236 ";
  color[10] = "194 255 172 ";
  color[11] = "  0   0   0 ";
}

void GetGrayScale(char** grayscale) {
  int grayValue;
  char grayElement[13];

  for (size_t ix = 0; ix < 129; ix++) {
    grayValue = ix;
    snprintf(grayElement, sizeof(grayElement), "%3d %3d %3d ", grayValue, grayValue, grayValue);
    grayscale[ix] = (char*)malloc(12*sizeof(char));
    memcpy(grayscale[ix], grayElement, 12);
  }
}









