#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <threads.h>
#include <complex.h>
#include <stdint.h>

#define PI 3.14159 //265358979323846; NOTE: Hur många värdesiffror?
typedef char TYPE_ATTR; //This should be char
typedef unsigned char TYPE_CONV; //This should be char

TYPE_ATTR ** attractors;
TYPE_CONV ** convergences;

void FunctionDerivate(double complex z,double complex derivate, int d);
void GetRoots( float ** roots,  int d);
void GetColors(char**color, int d);
void GetGrayScale(char**grayscale, int d);

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
  int numThreads = thrd_info->numThreads;
  int d = thrd_info->d;
  const int lines = thrd_info->lines;
  const int tx = thrd_info->tx;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;

  int i, j, k;
  double ix, jx;
  double complex z,derivate,derivateConj, functionValue, denom;

  //Warn: Vi tror koordinaterna funkar som de ska ej 100% säkra, kolla så att itereringen är rätt isf
  //Fel typ?? inget större än en character
  for (i = tx;  i < lines; i += numThreads ) {
    ix = 2.0 - (4.0 * i /(lines - 1));
    TYPE_ATTR*attractor = (TYPE_ATTR*) malloc(lines*sizeof(TYPE_ATTR)); //Hur blir minneshantering här?
    TYPE_CONV*convergence = (TYPE_CONV*) malloc(lines*sizeof(TYPE_CONV)); //Hur blir minneshantering här?

    for ( size_t cx = 0; cx < lines; ++cx ) {
        attractor[cx] = d;
        convergence[cx] = 128;
        }

    for (j = 0; j < lines; ++j){
       jx = -2.0 + (4.0 * j / (lines - 1));
       z = ix + jx * I;
       //printf("Complex number z: %lf + %lfi\n", creal(z), cimag(z));
       for (k = 0; k < 128; k++){
            if (creal(z) < 0.001 && creal(z) > -0.001 || cimag(z) < 0.001  && cimag(z) > -0.001 )
                break;
            if (creal(z) > 1000000000 && creal(z) < -1000000000  || cimag(z) > 1000000000 && cimag(z) < -1000000000)
                break;

            FunctionDerivate(z,derivate, d);
            derivateConj = conj(derivate);
            functionValue = derivate * z - 1;

            //För dyrt med absolutbelopp?
            if (cabs(functionValue) < 0.001*d && cabs(functionValue) > -0.001*d){
                convergence[j] = k;
                for(int ixd = 0; ixd < d; ixd++){
                    if((creal(z) <=  (roots[ixd][0] + 0.001) || creal(z) >=  (roots[ixd][0] - 0.001)) &&
                     (cimag(z) <=  (roots[ixd][1] + 0.001) || cimag(z) >=  (roots[ixd][1] - 0.001))){
                           attractor [j] = ixd;
                           break;}
                }
             break;
            }
            denom = derivateConj * derivate;
            z = z - (functionValue*derivateConj)/denom; //divison med komplext tal undvika
            //printf("lf",z);
       }
       //printf("%d", convergence[j]);
       //fflush(stdout);
      }

    //WARN: Är data transfer korrekt?
    mtx_lock(mtx);
    attractors[i] = attractor;
    convergences[i] = convergence;
    status[tx].val = i + tx;
    mtx_unlock(mtx);
    cnd_signal(cnd);
    //thrd_sleep(&(struct timespec){.tv_sec=0, .tv_nsec=10}, NULL);
   }
  return 0;
}

int
main_thrd_check(
    void *args
    )
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

  char*colorEntries = (char*) malloc(sizeof(char)* 12 * 11);
  char**color = (char**) malloc(sizeof(char*)*11);

  for ( size_t ix = 0, jx = 0; ix < 11; ++ix, jx+=21)
        color[ix] = colorEntries + jx;

  GetColors(color, d);
  printf("color = %d\n",color[0][0]);
 
  char*grayscaleEntries = (char*) malloc(sizeof(char)* 12 * 128);
  char**grayscale = (char**) malloc(sizeof(char*)*128);

  for ( size_t ix = 0, jx = 0; ix < 128; ++ix, jx+=12)
        grayscale[ix] = grayscaleEntries + jx;

  GetGrayScale(grayscale, d);
  printf("gray = %d\n",grayscale[0][0]);

  char*pixelBufferGrayEntries = (char*) malloc(sizeof(char)*12*lines);
  char**pixelBufferGray = (char**) malloc(sizeof(char*)* lines);

  for ( size_t ix = 0, jx = 0; ix < lines; ++ix, jx+= 12)
    pixelBufferGray[ix] = pixelBufferGrayEntries + jx;

  char*pixelBufferColorEntries = (char*) malloc(sizeof(char)*12*lines);
  char**pixelBufferColor = (char**) malloc(sizeof(char*)*12);

  for ( size_t ix = 0, jx = 0; ix < lines; ++ix, jx+= 12)
    pixelBufferColor[ix] = pixelBufferColorEntries + jx;

  int grayIndex, colorIndex;
  char elementColor[12];
  char elementGrayScale[12];

  //Vart vill vi öppna filen?
  FILE *fileColor = fopen("output.ppm", "w");
  if (fileColor == NULL) {
    perror("Error opening the file");
    exit(EXIT_FAILURE);
    }

  fprintf(fileColor, "P3\n");
  fprintf(fileColor, "%d %d \n", lines, lines);
  fprintf(fileColor, "255\n");

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

    char stringGray[12*lines*sizeof(char)];
    char stringColor[12*lines*sizeof(char)];

    //fprintf(stderr, "checking until %i\n", ibnd);

    for ( ; ix < ibnd; ++ix ) {
        for(int jx = 0; jx < lines; ++jx){
           //grayIndex = attractors[ix][jx];
           colorIndex = (int)attractors[ix][jx];
           memcpy(stringColor,color[colorIndex], sizeof(char)*12);
           //memcpy(stringGray,grayscale[grayIndex], sizeof(char)*12);
        }

    //fprintf(file_gray, "P3\n");
    //fprintf(file_gray, "%d %d \n", lines, lines);
    //fprintf(file_gray, "255\n");

    //fwrite(string_gray, sizeof(TYPE_CONV), 4*3*lines, file); 

    fwrite(stringColor, sizeof(char), 4*3*lines, fileColor);
     //memcpy(stringGray, pixelBufferGray[ix], sizeof(char)*12);

   // fwrite(stringColor, sizeof(TYPE_CONV), lines, file);

    //fwrite(pixelBuffer, sizeof(char), strlen(pixelBuffer), file);
    //Måste vi tömma pixelbuffer efter varje iteration?

    free(attractors[ix]);
    free(convergences[ix]);
  }
 }

 fclose(fileColor);
 free(colorEntries);
 free(color);
 free(pixelBufferGray);
 free(pixelBufferGrayEntries);
 free(grayscale);
 free(grayscaleEntries);
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
        if (strncmp(argv[i], "-t", 2) == 0) { //strncmp checks the argument 1 to 3
            numThreads = atoi(argv[i] + 2); //add 2 to get the number after "-t"
        } else if (strncmp(argv[i], "-l", 2) == 0) {
            lines = atoi(argv[i] + 2);
        } else if (strncmp(argv[i], "", 0) == 0) {
            d = atoi(argv[i] + 0);
        } else {
            printf("Wrong format\n");
            return 1;
        }
    }

  printf("numThreads = %d \n", numThreads);
  printf("lines = %d \n", lines);
  printf("d = %d \n", d);

float *rootsEntries = (float*) malloc(sizeof(float)* d * 2);
 float ** roots = (float**) malloc(sizeof(float*)*d);

 for ( size_t ix = 0, jx = 0; ix < d; ++ix, jx+=2)
    roots[ix] = rootsEntries + jx;

 for ( size_t ix = 0; ix < 2; ++ix )
    for ( size_t jx = 0; jx < d; ++jx )
        roots[ix][jx] = 0;

GetRoots(roots, d);
printf("roots = %f + %fi \n", roots[0][0], roots[0][1]);

double stepSize = 1/lines;
const double sz = lines/numThreads; //Används ej!

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
    double ib = 2.0 - (tx * stepSize);
    double ie = -2.0 + (numThreads-tx)*stepSize;
    thrds_info[tx].attractors = attractors;
    thrds_info[tx].convergences = convergences;
    thrds_info[tx].roots=roots;
    thrds_info[tx].numThreads = numThreads; //Casta till const?
    thrds_info[tx].d = d;
    thrds_info[tx].lines = lines; //Size of vectors in v
    thrds_info[tx].tx = tx; //Thread index
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

 mtx_destroy(&mtx);
 cnd_destroy(&cnd);

return 0;

}

void GetRoots( float ** roots, int d) {
    for( size_t ix = 0; ix < d; ix++) {
            float theta = (ix*2* PI) / d ;
            roots[ix][0] = cos(theta);
            roots[ix][1] = sin(theta);
    }
}

void FunctionDerivate(double complex z,double complex derivate, int d){
    switch (d) {
    case 1:
        derivate = 1;
        break;
    case 2:
        derivate = 2*z;
        break;
    case 3:
        derivate = 3*z*z;
        break;
    case 4:
        derivate = 4*z*z*z;
        break;
    case 5:
        derivate = 5*z*z*z*z;
        break;
    case 6:
        derivate = 6*z*z*z*z*z;
        break;
    case 7:
        derivate = 7*z*z*z*z*z*z;
        break;
    case 8:
        derivate = 8*z*z*z*z*z*z*z;
        break;
    case 9:
        derivate = 9*z*z*z*z*z*z*z*z;
        break;
    case 10:
        derivate = 10*z*z*z*z*z*z*z*z*z;
        break;

    default:
        fprintf(stderr, "unexpected degree\n");
        exit(1);
    }
}


void GetColors(char**color, int d){
    color[0] = "255 0   0   ";
    color[1] = "0   255 0   ";
    color[2] = "0   0   255 ";
    color[3] = "255 255 0   ";
    color[4] = "255 0   255 ";
    color[5] = "0   255 255 ";
    color[6] = "0   0   0   ";
    color[7] = "255 255 255 ";
    color[8] = "100 0   0   ";
    color[9] = "0   100 0   ";
    color[10] = "0   0   100 ";
    color[11] = "100 100 0   ";
}

void GetGrayScale(char**grayscale, int d){
    int stepSize = 255;
    int greyValue;

    char greyElement[12];

    for( size_t ix = 0; ix < 128; ix++) {
            greyValue = stepSize - ix*2;
            snprintf(greyElement, sizeof(greyElement), "%3d %3d %3d", greyValue,greyValue,greyValue);
            grayscale[ix]= greyElement;
    }

}
