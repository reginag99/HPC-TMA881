#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <threads.h>

typedef struct {
  int val;
  char pad[60]; // cacheline - sizeof(int)
} int_padded;

typedef struct {
  const float **v;
  float **w;
  int ib;
  int istep;
  int sz;
  int tx;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_t;

typedef struct {
  const float **v;
  float **w;
  int sz;
  int nthrds;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_check_t;

int main(int argc, char *argv[])
{
// ./newton -t5 -l1000 7
// ./newton -l1000 -t5 7
    int lines, numThreads, exponent_d;

    if (argc != 3 && argc != 4) {
        printf("Two or three arguments expected");
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i < 3; i++) { 
        if (strncmp(argv[i], "-t", 2) == 0) { 
            numThreads = atoi(argv[i] + 2); //add 2 to get the number after "-t"
        } else if (strncmp(argv[i], "-l", 2) == 0) {
            lines = atoi(argv[i] + 2); 
        } else if (strncmp(argv[i], "", 0) == 0) {
            exponent_d = atoi(argv[i] + 0);
          
        } else {
            printf("Wrong format\n");
            return 1;
        }   
    }

    omp_set_num_threads(numThreads);

    printf("numThreads = %d \n", numThreads);
    printf("lines = %d \n", lines);
    printf("d = %d \n", exponent_d);

    //Note: Defines the "size" of the pixels
    istep = 0.1;

    //WARN: Code copied from martin
    const thrd_info_t *thrd_info = (thrd_info_t*) args;
    const float **v = thrd_info->v;
    float **w = thrd_info->w;
    const int ib = thrd_info->ib;
    const int istep = thrd_info->istep;
    //const int sz = thrd_info->sz;
    const int lines = thrd_info->lines;
    const int tx = thrd_info->tx;
    mtx_t *mtx = thrd_info->mtx;
    cnd_t *cnd = thrd_info->cnd;
    int_padded *status = thrd_info->status;

    //Är detta möjligt?
    for( float ix = ib; ix < lines; ix += istep){
        const float *vix = v[ix];
        float *realValues = (float*) malloc(lines*sizeof(float)*1/istep);
        //WARN: Vilken typ ska vi använda här?
        float imgValues = ix - lines/2;
        //TODO: Allocate every row with rows in the image så varje tråd får en rad av bilden
        // x axeln sammma real del samma, imaginär del konstant värde i raden 
        //NOTE: Bilden ska ha storleken lines, en pixel per integer vill vi ha fler?
        for ( int jx = 0; jx < lines; jx += istep )
            realValues[jx] = ix - lines/2;

        //TO DO: Plocka ut ett x-pixel värde gör netwtons metod
        for (kx = 0; kx < lines*1/istep; kx += istep){
                
                //for loop 128 steps if solution not found what do we do with pixel?
                for (jx = 0; jx < 128; ix++){
                    for (int i = 0; i < exponent_d; i++){
                        x *= x[i]; //WARN: change for imaginary
                    }
                    //Suggestion: create function and hardcode computation for different values of d
                    for (int i = 1; i < exponent_d; i++){
                        dx *= x[i];
                    }

                x1 = x0 - (x-1)/(dx*exponent_d);
                x0 = x1;
                //WARN: Double check termination critera
                if(x0 < 0.001 && x0 > -0.001){
                    break;
                }
            //if (jx ==)
            }
    }


    }
    

   




  return 0;

}
