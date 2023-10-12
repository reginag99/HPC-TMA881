#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <threads.h>

//WARN: Change before submission hur många decimaler vill vi ha?
#define PI 3.14159265358979323846 

 //WARN: Code copied from martin
typedef struct {
  int val;
  char pad[60]; // cacheline - sizeof(int)
} int_padded;

 //WARN: Code copied from martin
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

void GetRoots(
    float ** roots,
    int d //Rätt?
);

 //WARN: Code copied from martin
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
    int lines, numThreads, d;

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
            d = atoi(argv[i] + 0);
          
        } else {
            printf("Wrong format\n");
            return 1;
        }   
    }

    omp_set_num_threads(numThreads);

    printf("numThreads = %d \n", numThreads);
    printf("lines = %d \n", lines);
    printf("d = %d \n", d);

    //Note: Defines the "size" of the pixels
    istep = 0.1;

    //WARN: Code copied from martin
    //const thrd_info_t *thrd_info = (thrd_info_t*) args;
    //const float **v = thrd_info->v;
    //float **w = thrd_info->w;
    //const int ib = thrd_info->ib;
    const int istep = thrd_info->istep;
    //const int sz = thrd_info->sz;
    const int lines = thrd_info->lines;
    //const int tx = thrd_info->tx;
    mtx_t *mtx = thrd_info->mtx;
    cnd_t *cnd = thrd_info->cnd;
    int_padded *status = thrd_info->status;

    // TO DO: Vi behöver hitta de riktiga rötterna för att ta reda på vilken de konvergerar
    // till för bilden med färg! Till svartvit bild måste vi spara antalet iterationer

    //NOTE: Vi vet redan två rötter!! 1 ev. -1 om jämnt d och complexa konjugatet 
    //Hur hittar vi de andra komplexa rötterna?

    //För rötterna skapar matrs storlek = d * 2 
    float *rootsEntries = (float*) malloc(sizeof(float)* d * 2);
    float ** roots = (float**) malloc(sizeof(float)* d); //Rätt?

    for ( size_t ix = 0, jx = 0; ix < d; ++ix, jx+=2)
        roots[ix] = rootsEntries + jx;

    for ( size_t ix = 0; ix < d; ++ix )
        for ( size_t jx = 0; jx < 2; ++jx )
            roots[ix][jx] = 0;

    
    //Tar ut rötterna för polynomet
    GetRoots( float ** roots, int d);

    //NOTE: Vill vi dela upp bilden i delar innan? En tråd gör en viss del?
    for( float ix = ib; ix < lines; ix += istep){
        //TODO: Vad är v?
        const float *vix = v[ix];
        float *realValues = (float*) malloc(lines*sizeof(float)*1/istep);
        //WARN: Vilken typ ska vi använda här är float för litet?
        float imgValues = ix - lines/2;
        //TODO: Allocate every row with rows in the image så varje tråd får en rad av bilden
        // x axeln sammma real del samma, imaginär del konstant värde i raden 
        //NOTE: Bilden ska ha storleken lines, en pixel per integer vill vi ha fler?
        for ( int jx = 0; jx < lines; jx += istep )
            realValues[jx] = ix - lines/2;

        //TO DO: Plocka ut ett x-pixel värde för netwtons metod
        //Hitta vilken rot man är närmast
        //Spara antalet iterationer
        for (kx = 0; kx < lines*1/istep; kx += istep){
                float x_re = realValues[kx];
                float x_im = imgValues;
                //for loop 128 steps if solution not found what do we do with pixel?
                for (jx = 0; jx < 128; ix++){
                    float dividor;
                    float dx_re, dx_im;

                    //Beräkna x^d och xd i en funktion?
                    // Vill vi använda en switch istället? som i beskrivningen av assignment 3?
                    //Hur funkar det i trådar?
                    for (int i = 0; i < (d-1); i++){
                        x_re = x_re * realValues[jx] + x_im * imgValues; 
                        x_im = x_im * realValues[jx] + x_re * imgValues; 
                    }
                    //WARN: create if-statement for case when d = 1 or 2
                    //Should we use sin/cos for of imaginary numbers instead?
                    for (int i = 0; i < (d-2); i++){
                        dx_re = dx_re * realValues[jx] + dx_im * imgValues; 
                        dx_im = dx_im * realValues[jx] + dx_re * imgValues; 
                    }

                //WARN: Dela upp i imaginär och real del!
                // För att dividera måste vi ha komplexa konjugatet! 
                // Check for division by zero!

                dx_im_conjugate = -dx_im; 
                //WARN: Double check these calculations!! and define everything somewhere else
                dividor = dx_re * dx_re + (dx_im) * (dx_im);
                float re_x =  ((x_re - 1)* dx_re + dx_im_conjugate * dx_im)/ dividor; // Real del divison
                float img_x = ((x_re - 1)* dx_im_conjugate + x_im * dx_re )/ dividor; // Imaginär del division

                //WARN: Byt notation detta blev otydligt!
                //Beräkning av värdena för nästa iteration
                x_re = x_re - re_x;
                x_im = x_im - img_x;
                
                //WARN: Detta gäller ej för komplext tal!!
                // Double check termination critera
                //Vill vi använda komplexa tal på polär form?
                if(x0 < 0.001 && x0 > -0.001){
                    break;
                }
            //if (jx ==)
            }
        }


    }
    
  return 0;

}

//WARN: Är detta rätt? Vill man inlina funktionen?
void GetRoots( float ** roots, int d) {

     // Hårdkoda de riktiga rötterna alltså 1 och -1? 
     // Beror på om d är jämnt eller ej
    if (d % 2 == 0) {
        roots[0][0] = 1.;
        roots[0][1] = 0.;
        roots[1][0] = -1.;
        roots[1][1] = 0.;

    } else {
        roots[0][0] = 1.; 
        roots[0][1] = 0.;
    }

    //TODO: Hitta de andra rötterna! De Moivre's Theorem?
    if (d > 2){
        for( size_t ix = 2; ix < d; ix++) {
            float theta = (ix*2* PI) / d ;        
            roots[ix][0] = cos(theta);
            roots[ix][1] = sin(theta);
        }

    }

}
