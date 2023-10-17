#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
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

float calculate_new_x(float x, int d) {
    float new_x;

    switch (d) {
    case 1:
        new_x = 1.0;

        break;

    case 2:
        new_x = 0.5*x + 1/(2*x);
        break;

    case 3:
        new_x = 2*x/3 + 1/(3*x*x);
        break;

    case 4:
        new_x = 3*x/4 + 1/(4*x*x*x);
        break;

    case 5:
        new_x = 4*x/5 + 1/(5*x*x*x*x);
        break;

    case 6:
        new_x = 5*x/6 + 1/(6*x*x*x*x*x);
        break;

    case 7:
        new_x = 6*x/7 + 1/(7*x*x*x*x*x*x);

        break;
    
    case 8:
        new_x = 7*x/8 + 1/(8*x*x*x*x*x*x*x);
        break;

    case 9:
        new_x = 8*x/9 + 1/(9*x*x*x*x*x*x*x*x);
        break;
    
    case 10:
        new_x = 9*x/10 + 1/(10*x*x*x*x*x*x*x*x*x);
        break;
    
    default:
        fprintf(stderr, "unexpected degree\n");
        exit(1);
    }
    return new_x;
}



typedef struct {
    float real;
    float imaginary;
} ComplexNumber;



void GetRoots(ComplexNumber *roots, int d) {
    // The 'roots' array will store the complex roots
    // The 'colours' array will store the colors associated with each root

    if (d % 2 == 0) {
        roots[0].real = 1.0;
        roots[0].imaginary = 0.0;
        roots[1].real = -1.0;
        roots[1].imaginary = 0.0;

        if (d > 2) {
            for (size_t ix = 2; ix < d; ix++) {
                float theta = (ix * 2 * PI) / d;
                roots[ix].real = cos(theta);
                roots[ix].imaginary = sin(theta);
            }
        }
    } else {
        roots[0].real = 1.0;
        roots[0].imaginary = 0.0;

        if (d > 1) {
            for (size_t ix = 1; ix < d; ix++) {
                float theta = (ix * 2 * PI) / d;
                roots[ix].real = cos(theta);
                roots[ix].imaginary = sin(theta);                
            
            }
        }
    }
}


int main(int argc, char *argv[])
{
    // ./newton -t5 -l1000 7
    // ./newton -l1000 -t5 7
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

    omp_set_num_threads(numThreads);

    printf("numThreads = %d \n", numThreads);
    printf("lines = %d \n", lines);
    printf("d = %d \n", d);

    //Note: Defines the "size" of the pixels
    // istep = 0.1;

    // //WARN: Code copied from martin
    // //const thrd_info_t *thrd_info = (thrd_info_t*) args;
    // //const float **v = thrd_info->v;
    // //float **w = thrd_info->w;
    // //const int ib = thrd_info->ib;
    // const int istep = thrd_info->istep;
    // //const int sz = thrd_info->sz;
    // const int lines = thrd_info->lines;
    // //const int tx = thrd_info->tx;
    // mtx_t *mtx = thrd_info->mtx;
    // cnd_t *cnd = thrd_info->cnd;
    // int_padded *status = thrd_info->status;

    // TO DO: Vi behöver hitta de riktiga rötterna för att ta reda på vilken de konvergerar
    // till för bilden med färg! Till svartvit bild måste vi spara antalet iterationer

    //NOTE: Vi vet redan två rötter!! 1 ev. -1 om jämnt d och complexa konjugatet 
    //Hur hittar vi de andra komplexa rötterna?

    //För rötterna skapar matris storlek = d * 2 
    float *rootsEntries = (float*) malloc(sizeof(float) * d * 2);
    float **roots = (float**) malloc(sizeof(float) * d); //Rätt?

    for ( size_t ix = 0, jx = 0; ix < d; ++ix, jx+=2)
        roots[ix] = rootsEntries + jx;

    for ( size_t ix = 0; ix < d; ++ix )
        for ( size_t jx = 0; jx < 2; ++jx )
            roots[ix][jx] = 0;

for (size_t ix = 0; ix < d; ++ix) {
    // For each 'ix', make 'roots[ix]' point to the corresponding elements in 'rootsEntries'.
    roots[ix] = rootsEntries + ix * 2;
}

// Initialize the values for each element in 'roots'.
for (size_t ix = 0; ix < d; ++ix) {
    roots[ix][0] = 0.0; // Real part
    roots[ix][1] = 0.0; // Imaginary part
}
    
    //Tar ut rötterna för polynomet
    GetRoots(roots, d, colours);
    float ib;
    float istep = 0.01;
    //NOTE: Vill vi dela upp bilden i delar innan? En tråd gör en viss del?
    for( float ix = ib; ix < lines; ix += istep){
        //TODO: Vad är v?
        // const float *vix = v[ix];

        double *realValues = (double*) malloc(sizeof(double)*lines*1/istep*128);
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
        for (int kx = 0; kx < lines*1/istep; kx += istep){
                float x_re = realValues[kx];
                float x_im = imgValues;
                //for loop 128 steps if solution not found what do we do with pixel?
                for (int jx = 0; jx < 128; ix++){

                    float new_x_re = calculate_new_x(x_re, d);
                    float new_x_im = calculate_new_x(x_im, d);
    
                    new_x_re = new_x_re - x_re;
                    new_x_im = new_x_im - x_im;
                    printf("new_x_re %f \n", new_x_re);
                    float new_x = new_x_re + new_x_im;
                
                    //WARN: Detta gäller ej för komplext tal!!
                    // Double check termination critera
                    //Vill vi använda komplexa tal på polär form?
                    if(new_x < 0.001 && new_x > -0.001){
                        break;
                    }
                }
            }
            free(realValues);
        }
    free(roots);
    free(rootsEntries);

    





// create colour matrix
    int colour_index = 11;
    int colour_coordinates = 3;
    int colour_array[colour_index][colour_coordinates] = {
        {255, 0, 0}
        {0, 255, 0}
        {0, 0, 255}
        {255, 255, 0}
        {0, 255, 255}
        {255, 0, 255}
        {255, 255, 255}
        {0, 0, 0}
        {100, 0, 100}
        {100, 100, 0}
        {0, 100, 100}
    };

//create gray matrix
    int gray_index = 128;
    int gray_coordinates = 3;  
    int gray_array[gray_index][gray_coordinates];

    for (int i = 0; i < gray_index; i++) {
        int intensity = i * 2; 
        gray_array[i][0] = intensity;  
        gray_array[i][1] = intensity;  
        gray_array[i][2] = intensity; 

    FILE *colour = fopen("newton_attractors_xd.ppm", "w");
        if (colour == NULL){
            printf("Error opening file.\n");
            return -1;
        }

    for(int ix = 0; ix < lines; ix++){ // loopar över rader
        for ( size_t cx = 0; cx < lines; ++cx ) {
            int gray_index_number = attractor[cx]; //gråskala
            int colour_index_number = convergence[cx]; //färgbild


        
//talet i attractor = rad i gray_array

            //fwrite(input data, size in bytes, number of elements, file)
            fwrite(input_colours, lines*lines*sizeof(int), 1, colour);

        }
    }


    fclose(colour);




  return 0;

}

