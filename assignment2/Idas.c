
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <ctype.h>
#include <math.h>   


/*


Har lagt in ett förslag på hur man kanske kan beräkna distanserna
Får dock segfault så kan inte kolla om den funkar...

*/



//int compare (const void * a, const void * b);

int main(int argc, char *argv[])
{

int numberRow = 10;
int numberCol = 3; 
//int numThreads;
//int argLenght = strlen(argv[1]);

double distance;
int block_row = 5;
int block_characters = 24;
int block_size = block_row * numberCol;
int number_blocks = numberRow/block_row;


int size_distances_in_block = block_row*(block_row-1)/2;
int sizeDistanceMatrix = numberRow*(numberRow-1)/2;
int size_distances_between_blocks = sizeDistanceMatrix - size_distances_in_block*number_blocks;

printf("size_distance_in_block %d\n", size_distances_in_block);
printf("number of blocks %d \n", number_blocks);
printf("size_distance_between_blocks %d\n", size_distances_between_blocks);
printf("sizeDistanceMatrix %d\n", sizeDistanceMatrix);


// if (argv[1][0] == '-' && argv[1][1] == 't' && argLenght > 2) {
//     numThreads = atoi(&argv[1][2]);
//    } else  {
//     printf("Wrong format for argument use format -t[INTEGER], where [INTEGER] must be valid whole number larger than 0.");
//     exit(EXIT_FAILURE);
//    }

// if (numThreads == 0){
//     printf("INTEGER must be valid whole number larger than 0.");
//     exit(EXIT_FAILURE);
//    }

// printf("Number of threads %d\n", numThreads);



float * in_blockEntries = (float*) malloc(sizeof(float) * block_size);
float **in_block = (float**) malloc(sizeof(float*) * block_row);
float * between_blockEntries = (float*) malloc(sizeof(float) * block_size);
float **between_block = (float**) malloc(sizeof(float*) * block_row);

double * distances_in_block =  (double*) malloc(sizeof(double) * size_distances_in_block);
double * distances_between_blocks = (double*)malloc(sizeof(double) * size_distances_between_blocks);
for ( size_t ix = 0, jx = 0; ix < block_row; ++ix, jx+=block_characters)
   in_block[ix] = in_blockEntries + jx;
for ( size_t ix = 0, jx = 0; ix < block_row; ++ix, jx+=block_characters)
   between_block[ix] = between_blockEntries + jx;



FILE *file = fopen("cells_10", "rb");
if (file == NULL){
    printf("Error opening file.\n");
    return -1;
   }


size_t block_bytes = block_size * sizeof(float);


// Calculate distances within each block
for (int i = 0; i < numberRow; i += block_row) {

        for (int ix = 0; ix < block_row; ix++) {
            fread(in_block[ix], sizeof(float), block_bytes, file);
            int element = 0;

            for (int j = 0; j < (block_row - 1); j++) {
               float x1 = in_block[j][0];
               float x2 = in_block[j][1];
               float x3 = in_block[j][2];

            for (int k = j + 1; k < block_row; k++){
               float y1 = in_block[k][0];
               float y2 = in_block[k][1];
               float y3 = in_block[k][2];

               distance = sqrt((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2)+(x3-y3)*(x3-y3));
               distances_in_block[element] = (int)(distance * 100 + 0.5) / 100.0;
               element += 1;
                printf("distances_in_block[%d] = %.2f\n", element - 1, distances_in_block[element - 1]);
               
            
            }
         }
    }
}



// Calculate distances between blocks
for (int i = 0; i < numberRow; i += block_row) {

    for (int ix = 0; ix < block_row; ix++)
        fread(between_block[ix], sizeof(float), block_bytes, file);

        for (int ix1 = i; ix1 < i + block_row; ix1++) 
        {   
            fread(between_block[ix1], sizeof(float), block_bytes, file);
            int element_between = 0;
            float x1 = between_block[ix1][0];
            float x2 = between_block[ix1][1];
            float x3 = between_block[ix1][2];

            for (int ix2 = i+1; ix2 < i + block_row; ix2++) {

                fread(between_block[ix2 + block_row], sizeof(float), block_bytes, file);
                float y1 = between_block[ix2 + block_row][0];
                float y2 = between_block[ix2 + block_row][1];
                float y3 = between_block[ix2 + block_row][2];

                distance = sqrt((x1 - y1) * (x1 - y1) + (x2 - y2) * (x2 - y2) + (x3 - y3) * (x3 - y3));
                distances_between_blocks[element_between] = (int)(distance * 100 + 0.5) / 100.0;
                element_between += 1;
                
                printf("distances_between_blocks[%d] = %.2f\n", element_between - 1, distances_between_blocks[element_between - 1]);
                
            }
        }
    }




double *distanceMatrix = (double*) malloc(sizeof(double) * sizeDistanceMatrix);

// Copy elements from distances_in_block to distanceMatrix
for (int i = 0; i < size_distances_in_block; i++) {
    distanceMatrix[i] = distances_in_block[i];
}

// Copy elements from distances_between_blocks to distanceMatrix
for (int i = 0; i < size_distances_between_blocks; i++) {
    distanceMatrix[size_distances_in_block + i] = distances_between_blocks[i];
}



// for (int ix = 0 ; ix < sizeDistanceMatrix; ++ix)                                                                                                                                       
//     printf("%f\n",distanceMatrix[ix]);     

fclose(file);

for (int i = 0; i < block_row; i++) {
    free(in_block[i]);
    free(between_block[i]);
}
free(in_block);
free(between_block);


free(in_blockEntries);
free(between_blockEntries);
free(distances_in_block);
free(distances_between_blocks);






for( int ix = 0; ix < sizeDistanceMatrix; ++ix)
    {
        for(int jx = ix + 1; jx < sizeDistanceMatrix; jx++)
        { if(distanceMatrix[ix]> distanceMatrix[jx])
            {
                double temp = distanceMatrix[ix];
                distanceMatrix[ix] = distanceMatrix[jx];
                distanceMatrix[jx] = temp;
            }
        }
    }

// for (int ix = 0 ; ix < sizeDistanceMatrix; ++ix)                                                                                                                                       
//     printf("%f\n",distanceMatrix[ix]);                                                                                                                                                 

float largestElement = distanceMatrix[sizeDistanceMatrix-1];

//printf("%f", largestElement);                                                                                                                                                          

int sizeFrequencyMatrix = largestElement * 100; //Antag att vi har två decimaler                                                                                                         

float *frequencyMatrixEntries = (float*) malloc(sizeof(float) * sizeFrequencyMatrix*2);
float **frequencyMatrix = (float**) malloc(sizeof(float*) * 2);

for ( size_t ix = 0, jx = 0; ix < 2; ++ix, jx+= sizeFrequencyMatrix) //Vilken ordning är bäst??                                                                                          
   frequencyMatrix[ix] = frequencyMatrixEntries + jx;


//Initiera 0or i frequencymatrix                                                                                                                                                         
for (size_t  ix = 0; ix < 2; ++ix)
    for (size_t jx = 0; jx <= sizeFrequencyMatrix; ++jx)
       frequencyMatrix[ix][jx] = 0.;


for (size_t ix = 0; ix < sizeDistanceMatrix; ++ix)//Vilken ordning är bäst??                                                                                                             
   {
    int jx = (int)(distanceMatrix[ix]*100);
    frequencyMatrix[0][jx] = distanceMatrix[ix];
    frequencyMatrix[1][jx] += 1;
   }



free(distanceMatrix);

printf("Before printing distanceMatrix\n");
for (int ix = 0; ix < sizeFrequencyMatrix; ++ix)
{
   if (frequencyMatrix[0][ix] != 0.00)
   {
    printf("%.2f %.2f \n",frequencyMatrix[0][ix],frequencyMatrix[1][ix]);
   }
}
printf("After printing distanceMatrix\n");

free(frequencyMatrix);
free(frequencyMatrixEntries);

 return 0;

}
