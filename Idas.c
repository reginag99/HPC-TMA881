
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <ctype.h>
#include <math.h>   




/*


Har lagt in ett förslag på hur man kanske kan beräkna distanserna
Får dock inte till allokeringarna när jag försöker lägga ihop alla distanser i distanceMatrix
Har därför inte kunnat kolla om beräkningarna funkar...


*/



int compare (const void * a, const void * b);

int main(int argc, char *argv[])
{
int numberRow = 17;
int numberCol = 3, numThreads;
int argLenght = strlen(argv[1]);

double distance, x1, x2, x3, y1, y2, y3;
int block_row = 1000;
int block_characters = 24;
int block_size = block_row * block_characters;


int size_distances_in_block = block_row*(block_row-1)/2;
int sizeDistanceMatrix = numberRow*(numberRow-1)/2;
int size_distances_between_blocks = sizeDistanceMatrix -size_distances_in_block;


if (argv[1][0] == '-' && argv[1][1] == 't' && argLenght > 2) {
    numThreads = atoi(&argv[1][2]);
   } else  {
    printf("Wrong format for argument use format -t[INTEGER], where [INTEGER] must be valid whole number larger than 0.");
    exit(EXIT_FAILURE);
   }

if (numThreads == 0){
    printf("INTEGER must be valid whole number larger than 0.");
    exit(EXIT_FAILURE);
   }

printf("Number of threads %d\n", numThreads);



float * blockEntries = (float*) malloc(sizeof(float) * block_size);
float **block = (float**) malloc(sizeof(float*) * block_row);

double * distances_in_block =  (double*) malloc(sizeof(double) * size_distances_in_block);
double * distances_between_blocks = (double*)malloc(sizeof(double) * size_distances_between_blocks);
for ( size_t ix = 0, jx = 0; ix < block_row; ++ix, jx+=block_characters)
   block[ix] = blockEntries + jx;




FILE *file = fopen("test_data/cells_17", "r");
if (file == NULL){
    printf("Error opening file.\n");
    return -1;
   }



for (int i = 0; i < numberRow; i += block_row) {
    for (int ix = 0; ix < block_row; ix++)
        fread(block[ix], sizeof(float), block_size, file);

    for (int w = i + 1000; w < numberRow; w += 1000) {
        for (int ix = 0; ix < block_row; ix++) {
            int element = 0;
            for (int j = 0; j < (block_row - 1); ++j) {
               x1 = block[j][0];
               x2 = block[j][1];
               x3 = block[j][2];
            for (int k = j + 1; k < block_row; ++k){
               y1 = block[k][0];
               y2 = block[k][1];
               y3 = block[k][2];

               distance = sqrt((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2)+(x3-y3)*(x3-y3));
               distances_in_block[element] = (int)(distance * 100 + 0.5) / 100.0;
               element += 1;
               
               }
            }
         }
        for (int ix1 = i; ix1 < i + block_row; ix1++) 
        {   int element_between = 0;
            for (int ix2 = w; ix2 < w + block_row; ix2++) 
            {
                x1 = block[ix1 - i][0];
                x2 = block[ix1 - i][1];
                x3 = block[ix1 - i][2];

                y1 = block[ix2 - w][0];
                y2 = block[ix2 - w][1];
                y3 = block[ix2 - w][2];

                distance = sqrt((x1 - y1) * (x1 - y1) + (x2 - y2) * (x2 - y2) + (x3 - y3) * (x3 - y3));
                distances_between_blocks[element_between] = (int)(distance * 100 + 0.5) / 100.0;
                element_between += 1;
                
            }
        }
    }
}



double *distanceMatrix = (double*) malloc(sizeof(double) * sizeDistanceMatrix);

for (int i = 0; i < size_distances_in_block; i++) {
    distanceMatrix[i] = distances_in_block[i];
}

for (int i = 0; i < size_distances_between_blocks; i++) {
    distanceMatrix[size_distances_in_block + i] = distances_between_blocks[i];
}



for (int ix = 0 ; ix < sizeDistanceMatrix; ++ix)                                                                                                                                       
    printf("%f\n",distanceMatrix[ix]);     

fclose(file);

free(block);
free(blockEntries);
free(distances_in_block);
free(distances_between_blocks);






for( int ix = 0; ix < sizeDistanceMatrix; ++ix)
    {
        for(int jx = ix + 1; jx < sizeDistanceMatrix; ++jx)
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