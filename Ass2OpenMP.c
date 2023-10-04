#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <ctype.h>
#include <math.h>

//Plocka ut siffran i -t5 och                                                                                                                              
//Måste lösningen kunna ta hänsyn till mellanslag?                                                                                                         
//Om ja leta efter första siffran sedan atoi?                                                                                                              

int compare (const void * a, const void * b);

int main(int argc, char *argv[])
{
int numberRow = 10;
int numberCol = 3, numThreads;
int argLenght = strlen(argv[1]);

if (argv[1][0] == '-' && argv[1][1] == 't' && argLenght > 2)
   {
    numThreads = atoi(&argv[1][2]);
   }
else
   {
    printf("Wrong format for argument use format -t[INTEGER], where [INTEGER] must be valid whole number larger than 0.");
    exit(EXIT_FAILURE);
   }

if (numThreads == 0)
   {
    printf("INTEGER must be valid whole number larger than 0.");
    exit(EXIT_FAILURE);
   }

printf("Number of threads %d\n", numThreads);

float * matrixEntries = (float*) malloc(sizeof(float) * numberRow*numberCol);
float ** matrix = (float**) malloc(sizeof(float*) * numberRow);

for ( size_t ix = 0, jx = 0; ix < numberRow; ++ix, jx+=numberCol)
   matrix[ix] = matrixEntries + jx;

FILE *file = fopen("test_data/cells_10.txt", "r");

if (file == NULL)
   {
    printf("Error opening file.\n");
    return -1;
   }

for (int ix = 0; ix < numberRow; ix++)
     for (int jx = 0; jx < numberCol; jx++)
     {
      fscanf(file, "%f", &matrix[ix][jx]);
     }

fclose(file);

//Räkna distanserna mellan koordinaterna                                                                                                                   

int sizeDistanceMatrix = (numberRow*(numberRow-1)/2);
double distance, x1, x2, x3, y1, y2, y3, temp;
int element = 0;

double * distanceMatrix =  (double*) malloc(sizeof(double) * sizeDistanceMatrix);

for (int ix = 0 ; ix < (numberRow - 1) ; ++ix)
   {
        x1 = matrix[ix][0];
        x2 = matrix[ix][1];
        x3 = matrix[ix][2];
    for (int jx = ix + 1; jx < numberRow; ++jx )
        {
        y1 = matrix[jx][0];
        y2 = matrix[jx][1];
        y3 = matrix[jx][2];

        distance = sqrt((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2)+(x3-y3)*(x3-y3));
        distanceMatrix[element] = round(distance*100)/100;
        element += 1;
        }
   }

free(matrix);
free(matrixEntries);

for( int ix = 0; ix < sizeDistanceMatrix; ++ix)
    {
        for(int jx = ix + 1; jx < sizeDistanceMatrix; ++jx)
        { if(distanceMatrix[ix]> distanceMatrix[jx])
            {
                temp = distanceMatrix[ix];
                distanceMatrix[ix] = distanceMatrix[jx];
                distanceMatrix[jx] = temp;
            }
        }
    }

//for (int ix = 0 ; ix < sizeDistanceMatrix; ++ix)                                                                                                         
//    printf("%f\n",distanceMatrix[ix]);                                                                                                                   

int largestElement = distanceMatrix[sizeDistanceMatrix-1];
int sizeFrequencyMatrix = largestElement * 100; //Antag att vi har två decimaler                                                                           

float * frequencyMatrixEntries = (float*) malloc(sizeof(float) * sizeFrequencyMatrix*2);
float ** frequencyMatrix = (float**) malloc(sizeof(float*) * sizeFrequencyMatrix);

for ( size_t ix = 0, jx = 0; ix <= sizeFrequencyMatrix; ++ix, jx+= 2) //Vilken ordning är bäst??                                                           
   frequencyMatrix[ix] = frequencyMatrixEntries + jx;


//Initiera 0or i frequencymatrix                                                                                                                           
for (size_t  ix = 0; ix <= sizeFrequencyMatrix; ++ix)
    for (size_t jx = 0; jx < 2; ++jx)
       frequencyMatrix[ix][jx] = 0;


for (size_t ix = 0; ix < sizeDistanceMatrix; ++ix)//Vilken ordning är bäst??                                                                               
   {
    int jx = (int)(distanceMatrix[ix]*100);
    frequencyMatrix[0][jx] = distanceMatrix[ix];
    frequencyMatrix[1][jx] =+ 1;
   }


printf("Distance: %.2f Frequency: %.2f",frequencyMatrix[0][300],frequencyMatrix[1][300]);


free(distanceMatrix);
free(frequencyMatrix);
free(frequencyMatrixEntries);

 return 0;

}
