#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <ctype.h>
#include <math.h>


int main(int argc, char *argv[])
{
  const int dim = 3;
  const size_t bytesPerCell = 24;
  // TODO: this will have to be much more
  const int blockSize = 1000 * 21;
  


  int numThreads;
  if (argv[1][0] == '-' && argv[1][1] == 't' && strlen(argv[1]) > 2)
  {
    numThreads = atoi(&argv[1][2]);
  }
  else
  {
    printf("Wrong format for argument use format -t[INTEGER], where [INTEGER] must be valid whole number larger than 0.");
    exit(EXIT_FAILURE);
  }
  
   omp_set_num_threads(numThreads);

  if (numThreads == 0)
  {
    printf("INTEGER must be valid whole number larger than 0.");
    exit(EXIT_FAILURE);
  }


//TODO:räkna igenom så vi tillräckligt med minne
  unsigned long int distCount[3465*numThreads];
  for (size_t ix = 0; ix < 3465*numThreads; ++ix )
    distCount[ix] = 0;

  float * block0Entries =  (float*) malloc(sizeof(float) * blockSize * dim);
  float ** block0 = (float**) malloc(sizeof(float*) * blockSize);
  for (size_t ix = 0, jx = 0; ix < blockSize; ++ix, jx += dim)
  {
    block0[ix] = block0Entries + jx;
  }

  float * block1Entries =  (float*) malloc(sizeof(float) * blockSize * dim);
  float ** block1 = (float**) malloc(sizeof(float*) * blockSize);
  for ( size_t ix = 0, jx = 0; ix < blockSize; ++ix, jx += dim)
  {
    block1[ix] = block1Entries + jx;
  }


  // NOTE: change before submission
  FILE *file = fopen("test_data/cells_1e5", "r");
  if (file == NULL)
  {
    printf("Error opening file.\n");
    return -1;
  }

  fseek(file, 0, SEEK_END);
  const size_t fileSize = ftell(file);
  const size_t nmbCells = fileSize / bytesPerCell;

  
  for ( size_t ib = 0; ib < nmbCells; ib += blockSize ) {
    size_t ie = ib + blockSize < nmbCells ? ib + blockSize : nmbCells;

    fseek(file, ib*bytesPerCell, SEEK_SET);
    for (size_t ix = 0; ix < ie-ib; ++ix)
    {
      fscanf(file, "%f %f %f", block0[ix], block0[ix] + 1,  block0[ix] + 2);
    }  
    #pragma omp parallel for schedule(dynamic) shared(distCount)                                                          
    for ( size_t ix = 0; ix < ie-ib; ++ix ) {
      int threadNr = omp_get_thread_num();
      float x1 = block0[ix][0];
      float y1 = block0[ix][1];
      float z1 = block0[ix][2];

      for ( size_t jx = ix+1; jx < ie-ib; ++jx ) {
        float x2 = block0[jx][0];
        float y2 = block0[jx][1];
        float z2 = block0[jx][2];

        float xd = x1 - x2;
        float yd = y1 - y2;
        float zd = z1 - z2;
  //Warning check number of threads
        float dist = sqrtf(xd*xd + yd*yd + zd*zd);
        distCount[(int) (dist * 100.f)+3465*(threadNr)] += 1;
      }
    }
     
    for ( size_t jb = ie; jb < nmbCells; jb += blockSize ) {
      size_t je = jb + blockSize < nmbCells ? jb + blockSize : nmbCells;

      fseek(file, jb*bytesPerCell, SEEK_SET);
      for ( size_t jx = 0; jx < je-jb; ++jx)
      {
        fscanf(file, "%f %f %f", block1[jx], block1[jx] + 1,  block1[jx] + 2);
      }
      #pragma omp parallel for schedule(dynamic)  shared(distCount) 
      for ( size_t ix = 0; ix < ie-ib; ++ix ) {
	int threadNr = omp_get_thread_num();
        float x1 = block0[ix][0];
        float y1 = block0[ix][1];
        float z1 = block0[ix][2];
        for ( size_t jx = 0; jx < je-jb; ++jx ) {
          float x2 = block1[jx][0];
          float y2 = block1[jx][1];
          float z2 = block1[jx][2];

          float xd = x1 - x2;
          float yd = y1 - y2;
          float zd = z1 - z2;

          float dist = sqrtf(xd*xd + yd*yd + zd*zd);
          distCount[(int) (dist * 100.f)+3465*(threadNr)] += 1;        
	}
      }
    }
  }
   fclose(file);

   //TODO: Summera igen för att se att det blir rätt

   for (size_t ix = 0; ix < 3465; ++ix )
	   for(size_t jx = 1; jx < numThreads ; ++jx)
	   distCount[ix] += distCount[ix + jx*3465];

   long int sum = 0;
   for ( int ix = 0; ix < 3465; ++ix ){
    if (distCount[ix] != 0)
      printf("%.2f %d\n", ix/100.0, distCount[ix]);
      sum += distCount[ix];
   }
   printf("%ld", sum);

 
  free(block0Entries);
  free(block0);
  free(block1Entries);
  free(block1);

  return 0;

}
