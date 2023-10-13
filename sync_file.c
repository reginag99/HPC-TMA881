// Synchronization

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <threads.h>

//WARN: Code copied from martin
typedef struct {
  int val;
  char pad[60]; // cacheline - sizeof(int)
} int_padded; //

//WARN: Code copied from martin
//Defines a thread structure
//Outpuds needs to be somewhat continuous to be able to send to 
// the other threads in the other thread structure

//We have a continous requirement for the work process

//TODO: Create struct for threads for creating "coordinates" and doing computation
// Lämna varje pixel till reading thread
//Definera allt vi vill skicka in i funktionen för trådarna!!
typedef struct {
  const float **vReal; //Double pointer "collection of lines" ska denna vara const??
  float **vImg; //För imaginär delen endast en float EJ array
  float **wReal; // Vill vi använda denna för att utföra operationer? 
  float **wImg; // Vill vi använda denna för att utföra operationer?
  int ib; //increment ib with istep until size sz is reached, ska variera med trådar
  int stepSize;
  int lines; //Stoppa iterationen när vi har itererat igenom alla y värden!!
  int tx; 
  mtx_t *mtx; //Mutex, used to syncronize work need to be destoryed at the end of program
  cnd_t *cnd; // conditional variable, used to syncronize work need to be destoryed at the end of program
  int_padded *status;// padded used for treating false sharing, used to communicate status/progress of threads in processing its lines
} thrd_info_t;


//TODO: Create struct for threads for writing to file needs to work from top to the bottom?
typedef struct {
  const float **v; //Double pointer "collection of lines" "pixel"??
  float **w;
  int ib; //increment ib with istep until size sz is reached
  int stepSize;
  int sz; //Vill vi ha en annan storlek?
  int tx; 
  mtx_t *mtx; //Mutex, used to syncronize work need to be destoryed at the end of program
  cnd_t *cnd; // conditional variable, used to syncronize work need to be destoryed at the end of program
  int_padded *status;// padded used for treating false sharing, used to communicate status/progress of threads in processing its lines
} thrd_info_t;

//This should do the actual computation!!
int
main_thrd(
    void *args
    )
{ 
  //Copying things to local variables
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  const float **vReal = thrd_info->vReal; //Kan hårdkodas
  const float **vImg = thrd_info->vImg; //Kan hårdkodas
  float **wRel = thrd_info->wRel;
  float **wImg = thrd_info->wImg;
  const int ib = thrd_info->ib;//Vad är ib skickar vi in det?
  const int stepSize = thrd_info->stepSize;
  const int lines = thrd_info->lines;
  const int tx = thrd_info->tx;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;

  //const float *vix_Img = vImg[]; //Hämta ut en konstant float

  for ( int ix = ib; ix < sz; ix += stepSize ) {
    const float *vix_real = vReal[ix]; 
    // We allocate the rows of the result before computing, and free them in another thread.
    float *wix_real = (float*) malloc(sz*sizeof(float)); //Beräknar vi en rad i taget eller 1 element?
    float *wix_img = (float*) malloc(sz*sizeof(float)); //

    for ( int jx = 0; jx < sz; ++jx )
      wix_real[jx] = vix_real[jx];
      //wix_img[jx] = vix[jx];

    //Check thread to see if further progress has been made?
    mtx_lock(mtx);
    wix_real[jx] = wix_real;
    status[tx].val = ix + stepSize;
    mtx_unlock(mtx);
    cnd_signal(cnd);

    // In order to illustrate thrd_sleep and to force more synchronization
    // points, we sleep after each line for one micro seconds.
    thrd_sleep(&(struct timespec){.tv_sec=0, .tv_nsec=1000}, NULL);
  }

  return 0;
}


int
main()
{ 
  //Talen kommer variera mellan -2 och 2!!!
  //WARN: Remove declaration of lines when merging codes!!
  const int lines = 1000; 
  float stepSize = 1/lines; //WARN: Är detta rätt??
  int numThreads = 4; 
  //Är detta storleken på ett block eller en rad eller hela bilden?
  const int sz = lines/numThreads;  //Storleken i y-led! storlek i x-led är konstant!

  //Plocka ut ett tal från en array behöver vi spara den eller inte 

  //Deklarera matrix ska detta vara storleken på ett block en rad eller hela bilden? Just nu hela bilden
  //WARN: Måste vi deklarera något minne utanför trådarna??
  //float **vReal = (float**) malloc(sz*sizeof(float*)); 

  //Allokera reella talen mellan -2 och 2 anpassat till storleken på bilden!
  //Detta kan eventuellt hårdkodas till main thread function
  float *vReal = (float*) malloc(lines*sizeof(float));

  for ( int ix = -2; ix <= 2; ix += stepSize )
    vReal[ix] = ix;

  //Allokera imaginära talen mellan -2 och 2 anpassat till storleken på bilden!
  //Detta kan eventuellt hårdkodas till main thread function
  float *vImg = (float*) malloc(lines*sizeof(float));

  for ( int ix = -2; ix <= 2; ix += stepSize )
    vImg[ix] = ix;

  //Allt i main thread funktionen kommer göras av alla trådar! 

  //Fel dessa vill vi ska vara det vi sparar under newton iterationerna
  // ska innehålla 1 real del + imaginär del + antalet iterationer och vilken rot vi konvergerar till??
  float **colors = (float**) malloc(lines*sizeof(float*));
  float **nIterates = (float**) malloc(lines*sizeof(float*));
  
  
  //Vill vi spara koordinaten med??
  //float **wReal = (float**) malloc(lines*sizeof(float*));
  //float **wImg = (float**) malloc(lines*sizeof(float*));
  // The entries of w will be allocated in the computation threads are freed in
  // the check thread.

  //Setting up thread info
  thrd_t thrds[numThreads];
  thrd_info_t thrds_info[numThreads];

  //Creating object for thread checks
  thrd_t thrd_check;
  thrd_info_check_t thrd_info_check;
  
  //Initialize mutex
  mtx_t mtx;
  mtx_init(&mtx, mtx_plain);

  //Initialize contitional variables
  cnd_t cnd;
  cnd_init(&cnd);

  int_padded status[numThreads]; //Each thread needs its own status variable

//Martin gör for loop för att skapa numThreads antal trådar och kalla på main_thrd function iterera mellan antalet trådar, 
//skicka in vReal, vImg, trueRoots  vad mer??
   for ( int tx = 0; tx < numThreads; ++tx ) {
    //Initializes local variables for threads

    //Kalla på main_thrd funktionen!
    int r = thrd_create(thrds+tx, main_thrd, (void*) (thrds_info+tx));//thrds_info+tx <- start punkten, vad är ib?
    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
   }

   return 0;

}

