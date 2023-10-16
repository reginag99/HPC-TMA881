#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <threads.h>


typedef struct {
  int val;
  char pad[60];
} int_padded;


typedef struct {
<<<<<<< HEAD
  const float **vReal; //Double pointer "collection of lines" ska denna vara const?? nej tänker jag spontant, men imaginära kan väl det?
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
=======
  double *wReal;
  double *wImg;
  double ib;//Ska denna vara const?                                                                                                                                   
  float stepSize;
  int lines;
  int tx;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
>>>>>>> 3ea401b (Update sync_file.c)
} thrd_info_t;


typedef struct {
  double *wReal;
  double *wImg;
  int lines;
  int numThreads;
  //const float sz;                                                                                                                                                   
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_check_t;

int main_thrd(void *args)
{
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
<<<<<<< HEAD
  const float **vReal = thrd_info->vReal; //Kan hårdkodas. VAD MENAR DU
  const float **vImg = thrd_info->vImg; //Kan hårdkodas
  float **wRel = thrd_info->wRel;
  float **wImg = thrd_info->wImg;
  const int ib = thrd_info->ib;//Vad är ib skickar vi in det?
  const int stepSize = thrd_info->stepSize;
=======
  double *wReal = thrd_info->wReal; //Ska bytas till färg, koordinater för testning                                                                                   
  double *wImg = thrd_info->wImg; //Ska bytas till antalet, koordinater för testning                                                                                  
  double ib = thrd_info->ib;
  const double stepSize = thrd_info->stepSize;
>>>>>>> 3ea401b (Update sync_file.c)
  const int lines = thrd_info->lines;
  const int tx = thrd_info->tx;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;

  int i;
  double jx;

  //ib ska gå från 2 till -2                                                                                                                                          
  for ( double ix = ib; ix >= -2.; ix -= (stepSize*4.0)) { //Skickar in vart vi börjar i ib!                                                                          
    double *wix_real = (double*) malloc(lines*sizeof(double));
    double *wix_img = (double*) malloc(lines*sizeof(double));

    for ( i = 0, jx = -2.; i < lines, jx <= 2.; ++i , jx += stepSize){
      wix_real[i] = jx;
      wix_img[i] = ix;
      }

    //TODO: Coordinate not printing as expected and cannot see that threads loop through everything                                                                   
    //Also how is this information "sent forward??"                                                                                                                   
    printf("Coordinate: [%.4lf, %.4lf]",wix_real[0],wix_img[0]);
    fflush(stdout);

    mtx_lock(mtx);
    status[tx].val = ix - (stepSize*4);
    mtx_unlock(mtx);
    cnd_signal(cnd);

    thrd_sleep(&(struct timespec){.tv_sec=0, .tv_nsec=1000}, NULL);
  }
  return 0;
}

int
main_thrd_check(
    void *args
    )
{
  const thrd_info_check_t *thrd_info = (thrd_info_check_t*) args;
  double *wReal = thrd_info->wReal; //Ska bytas till färg eller antal iterationer                                                                                     
  double *wImg = thrd_info->wImg; //Ska bytas till färg eller antal iterationer!!                                                                                     
  const int lines = thrd_info->lines;
  const int numThreads = thrd_info->numThreads;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;

  for ( int ix = 0, ibnd; ix < lines; ) {
    for ( mtx_lock(mtx); ; ) {
      ibnd = lines;

      for ( int tx = 0; tx < numThreads; ++tx ){
        if ( ibnd > status[tx].val )
          ibnd = status[tx].val;
      if ( ibnd <= ix )
        cnd_wait(cnd,mtx);
      else {
        mtx_unlock(mtx);
        break;
      }
     }
    }

    fprintf(stderr, "checking until %i\n", ibnd);
  
    for ( ; ix < ibnd; ++ix ) {
      int jx = lines-1;
      printf("[%.4lf, %.4lf]", wReal[ix],wImg[ix]);

      free(wReal);
      free(wImg);
    }
  }

  return 0;
}

int
main()
{
  const float lines = 1000.;
  double stepSize = 1/lines;
  int numThreads = 4;
  const double sz = lines/numThreads; //Används ej!                                                                                                                   
  double *wReal = (double*) malloc(lines*sizeof(double*));
  double *wImg = (double*) malloc(lines*sizeof(double*));
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

  int_padded status[numThreads];

//Vi skickar iväg trådarna 1 gång!! Varje tråd 1 gång med 1 ib!!                                                                                                      
// tx 0 ib = 2, tx 1 ib = 2-1/stepSize ...                                                                                                                            
//Hur fungerar detta med behovet av en kontinuerlig minneshantering?     

 for ( int tx = 0; tx < numThreads; ++tx ) {
    double ib = 2.0 - (tx * stepSize);
    thrds_info[tx].wReal = wReal;
    thrds_info[tx].wImg = wImg;
    thrds_info[tx].ib = ib; // Är detta för att float ej räcker??                                                                                                     
    thrds_info[tx].stepSize = stepSize;
    thrds_info[tx].lines = lines; //Size of vectors in v                                                                                                              
    thrds_info[tx].tx = tx; //Thread index                                                                                                                            
    thrds_info[tx].mtx = &mtx;
    thrds_info[tx].cnd = &cnd;
    thrds_info[tx].status = status;
    status[tx].val = 0;

    printf("thrds_info[%d].ib = %.4lf\n", tx, thrds_info[tx].ib);
    fflush(stdout);

    int r = thrd_create(thrds+tx, main_thrd, (void*) (thrds_info+tx));

    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
   }

   return 0;

}









