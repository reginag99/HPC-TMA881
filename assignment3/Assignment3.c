#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <threads.h>

#define PI 3.14159265358979323846 

typedef struct {
  int val;
  char pad[60];
} int_padded;


typedef struct {
  double **wReal;
  double **wImg;
  int numThreads;
  double ib;//Ska denna vara const?                                                                                                                                                       
  double ie;
  float stepSize;
  int lines;
  int tx;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_t;


typedef struct {
  double **wReal;
  double **wImg;
  int lines;
  int numThreads;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_check_t;

int main_thrd(void *args)
{
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  double **wReal = thrd_info->wReal; //Bytt typ till int, Ska bytas till färg, koordingater för testning, ska bli globala                                                                 
  double **wImg = thrd_info->wImg; //Bytt typ till int, Ska bytas till antalet, koordinater för testning, ska bli globala                                                                 
  //double *dReal = thrd_info->wReal; //Bytt typ till int, Data transfer matris ska bytas till färg eller antal iterationer                                                               
  //double *dImg = thrd_info->wImg; //Bytt typ till int, Data transfer matris ska bytas till färg eller antal iterationer!!                                                               
  int numThreads = thrd_info->numThreads;
  double ib = thrd_info->ib;
  double ie = thrd_info->ie;
  const double stepSize = thrd_info->stepSize;
  const int lines = thrd_info->lines;
  const int tx = thrd_info->tx;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;

  int i;
  double ix;

  int j;
  double jx;

  //ib ska gå från 2 till -2                                                                                                                                                              
  // Vi måste ha en ie då alla trådar inte slutar på -2                                                                                                                                   
  // i måste börja på sitt tråd index                                                                                                                                                     
  for (i = tx, ix = ib; i < lines, ix >= ie; i += numThreads, ix -= (stepSize*numThreads)) { //Skickar in vart vi börjar i ib!                                                            
    double *wix_real = (double*) malloc(lines*sizeof(double));
    double *wix_img = (double*) malloc(lines*sizeof(double));

    //WARN: Beräkna så att indexeringen är rätt bör ej vara stepSize                                                                                                                      
    for ( j = 0, jx = -2.; j < lines, jx <= 2.; ++j , jx += (stepSize*numThreads)){
      wix_real[j] = jx;
      wix_img[j] = ix;
      }

    //Spara bara siffran för roten. Det är inte så många olika rötter. 8 rötter behöver                                                                                                   
    //inget större än en character                                                                                                                                                        

    mtx_lock(mtx); //Vi vet ej ordningen på trådar                                                                                                                                        
    wReal[i] = wix_real;
    wImg[i] = wix_img;
    status[tx].val = i + tx;
    //Matris med antalet trådar rader varje tråd lägger sin rad i sin kolumn?                                                                                                             
    //Varje tråd lägger sin pixel i sin column.                                                                                                                                           
    // Vill vi göra data transfer här??                                                                                                                                                   
    // Vi kan skicka iväg tråden :)                                                                                                                                                       
    mtx_unlock(mtx);
    cnd_signal(cnd);
    thrd_sleep(&(struct timespec){.tv_sec=0, .tv_nsec=10}, NULL);
   }
  return 0;
}

int
main_thrd_check(
    void *args
    )
{
  const thrd_info_check_t *thrd_info = (thrd_info_check_t*) args;
  double **wReal = thrd_info->wReal; //Ska bytas till färg eller antal iterationer                                                                                                        
  double **wImg = thrd_info->wImg; //Ska bytas till färg eller antal iterationer!!                                                                                                        
  const int lines = thrd_info->lines;
  const int numThreads = thrd_info->numThreads;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;

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

    fprintf(stderr, "checking until %i\n", ibnd);

    for ( ; ix < ibnd; ++ix ) {
    free(wReal[ix]);
    free(wImg[ix]);
  }
 }

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

    for (int i = 1; i < 3; i++) { 
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

//För rötterna skapar matrs storlek = d * 2. AHA DÅ VI HAR EN KOLUMN FÖR KOMPLEXA OCH EN FÖR REELA
float *rootsEntries = (float*) malloc(sizeof(float)* d * 2);
float ** roots = (float**) malloc(sizeof(float)* d); //Rätt?

for ( size_t ix = 0, jx = 0; ix < d; ++ix, jx+=2)
    roots[ix] = rootsEntries + jx;

for ( size_t ix = 0; ix < d; ++ix )
    for ( size_t jx = 0; jx < 2; ++jx )
        roots[ix][jx] = 0;

    
    //Tar ut rötterna för polynomet
GetRoots( float ** roots, int d);


const float lines = 100.;
double stepSize = 1/lines;
int numThreads = 10;
const double sz = lines/numThreads; //Används ej!                                                                                                                                       
double **wReal = (double**) malloc(lines*sizeof(double*));
double **wImg = (double**) malloc(lines*sizeof(double*));
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
    double ie = -2.0 + (numThreads-tx)*stepSize;
    thrds_info[tx].wReal = wReal;
    thrds_info[tx].wImg = wImg;
    thrds_info[tx].numThreads = numThreads; //Casta till const?                                                                                                                           
    thrds_info[tx].ib = ib; //Ska vi casta till const?                                                                                                                                    
    thrds_info[tx].ie = ie; //Ska vi casta till const?                                                                                                                                    
    thrds_info[tx].stepSize = stepSize;
    thrds_info[tx].lines = lines; //Size of vectors in v                                                                                                                                  
    thrds_info[tx].tx = tx; //Thread index                                                                                                                                                
    thrds_info[tx].mtx = &mtx;
    thrds_info[tx].cnd = &cnd;
    thrds_info[tx].status = status;
    status[tx].val = -1;
    //printf("thrds_info[%d].ib = %.4lf\n", tx, thrds_info[tx].ib);                                                                                                                       

    int r = thrd_create(thrds+tx, main_thrd, (void*) (thrds_info+tx)); //Varför +tx?                                                                                                      

    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
   }

  {
    thrd_info_check.wReal = wReal;
    thrd_info_check.wImg = wImg;
    thrd_info_check.lines = lines;
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
