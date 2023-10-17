#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <threads.h>
#include <complex.h>

#define PI 3.14159265358979323846 
typedef uint8_t TYPE_ATTR;  
typedef uint8_t TYPE_CONV; 

TYPE_ATTR ** attractors;
TYPE_CONV ** convergences;

TYPE_ATTR * attractor;
TYPE_CONV * convergence;

double complex StepLength(z, int d);
void GetRoots( float ** roots,  int d);

typedef struct {
  int val;
  char pad[60];
} int_padded;

typedef struct {
  int **attractors;
  int **convergences;
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
  int **attractors;
  int **convergences;
  int lines;
  int numThreads;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_check_t;

int main_thrd(void *args)
{
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  int **attractors = thrd_info->attractors; //Bytt typ till int, Ska bytas till färg, koordingater för testning, ska bli globala                                                                 
  int **convergences = thrd_info->convergences; //Bytt typ till int, Ska bytas till antalet, koordinater för testning, ska bli globala                                                                 
  //double *dReal = thrd_info->attractors; //Bytt typ till int, Data transfer matris ska bytas till färg eller antal iterationer                                                               
  //double *dImg = thrd_info->convergences; //Bytt typ till int, Data transfer matris ska bytas till färg eller antal iterationer!!                                                               
  int numThreads = thrd_info->numThreads;
  double ib = thrd_info->ib;
  double ie = thrd_info->ie;
  const double stepSize = thrd_info->stepSize;
  const int lines = thrd_info->lines;
  const int tx = thrd_info->tx;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;

  int i, j, k;
  double ix, jx;
  double complex z;
  struct TwoValues result;

  //Warn: Vi tror koordinaterna funkar som de ska ej 100% säkra
  //Kolla så att itereringen är rätt isf 
  //Spara bara siffran för roten. Det är inte så många olika rötter. 8 rötter behöver                                                                                                   
  //inget större än en character                                                                                                                                                        
                                                                                                                                                   
  for (i = tx, ix = ib; i < lines, ix >= ie; i += numThreads, ix -= (stepSize*numThreads)) { //Skickar in vart vi börjar i ib!                                                            
    *attractor = (TYPE_ATTR *) malloc(lines*sizeof(TYPE_ATTR)); //Vill vi initiera till -1
    *convergence = (TYPE_CONV *) malloc(lines*sizeof(TYPE_CONV));
                                                                                                                    
    for ( j = 0, jx = -2.; j < lines, jx <= 2.; ++j , jx += (stepSize*numThreads)){
       z = ix + jx * I;
       for (k = 0; k < 128; k++){ 
            
            //Ifsats kolla om x är nära origin
            if (0.001 < creal(z)  && creal(z) > 0.001 || 0.001 < cimag(z)  && cimag(z) > 0.001 ){
                convergence [j] = d + 1;
                break;
            }

            //Kolla om real or img part is bigger than 10^10
            if (1000000000 > creal(z)  && creal(z) < -1000000000  || 1000000000 > cimag(z)  && cimag(z) < -1000000000  ){
                convergence [j] = d + 1;
                break;
            }

            result = double complex StepLength(z, int d);
            //WARN: Blir det fler eller färre värdesiffror när vi kollar på funk värdet?
            //Vill vi verkligen kolla på funktionsvärdet
            //Breakas det correct och vill vi ha warning om det inte funkar?
            if (result.nom < 0.001*d && result.nom > -0.001*d){
                for(int ixd = 0; ixd < d; ixd++){
                    if((ix <=  (roots[0][ixd] + 0.001) || ix >=  (roots[0][ixd] - 0.001)) &&
                     (jx <=  (roots[1][ixd] + 0.001) || jx >=  (roots[1][ixd] - 0.001))){
                           convergence [j] = ixd;
                           break;
                        }
                }
                break;
            }

            if (k == 128){
                 convergence [j] = d + 1;
            }
            z = z - result.nom/result.denom;
       }
      }
        
    attractor[j] = k;

    mtx_lock(mtx); //Vi vet ej ordningen på trådar     

    attractors[i] = attractor;
    convergences[i] = convergence;

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
  int **attractors = thrd_info->attractors; //Ska bytas till färg eller antal iterationer                                                                                                        
  int **convergences = thrd_info->convergences; //Ska bytas till färg eller antal iterationer!!                                                                                                        
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
    free(attractors[ix]);
    free(convergences[ix]);
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
**attractors = (TYPE_ATTR **) malloc(lines*sizeof(TYPE_ATTR *));
**convergences = (TYPE_CONV **) malloc(lines*sizeof(TYPE_CONV *));
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
    thrds_info[tx].attractors = attractors;
    thrds_info[tx].convergences = convergences;
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
 thrd_info_check.attractors = attractors;
 thrd_info_check.convergences = convergences;
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

double complex StepLength(double ix, double jx, int d) {

    struct TwoValues result;

    result.value1 = 42;
    result.value2 = 99;
    return result;

    switch (d) {
    case 1:
        result.nom = x - 1;
        result.denom = 1;
        return result;
    case 2:
        result.nom = x*x - 1;
        result.denom = 2*x;
        return result;
    case 3:
        result.nom = x*x*x - 1;
        result.denom = 3*x*x;
        return result;
    case 4:
        result.nom = x*x*x*x - 1;
        result.denom = 4*x*x*x;
        return result;
    case 5:
        result.nom = x*x*x*x*x - 1;
        result.denom = 5*x*x*x*x;
        return result;
    case 6:
        result.nom = x*x*x*x*x*x - 1;
        result.denom = 6*x*x*x*x*x;
        return result;
    case 7:
        result.nom = x*x*x*x*x*x*x - 1;
        result.denom = 7*x*x*x*x*x*x;
        return result; 
    case 8:
        result.nom = x*x*x*x*x*x*x*x - 1;
        result.denom = 8*x*x*x*x*x*x*x;
        return result;
    case 9:
        result.nom = x*x*x*x*x*x*x*x*x - 1;
        result.denom = 9*x*x*x*x*x*x*x*x;
        return result;  
    case 10:
        result.nom = x*x*x*x*x*x*x*x*x*x - 1;
        result.denom = 10*x*x*x*x*x*x*x*x*x;
        return result;
    
    default:
        fprintf(stderr, "unexpected degree\n");
        exit(1);
    }
}


/*
    switch (d) {
    case 1:
        return 1.0;
    case 2:
        return 0.5*x + 1/(2*x);
    case 3:
        return 2*x/3 + 1/(3*x*x);
    case 4:
        return 3*x/4 + 1/(4*x*x*x);
    case 5:
        return 4*x/5 + 1/(5*x*x*x*x);
    case 6:
        return 5*x/6 + 1/(6*x*x*x*x*x);
    case 7:
        return 6*x/7 + 1/(7*x*x*x*x*x*x); 
    case 8:
        return 7*x/8 + 1/(8*x*x*x*x*x*x*x);
    case 9:
        return 8*x/9 + 1/(9*x*x*x*x*x*x*x*x);  
    case 10:
        return 9*x/10 + 1/(10*x*x*x*x*x*x*x*x*x);
    
    default:
        fprintf(stderr, "unexpected degree\n");
        exit(1);
    }

*/
