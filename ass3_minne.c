#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <threads.h>
#include <complex.h>
#include <stdint.h>

#define PI 3.14159//265358979323846;                                                                                                                                                    
typedef uint8_t TYPE_ATTR;
typedef uint8_t TYPE_CONV;

struct result{
    double complex nom;
    double complex denom;
    };

TYPE_ATTR ** attractors;
TYPE_CONV ** convergences;

TYPE_ATTR * attractor;
TYPE_CONV * convergence;

void StepLength(double complex z,double complex nom,double complex denom, int d);
void GetRoots( float ** roots,  int d);
//void GetColors(TYPE_CONV** color, int d);
//void GetGrayScale(TYPE_ATTR**grayscale, int d);

typedef struct {
  int val;
  char pad[60];
} int_padded;

typedef struct {
  TYPE_ATTR **attractors;
  TYPE_CONV **convergences;
  float**roots;
  int numThreads;
  double ib;//Ska denna vara const?                                                                                                                                                     
  double ie;
  int d;
  float stepSize;
  int lines;
  int tx;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_t;

typedef struct {
  TYPE_ATTR **attractors;
  TYPE_CONV **convergences;
  int lines;
  int d;
  int numThreads;
  mtx_t *mtx;
  cnd_t *cnd;
  int_padded *status;
} thrd_info_check_t;

int main_thrd(void *args)
{
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  TYPE_ATTR **attractors = thrd_info->attractors; //Bytt typ till int, Ska bytas till färg, koordingater för testning, ska bli globala                                                 \
                                                                                                                                                                                        
  TYPE_CONV **convergences = thrd_info->convergences; //Bytt typ till int, Ska bytas till antalet, koordinater för testning, ska bli globala                                           \
                                                                                                                                                                                        
  float**roots = thrd_info->roots;
  int numThreads = thrd_info->numThreads;
  double ib = thrd_info->ib;
  double ie = thrd_info->ie;
  int d = thrd_info->d;
  const double stepSize = thrd_info->stepSize;
  const int lines = thrd_info->lines;
  const int tx = thrd_info->tx;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;

  int i, j, k;
  double ix, jx;
  double complex z,nom,denom;

  //Warn: Vi tror koordinaterna funkar som de ska ej 100% säkra                                                                                                                         
  //Kolla så att itereringen är rätt isf                                                                                                                                                
  //Spara bara siffran för roten. Det är inte så många olika rötter. 8 rötter behöver                                                                                                  \
                                                                                                                                                                                        
  //inget större än en character  
  
  for (i = tx, ix = ib; i < lines, ix >= ie; i += numThreads, ix -= (stepSize*numThreads)) { //Skickar in vart vi börjar i ib!                                                         \
                                                                                                                                                                                        
    TYPE_ATTR*attractor = (TYPE_ATTR *) malloc(lines*sizeof(TYPE_ATTR)); //Vill vi initiera till -1                                                                                     
    TYPE_CONV*convergence = (TYPE_CONV *) malloc(lines*sizeof(TYPE_CONV));

    for ( size_t cx = 0; cx < lines; ++cx ) {
        attractor[cx] = 0;
        convergence[cx] = 0;
        }

    for ( j = 0, jx = -2.; j < lines, jx <= 2.; ++j , jx += (stepSize*numThreads)){
       z = ix + jx * I;
       //Print worked here!!                                                                                                                                                            
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
            StepLength(z,nom,denom, d);
            //WARN: Blir det fler eller färre värdesiffror när vi kollar på funk värdet?                                                                                                
            //Vill vi verkligen kolla på funktionsvärdet                                                                                                                                
            //Breakas det correct och vill vi ha warning om det inte funkar?                                                                                                            
            //För dyrt med absolutbelopp?
            if (cabs(nom) < 0.001*d && cabs(nom) > -0.001*d){
                for(int ixd = 0; ixd < d; ixd++){
                    //Blir det typfel här? Måste vi använda complext tal??                                                                                                              
                    if((creal(z) <=  (roots[ixd][0] + 0.001) || creal(z) >=  (roots[ixd][0] - 0.001)) &&
                     (cimag(z) <=  (roots[ixd][1] + 0.001) || cimag(z) >=  (roots[ixd][1] - 0.001))){
                           attractor [j] = ixd;
                           break;
                        }
                }
                break;
            }

            if (k == 128){
                 convergence [j] = d + 1;
            }
            z = z - nom/denom;
            //Print does not work here!!                                                                                                                                                
            printf("Real: %f\n", creal(z));     // Print the real part                                                                                                                  
            printf("Imaginary: %f\n", cimag(z)); //                                                                                                                                     
       }

      }

    mtx_lock(mtx); //Vi vet ej ordningen på trådar

    attractors[i] = attractor; //Detta är data transfer? Hur garanterar vi att den skriver på rätt rad?                                                                                 
    convergences[i] = convergence; //Detta är data transfer? Hur garanterar vi att den skriver på rätt rad?                                                                             

    status[tx].val = i + tx;
    //Matris med antalet trådar rader varje tråd lägger sin rad i sin kolumn?                                                                                                          \
                                                                                                                                                                                        
    //Varje tråd lägger sin pixel i sin column.                                                                                                                                        \
                                                                                                                                                                                        
    // Vill vi göra data transfer här??                                                                                                                                                \
                                                                                                                                                                                        
    // Vi kan skicka iväg tråden :)                                                                                                                                                    \
                                                                                                                                                                                        
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
  TYPE_ATTR **attractors = thrd_info->attractors;
  TYPE_CONV **convergences = thrd_info->convergences;
  const int lines = thrd_info->lines;
  const int d = thrd_info->d;
  const int numThreads = thrd_info->numThreads;
  mtx_t *mtx = thrd_info->mtx;
  cnd_t *cnd = thrd_info->cnd;
  int_padded *status = thrd_info->status;

  //Vart vill vi öppna filen?                                                                                                                                                           
  FILE *file = fopen("output.ppm", "w");
    if (file == NULL) {
    perror("Error opening the file");
    exit(EXIT_FAILURE);
    }

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
    int grey_step = 255/gray_index;

    int gray_coordinates = 3;  
    int gray_array[gray_index][gray_coordinates];

  int grayIndex, colorIndex;
/*
  TYPE_ATTR *colorEntries = (TYPE_ATTR*) malloc(sizeof(TYPE_ATTR)* (d + 1)* 3);
  TYPE_ATTR** color = (TYPE_ATTR**) malloc(sizeof(TYPE_ATTR*) * 3);

  for ( size_t ix = 0, jx = 0; ix < 3; ++ix, jx+=(d+1))
        color[ix] = colorEntries + jx;

  //GetColors(color, d);

  TYPE_CONV *grayscaleEntries = (TYPE_CONV*) malloc(sizeof(TYPE_CONV)* (d + 1)* 3);
  TYPE_CONV** grayscale = (TYPE_CONV**) malloc(sizeof(TYPE_CONV*) * 3);

  for ( size_t ix = 0, jx = 0; ix < 3; ++ix, jx+=(d+1))
        grayscale[ix] = grayscaleEntries + jx;

  //GetGrayScale(grayscale, d);
  //printf("%d", grayscale);

  TYPE_CONV* pixelBufferGrayEntries = (TYPE_CONV*) malloc(sizeof(TYPE_CONV)* lines * 3);
  TYPE_CONV**pixelBufferGray= (TYPE_CONV**) malloc(sizeof(TYPE_CONV*) * lines);

  for ( size_t ix = 0, jx = 0; ix < lines; ++ix, jx+= 3)
    pixelBufferGray[ix] = pixelBufferGrayEntries + jx;


  TYPE_ATTR*pixelBufferColor = (TYPE_ATTR*) malloc(sizeof(TYPE_ATTR)* lines);

  int grayIndex, colorIndex;
  TYPE_CONV Rg,Gg,Bg;
  TYPE_ATTR Rc,Gc,Bc;

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
*/
    //fprintf(stderr, "checking until %i\n", ibnd); 

    for ( ; ix < ibnd; ++ix ) {
        for(int jx = 0; jx < lines; ++jx){
            grayIndex = attractors[ix][jx];
            colorIndex = attractors[ix][jx];

            Rc = color[colorIndex][0];
            Gc = color[colorIndex][1];
            Bc = color[colorIndex][2];

            Rg = gray_index*grey_step;
            Gg = gray_index*grey_step;
            Bg = gray_index*grey_step;
            //Skriv till en buffer                                                                                                                                                      

            pixelBufferGray[jx][0] = Rg;
            pixelBufferGray[jx][1] = Gg;
            pixelBufferGray[jx][2] = Bg;
        }
    //Här skriver vi hela buffer till filen!!                                                                                                                                           
    //Pixel buffer har samma storlek som lines!                                                                                                                                         
    //fwrite(pixelBufferGray, sizeof(TYPE_CONV), lines, file);                                                                                                                          
    //fwrite(pixelBuffer, sizeof(char), strlen(pixelBuffer), file);                                                                                                                     
    //Måste vi tömma pixelbuffer efter varje iteration?                                                                                                                                 

    free(attractors[ix]);
    free(convergences[ix]);
  }
 }

 //Vill vi stänga filen här?                                                                                                                                                            
 fclose(file);
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

  printf("numThreads = %d \n", numThreads);
  printf("lines = %d \n", lines);
  printf("d = %d \n", d);

float *rootsEntries = (float*) malloc(sizeof(float)* d * 2);
float ** roots = (float**) malloc(sizeof(float*)* 2); //Rätt?                                                                                                                           

for ( size_t ix = 0, jx = 0; ix < 2; ++ix, jx+=d)
    roots[ix] = rootsEntries + jx;

for ( size_t ix = 0; ix < 2; ++ix )
    for ( size_t jx = 0; jx < d; ++jx )
        roots[ix][jx] = 0; 

//Tar ut rötterna för polynomet                                                                                                                                                         
GetRoots(roots, d);
printf("roots = %f + %fi \n", roots[0][0], roots[0][1]);

double stepSize = 1/lines;
const double sz = lines/numThreads; //Används ej!
                                                                                                                                                                                        
TYPE_ATTR**attractors = (TYPE_ATTR **) malloc(lines*sizeof(TYPE_ATTR *));
TYPE_CONV**convergences = (TYPE_CONV **) malloc(lines*sizeof(TYPE_CONV *));
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
    thrds_info[tx].roots=roots;
    thrds_info[tx].numThreads = numThreads; //Casta till const? 
                                                                                                                                                                                        
    thrds_info[tx].ib = ib; //Ska vi casta till const? 
                                                                                                                                                                                        
    thrds_info[tx].ie = ie; //Ska vi casta till const?                                                                                                                                  
    thrds_info[tx].d = d;
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
 thrd_info_check.d = d;
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


double complex StepLength(double complex z, int d) {
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
void GetRoots( float ** roots, int d) {
     // Hårdkoda de riktiga rötterna alltså 1 och -1? 
     // Beror på om d är jämnt eller ej
    if (d % 2 == 0) {
        roots[0][0] = 1.;
        roots[0][1] = 0.;
        roots[1][0] = -1.;
        roots[1][1] = 0.;

     if (d > 2){
        for( size_t ix = 2; ix < d; ix++) {
            float theta = (ix*2* PI) / d ;        
            roots[ix][0] = cos(theta);
            roots[ix][1] = sin(theta);
        }
    }

    } else {
        roots[0][0] = 1.; 
        roots[0][1] = 0.;

        if (d > 1){
        for( size_t ix = 2; ix < d; ix++) {
            float theta = (ix*2* PI) / d ;        
            roots[ix][0] = cos(theta);
            roots[ix][1] = sin(theta);
            }
        } 
    }

}
/*
void GetColors(TYPE_CONV** color, int d);{
TYPE_ATTR *colorEntries = (TYPE_ATTR*) malloc(sizeof(TYPE_ATTR)* (d + 1)* 3);
TYPE_ATTR** color = (TYPE_ATTR**) malloc(sizeof(TYPE_ATTR*) * 3);

for ( size_t ix = 0, jx = 0; ix < 3; ++ix, jx+=(d+1))
  color[ix] = colorEntries + jx;


    color = {
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
}


 



  //GetColors(color, d);

  TYPE_CONV *grayscaleEntries = (TYPE_CONV*) malloc(sizeof(TYPE_CONV)* (d + 1)* 3);
  TYPE_CONV** grayscale = (TYPE_CONV**) malloc(sizeof(TYPE_CONV*) * 3);

  for ( size_t ix = 0, jx = 0; ix < 3; ++ix, jx+=(d+1))
        grayscale[ix] = grayscaleEntries + jx;


void GetGrayScale(TYPE_ATTR**grayscale, int d);
*/
