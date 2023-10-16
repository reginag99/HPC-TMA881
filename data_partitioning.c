//Recall that the result of a single computation in the video was represented by an array

//int * result;
//Since items correspond to rows in this implementation, each result has to include the attractor index and the convergence 
//index for each pixel in a row.

//Given the assumption that the degree and hence number of attractors is small, what is the smallest data type that can 
//encode the attractor index associated with a pixel? Call this type TYPE_ATTR.

//Given the assumption that the number of iterations may be capped at a small value, what is the smallest data type that 
//can encode the (capped) number of iterations till convergence associated with a pixel? Call this type TYPE_CONV.

//Implement, allocate, and free two global arrays

//TYPE_ATTR ** attractors;
//TYPE_CONV ** convergences;
//and corresponding local arrays in the computation and writing threads:

//TYPE_ATTR * attractor;
//TYPE_CONV * convergence;

#include <stdatomic.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <threads.h>


double
rand_double()
{
    return (double)rand()/RAND_MAX;
}


typedef struct {
  atomic_int *nthrds_completed;
  atomic_char *thrds_block_even;
  atomic_char *thrds_block_odd;
  int nthrds;
  int niter;
  int sz;
  const double **m;
  double *v;
  double *w;
  double *norm;
  double *norms;
  int ib;
  int ie;
} thrd_info_t;

int
main_thrd(
    void *args
    )
{
  const thrd_info_t *thrd_info = (thrd_info_t*) args;
  const int nthrds = thrd_info->nthrds;
  const int niter = thrd_info->niter;
  const int sz = thrd_info->sz;
  const double **m = thrd_info->m;
  double *v = thrd_info->v;
  double *w = thrd_info->w;
  const int ib = thrd_info->ib;
  const int ie = thrd_info->ie;


  for ( int lx = 0; lx < niter; ++lx ) {
    double norm = 0.;
    for ( int ix = ib; ix < ie; ++ix ) {
      const double *mix = m[ix];
      double wix = 0.;
      for ( int kx = 0; kx < sz; ++kx )
        wix += mix[kx] * v[kx];
      w[ix] = wix;
      norm += wix > 0 ? wix : -wix;
    }
    *thrd_info->norm = norm;

    // Swap pointers v and w.
    { double *tmp = v;
      v = w;
      w = tmp;
    }

    // Synchronize threads
    int nthrds_completed = 1 + atomic_fetch_add_explicit(
        thrd_info->nthrds_completed, 1, memory_order_relaxed);
    if ( nthrds_completed == nthrds ) {
      // Synchronize norm computation
      double norm = thrd_info->norms[0];
      for ( int tx = 1; tx < nthrds; ++tx )
        norm += thrd_info->norms[tx];
      for ( int tx = 0; tx < nthrds; ++tx )
        thrd_info->norms[tx] = norm;

      // We can reset the number of completed threads, since we will enforce
      // consistentcy when setting the other flag.
      atomic_store_explicit(thrd_info->nthrds_completed, 0,
          memory_order_relaxed);

      // We can reset this flag in relaxed order, since we will enforce
      // consistentcy when setting the other flag.
      atomic_store_explicit(thrd_info->thrds_block_odd, 1,
          memory_order_relaxed);
      atomic_store_explicit(thrd_info->thrds_block_even, 0,
          memory_order_release);
    } else {
      // We use a spin/busy lock, since we expect that threads take about the
      // same time.
      while ( atomic_load_explicit(thrd_info->thrds_block_even,
                                   memory_order_acquire) );
    }

    // Do not normalize in the last step so that we can extract the eigenvalue.
    if ( lx + 1 == niter ) break;

    // Normalize the vector v to avoid numerical instability.
    norm = *thrd_info->norm;
    for ( int ix = ib; ix < ie; ++ix )
      v[ix] /= norm;

    // Synchronize threads
    nthrds_completed = 1 + atomic_fetch_add_explicit(
        thrd_info->nthrds_completed, 1, memory_order_relaxed);
    if ( nthrds_completed == nthrds ) {

      // We can reset the number of completed threads, since we will enforce
      // consistentcy when setting the other flag.
      atomic_store_explicit(thrd_info->nthrds_completed, 0,
          memory_order_relaxed);

      // We can reset this flag in relaxed order, since we will enforce
      // consistentcy when setting the other flag.
      atomic_store_explicit(thrd_info->thrds_block_even, 1,
          memory_order_relaxed);
      atomic_store_explicit(thrd_info->thrds_block_odd, 0,
          memory_order_release);
    } else {
      // We use a spin/busy lock, since we expect that threads take about the
      // same time.
      while ( atomic_load_explicit(thrd_info->thrds_block_odd,
                                   memory_order_acquire) );
    }
    
  }
}


int
main()
{
  const int nthrds = 4;
  thrd_t thrds[nthrds];
  thrd_info_t thrds_info[nthrds];

  const int sz = 1000;
  const int niter = 1000;

  double *mentries;
  double **m;
  double *v;
  double *w;

  mentries = (double*) malloc(sz*sz*sizeof(double));
  m = (double**) malloc(sz*sizeof(double*));
  v = (double*) malloc(sz*sizeof(double));
  w = (double*) malloc(sz*sizeof(double));

  for ( int ix = 0, jx = 0; ix < sz; ++ix, jx += sz )
    m[ix] = mentries + jx;

  srand(clock());
  for ( double *e = mentries, *end = mentries + sz*sz; e != end; ++e )
    *e = rand_double();
  for ( double *e = v, *end = v + sz; e != end; ++e )
    *e = rand_double();

  atomic_int nthrds_completed;
  atomic_char thrds_block_even;
  atomic_char thrds_block_odd;
  atomic_init(&nthrds_completed, 0);
  atomic_init(&thrds_block_even, 1);

  double norms[nthrds];

  const int szloc = sz / nthrds;
  for ( int tx = 0, ib = 0; tx < nthrds; ++tx, ib += szloc ) {
    thrds_info[tx].nthrds_completed = &nthrds_completed;
    thrds_info[tx].thrds_block_even = &thrds_block_even;
    thrds_info[tx].thrds_block_odd = &thrds_block_odd;
    thrds_info[tx].nthrds = nthrds;
    thrds_info[tx].niter = niter;
    thrds_info[tx].sz = sz;
    thrds_info[tx].m = (const double**) m;
    thrds_info[tx].v = v;
    thrds_info[tx].w = w;
    norms[tx] = 1.;
    thrds_info[tx].norm = norms+tx;
    thrds_info[tx].norms = norms;
    thrds_info[tx].ib = ib;
    thrds_info[tx].ie = tx != nthrds - 1 ? ib + szloc : sz;

    int r = thrd_create(thrds+tx, main_thrd, (void*) (thrds_info+tx));
    if ( r != thrd_success ) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
  }

  for ( int tx = 0; tx < nthrds; ++tx ) {
    int r;
    thrd_join(thrds[tx], &r);
  }


  for ( int ix = 0; ix < 10; ++ix )
    printf("%f ", v[ix] / w[ix]);
  printf("\n");

  free(m);
  free(mentries);
  free(v);
  free(w);

  return 0;
}