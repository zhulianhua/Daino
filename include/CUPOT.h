
#ifndef __CUPOT_H__
#define __CUPOT_H__



#ifdef GRAVITY



// *****************************************************************
// ** This header will be included by all CPU/GPU solvers **
// *****************************************************************


// include "Typedef" here since the header "DAINO.h" is NOT included in GPU solvers
#include "Typedef.h"


// experiments show that the following macros give lower performance even in Fermi GPUs
/* 
// faster integer multiplication in GPU
#if ( defined __CUDACC__  &&  __CUDA_ARCH__ >= 200 )
   #define __umul24( a, b )   ( (a)*(b) )
   #define  __mul24( a, b )   ( (a)*(b) )
#endif
*/


// ####################
// ## macros for SOR ##
// ####################
#if   ( POT_SCHEME == SOR )

// determine the version of the GPU Poisson solver (10to14cube vs. 16to18cube)
#if ( POT_GHOST_SIZE <= 3  ||  defined FERMI )
   #define USE_PSOLVER_10TO14
#endif

// blockDim.z for the GPU Poisson solver
#if   ( POT_GHOST_SIZE == 1 )
      #define POT_BLOCK_SIZE_Z      4
#elif ( POT_GHOST_SIZE == 2 )
      #define POT_BLOCK_SIZE_Z      2
#elif ( POT_GHOST_SIZE == 3 )
      #define POT_BLOCK_SIZE_Z      2
#elif ( POT_GHOST_SIZE == 4 )
   #ifdef USE_PSOLVER_10TO14
      #define POT_BLOCK_SIZE_Z      5
   #else
      #define POT_BLOCK_SIZE_Z      1
   #endif
#elif ( POT_GHOST_SIZE == 5 )
   #ifdef USE_PSOLVER_10TO14
      #ifdef FLOAT8
      #define POT_BLOCK_SIZE_Z      2
      #else
      #define POT_BLOCK_SIZE_Z      4
      #endif
   #else
      #define POT_BLOCK_SIZE_Z      1
   #endif
#endif


// ###################
// ## macros for MG ##
// ###################
#elif ( POT_SCHEME == MG  )

// blockDim.x for the GPU Poisson solver
#if   ( POT_GHOST_SIZE == 1 )
      #define POT_BLOCK_SIZE_X      256
#elif ( POT_GHOST_SIZE == 2 )
      #define POT_BLOCK_SIZE_X      256
#elif ( POT_GHOST_SIZE == 3 )
      #define POT_BLOCK_SIZE_X      256
#elif ( POT_GHOST_SIZE == 4 )
      #define POT_BLOCK_SIZE_X      256
#elif ( POT_GHOST_SIZE == 5 )
      #define POT_BLOCK_SIZE_X      256
#endif

#endif // POT_SCHEME


// blockDim.z for the GPU Gravity solver
#define GRA_BLOCK_SIZE_Z            2



#endif // #ifdef GRAVITY



#endif // #ifndef __CUPOT_H__
