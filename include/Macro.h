
#ifndef __MACRO_H__
#define __MACRO_H__



// ****************************************************************************
// ** This header defines the symbolic constants and macros used in DAINO.   **
// ** For clarity, useless options defined in the makefile will be "undef"   **
// ** in the end of this file.                                               **
// ****************************************************************************


// ########################
// ## Symbolic Constants ##
// ########################
   
// option == NONE --> the option is turned off
#define NONE      0


// models
#define HYDRO     1
#define MHD       2
#define ELBDM     3


// hydrodynamic schemes
#define RTVD      1
#define WAF       2
#define MHM       3
#define MHM_RP    4
#define CTU       5


// data reconstruction schemes
#define PLM       1
#define PPM       2


// Riemann solvers
#define EXACT     1
#define ROE       2
#define HLLE      3
#define HLLC      4


// Poisson solvers
#define SOR       1
#define MG        2


// load-balance parallelization
#define HILBERT   1


// number of components in each cell (FLU_NIN/NOUT : number of input/output variables in the fluid solver)
#if   ( MODEL == HYDRO )
#  define NCOMP              5
#  define FLU_NIN        NCOMP
#  define FLU_NOUT       NCOMP

#elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  define NCOMP              8
#  define FLU_NIN        NCOMP
#  define FLU_NOUT       NCOMP

// for ELBDM, we do not need to transfer the density component into GPU here
#elif ( MODEL == ELBDM )
#  define NCOMP              3
#  define FLU_NIN    ( NCOMP-1 )
#  define FLU_NOUT   ( NCOMP-0 )

#else
#  error : ERROR : unsupported MODEL (please edit NCOMP, FLU_NIN, and FLU_NOUT in the new MODEL) !!
#endif // MODEL


// main variables in different models
#if   ( MODEL == HYDRO )

// **variable indices in the array "fluid"**
#  define  DENS   0
#  define  MOMX   1
#  define  MOMY   2
#  define  MOMZ   3
#  define  ENGY   4

// **symbolic constants used as function parameters (e.g., Prepare_PatchGroupData)**
#  define _DENS   ( 1 << (DENS) ) 
#  define _MOMX   ( 1 << (MOMX) )
#  define _MOMY   ( 1 << (MOMY) )
#  define _MOMZ   ( 1 << (MOMZ) )
#  define _ENGY   ( 1 << (ENGY) )
#  define _FLU    ( _DENS | _MOMX | _MOMY | _MOMZ | _ENGY )

#  ifdef GRAVITY
#  define _POTE   ( 1 << (NCOMP+0) )
#  endif

#  ifdef OOC
#  define _FLUX   ( 1 << (NCOMP+1) )
#  define _STRC   ( 1 << (NCOMP+2) )
#  define _REAL   ( 1 << (NCOMP+3) )
#  define _BUFF   ( 1 << (NCOMP+4) )
#  endif


#elif ( MODEL == MHD )
#  warning : WAIT MHD !!!


#elif ( MODEL == ELBDM )
// **variable indices in the array "fluid"**
#  define  DENS   0
#  define  REAL   1
#  define  IMAG   2

// **symbolic constants used as function parameters (e.g., Prepare_PatchGroupData)**
#  define _DENS   ( 1 << (DENS) ) 
#  define _REAL   ( 1 << (REAL) )
#  define _IMAG   ( 1 << (IMAG) )
#  define _FLU    ( _DENS | _REAL | _IMAG )

#  ifdef GRAVITY
#  define _POTE   ( 1 << (NCOMP) )
#  endif

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL


// number of fluid ghost zones for the fluid solver
#if   ( MODEL == HYDRO )   // hydro
#  if   ( FLU_SCHEME == RTVD )
#        define FLU_GHOST_SIZE      3
#  elif ( FLU_SCHEME == WAF )
#        define FLU_GHOST_SIZE      2
#  elif ( FLU_SCHEME == MHM )
#     if ( LR_SCHEME == PLM )
#        define FLU_GHOST_SIZE      2
#     else // PPM
#        define FLU_GHOST_SIZE      3
#     endif
#  elif ( FLU_SCHEME == MHM_RP )
#     if ( LR_SCHEME == PLM )
#        define FLU_GHOST_SIZE      3
#     else // PPM
#        define FLU_GHOST_SIZE      4
#     endif
#  elif ( FLU_SCHEME == CTU )
#     if ( LR_SCHEME == PLM )
#        define FLU_GHOST_SIZE      2
#     else // PPM
#        define FLU_GHOST_SIZE      3
#     endif
#  endif

#elif ( MODEL == MHD )     // MHD
#        warning : WAIT MHD !!!
#        define FLU_GHOST_SIZE      ?

#elif ( MODEL == ELBDM )   // ELBDM
#        define FLU_GHOST_SIZE      3

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL


// self-gravity constants
#ifdef GRAVITY

// number of input and output variables in the gravity solver
#  if   ( MODEL == HYDRO )
#     define GRA_NIN               NCOMP

#  elif ( MODEL == MHD )
#     warning : WAIT MHD !!!
#     define GRA_NIN               NCOMP

// for ELBDM, we do not need to transfer the density component
#  elif ( MODEL == ELBDM )
#     define GRA_NIN             ( NCOMP-1 )

#  else
#     error Error : unsupported MODEL (please edit GRA_NIN in the new MODEL) !!
#  endif // MODEL


// number of potential ghost zones for evaluating potential (maximum=5) ~ Poisson solver
#     define POT_GHOST_SIZE      5


// number of potential ghost zones for advancing fluid by gravity ~ Gavity solver
#  if   ( MODEL == HYDRO )
#     define GRA_GHOST_SIZE      1

#  elif ( MODEL == MHD )
#     warning : WAIT MHD !!!
#     define GRA_GHOST_SIZE      ?

#  elif ( MODEL == ELBDM )
#     define GRA_GHOST_SIZE      0

#  else
#     error : ERROR : unsupported MODEL !!
#  endif // MODEL


// number of density ghost zones for the Poisson solver
#     define RHO_GHOST_SIZE      ( POT_GHOST_SIZE-1 )

#else
#     define POT_GHOST_SIZE      0
#     define GRA_GHOST_SIZE      0
#     define RHO_GHOST_SIZE      0
#endif


// patch size (number of cells of a single patch in the x/y/z directions)
#define PATCH_SIZE                   8
#define PS1             ( 1*PATCH_SIZE )
#define PS2             ( 2*PATCH_SIZE )


// the size of arrays (in one dimension) sending into GPU
//###REVISE: support interpolation schemes requiring 2 ghost cells in each side for POT_NXT
#define FLU_NXT         ( 2*(PATCH_SIZE+FLU_GHOST_SIZE)   )
#define POT_NXT         ( PATCH_SIZE/2 + 2*( (POT_GHOST_SIZE+3)/2 ) ) // assuming interpolation ghost-zone == 1
#define RHO_NXT         ( PATCH_SIZE   + 2*RHO_GHOST_SIZE )
#define GRA_NXT         ( PATCH_SIZE   + 2*GRA_GHOST_SIZE )


// constant to ensure the positive pressure
#ifdef FLOAT8
#  define MIN_VALUE        1.e-15
#else
#  ifdef COMOVING
#  define MIN_VALUE        1.e-10f
#  else
#  define MIN_VALUE        1.e-06f
#  endif
#endif


// extreme values
#ifndef __INT_MAX__
#  define __INT_MAX__      2147483647
#endif

#ifndef __UINT_MAX__
#  define __UINT_MAX__     ( __INT_MAX__ * 2U + 1U )
#endif

#ifndef __FLT_MAX__
#  define __FLT_MAX__      3.40282347e+38F
#endif

#ifndef __FLT_MIN__
#  define __FLT_MIN__      1.17549435e-38F
#endif


// NULL values
#ifndef NULL
#  define NULL             0
#endif

#ifndef NULL_INT
#  define NULL_INT         __INT_MAX__
#endif

#ifndef NULL_REAL
#  define NULL_REAL        __FLT_MAX__ 
#endif

#ifndef NULL_BOOL
#  define NULL_BOOL        false 
#endif


// macro for the function "Aux_Error"
#define ERROR_INFO         __FILE__, __LINE__, __FUNCTION__


// rank information of the entire domain decomposition ( MPI + OOC )
#ifdef OOC
#  define DAINO_NRANK            ( MPI_NRank*ooc.NRank                                          )
#  define DAINO_NRANK_X( dim )   ( MPI_NRank_X[dim]*ooc.NRank_X[dim]                            )
#  define DAINO_RANK             ( MPI_Rank*ooc.NRank + ooc.Rank                                )
#  define DAINO_RANK_X( dim )    ( MPI_Rank_X[dim]*ooc.NRank_X[dim] + ooc.Rank_X[ooc.Rank][dim] )
#else
#  define DAINO_NRANK            ( MPI_NRank        )
#  define DAINO_NRANK_X( dim )   ( MPI_NRank_X[dim] )
#  define DAINO_RANK             ( MPI_Rank         )
#  define DAINO_RANK_X( dim )    ( MPI_Rank_X[dim]  )
#endif



// ############
// ## Macros ##
// ############

// single/double-precision mathematic functions
#ifdef FLOAT8
#  define  FABS( a )        fabs( a )
#  define  SQRT( a )        sqrt( a )
#  define  FMAX( a, b )     fmax( a, b )
#  define  FMIN( a, b )     fmin( a, b )
#  define   POW( a, b )      pow( a, b )
#  define   SIN( a )         sin( a )
#  define   COS( a )         cos( a )
#  define ATAN2( a, b )    atan2( a, b )
#else
#  define  FABS( a )        fabsf( a )
#  define  SQRT( a )        sqrtf( a )
#  define  FMAX( a, b )     fmaxf( a, b )
#  define  FMIN( a, b )     fminf( a, b )
#  define   POW( a, b )      powf( a, b )
#  define   SIN( a )         sinf( a )
#  define   COS( a )         cosf( a )
#  define ATAN2( a, b )    atan2f( a, b )
#endif


// sign function
#define SIGN( a )       (  ( (a) < (real)0.0 ) ? (real)-1.0 : (real)+1.0  )


// max/min functions
#define MAX( a, b )     (  ( (a) > (b) ) ? (a) : (b)  )
#define MIN( a, b )     (  ( (a) < (b) ) ? (a) : (b)  )


// square function
#define SQR( a )        ( (a)*(a) )
#define CUBE( a )       ( (a)*(a)*(a) )


// 3D to 1D array indices transformation
#define IDX321( i, j, k, Ni, Nj )   (  ( (k)*(Nj) + (j) )*(Ni) + (i)  )



// ################################
// ## Remove useless definitions ##
// ################################
#if ( MODEL == HYDRO )
#  if ( FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU )
#  undef LR_SCHEME
#  endif

#  if ( FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU  &&  FLU_SCHEME != WAF )
#  undef RSOLVER
#  endif

#elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#endif // MODEL

#if ( MODEL != HYDRO  &&  MODEL != MHD )
#  undef FLU_SCHEME
#  undef LR_SCHEME
#  undef RSOLVER
#endif

#ifndef GRAVITY
#  undef POT_SCHEME
#endif



#endif  // #ifndef __MACRO_H__
