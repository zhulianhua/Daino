
#ifndef __TYPEDEF_H__
#define __TYPEDEF_H__



// ****************************************************************************
// ** This header defines the "typedef and enum" in DAINO.                   **
// ** Please DO NOT modify the number assigned to any enumerator constant !! ** 
// ****************************************************************************

#include "Macro.h"


// single/double precision
#ifdef FLOAT8
typedef double real;
#else
typedef float  real;
#endif


// short names for unsigned type
typedef unsigned short     ushort;
typedef unsigned int       uint;
typedef unsigned long int  ulong;


// options of program initialization
enum OptInit_t { INIT_STARTOVER=1, INIT_RESTART=2, INIT_UM=3 };


// options of program restart
enum OptRestartH_t { RESTART_HEADER_SKIP=0, RESTART_HEADER_CHECK=1 };


// interpolation scheme
enum IntScheme_t { INT_DEFAULT=-1, INT_CENTRAL=1, INT_MINMOD=2, INT_VANLEER=3, INT_CQUAD=4, INT_QUAD=5,
                   INT_CQUAR=6, INT_QUAR=7 };


// data reconstruction TVD limiters
enum LR_Limiter_t  { LR_LIMITER_NONE=0, VANLEER=1, GMINMOD=2, ALBADA=3, VL_GMINMOD=4, EXTPRE=5 };
enum WAF_Limiter_t { WAF_LIMITER_NONE=0, WAF_SUPERBEE=1, WAF_VANLEER=2, WAF_ALBADA=3, WAF_MINBEE=4 };


// options of data output
enum OptOutputMode_t { OUTPUT_CONST_STEP=1, OUTPUT_CONST_DT=2, OUTPUT_USE_TABLE=3 };


// options of output part
enum OptOutputPart_t { OUTPUT_NONE=0, OUTPUT_XY=1, OUTPUT_YZ=2, OUTPUT_XZ=3, OUTPUT_X=4, OUTPUT_Y=5, OUTPUT_Z=6,
                       OUTPUT_DIAG=7 };


// options in "Prepare_PatchGroupData"
enum PrepUnit_t { UNIT_PATCH=1, UNIT_PATCHGROUP=2 };
enum NSide_t    { NSIDE_06=6, NSIDE_26=26 };


// use the load-balance alternative function in "Buf_GetBufferData" and "Flag_Real"
enum UseLBFunc_t { USELB_NO=1, USELB_YES=2 };


// enable check or not
enum Check_t { CHECK_ON=1, CHECK_OFF=2 };


// targeted solver in "InvokeSolvers"
#ifdef GRAVITY
enum Solver_t { FLUID_SOLVER=0, POISSON_SOLVER=1, GRAVITY_SOLVER=2, POISSON_AND_GRAVITY_SOLVER=3 };
#else
enum Solver_t { FLUID_SOLVER=0 };
#endif


// targeted mode in "Buf_GetBufferData and LB_GetBufferData"
#ifdef GRAVITY
enum GetBufMode_t { DATA_GENERAL=1, DATA_AFTER_FIXUP=2, DATA_AFTER_REFINE=3, DATA_RESTRICT=4, COARSE_FINE_FLUX=5,
                    POT_FOR_POISSON=6, POT_AFTER_REFINE=7 };
#else
enum GetBufMode_t { DATA_GENERAL=1, DATA_AFTER_FIXUP=2, DATA_AFTER_REFINE=3, DATA_RESTRICT=4, COARSE_FINE_FLUX=5};
#endif // #ifdef GRAVITY ... else ...




#endif  // #ifndef __TYPEDEF_H__
