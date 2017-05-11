
#include "Macro.h"
#include "CUFLU.h"

#if ( defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == CTU )



#include "CUFLU_Shared_FluUtility.cu"
#include "CUFLU_Shared_DataReconstruction.cu"
#include "CUFLU_Shared_ComputeFlux.cu"
#include "CUFLU_Shared_FullStepUpdate.cu"

static __device__ void CUFLU_TGradient_Correction( real g_FC_Var_xL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                                   real g_FC_Var_xR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                                   real g_FC_Var_yL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                                   real g_FC_Var_yR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                   real g_FC_Var_zL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                                   real g_FC_Var_zR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                   const real g_FC_Flux_x[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                   const real g_FC_Flux_y[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                   const real g_FC_Flux_z[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                   const real dt, const real _dh );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FluidSolver_CTU
// Description :  GPU fluid solver based on the Corner-Transport-Upwind (CTU) scheme
//
// Note        :  1. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                2. Ref : Stone et al., ApJS, 178, 137 (2008)
//                3. Each patch group requires about 2.9*10^7 flops with CTU + PPM + Roe solver
//                   --> 206 GFLOPS is achieved in one C2050 GPU
//
// Parameter   :  g_Fluid_In     : Global memory array storing the input fluid variables
//                g_Fluid_Out    : Global memory array to store the output fluid variables
//                g_Flux         : Global memory array to store the output fluxes
//                g_PriVar       : Global memory array to store the primitive variables
//                g_Slope_PPM_x  : Global memory array to store the x-slope for the PPM reconstruction
//                g_Slope_PPM_y  : Global memory array to store the y-slope for the PPM reconstruction
//                g_Slope_PPM_z  : Global memory array to store the z-slope for the PPM reconstruction
//                g_FC_Var_xL    : Global memory array to store the half-step variables on the -x surface
//                g_FC_Var_xR    : Global memory array to store the half-step variables on the +x surface
//                g_FC_Var_yL    : Global memory array to store the half-step variables on the -y surface
//                g_FC_Var_yR    : Global memory array to store the half-step variables on the +y surface
//                g_FC_Var_zL    : Global memory array to store the half-step variables on the -z surface
//                g_FC_Var_zR    : Global memory array to store the half-step variables on the +z surface
//                g_FC_Flux_x    : Global memory array to store the face-centered fluxes in the x direction
//                g_FC_Flux_y    : Global memory array to store the face-centered fluxes in the y direction
//                g_FC_Flux_z    : Global memory array to store the face-centered fluxes in the z direction
//                dt             : Time interval to advance solution
//                _dh            : 1 / grid size
//                Gamma          : Ratio of specific heats
//                StoreFlux      : true --> store the coarse-fine fluxes
//                LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                 (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//                EP_Coeff       : Coefficient of the extrema-preserving limiter
//-------------------------------------------------------------------------------------------------------
__global__ void CUFLU_FluidSolver_CTU( const real g_Fluid_In[][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                       real g_Fluid_Out  [][5][ PS2*PS2*PS2 ],
                                       real g_Flux    [][9][5][ PS2*PS2 ], 
                                       real g_PriVar     [][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Slope_PPM_x[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_y[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_z[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_FC_Var_xL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_xR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_zR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Flux_x  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_y  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_z  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, 
                                       const real EP_Coeff )
{

// 1. conserved variables --> primitive variables
   CUFLU_Con2Pri_AllGrids( g_Fluid_In, g_PriVar, Gamma );
   __syncthreads();


// 2. evaluate the half-step face-centered solution
   CUFLU_DataReconstruction( g_PriVar, g_Slope_PPM_x, g_Slope_PPM_y, g_Slope_PPM_z, g_FC_Var_xL, g_FC_Var_xR, 
                             g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR, FLU_NXT, FLU_GHOST_SIZE-1, 
                             Gamma, LR_Limiter, MinMod_Coeff, EP_Coeff, dt, _dh );
   __syncthreads();


// 3. evaluate the face-centered half-step fluxes by solving the Riemann problem
   CUFLU_ComputeFlux( g_FC_Var_xL, g_FC_Var_xR, g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR, 
                      g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z, NULL, false, 0, Gamma );
   __syncthreads();


// 4. correct the face-centered variables by the transverse flux gradients
   CUFLU_TGradient_Correction( g_FC_Var_xL, g_FC_Var_xR, g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR, 
                               g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z, dt, _dh );
   __syncthreads();


// 5. evaluate the face-centered full-step fluxes by solving the Riemann problem with the corrected data
   CUFLU_ComputeFlux( g_FC_Var_xL, g_FC_Var_xR, g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR, 
                      g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z, g_Flux, StoreFlux, 1, Gamma );
   __syncthreads();


// 6. evaluate the full-step solution
   CUFLU_FullStepUpdate( g_Fluid_In, g_Fluid_Out, g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z, dt, _dh, Gamma );

} // FUNCTION : CUFLU_FluidSolver_CTU



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_TGradient_Correction
// Description :  Correct the face-centered variables by the transverse flux gradients 
//
// Note        :  1. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                2. The function is asynchronous
//                   --> "__syncthreads" must be called before using the output data
//                3. The sizes of the arrays g_FC_Var_x/y/z and g_FC_Flux_XX in each direction are assumed
//                   to be "N_FC_VAR" (which is always PS2+2 in the currennt cases)
//
// Parameter   :  g_FC_Var_xL : Global memory array to store the input and output -x face-centered variables 
//                g_FC_Var_xR : Global memory array to store the input and output +x face-centered variables 
//                g_FC_Var_yL : Global memory array to store the input and output -y face-centered variables 
//                g_FC_Var_yR : Global memory array to store the input and output +y face-centered variables 
//                g_FC_Var_zL : Global memory array to store the input and output -z face-centered variables 
//                g_FC_Var_zR : Global memory array to store the input and output +z face-centered variables 
//                g_FC_Flux_x : Global memory array storing the face-centered fluxes in the x direction
//                g_FC_Flux_y : Global memory array storing the face-centered fluxes in the y direction
//                g_FC_Flux_z : Global memory array storing the face-centered fluxes in the z direction
//                dt          : Time interval to advance solution
//                _dh         : 1 / grid size
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_TGradient_Correction( real g_FC_Var_xL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                            real g_FC_Var_xR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                            real g_FC_Var_yL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                            real g_FC_Var_yR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                            real g_FC_Var_zL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                            real g_FC_Var_zR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                            const real g_FC_Flux_x[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                            const real g_FC_Flux_y[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                            const real g_FC_Flux_z[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                            const real dt, const real _dh )
{

   const uint  bx         = blockIdx.x;
   const uint  tx         = threadIdx.x;
   const uint  dID_FC_Var = blockDim.x;
   const uint3 dID_In     = make_uint3( 1, N_FC_VAR, N_FC_VAR*N_FC_VAR );
   const real  dt_dh2     = (real)0.5*dt*_dh;

   uint  ID_In_xL, ID_In_yL, ID_In_zL, ID_In_R, ID_FC_Var;
   uint3 ID3d;
   real  FC_Var_xL, FC_Var_xR, FC_Var_yL, FC_Var_yR, FC_Var_zL, FC_Var_zR;
   real  FC_Flux_xL, FC_Flux_xR, FC_Flux_yL, FC_Flux_yR, FC_Flux_zL, FC_Flux_zR;
   real  TGrad1, TGrad2, Corr;
   bool  Inner_x, Inner_y, Inner_z;
   bool  Inner_xy, Inner_yz, Inner_xz;


#  define Load( Input, Output, ID, v )    (  Output = Input[bx][v][ID]  )
#  define Dump( Input, Output, ID, v )    (  Output[bx][v][ID] = Input  )

#  define Correct( FC_Var_L, FC_Var_R, FC_Flux_T1_L, FC_Flux_T1_R, FC_Flux_T2_L, FC_Flux_T2_R )    \
   {                                                                                               \
      TGrad1 = FC_Flux_T1_R - FC_Flux_T1_L;                                                        \
      TGrad2 = FC_Flux_T2_R - FC_Flux_T2_L;                                                        \
      Corr   = -dt_dh2*( TGrad1 + TGrad2 );                                                        \
                                                                                                   \
      FC_Var_L += Corr;                                                                            \
      FC_Var_R += Corr;                                                                            \
   } // Correct

#  define Correct_1v( v )                                                                                   \
   {                                                                                                        \
      /* load the face-centered variables */                                                                \
      if ( Inner_yz )                                                                                       \
      {                                                                                                     \
         Load( g_FC_Var_xL, FC_Var_xL, ID_FC_Var, v );                                                      \
         Load( g_FC_Var_xR, FC_Var_xR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_xz )                                                                                       \
      {                                                                                                     \
         Load( g_FC_Var_yL, FC_Var_yL, ID_FC_Var, v );                                                      \
         Load( g_FC_Var_yR, FC_Var_yR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_xy )                                                                                       \
      {                                                                                                     \
         Load( g_FC_Var_zL, FC_Var_zL, ID_FC_Var, v );                                                      \
         Load( g_FC_Var_zR, FC_Var_zR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
                                                                                                            \
                                                                                                            \
      /* load the face-ceneterd fluxes */                                                                   \
      if ( Inner_x )                                                                                        \
      {                                                                                                     \
         Load( g_FC_Flux_x, FC_Flux_xL, ID_In_xL, v );                                                      \
         Load( g_FC_Flux_x, FC_Flux_xR, ID_In_R,  v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_y )                                                                                        \
      {                                                                                                     \
         Load( g_FC_Flux_y, FC_Flux_yL, ID_In_yL, v );                                                      \
         Load( g_FC_Flux_y, FC_Flux_yR, ID_In_R,  v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_z )                                                                                        \
      {                                                                                                     \
         Load( g_FC_Flux_z, FC_Flux_zL, ID_In_zL, v );                                                      \
         Load( g_FC_Flux_z, FC_Flux_zR, ID_In_R,  v );                                                      \
      }                                                                                                     \
                                                                                                            \
                                                                                                            \
      /* compute the transverse gradient and correct the face-centered variables */                         \
      if ( Inner_yz )   Correct( FC_Var_xL, FC_Var_xR, FC_Flux_yL, FC_Flux_yR, FC_Flux_zL, FC_Flux_zR );    \
      if ( Inner_xz )   Correct( FC_Var_yL, FC_Var_yR, FC_Flux_xL, FC_Flux_xR, FC_Flux_zL, FC_Flux_zR );    \
      if ( Inner_xy )   Correct( FC_Var_zL, FC_Var_zR, FC_Flux_xL, FC_Flux_xR, FC_Flux_yL, FC_Flux_yR );    \
                                                                                                            \
                                                                                                            \
      /* store the corrected face-centered variables back to the global arrays */                           \
      if ( Inner_yz )                                                                                       \
      {                                                                                                     \
         Dump( FC_Var_xL, g_FC_Var_xL, ID_FC_Var, v );                                                      \
         Dump( FC_Var_xR, g_FC_Var_xR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_xz )                                                                                       \
      {                                                                                                     \
         Dump( FC_Var_yL, g_FC_Var_yL, ID_FC_Var, v );                                                      \
         Dump( FC_Var_yR, g_FC_Var_yR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_xy )                                                                                       \
      {                                                                                                     \
         Dump( FC_Var_zL, g_FC_Var_zL, ID_FC_Var, v );                                                      \
         Dump( FC_Var_zR, g_FC_Var_zR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
   } // Correct_1v


   ID_FC_Var = tx;

// loop over all cells
   while ( ID_FC_Var < N_FC_VAR*N_FC_VAR*N_FC_VAR )
   {
      ID_In_R  = ID_FC_Var;
      ID_In_xL = ID_In_R - dID_In.x;
      ID_In_yL = ID_In_R - dID_In.y;
      ID_In_zL = ID_In_R - dID_In.z;

      ID3d.x   = ID_FC_Var%N_FC_VAR;
      ID3d.y   = ID_FC_Var%(N_FC_VAR*N_FC_VAR)/N_FC_VAR;
      ID3d.z   = ID_FC_Var/(N_FC_VAR*N_FC_VAR);

      Inner_x  = ( ID3d.x != 0U  &&  ID3d.x != N_FC_VAR-1 );
      Inner_y  = ( ID3d.y != 0U  &&  ID3d.y != N_FC_VAR-1 );
      Inner_z  = ( ID3d.z != 0U  &&  ID3d.z != N_FC_VAR-1 );
      Inner_xy = Inner_x && Inner_y;
      Inner_yz = Inner_y && Inner_z;
      Inner_xz = Inner_x && Inner_z;


      Correct_1v( 0 );
      Correct_1v( 1 );
      Correct_1v( 2 );
      Correct_1v( 3 );
      Correct_1v( 4 );


      ID_FC_Var += dID_FC_Var;

   } // while ( ID_FC_Var < N_FC_VAR*N_FC_VAR*N_FC_VAR )

#  undef Load
#  undef Dump
#  undef Correct
#  undef Correct_1v

} // FUNCTION : CUFLU_TGradient_Correction



#endif // #if ( defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == CTU )
