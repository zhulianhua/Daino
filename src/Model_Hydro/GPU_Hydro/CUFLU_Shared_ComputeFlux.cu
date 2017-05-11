
#ifndef __CUFLU_COMPUTEFLUX_CU__
#define __CUFLU_COMPUTEFLUX_CU__



#include "Macro.h"
#include "CUFLU.h"
#if   ( RSOLVER == EXACT )
#include "CUFLU_Shared_RiemannSolver_Exact.cu"
#elif ( RSOLVER == ROE )
#include "CUFLU_Shared_RiemannSolver_Roe.cu"
#elif ( RSOLVER == HLLE )
#include "CUFLU_Shared_RiemannSolver_HLLE.cu"
#elif ( RSOLVER == HLLC )
#include "CUFLU_Shared_RiemannSolver_HLLC.cu"
#endif

static __device__ void CUFLU_ComputeFlux( const real g_FC_Var_xL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                          const real g_FC_Var_xR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                          const real g_FC_Var_yL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                          const real g_FC_Var_yR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          const real g_FC_Var_zL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                          const real g_FC_Var_zR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Flux_x[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                          real g_FC_Flux_y[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                          real g_FC_Flux_z[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                          real g_Flux[][9][5][ PS2*PS2 ], const bool DumpFlux,
                                          const uint Gap, const real Gamma );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_ComputeFlux
// Description :  Compute the face-centered fluxes by Riemann solver 
//
// Note        :  1. Currently support the exact, Roe, HLLE, and HLLC solvers
//                2. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                3. The function is asynchronous
//                   --> "__syncthreads" must be called before using the output data
//                4. The sizes of the array g_FC_Var_x/y/z are assumed to be "N_FC_VAR" 
//                   --> "N_FC_VAR-1" fluxes will be computed along each normal direction 
//                5. The sizes of the arrays g_FC_Flux_XX in each direction are assumed to be "N_FL_FLUX" 
//                   --> The (i,j,k) flux will be stored in the array g_FC_Flux_XX with 
//                       the index "(k*N_FL_FLUX+j)*N_FL_FLUX+i"
//                   --> The (i,j,k) FC_Flux_x is defined at the +x surfaces of the cell (i,     j-Gap, k-Gap)
//                       The (i,j,k) FC_Flux_y is defined at the +y surfaces of the cell (i-Gap, j,     k-Gap)
//                       The (i,j,k) FC_Flux_z is defined at the +z surfaces of the cell (i-Gap, j-Gap, k    )
//                6. This function is shared by MHM, MHM_RP, and CTU schemes
//                7. For the performance consideration, this function will also be responsible for store the
//                   inter-patch fluxes 
//                   --> Option "DumpFlux"
//                8. The "__forceinline__" qualifier is added for higher performance
//
// Parameter   :  g_FC_Var_xL : Global memory array storing the face-centered variables on the -x surface
//                g_FC_Var_xR : Global memory array storing the face-centered variables on the +x surface
//                g_FC_Var_yL : Global memory array storing the face-centered variables on the -y surface
//                g_FC_Var_yR : Global memory array storing the face-centered variables on the +y surface
//                g_FC_Var_zL : Global memory array storing the face-centered variables on the -z surface
//                g_FC_Var_zR : Global memory array storing the face-centered variables on the +z surface
//                g_FC_Flux_x : Global memory array to store the face-centered fluxes in the x direction
//                g_FC_Flux_y : Global memory array to store the face-centered fluxes in the y direction
//                g_FC_Flux_z : Global memory array to store the face-centered fluxes in the z direction
//                g_Flux      : Global memory array to store the output fluxes
//                DumpFlux    : true --> store the inter-patch fluxes to "g_Flux" for the AMR correction
//                Gap         : Number of grids to be skipped in the transverse direction
//                              --> "(N_FC_VAR-2*Gap)^2" fluxes will be computed in each surface
//                Gamma       : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
__forceinline__
__device__ void CUFLU_ComputeFlux( const real g_FC_Var_xL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                   const real g_FC_Var_xR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                   const real g_FC_Var_yL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                   const real g_FC_Var_yR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                   const real g_FC_Var_zL[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                   const real g_FC_Var_zR[][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                   real g_FC_Flux_x[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                   real g_FC_Flux_y[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                   real g_FC_Flux_z[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                   real g_Flux[][9][5][ PS2*PS2 ], const bool DumpFlux,
                                   const uint Gap, const real Gamma )
                                   
{

   const uint  bx       = blockIdx.x;
   const uint  tx       = threadIdx.x;
   const uint  dID      = blockDim.x;
   const uint3 dID_In   = make_uint3( 1, N_FC_VAR, N_FC_VAR*N_FC_VAR );
   const uint  N_N      = N_FC_VAR - 1;      // number of fluxes to be calculated in the normal direction
   const uint  N_T      = N_FC_VAR - 2*Gap;  // number of fluxes to be calculated in the transverse direction

#  if ( RSOLVER == EXACT )
   const real  Gamma_m1 = Gamma - (real)1.0;
#  endif

   uint   ID, ID_In, ID_Out, Nxy, SurfID;
   uint3  ID3d;
   FluVar VarL, VarR, FC_Flux;

#  if ( RSOLVER == EXACT )
   FluVar *Useless = NULL;
#  endif


#  define Load( Input, Output, ID )    \
   {                                   \
      Output.Rho = Input[bx][0][ID];   \
      Output.Px  = Input[bx][1][ID];   \
      Output.Py  = Input[bx][2][ID];   \
      Output.Pz  = Input[bx][3][ID];   \
      Output.Egy = Input[bx][4][ID];   \
   } // Load

#  define Dump_FC_Flux( Input, Output, ID )     \
   {                                            \
      Output[bx][0][ID] = Input.Rho;            \
      Output[bx][1][ID] = Input.Px;             \
      Output[bx][2][ID] = Input.Py;             \
      Output[bx][3][ID] = Input.Pz;             \
      Output[bx][4][ID] = Input.Egy;            \
   } // Dump_FC_Flux

#  define Dump_InterPatch_Flux( FC_Flux, SurfID, ID_Flux )     \
   {                                                           \
      g_Flux[bx][SurfID][0][ID_Flux] = FC_Flux.Rho;            \
      g_Flux[bx][SurfID][1][ID_Flux] = FC_Flux.Px;             \
      g_Flux[bx][SurfID][2][ID_Flux] = FC_Flux.Py;             \
      g_Flux[bx][SurfID][3][ID_Flux] = FC_Flux.Pz;             \
      g_Flux[bx][SurfID][4][ID_Flux] = FC_Flux.Egy;            \
   } // Dump_InterPatch_Flux

#  if   ( RSOLVER == EXACT )

      /* exact solver */
      #define RiemannSolver( Dir, VarL, VarR )                                                           \
      {                                                                                                  \
         VarL = CUFLU_Con2Pri( VarL, Gamma_m1 );                                                         \
         VarR = CUFLU_Con2Pri( VarR, Gamma_m1 );                                                         \
                                                                                                         \
         FC_Flux = CUFLU_RiemannSolver_Exact( Dir, *Useless, *Useless, *Useless, VarL, VarR, Gamma );    \
      } // RiemannSolver

#  elif ( RSOLVER == ROE )

      /* Roe solver */
      #define RiemannSolver( Dir, VarL, VarR )                                                           \
      {                                                                                                  \
         FC_Flux = CUFLU_RiemannSolver_Roe( Dir, VarL, VarR, Gamma );                                    \
      } // RiemannSolver

#  elif ( RSOLVER == HLLE )

      /* HLLE solver */
      #define RiemannSolver( Dir, VarL, VarR )                                                           \
      {                                                                                                  \
         FC_Flux = CUFLU_RiemannSolver_HLLE( Dir, VarL, VarR, Gamma );                                   \
      } // RiemannSolver

#  elif ( RSOLVER == HLLC )

      /* HLLC solver */
      #define RiemannSolver( Dir, VarL, VarR )                                                           \
      {                                                                                                  \
         FC_Flux = CUFLU_RiemannSolver_HLLC( Dir, VarL, VarR, Gamma );                                   \
      } // RiemannSolver

#  else

#     error : ERROR : unsupported Riemann solver (EXACT/ROE/HLLE/HLLC) !!

#  endif

#  define GetFlux( Dir, Nx, Ny, Gap_x, Gap_y, Gap_z, dID_In, g_FC_Var_L, g_FC_Var_R, g_FC_Flux )               \
   {                                                                                                           \
      ID  = tx;                                                                                                \
      Nxy = (Nx)*(Ny);                                                                                         \
                                                                                                               \
      while ( ID < N_N*N_T*N_T )                                                                               \
      {                                                                                                        \
         ID3d.x = ID%(Nx);                                                                                     \
         ID3d.y = ID%Nxy/(Nx);                                                                                 \
         ID3d.z = ID/Nxy;                                                                                      \
         ID_In  = __umul24( __umul24( ID3d.z+Gap_z, N_FC_VAR  ) + ID3d.y+Gap_y, N_FC_VAR  ) + ID3d.x+Gap_x;    \
         ID_Out = __umul24( __umul24( ID3d.z,       N_FL_FLUX ) + ID3d.y,       N_FL_FLUX ) + ID3d.x;          \
                                                                                                               \
         Load( g_FC_Var_R, VarL, ID_In        );                                                               \
         Load( g_FC_Var_L, VarR, ID_In+dID_In );                                                               \
                                                                                                               \
         RiemannSolver( Dir, VarL, VarR );                                                                     \
                                                                                                               \
         Dump_FC_Flux( FC_Flux, g_FC_Flux, ID_Out );                                                           \
                                                                                                               \
         /* store the inter-patch fluxes */                                                                    \
         if ( DumpFlux )                                                                                       \
         {                                                                                                     \
            if      (  Dir == 0  &&  ( ID3d.x == 0 || ID3d.x == PS1 || ID3d.x == PS2 )  )                      \
            {                                                                                                  \
               SurfID = 3*Dir + ID3d.x/PS1;                                                                    \
               Dump_InterPatch_Flux( FC_Flux, SurfID, __umul24( ID3d.z, PS2 ) + ID3d.y );                      \
            }                                                                                                  \
                                                                                                               \
            else if (  Dir == 1  &&  ( ID3d.y == 0 || ID3d.y == PS1 || ID3d.y == PS2 )  )                      \
            {                                                                                                  \
               SurfID = 3*Dir + ID3d.y/PS1;                                                                    \
               Dump_InterPatch_Flux( FC_Flux, SurfID, __umul24( ID3d.z, PS2 ) + ID3d.x );                      \
            }                                                                                                  \
                                                                                                               \
            else if (  Dir == 2  &&  ( ID3d.z == 0 || ID3d.z == PS1 || ID3d.z == PS2 )  )                      \
            {                                                                                                  \
               SurfID = 3*Dir + ID3d.z/PS1;                                                                    \
               Dump_InterPatch_Flux( FC_Flux, SurfID, __umul24( ID3d.y, PS2 ) + ID3d.x );                      \
            }                                                                                                  \
         }                                                                                                     \
                                                                                                               \
         ID += dID;                                                                                            \
                                                                                                               \
      } /* while ( ID < N_N*N_T*N_T ) */                                                                       \
   } // GetFlux

   GetFlux( 0, N_N, N_T,   0, Gap, Gap, dID_In.x, g_FC_Var_xL, g_FC_Var_xR, g_FC_Flux_x );
   GetFlux( 1, N_T, N_N, Gap,   0, Gap, dID_In.y, g_FC_Var_yL, g_FC_Var_yR, g_FC_Flux_y );
   GetFlux( 2, N_T, N_T, Gap, Gap,   0, dID_In.z, g_FC_Var_zL, g_FC_Var_zR, g_FC_Flux_z );

#  undef Load
#  undef Dump_FC_Flux
#  undef Dump_InterPatch_Flux
#  undef RiemannSolver
#  undef GetFlux

} // FUNCTION : CUFLU_ComputeFlux



#endif // #ifndef __CUFLU_COMPUTEFLUX_CU__
