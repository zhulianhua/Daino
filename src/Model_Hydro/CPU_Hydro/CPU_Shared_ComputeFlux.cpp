
#include "DAINO.h"
#include "CUFLU.h"

#if ( !defined GPU  &&  MODEL == HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU) )



#if   ( RSOLVER == EXACT )
extern void CPU_Con2Pri( const real In[], real Out[], const real  Gamma_m1 );
extern void CPU_RiemannSolver_Exact( const int XYZ, real eival_out[], real L_star_out[], real R_star_out[], 
                                     real Flux_Out[], const real L_In[], const real R_In[], const real Gamma ); 
#elif ( RSOLVER == ROE )
extern void CPU_RiemannSolver_Roe( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], 
                                   const real Gamma );
#elif ( RSOLVER == HLLE )
extern void CPU_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], 
                                    const real Gamma );
#elif ( RSOLVER == HLLC )
extern void CPU_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], 
                                    const real Gamma );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_ComputeFlux
// Description :  Compute the face-centered fluxes by Riemann solver 
//
// Note        :  1. Currently support the exact and Roe solvers
//                2. The size of the input array "FC_Var" is assumed to be N_FC_VAR^3
//                   --> "N_FC_VAR-1" fluxes will be computed along each direction 
//
// Parameter   :  FC_Var   : Array storing the input face-centered conserved variables
//                FC_Flux  : Array to store the output face-centered flux
//                NFlux    : Size of the array FC_Flux in each direction (must be >= N_FC_VAR-1)
//                           --> The (i,j,k) flux will be stored in the array FC_Flux with 
//                               the index "(k*NFlux+j)*NFlux+i"
//                           --> The (i,j,k) FC_Flux_x/y/z are defined at the +x/+y/+z surfaces of the 
//                               cell (i,j,k)
//                Gap      : Number of grids to be skipped in the transverse direction
//                           --> "(N_FC_VAR-2*Gap)^2" fluxes will be computed in each surface
//                Gamma    : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
void CPU_ComputeFlux( const real FC_Var[][6][5], real FC_Flux[][3][5], const int NFlux, const int Gap,
                      const real Gamma )
{

   const int dID2[3] = { 1, N_FC_VAR, N_FC_VAR*N_FC_VAR };

   const real *ConVar_L=NULL, *ConVar_R=NULL;
   int ID1, ID2, dL, dR, start2[3]={0}, end1[3]={0};

#  if ( RSOLVER == EXACT )
   const real Gamma_m1 = Gamma - (real)1.0;
   real PriVar_L[5], PriVar_R[5];
#  endif


// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
      dL = 2*d;
      dR = dL+1;

      switch ( d )
      {
         case 0 : start2[0] = 0;                start2[1] = Gap;              start2[2] = Gap;  
                  end1  [0] = N_FC_VAR-1;       end1  [1] = N_FC_VAR-2*Gap;   end1  [2] = N_FC_VAR-2*Gap;   break;

         case 1 : start2[0] = Gap;              start2[1] = 0;                start2[2] = Gap;  
                  end1  [0] = N_FC_VAR-2*Gap;   end1  [1] = N_FC_VAR-1;       end1  [2] = N_FC_VAR-2*Gap;   break;

         case 2 : start2[0] = Gap;              start2[1] = Gap;              start2[2] = 0;  
                  end1  [0] = N_FC_VAR-2*Gap;   end1  [1] = N_FC_VAR-2*Gap;   end1  [2] = N_FC_VAR-1;       break;
      }

      for (int k1=0, k2=start2[2];  k1<end1[2];  k1++, k2++)
      for (int j1=0, j2=start2[1];  j1<end1[1];  j1++, j2++)
      for (int i1=0, i2=start2[0];  i1<end1[0];  i1++, i2++)
      {
         ID1 = (k1*NFlux    + j1)*NFlux    + i1;
         ID2 = (k2*N_FC_VAR + j2)*N_FC_VAR + i2;

         ConVar_L = FC_Var[ ID2         ][dR];
         ConVar_R = FC_Var[ ID2+dID2[d] ][dL];

#        if   ( RSOLVER == EXACT )
         CPU_Con2Pri( ConVar_L, PriVar_L, Gamma_m1 );
         CPU_Con2Pri( ConVar_R, PriVar_R, Gamma_m1 );

         CPU_RiemannSolver_Exact( d, NULL, NULL, NULL, FC_Flux[ID1][d], PriVar_L, PriVar_R, Gamma );
#        elif ( RSOLVER == ROE )
         CPU_RiemannSolver_Roe ( d, FC_Flux[ID1][d], ConVar_L, ConVar_R, Gamma );
#        elif ( RSOLVER == HLLE )
         CPU_RiemannSolver_HLLE( d, FC_Flux[ID1][d], ConVar_L, ConVar_R, Gamma );
#        elif ( RSOLVER == HLLC )
         CPU_RiemannSolver_HLLC( d, FC_Flux[ID1][d], ConVar_L, ConVar_R, Gamma );
#        else
#        error : ERROR : unsupported Riemann solver (EXACT/ROE) !!
#        endif
      }

   } // for (int d=0; d<3; d++)

} // FUNCTION : CPU_ComputeFlux



#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  (FLU_SCHEME == MHM || MHM_RP || CTU) )
