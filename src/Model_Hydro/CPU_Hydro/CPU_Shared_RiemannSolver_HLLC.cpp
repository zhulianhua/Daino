
#include "DAINO.h"
#include "CUFLU.h"

#if (  !defined GPU  &&  MODEL == HYDRO  &&  \
       ( RSOLVER == HLLC || CHECK_INTERMEDIATE == HLLC )  &&  \
       ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )



extern void CPU_Rotate3D( real InOut[], const int XYZ, const bool Forward );
extern void CPU_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma );




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_RiemannSolver_HLLC
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer.
//                The wave speed is estimated by the same formula in HLLE solver
//
// Note        :  1. The input data should be conserved variables 
//                2. Ref : a. Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro 
//                         b. Batten, P., Clarke, N., Lambert, C., & Causon, D. M. 1997, SIAM J. Sci. Comput.,
//                            18, 1553
//                3. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  XYZ      : Targeted spatial direction : (0/1/2) --> (x/y/z)  
//                Flux_Out : Array to store the output flux
//                L_In     : Input left  state (conserved variables)
//                R_In     : Input right state (conserved variables)
//                Gamma    : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
void CPU_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], 
                             const real Gamma )
{

// 1. reorder the input variables for different spatial directions
   real L[5], R[5];

   for (int v=0; v<5; v++)    
   {
      L[v] = L_In[v];
      R[v] = R_In[v];
   }

   CPU_Rotate3D( L, XYZ, true );
   CPU_Rotate3D( R, XYZ, true );


// 2. evaluate the Roe's average values
   const real Gamma_m1 = Gamma - (real)1.0;
   real _RhoL, _RhoR, P_L, P_R, H_L, H_R, u, v, w, V2, H, Cs; 
   real RhoL_sqrt, RhoR_sqrt, _RhoL_sqrt, _RhoR_sqrt, _RhoLR_sqrt_sum, TempP_Rho;
   
   _RhoL = (real)1.0 / L[0];
   _RhoR = (real)1.0 / R[0];
   P_L   = Gamma_m1*(  L[4] - (real)0.5*( L[1]*L[1] + L[2]*L[2] + L[3]*L[3] )*_RhoL  );
   P_R   = Gamma_m1*(  R[4] - (real)0.5*( R[1]*R[1] + R[2]*R[2] + R[3]*R[3] )*_RhoR  );
#  ifdef ENFORCE_POSITIVE
   P_L   = FMAX( P_L, MIN_VALUE );
   P_R   = FMAX( P_R, MIN_VALUE );
#  endif
   H_L   = ( L[4] + P_L )*_RhoL;  
   H_R   = ( R[4] + P_R )*_RhoR;  

   RhoL_sqrt       = SQRT( L[0] );
   RhoR_sqrt       = SQRT( R[0] );
   _RhoL_sqrt      = (real)1.0 / RhoL_sqrt;
   _RhoR_sqrt      = (real)1.0 / RhoR_sqrt;
   _RhoLR_sqrt_sum = (real)1.0 / (RhoL_sqrt + RhoR_sqrt); 

   u  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[1] + _RhoR_sqrt*R[1] );
   v  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[2] + _RhoR_sqrt*R[2] );
   w  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[3] + _RhoR_sqrt*R[3] );
   V2 = u*u + v*v + w*w;
   H  = _RhoLR_sqrt_sum*(  RhoL_sqrt*H_L  +  RhoR_sqrt*H_R  );

   TempP_Rho = H - (real)0.5*V2;
#  ifdef ENFORCE_POSITIVE
   TempP_Rho = FMAX( TempP_Rho, MIN_VALUE );
#  endif
   Cs        = SQRT( Gamma_m1*TempP_Rho );


// 3. estimate the maximum wave speeds
   const real EVal[NCOMP] = { u-Cs, u, u, u, u+Cs };
   real u_L, u_R, Cs_L, Cs_R, W_L, W_R, MaxV_L, MaxV_R;

   u_L    = _RhoL*L[1];
   u_R    = _RhoR*R[1];
   Cs_L   = SQRT( Gamma*P_L*_RhoL );
   Cs_R   = SQRT( Gamma*P_R*_RhoR );
   W_L    = FMIN( EVal[      0], u_L-Cs_L );
   W_R    = FMAX( EVal[NCOMP-1], u_R+Cs_R );
   MaxV_L = FMIN( W_L, (real)0.0 );
   MaxV_R = FMAX( W_R, (real)0.0 );


// 4. evaluate the star-region velocity (V_S) and pressure (P_S)
   real V_S, P_S, temp1_L, temp1_R, temp2_L, temp2_R, temp3;

   temp1_L = L[0]*( u_L - W_L );
   temp1_R = R[0]*( u_R - W_R );
   temp2_L = P_L + temp1_L*u_L;
   temp2_R = P_R + temp1_R*u_R;
   temp3   = real(1.0) / ( temp1_L - temp1_R );

   V_S = temp3*( P_L - P_R + temp1_L*u_L - temp1_R*u_R );
   P_S = temp3*( temp1_L*temp2_R - temp1_R*temp2_L );

#  ifdef ENFORCE_POSITIVE
   P_S = FMAX( P_S, MIN_VALUE );
#  endif


// 5. evaluate the weightings of the left(right) fluxes and contact wave
   real Flux_LR[5], temp4, Coeff_LR, Coeff_S;

   if ( V_S >= (real)0.0 )
   {
      CPU_Con2Flux( 0, Flux_LR, L, Gamma );

      for (int v=0; v<5; v++)    Flux_LR[v] -= MaxV_L*L[v];    // fluxes along the maximum wave speed

      temp4    = (real)1.0 / ( V_S - MaxV_L );
      Coeff_LR = temp4*V_S;
      Coeff_S  = -temp4*MaxV_L*P_S;
   }

   else // V_S < 0.0
   {
      CPU_Con2Flux( 0, Flux_LR, R, Gamma );

      for (int v=0; v<5; v++)    Flux_LR[v] -= MaxV_R*R[v];    // fluxes along the maximum wave speed

      temp4    = (real)1.0 / ( V_S - MaxV_R );
      Coeff_LR = temp4*V_S;
      Coeff_S  = -temp4*MaxV_R*P_S;
   }


// 6. evaluate the HLLC fluxes
   for (int v=0; v<5; v++)    Flux_Out[v] = Coeff_LR*Flux_LR[v];

   Flux_Out[1] += Coeff_S;
   Flux_Out[4] += Coeff_S*V_S;


// 7. restore the correct order
   CPU_Rotate3D( Flux_Out, XYZ, false );

} // FUNCTION : CPU_RiemannSolver_HLLC



#endif // #if ( !GPU && HYDRO && ( RSOLVER == HLLC || CHECK_INTER == HLLC ) && ( FLU_SCHEME==MHM/MHM_RP/CTU ) )
