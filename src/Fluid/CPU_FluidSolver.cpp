
#ifndef GPU



#include "DAINO.h"

#if   ( MODEL == HYDRO )
#if   ( FLU_SCHEME == RTVD )
void CPU_FluidSolver_RTVD( real Flu_Array_In [][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                           real Flu_Array_Out[][5][ PS2*PS2*PS2 ], 
                           real Flux_Array[][9][5][ PS2*PS2 ], 
                           const int NPatchGroup, const real dt, const real dh, const real Gamma,
                           const bool StoreFlux, const bool XYZ );
#elif ( FLU_SCHEME == WAF )
void CPU_FluidSolver_WAF( real Flu_Array_In [][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                          real Flu_Array_Out[][5][ PS2*PS2*PS2 ], 
                          real Flux_Array[][9][5][ PS2*PS2 ], 
                          const int NPatchGroup, const real dt, const real dh, const real Gamma, 
                          const bool StoreFlux, const bool XYZ, const WAF_Limiter_t WAF_Limiter );
#elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
void CPU_FluidSolver_MHM( const real Flu_Array_In[][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                          real Flu_Array_Out[][5][ PS2*PS2*PS2 ], 
                          real Flux_Array[][9][5][ PS2*PS2 ], 
                          const int NPatchGroup, const real dt, const real dh, const real Gamma, 
                          const bool StoreFlux, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, 
                          const real EP_Coeff );
#elif ( FLU_SCHEME == CTU )
void CPU_FluidSolver_CTU( const real Flu_Array_In[][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                          real Flu_Array_Out[][5][ PS2*PS2*PS2 ], 
                          real Flux_Array[][9][5][ PS2*PS2 ], 
                          const int NPatchGroup, const real dt, const real dh, const real Gamma, 
                          const bool StoreFlux, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, 
                          const real EP_Coeff );
#endif // FLU_SCHEME

#elif ( MODEL == MHD )
#warning : WAIT MHD !!!

#elif ( MODEL == ELBDM )
void CPU_ELBDMSolver( real Flu_Array_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                      real Flu_Array_Out[][FLU_NOUT][ PS2*PS2*PS2 ], 
                      const int NPatchGroup, const real dt, const real dh, const real Eta, const bool XYZ );

#else 
#error : ERROR : unsupported MODEL !!
#endif // MODEL




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_FluidSolver
// Description :  1. MODEL == HYDRO : use CPU to solve the Euler equations by different schemes
//                2. MODEL == ELBDM : use CPU to solve the kinematic operator in the Schrodinger's equation
//
// Note        :  Currently five hydro schemes are supported in HYDRO : 
//                   1. Relaxing TVD scheme                            (RTVD  ) -->   split
//                   2. Weighted-Average-Flux scheme                   (WAF   ) -->   split
//                   3. MUSCL-Hancock scheme                           (MHM   ) --> unsplit
//                   4. MUSCL-Hancock scheme with Riemann prediction   (MHM_RP) --> unsplit
//                   5. Corner-Transport-Upwind scheme                 (CTU   ) --> unsplit
//
//
// Parameter   :  h_Flu_Array_In    : Host array storing the input variables
//                h_Flu_Array_Out   : Host array to store the output variables
//                h_Flux_Array      : Host array to store the output fluxes
//                h_MinDtInfo_Array : Host array to store the minimum time-step information in each patch group
//                                    --> useful only if "GetMinDtInfo == true"
//                                    --> NOT supported yet
//                NPatchGroup       : Number of patch groups to be evaluated
//                dt                : Time interval to advance solution
//                dh                : Grid size
//                Gamma             : Ratio of specific heats
//                StoreFlux         : true --> store the coarse-fine fluxes
//                XYZ               : true   : x->y->z ( forward sweep)
//                                    false1 : z->y->x (backward sweep)
//                                    --> only useful for the RTVD and WAF schemes
//                LR_Limiter        : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                    (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                   vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff      : Coefficient of the generalized MinMod limiter
//                EP_Coeff          : Coefficient of the extrema-preserving limiter
//                WAF_Limiter       : Flux limiter for the WAF scheme
//                                    (0/1/2/3) = (SuperBee/vanLeer/vanAlbada/MinBee)
//                Eta               : Particle mass / Planck constant
//                GetMinDtInfo      : true --> Gather the minimum time-step information (the CFL condition in 
//                                             HYDRO) in each patch group
//                                         --> NOT supported yet
//
// Useless parameters in HYDRO : Eta
// Useless parameters in ELBDM : h_Flux_Array, Gamma, StoreFlux, LR_Limiter, MinMod_Coeff, EP_Coeff, WAF_Limiter
//-------------------------------------------------------------------------------------------------------
void CPU_FluidSolver( real h_Flu_Array_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                      real h_Flu_Array_Out[][FLU_NOUT][ PS2*PS2*PS2 ], 
                      real h_Flux_Array[][9][NCOMP   ][ PS2*PS2 ], 
                      real h_MinDtInfo_Array[],
                      const int NPatchGroup, const real dt, const real dh, const real Gamma, const bool StoreFlux,
                      const bool XYZ, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const real EP_Coeff,
                      const WAF_Limiter_t WAF_Limiter, const real Eta, const bool GetMinDtInfo )
{

#  if   ( MODEL == HYDRO )

#     if   ( FLU_SCHEME == RTVD )

      CPU_FluidSolver_RTVD( h_Flu_Array_In, h_Flu_Array_Out, h_Flux_Array, NPatchGroup, dt, dh, Gamma, StoreFlux,
                            XYZ );

#     elif ( FLU_SCHEME == WAF )

      CPU_FluidSolver_WAF ( h_Flu_Array_In, h_Flu_Array_Out, h_Flux_Array, NPatchGroup, dt, dh, Gamma, StoreFlux,
                            XYZ, WAF_Limiter );

#     elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )

      CPU_FluidSolver_MHM ( h_Flu_Array_In, h_Flu_Array_Out, h_Flux_Array, NPatchGroup, dt, dh, Gamma, StoreFlux,
                            LR_Limiter, MinMod_Coeff, EP_Coeff );

#     elif ( FLU_SCHEME == CTU )

      CPU_FluidSolver_CTU ( h_Flu_Array_In, h_Flu_Array_Out, h_Flux_Array, NPatchGroup, dt, dh, Gamma, StoreFlux,
                            LR_Limiter, MinMod_Coeff, EP_Coeff );

#     else

#     error : unsupported CPU hydro scheme

#     endif


#  elif ( MODEL == MHD )
#     warning : WAIT MHD !!!


#  elif ( MODEL == ELBDM )
      CPU_ELBDMSolver( h_Flu_Array_In, h_Flu_Array_Out, NPatchGroup, dt, dh, Eta, XYZ );

#  else 
#     error : ERROR : unsupported MODEL !!
#  endif // MODEL

} // FUNCTION : CPU_FluidSolver



#endif // #ifndef GPU
