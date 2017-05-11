
#include "DAINO.h"

#if ( defined GRAVITY  &&  !defined GPU )



// Poisson solver prototypes
#if   ( POT_SCHEME == SOR )
void CPU_PoissonSolver_SOR( const real Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT], 
                            const real Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT], 
                                  real Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT], 
                            const int NPatchGroup, const real dh, const int Min_Iter, const int Max_Iter, 
                            const real Omega, const real Poi_Coeff, const IntScheme_t IntScheme );

#elif ( POT_SCHEME == MG  )
void CPU_PoissonSolver_MG( const real Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT],
                           const real Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                 real Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                           const int NPatchGroup, const real dh_Min, const int Max_Iter, const int NPre_Smooth,
                           const int NPost_Smooth, const real Tolerated_Error, const real Poi_Coeff, 
                           const IntScheme_t IntScheme );
#endif // POT_SCHEME


// Gravity solver prototypes
#if   ( MODEL == HYDRO )
void CPU_HydroGravitySolver(       real Flu_Array[][GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE], 
                             const real Pot_Array[][GRA_NXT][GRA_NXT][GRA_NXT],
                             const int NPatchGroup, const real dt, const real dh, const bool P5_Gradient );

#elif ( MODEL == MHD )
#warning : WAIT MHD !!!

#elif ( MODEL == ELBDM )
void CPU_ELBDMGravitySolver(       real Flu_Array[][GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE],
                             const real Pot_Array[][GRA_NXT][GRA_NXT][GRA_NXT],
                             const int NPatchGroup, const real EtaDt );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_PoissonGravitySolver
// Description :  Invoke the CPU_PoissonSolver and/or CPU_GravitySolver to evaluate the potential and/or 
//                advance the fluid variables by the gravitational acceleration for a group of patches
//
// Parameter   :  h_Rho_Array          : Host array to store the input density 
//                h_Pot_Array_In       : Host array to store the input "coarse-grid" potential for interpolation
//                h_Pot_Array_Out      : Host array to store the output potential
//                h_Flu_Array          : Host array to store the fluid variables for the Gravity solver
//                NPatchGroup          : Number of patch groups evaluated simultaneously by GPU 
//                dt                   : Time interval to advance solution
//                dh                   : Grid size
//                SOR_Min_Iter         : Minimum number of iterations for SOR
//                SOR_Max_Iter         : Maximum number of iterations for SOR
//                SOR_Omega            : Over-relaxation parameter
//                MG_Max_Iter          : Maximum number of iterations for multigrid
//                MG_NPre_Smooth       : Number of pre-smoothing steps for multigrid
//                MG_NPos_tSmooth      : Number of post-smoothing steps for multigrid
//                MG_Tolerated_Error   : Maximum tolerated error for multigrid
//                Poi_Coeff            : Coefficient in front of the RHS in the Poisson eq.
//                IntScheme            : Interpolation scheme for potential
//                                       --> currently supported schemes include
//                                           INT_CENTRAL : central interpolation
//                                           INT_CQUAD   : conservative quadratic interpolation 
//                                           INT_QUAD    : quadratic interpolation 
//                P5_Gradient          : Use 5-points stencil to evaluate the potential gradient
//                Eta                  : Particle mass / Planck constant
//                Poisson              : true --> invoke the Poisson solver
//                GraAcc               : true --> invoke the Gravity solver
//
// Useless parameters in HYDRO : Eta
// Useless parameters in ELBDM : P5_Gradient
//-------------------------------------------------------------------------------------------------------
void CPU_PoissonGravitySolver( const real h_Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT], 
                                     real h_Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                     real h_Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                                     real h_Flu_Array    [][GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE], 
                               const int NPatchGroup, const real dt, const real dh, const int SOR_Min_Iter, 
                               const int SOR_Max_Iter, const real SOR_Omega, const int MG_Max_Iter,
                               const int MG_NPre_Smooth, const int MG_NPost_Smooth, const real MG_Tolerated_Error,
                               const real Poi_Coeff, const IntScheme_t IntScheme, const bool P5_Gradient,
                               const real Eta, const bool Poisson, const bool GraAcc )
{

// check
   if ( Poisson  &&  ( IntScheme != INT_CENTRAL  &&  IntScheme != INT_CQUAD  &&  IntScheme != INT_QUAD )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IntScheme", IntScheme );


// Poisson solver
   if ( Poisson )
   {
#     if   ( POT_SCHEME == SOR )

      CPU_PoissonSolver_SOR( h_Rho_Array, h_Pot_Array_In, h_Pot_Array_Out, NPatchGroup, dh, 
                             SOR_Min_Iter, SOR_Max_Iter, SOR_Omega, 
                             Poi_Coeff, IntScheme );

#     elif ( POT_SCHEME == MG  )

      CPU_PoissonSolver_MG ( h_Rho_Array, h_Pot_Array_In, h_Pot_Array_Out, NPatchGroup, dh, 
                             MG_Max_Iter, MG_NPre_Smooth, MG_NPost_Smooth, MG_Tolerated_Error, 
                             Poi_Coeff, IntScheme );

#     else

#     error : ERROR : unsupported CPU Poisson solver !!

#     endif
   } // if ( Poisson )


// Gravity solver
   if ( GraAcc )
   {
#     if   ( MODEL == HYDRO )
      CPU_HydroGravitySolver( h_Flu_Array, h_Pot_Array_Out, NPatchGroup, dt, dh, P5_Gradient );

#     elif ( MODEL == MHD )
#     error : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      CPU_ELBDMGravitySolver( h_Flu_Array, h_Pot_Array_Out, NPatchGroup, Eta*dt );

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL
   } // if ( GraAcc )

} // FUNCTION : CPU_PoissonGravitySolver



#endif // #if ( defined GRAVITY  &&  !defined GPU )
