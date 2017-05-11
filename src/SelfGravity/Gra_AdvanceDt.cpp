
#include "DAINO.h"

#ifdef GRAVITY



extern Timer_t *Timer_Gra_Advance[NLEVEL];
extern Timer_t *Timer_GetBuf     [NLEVEL][4];




//-------------------------------------------------------------------------------------------------------
// Function    :  Gra_AdvanceDt
// Description :  Solve the Poisson equation and advance the fluid variables by the gravitational acceleration
//
// Note        :  a. Poisson solver : lv = 0 : invoke the function "CPU_PoissonSolver_FFT"
//                                    lv > 0 : invoke the function "InvokeSolver"
// Note        :  b. Gravity solver : invoke the function "InvokeSolver"
//             :  c. The updated potential and fluid variables will be stored in the same sandglass 
//
// Parameter   :  lv             : Targeted refinement level 
//                PrepTime       : Targeted physical time to prepare the coarse-grid data
//                dt             : Time interval to advance solution
//                SaveSg         : Sandglass to store the updated data 
//                GraAcc         : true/false --> invoke the (Poisson+Gravity)/(Poisson only) solver(s)
//                OverlapMPI     : true --> Overlap MPI time with CPU/GPU computation
//                Overlap_Sync   : true  --> Advance the patches which cannot be overlapped with MPI communication
//                                 false --> Advance the patches which can    be overlapped with MPI communication
//                                 (useful only if "OverlapMPI == true")
//-------------------------------------------------------------------------------------------------------
void Gra_AdvanceDt( const int lv, const double PrepTime, const double dt, const int SaveSg, const bool GraAcc,
                    const bool OverlapMPI, const bool Overlap_Sync )
{

// coefficient in front of the RHS in the Poisson eq.
#  ifdef COMOVING
   const real Poi_Coeff = 4.0*M_PI*NEWTON_G*Time[lv];
#  else
   const real Poi_Coeff = 4.0*M_PI*NEWTON_G;
#  endif


// the base-level Poisson solver is implemented by using the FFTW library (with CPUs only)
   if ( lv == 0 )    
   {
#     ifdef OOC
                     CPU_PoissonSolver_FFT( Poi_Coeff, SaveSg );
#     else
      TIMING_FUNC(   CPU_PoissonSolver_FFT( Poi_Coeff, SaveSg ),   Timer_Gra_Advance[lv],   false   );   
#     endif

      patch    ->PotSg[0] = SaveSg;
#     ifdef OOC
      OOC_patch->PotSg[0] = patch->PotSg[0];
#     endif

#     ifndef OOC
      TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, SaveSg, POT_FOR_POISSON, _POTE, Pot_ParaBuf, USELB_YES ),
                     Timer_GetBuf[lv][1],   true   );
#     endif

      if ( GraAcc )  
      {
#ifndef OOC

//       TIMING_FUNC(   InvokeSolver( GRAVITY_SOLVER, lv, PrepTime, dt, NULL_REAL, SaveSg, OverlapMPI, 
//                                    Overlap_Sync ),
         TIMING_FUNC(   InvokeSolver( GRAVITY_SOLVER, lv, PrepTime, dt, NULL_REAL, SaveSg, false, false ),
                        Timer_Gra_Advance[lv],   true   );

         TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg, NULL_INT, DATA_GENERAL, _FLU, Flu_ParaBuf, USELB_YES ),
                        Timer_GetBuf[lv][2],   true   );

#else // OOC

         OOC_Gra_AdvanceDt_InvokeSolver( lv, PrepTime, dt, SaveSg );

#endif 

         patch    ->FluSg[0] = SaveSg;
#        ifdef OOC
         OOC_patch->FluSg[0] = patch->FluSg[0];
#        endif
      } // if ( GraAcc )
   } // if ( lv == 0 )

   else // lv > 0 
   {
      if ( GraAcc )  
         InvokeSolver( POISSON_AND_GRAVITY_SOLVER, lv, PrepTime, dt,        Poi_Coeff, SaveSg,
                       OverlapMPI, Overlap_Sync );
      else           
         InvokeSolver( POISSON_SOLVER,             lv, PrepTime, NULL_REAL, Poi_Coeff, SaveSg,
                       OverlapMPI, Overlap_Sync );
   }

} // FUNCTION : Gra_AdvanceDt



#endif // #ifdef GRAVITY
