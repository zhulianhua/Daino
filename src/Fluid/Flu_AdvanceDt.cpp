
#include "DAINO.h"
#include <climits>




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_AdvanceDt
// Description :  Advance the fluid attributes by the flux gradients
//
// Note        :  a. Invoke the function "InvokeSolver"
//                b. Currently the updated data can only be stored in the different sandglass from the 
//                   input fluid data
//
// Parameter   :  lv             : Targeted refinement level 
//                PrepTime       : Targeted physical time to prepare the coarse-grid data
//                dt             : Time interval to advance solution
//                SaveSg         : Sandglass to store the updated data 
//                OverlapMPI     : true --> Overlap MPI time with CPU/GPU computation
//                Overlap_Sync   : true  --> Advance the patches which cannot be overlapped with MPI communication
//                                 false --> Advance the patches which can    be overlapped with MPI communication
//                                 (useful only if "OverlapMPI == true")
//-------------------------------------------------------------------------------------------------------
void Flu_AdvanceDt( const int lv, const double PrepTime, const double dt, const int SaveSg,
                    const bool OverlapMPI, const bool Overlap_Sync )
{

   InvokeSolver( FLUID_SOLVER, lv, PrepTime, dt, NULL_REAL, SaveSg, OverlapMPI, Overlap_Sync );

   if ( OPT__FIXUP_FLUX )  Buf_ResetBufferFlux( lv );

} // FUNCTION : Flu_AdvanceDt
