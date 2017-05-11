
#include "DAINO.h"

#ifdef GRAVITY



#ifdef SERIAL
rfftwnd_plan     FFTW_Plan, FFTW_Plan_Inv;
#else
rfftwnd_mpi_plan FFTW_Plan, FFTW_Plan_Inv;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_FFTW
// Description :  Create the FFTW plans 
//-------------------------------------------------------------------------------------------------------
void Init_FFTW()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );


#  ifdef SERIAL
   FFTW_Plan     = rfftw3d_create_plan( NX0_TOT[2], NX0_TOT[1], NX0_TOT[0], FFTW_REAL_TO_COMPLEX, 
                                        FFTW_ESTIMATE | FFTW_IN_PLACE );

   FFTW_Plan_Inv = rfftw3d_create_plan( NX0_TOT[2], NX0_TOT[1], NX0_TOT[0], FFTW_COMPLEX_TO_REAL, 
                                        FFTW_ESTIMATE | FFTW_IN_PLACE );

#  else

   FFTW_Plan     = rfftw3d_mpi_create_plan( MPI_COMM_WORLD, NX0_TOT[2], NX0_TOT[1], NX0_TOT[0], 
                                            FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE );

   FFTW_Plan_Inv = rfftw3d_mpi_create_plan( MPI_COMM_WORLD, NX0_TOT[2], NX0_TOT[1], NX0_TOT[0], 
                                            FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" ); 

} // FUNCTION : Init_FFTW



//-------------------------------------------------------------------------------------------------------
// Function    :  End_FFTW
// Description :  Delete the FFTW plans 
//-------------------------------------------------------------------------------------------------------
void End_FFTW()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );

#  ifdef SERIAL
   rfftwnd_destroy_plan    ( FFTW_Plan );
   rfftwnd_destroy_plan    ( FFTW_Plan_Inv );
#  else
   rfftwnd_mpi_destroy_plan( FFTW_Plan );
   rfftwnd_mpi_destroy_plan( FFTW_Plan_Inv );
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : End_FFTW



#endif // #ifdef GRAVITY
