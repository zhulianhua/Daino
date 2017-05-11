
#include "DAINO.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_GetAverageDensity
// Description :  Evaluate the average density for the Poisson solver
//
// Note        :  1. The average density will be subtracted from the total density at each cell when 
//                   solving the Poisson equation at refined level in order to be consistent with the 
//                   base-level FFT solver
//                2. For the out-of-core computing, we must record the average density of each OOC rank
//                   in "ooc.AveRho" in advance
//                3. For the debug mode, we perform summation in a specific order in order to ensure that
//                   the round-off errors will be the same in runs with different numbers of MPI ranks
//-------------------------------------------------------------------------------------------------------
void Poi_GetAverageDensity()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... \n", __FUNCTION__ );


// check
#  ifndef DENS
#  error : ERROR : VARIABLE "DENS" IS NOT DEFINED IN THE FUNCTION "Poi_GetAverageDensity" !!
#  endif

   if ( !OPT__INIT_RESTRICT  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : option \"%s\" is NOT turned on when evaluating the average density !!\n",
                   "OPT__INIT_RESTRICT" );

// initialize it to zero
   AveDensity = 0.0;


// 1. for OOC computing
// ==================================================================================================
#  ifdef OOC

   double AveDensity_local = 0.0;

// sum over all OOC ranks
   for (int r=0; r<ooc.NRank; r++)     AveDensity_local += ooc.AveRho[r];

// sum over all MPI ranks
   MPI_Allreduce( &AveDensity_local, &AveDensity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

// average
   AveDensity /= double( DAINO_NRANK );


#  else // #ifdef OOC


// 2. for debug mode
// ==================================================================================================
#  ifdef DAINO_DEBUG

// only use rank 0 for summation in the debug mode
   const int NP[3]  = { NX0_TOT[0]/PATCH_SIZE, NX0_TOT[1]/PATCH_SIZE, NX0_TOT[2]/PATCH_SIZE };
   const int PScale = PATCH_SIZE * patch->scale[0];

   int     Cr3D[3], NPatch_All[MPI_NRank], Disp[MPI_NRank];
   int     NPatch_Sum    = 0;
   double *Rho_All       = NULL;
   long   *Cr1D_All      = NULL;
   int    *Cr1D_IdxTable = NULL;
   double *Rho_Local     = new double [ patch->NPatchComma[0][1] ];
   long   *Cr1D_Local    = new long   [ patch->NPatchComma[0][1] ];

   MPI_Gather( &patch->NPatchComma[0][1], 1, MPI_INT, NPatch_All, 1, MPI_INT, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      for (int r=0; r<MPI_NRank; r++)  NPatch_Sum += NPatch_All[r];

      Disp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  Disp[r] = Disp[r-1] + NPatch_All[r-1];

      Rho_All       = new double [NPatch_Sum];
      Cr1D_All      = new long   [NPatch_Sum];
      Cr1D_IdxTable = new int    [NPatch_Sum];
   }

// prepare density and 1D-corner arrays
   for (int PID=0; PID<patch->NPatchComma[0][1]; PID++)
   {
      for (int d=0; d<3; d++)    Cr3D[d] = patch->ptr[0][0][PID]->corner[d]/PScale;

      Cr1D_Local[PID] = ( Cr3D[2]*NP[1] + Cr3D[1] )*NP[0] + Cr3D[0];
      Rho_Local [PID] = 0.0;

      for (int k=0; k<PATCH_SIZE; k++)
      for (int j=0; j<PATCH_SIZE; j++)
      for (int i=0; i<PATCH_SIZE; i++)
         Rho_Local[PID] += patch->ptr[ patch->FluSg[0] ][0][PID]->fluid[DENS][k][j][i];
   }

// gather data
   MPI_Gatherv( Cr1D_Local, patch->NPatchComma[0][1], MPI_LONG,   Cr1D_All, NPatch_All, Disp, MPI_LONG,
                0, MPI_COMM_WORLD );

   MPI_Gatherv( Rho_Local,  patch->NPatchComma[0][1], MPI_DOUBLE, Rho_All,  NPatch_All, Disp, MPI_DOUBLE,
                0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
//    sort
      Mis_Heapsort( NPatch_Sum, Cr1D_All, Cr1D_IdxTable );
   
//    get averaged density
      for (int t=0; t<NPatch_Sum; t++)    AveDensity += Rho_All[ Cr1D_IdxTable[t] ];
      AveDensity /= double( NX0_TOT[0]*NX0_TOT[1]*NX0_TOT[2] );
   }

// broadcast
   MPI_Bcast( &AveDensity, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

// free memory
   delete [] Rho_Local;
   delete [] Cr1D_Local;
   if ( MPI_Rank == 0 )
   {
      delete [] Rho_All;
      delete [] Cr1D_All;
      delete [] Cr1D_IdxTable;
   }


// 3. for general cases
// ==================================================================================================
#  else // #ifdef DAINO_DEBUG

// evaluate the sum of density (we only use the base-level data because the restriction condition is assumed 
// to be fulfilled
   double AveDensity_local = 0.0;

   for (int PID=0; PID<patch->NPatchComma[0][1]; PID++)
   for (int k=0; k<PATCH_SIZE; k++)
   for (int j=0; j<PATCH_SIZE; j++)
   for (int i=0; i<PATCH_SIZE; i++)
      AveDensity_local += patch->ptr[ patch->FluSg[0] ][0][PID]->fluid[DENS][k][j][i];

// sum over all MPI ranks
   MPI_Allreduce( &AveDensity_local, &AveDensity, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

// average
   AveDensity /= double( NX0_TOT[0]*NX0_TOT[1]*NX0_TOT[2] );

#  endif // #ifdef DAINO_DEBUG ... else ...
#  endif // #ifdef OOC ... else ...


// 4. output results
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "NOTE : background density = %20.14e\n", AveDensity );

#     ifdef COMOVING
      const double Deviation = fabs( AveDensity - 1.0 );
      if ( Deviation > 1.0e-5 )
      {
         Aux_Message( stderr, "WARNING : background density deviates from unity by %20.14e ", Deviation );
         Aux_Message( stderr,           "(UNITY is assumed in COMOVING) !!\n" );
      }
#     endif
   }


// check
   if ( AveDensity <= 0.0 )   Aux_Error( ERROR_INFO, "average density (%14.7e) <= 0.0 !!\n", AveDensity );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Poi_GetAverageDensity



#endif // #ifdef GRAVITY
