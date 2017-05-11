
#include "DAINO.h"

#if ( MODEL == HYDRO  &&  defined GRAVITY )




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetTimeStep_Gravity
// Description :  Estimate the evolution time-step and physical time interval by gravity
//
// Note        :  Physical coordinates : dTime == dt
//                Comoving coordinates : dTime == dt*(Hubble parameter)*(scale factor)^3 == delta(scale factor)
// 
// Parameter   :  dt       : Time interval to advance solution
//                dTime    : Time interval to update physical time 
//                MinDtLv  : Refinement level determining the smallest time-step
//                MinDtVar : Maximum acceleration determining the minimum time-step
//                dt_dTime : dt/dTime (== 1.0 if COMOVING is off)
//
// Return      :  dt, dTime, MinDtLv, MinDtVar
//-------------------------------------------------------------------------------------------------------
void Hydro_GetTimeStep_Gravity( double &dt, double &dTime, int &MinDtLv, real &MinDtVar, const double dt_dTime )
{

   real  *MaxAcc   = MinDtInfo_Gravity;   // "MinDtInfo_Gravity" is a global variable
   double dt_local = __FLT_MAX__;         // initialize it as an extremely large number
   double dt_min, dt_tmp;


// get the maximum gravitational acceleration
   if ( !OPT__ADAPTIVE_DT )   Hydro_GetMaxAcc( MaxAcc );


// get the time-step in one rank
   for (int lv=0; lv<NLEVEL; lv++)
   {
#     ifdef INDIVIDUAL_TIMESTEP
      dt_tmp = sqrt( 2.0*patch->dh[lv]/MaxAcc[lv] ) * (1<<lv);
#     else
      dt_tmp = sqrt( 2.0*patch->dh[lv]/MaxAcc[lv] );
#     endif

      if ( dt_tmp <= 0.0 )
         Aux_Error( ERROR_INFO, "dt_tmp = %14.7e <= 0.0 (Rank %d, Lv %d) !!\n", dt_tmp, MPI_Rank, lv );
      
      if ( dt_tmp < dt_local )    
      {
         dt_local = dt_tmp;
         MinDtLv  = lv;
         MinDtVar = MaxAcc[lv];
      }
   }


// get the minimum time-step from all ranks
   MPI_Allreduce( &dt_local, &dt_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );


// verify the minimum time-step
   if ( dt_min == __FLT_MAX__ )
      Aux_Error( ERROR_INFO, "time-step estimation by gravity is incorrect (dt_min = %13.7e) !!\n", dt_min );


// gather the minimum time-step information from all ranks
#  ifndef SERIAL
   double *dt_AllRank     = new double [MPI_NRank];
   int    *MinDtLv_AllRank  = new int    [MPI_NRank];
   real   *MinDtVar_AllRank = new real   [MPI_NRank];
   
   MPI_Gather( &dt_local, 1, MPI_DOUBLE, dt_AllRank,       1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Gather( &MinDtLv,  1, MPI_INT,    MinDtLv_AllRank,  1, MPI_INT,    0, MPI_COMM_WORLD );
#  ifdef FLOAT8
   MPI_Gather( &MinDtVar, 1, MPI_DOUBLE, MinDtVar_AllRank, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#  else
   MPI_Gather( &MinDtVar, 1, MPI_FLOAT,  MinDtVar_AllRank, 1, MPI_FLOAT,  0, MPI_COMM_WORLD );
#  endif

   if ( MPI_Rank == 0 )
   {
      for (int Rank=0; Rank<MPI_NRank; Rank++)
      {
         if ( dt_AllRank[Rank] == dt_min )
         {
            MinDtLv  = MinDtLv_AllRank [Rank];
            MinDtVar = MinDtVar_AllRank[Rank];
            break;
         }

         if ( Rank == MPI_NRank-1 )    Aux_Message( stderr, "WARNING : no match of \"dt_min\" was found !!\n" );
      }
   }

   delete [] dt_AllRank;
   delete [] MinDtLv_AllRank;
   delete [] MinDtVar_AllRank;
#  endif // #ifndef SERIAL 


// return "2*time-step" since at the base level each step actually includes two sub-steps
   dt    = DT__GRAVITY * 2.0 * dt_min;
   dTime = dt / dt_dTime;

} // FUNCTION : Hydro_GetTimeStep_Gravity



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetMaxAcc
// Description :  Evaluate the maximum gravitational acceleration for the time-step estimation
//
// Note        :  This function is also invoked in "Init_DAINO"
// 
// Parameter   :  MaxAcc   : Array to store the maximum gravitational acceleration at each level
//-------------------------------------------------------------------------------------------------------
void Hydro_GetMaxAcc( real MaxAcc[] )
{

   const bool IntPhase_No = false;
   const int NPG = 1;               // number of patch groups
   const int NP  = NPG*8;           // number of patches
#  ifdef OPENMP
   const int NT  = OMP_NTHREAD;     // number of OpenMP threads
#  else
   const int NT  = 1;
#  endif

   real (*Acc_Array)[GRA_NXT][GRA_NXT][GRA_NXT] = NULL;
   real *MaxAcc_OMP = new real [NT];
   real Acc[3], Coeff;
   int  TID = 0;  // thread ID


// loop over all levels   
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    initialize the arrays MaxAcc and MaxAcc_OMP as an extremely small number
      MaxAcc[lv] = __FLT_MIN__;
      for (int t=0; t<NT; t++)   MaxAcc_OMP[t] = __FLT_MIN__;

//    set the constant coefficient
      if ( OPT__GRA_P5_GRADIENT )   Coeff = 1.0/(12.0*patch->dh[lv]);
      else                          Coeff = 1.0/( 2.0*patch->dh[lv]);

#ifndef OOC

#     pragma omp parallel private( Acc_Array, Acc, TID )
      {
         Acc_Array = new real [NP][GRA_NXT][GRA_NXT][GRA_NXT];

#        ifdef OPENMP
         TID = omp_get_thread_num();
#        endif

//       loop over all patches
#        pragma omp for
         for (int PID0=0; PID0<patch->NPatchComma[lv][1]; PID0+=NP)
         {
//          prepare the potential data with ghost zone
            Prepare_PatchGroupData( lv, Time[lv], &Acc_Array[0][0][0][0], GRA_GHOST_SIZE, NPG, &PID0, _POTE,
                                    OPT__GRA_INT_SCHEME, UNIT_PATCH, NSIDE_06, IntPhase_No );


//          loop over eight patches within the same patch group
            for (int ID=0; ID<NP; ID++)
            {
               for (int k=GRA_GHOST_SIZE; k<GRA_NXT-GRA_GHOST_SIZE; k++)
               for (int j=GRA_GHOST_SIZE; j<GRA_NXT-GRA_GHOST_SIZE; j++)
               for (int i=GRA_GHOST_SIZE; i<GRA_NXT-GRA_GHOST_SIZE; i++)
               {
//                evaluate the gravitational acceleration
                  if ( OPT__GRA_P5_GRADIENT )
                  {
                     Acc[0] = Coeff * ( -     Acc_Array[ID][k  ][j  ][i+2] +     Acc_Array[ID][k  ][j  ][i-2]
                                        + 8.0*Acc_Array[ID][k  ][j  ][i+1] - 8.0*Acc_Array[ID][k  ][j  ][i-1] );
                     Acc[1] = Coeff * ( -     Acc_Array[ID][k  ][j+2][i  ] +     Acc_Array[ID][k  ][j-2][i  ]
                                        + 8.0*Acc_Array[ID][k  ][j+1][i  ] - 8.0*Acc_Array[ID][k  ][j-1][i  ] );
                     Acc[2] = Coeff * ( -     Acc_Array[ID][k+2][j  ][i  ] +     Acc_Array[ID][k-2][j  ][i  ]
                                        + 8.0*Acc_Array[ID][k+1][j  ][i  ] - 8.0*Acc_Array[ID][k-1][j  ][i  ] );
                  }

                  else
                  {
                     Acc[0] = Coeff * ( Acc_Array[ID][k  ][j  ][i+1] - Acc_Array[ID][k  ][j  ][i-1] );
                     Acc[1] = Coeff * ( Acc_Array[ID][k  ][j+1][i  ] - Acc_Array[ID][k  ][j-1][i  ] );
                     Acc[2] = Coeff * ( Acc_Array[ID][k+1][j  ][i  ] - Acc_Array[ID][k-1][j  ][i  ] );
                  }


//                record the maximum acceleration
                  for (int d=0; d<3; d++)
                  {
                     Acc[d]          = fabs( Acc[d] );
                     MaxAcc_OMP[TID] = ( Acc[d] > MaxAcc_OMP[TID] ) ? Acc[d] : MaxAcc_OMP[TID];
                  }

               } // i,j,k
            } // for (int ID=0; ID<NP; ID++)
         } // for (int PID0=0; PID0<patch->NPatchComma[lv][1]; PID0+=NP)

         delete [] Acc_Array;

      } // OpenMP parallel region


//    compare the maximum acceleration evaluated by different OMP threads
      for (int t=0; t<NT; t++)   MaxAcc[lv] = ( MaxAcc_OMP[t] > MaxAcc[lv] ) ? MaxAcc_OMP[t] : MaxAcc[lv];


#else // OOC

      Acc_Array = new real [NP][GRA_NXT][GRA_NXT][GRA_NXT];

      OOC_Mis_GetMaxAcc( lv, Coeff, Acc_Array, MaxAcc );

      delete [] Acc_Array;

#endif
   } // for (int lv=0; lv<NLEVEL; lv++)


   delete [] MaxAcc_OMP;

} // FUNCTION : Hydro_GetMaxAcc



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )
