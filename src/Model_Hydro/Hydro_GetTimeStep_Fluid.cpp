
#include "DAINO.h"
#include "CUFLU.h"
#ifdef INTEL
#include <mathimf.h>
#endif

#if ( MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetTimeStep_Fluid
// Description :  Estimate the evolution time-step and physical time interval from the fluid solver condition
//
// Note        :  1. Physical coordinates : dTime == dt
//                   Comoving coordinates : dTime == dt*(Hubble parameter)*(scale factor)^3 == delta(scale factor)
//                2. time-step is estimated by the stability criterion from the von Neumann stability analysis
// 
// Parameter   :  dt       : Time interval to advance solution
//                dTime    : Time interval to update physical time 
//                MinDtLv  : Refinement level determining the smallest time-step
//                MinDtVar : Array to store the variables with the maximum speed (minimum time-step) at each level
//                dt_dTime : dt/dTime (== 1.0 if COMOVING is off)
//
// Return      :  dt, dTime, MinDtLv, MinDtVar
//-------------------------------------------------------------------------------------------------------
void Hydro_GetTimeStep_Fluid( double &dt, double &dTime, int &MinDtLv, real MinDtVar[], const double dt_dTime )
{

   real  *MaxCFL   = MinDtInfo_Fluid;  // "MinDtInfo_Fluid" is a global variable
   double dt_local = __FLT_MAX__;      // initialize it as an extremely large number
   double dt_min, dt_tmp; 
   real   MinDtVar_AllLv[NLEVEL][NCOMP];


// get the maximum CFL velocity ( sound speed + fluid velocity )
   if ( !OPT__ADAPTIVE_DT )   Hydro_GetMaxCFL( MaxCFL, MinDtVar_AllLv );


// get the time-step at the base level per sub-step in one rank
   for (int lv=0; lv<NLEVEL; lv++)
   {
#     ifdef INDIVIDUAL_TIMESTEP
      dt_tmp = (double)patch->dh[lv] / (double)MaxCFL[lv] * (double)(1<<lv);
#     else
      dt_tmp = (double)patch->dh[lv] / (double)MaxCFL[lv];
#     endif

      if ( dt_tmp <= 0.0 )
         Aux_Error( ERROR_INFO, "dt_tmp = %14.7e <= 0 (Rank %d, Lv %d) !!\n", dt_tmp, MPI_Rank, lv );
      
      if ( dt_tmp < dt_local )    
      {
         dt_local = dt_tmp;
         MinDtLv  = lv;

         for (int v=0; v<NCOMP; v++)   MinDtVar[v] = MinDtVar_AllLv[lv][v];
      }
   } // for (int lv=0; lv<NLEVEL; lv++)


// get the minimum time-step from all ranks
   MPI_Allreduce( &dt_local, &dt_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );


// verify the minimum time-step
   if ( dt_min == __FLT_MAX__ )
      Aux_Error( ERROR_INFO, "time-step estimation by %s is incorrect (dt_min = %13.7e) !!\n",
                 __FUNCTION__, dt_min );


// gather the minimum time-step information from all ranks
#  ifndef SERIAL
   double *dt_AllRank                = new double [MPI_NRank];
   int    *MinDtLv_AllRank           = new int    [MPI_NRank];
   real   (*MinDtVar_AllRank)[NCOMP] = new real   [MPI_NRank][NCOMP];
   
   MPI_Gather( &dt_local, 1,     MPI_DOUBLE, dt_AllRank,       1,     MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Gather( &MinDtLv,  1,     MPI_INT,    MinDtLv_AllRank,  1,     MPI_INT,    0, MPI_COMM_WORLD );
#  ifdef FLOAT8
   MPI_Gather( MinDtVar,  NCOMP, MPI_DOUBLE, MinDtVar_AllRank, NCOMP, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#  else
   MPI_Gather( MinDtVar,  NCOMP, MPI_FLOAT,  MinDtVar_AllRank, NCOMP, MPI_FLOAT,  0, MPI_COMM_WORLD );
#  endif

   if ( MPI_Rank == 0 )
   {
      for (int Rank=0; Rank<MPI_NRank; Rank++)
      {
         if ( dt_AllRank[Rank] == dt_min )
         {
            MinDtLv = MinDtLv_AllRank[Rank];
            for (int v=0; v<NCOMP; v++)   MinDtVar[v] = MinDtVar_AllRank[Rank][v];
            break;
         }

         if ( Rank == MPI_NRank-1 )    Aux_Message( stderr, "WARNING : no match of \"dt_min\" was found !!\n" );
      }
   }

   delete [] dt_AllRank;
   delete [] MinDtVar_AllRank;
   delete [] MinDtLv_AllRank;
#  endif // #ifndef SERIAL 


// return "2*time-step" since at the base level each step actually includes two sub-steps
   dt    = DT__FLUID * 2.0 * dt_min;
   dTime = dt / dt_dTime;

} // FUNCTION : Hydro_GetTimeStep_Fluid



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetMaxCFL
// Description :  Evaluate the maximum (fluid velocity + sound speed) for the time-step estimation
//
// Note        :  This function is also invoked in "Init_DAINO"
// 
// Parameter   :  MaxCFL         : Array to store the maximum speed at each level
//                MinDtVar_AllLv : Array to store the variables with the maximum speed at each level
//-------------------------------------------------------------------------------------------------------
void Hydro_GetMaxCFL( real MaxCFL[], real MinDtVar_AllLv[][NCOMP] )
{

   const real Gamma_m1 = GAMMA - 1.0;
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

   real  *MaxCFL_OMP           = new real [NT];
   real (*MinDtVar_OMP)[NCOMP] = new real [NT][NCOMP];
   real Rho, _Rho, Vx, Vy, Vz, Egy, Ek, Pres, Cs, MaxV, MaxCFL_candidate;
   int  TID = 0;  // thread ID


// loop over all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    initialize the arrays MaxCFL and MaxCFL_OMP as extremely small numbers
      MaxCFL[lv] = __FLT_MIN__;
      for (int t=0; t<NT; t++)   MaxCFL_OMP[t] = __FLT_MIN__;

#ifndef OOC

#     pragma omp parallel private( Rho, _Rho, Vx, Vy, Vz, Egy, Ek, Pres, Cs, MaxV, MaxCFL_candidate, TID )
      {
#        ifdef OPENMP
         TID = omp_get_thread_num();
#        endif

//       loop over all patches
#        pragma omp for
         for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
         {
//          if ( patch->ptr[0][lv][PID]->son == -1 )  // only check the patches without son
            if ( true )
            {
               for (int k=0; k<PATCH_SIZE; k++)
               for (int j=0; j<PATCH_SIZE; j++)
               for (int i=0; i<PATCH_SIZE; i++)
               {
                  Rho =       patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
                 _Rho = 1.0 / Rho;
                  Vx  = fabs( patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[MOMX][k][j][i] ) * _Rho;
                  Vy  = fabs( patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[MOMY][k][j][i] ) * _Rho;
                  Vz  = fabs( patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[MOMZ][k][j][i] ) * _Rho;
                  Egy =       patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[ENGY][k][j][i];

                  Ek   = 0.5*Rho*( Vx*Vx + Vy*Vy + Vz*Vz );
                  Pres = Gamma_m1*( Egy - Ek );
#                 ifdef ENFORCE_POSITIVE
//                Pres = fmax( Pres, MIN_VALUE );  // replaced by the following line for detecting NaN
                  Pres = ( Pres < MIN_VALUE ) ? MIN_VALUE : Pres;
#                 endif               
                  Cs   = sqrt( GAMMA*Pres*_Rho );

#                 if   ( FLU_SCHEME == RTVD  ||  FLU_SCHEME == CTU  ||  FLU_SCHEME == WAF )
                  MaxV             = ( Vx > Vy   ) ? Vx : Vy;
                  MaxV             = ( Vz > MaxV ) ? Vz : MaxV;
                  MaxCFL_candidate = MaxV + Cs;

#                 elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
                  MaxV             = Vx + Vy + Vz;
                  MaxCFL_candidate = MaxV + 3.0*Cs;

#                 else
                  MaxV             = ( Vx > Vy   ) ? Vx : Vy;
                  MaxV             = ( Vz > MaxV ) ? Vz : MaxV;
                  MaxCFL_candidate = MaxV + Cs;
#                 endif

                  if ( MaxCFL_candidate > MaxCFL_OMP[TID] )
                  {
                     MaxCFL_OMP[TID] = MaxCFL_candidate;

                     MinDtVar_OMP[TID][0] = Rho;
                     MinDtVar_OMP[TID][1] = Vx;
                     MinDtVar_OMP[TID][2] = Vy;
                     MinDtVar_OMP[TID][3] = Vz;
                     MinDtVar_OMP[TID][4] = Cs;
                  }


//                verify the CFL number
                  if (  ! isfinite( MaxCFL_candidate )  )
                  {
                     Aux_Message( stderr, "ERROR : \"incorrect CFL evaluation (%14.7e)\" !!\n", MaxCFL_candidate);
                     Aux_Message( stderr, "        t %14.7e, Step %ld, Cs %14.7e, MaxV %14.7e, Rho %14.7e\n",
                                  Time[0], Step, Cs, MaxV, Rho );
                     Aux_Message( stderr, "        Rank %d, Lv %d, PID  %d, Coordinates (%5d,%5d,%5d)\n",
                                  DAINO_RANK, lv, PID, patch->ptr[0][lv][PID]->corner[0] + i*patch->scale[lv],  
                                                       patch->ptr[0][lv][PID]->corner[1] + j*patch->scale[lv], 
                                                       patch->ptr[0][lv][PID]->corner[2] + k*patch->scale[lv]  );
                     Aux_Message( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  
                                                                                           __FUNCTION__  );
                     Aux_Message( stderr, "        **********************************************************\n");
                     Aux_Message( stderr, "        ** Usually this error is caused by the negative density **\n");
                     Aux_Message( stderr, "        ** or pressure evaluated in the fluid solver            **\n");
                     Aux_Message( stderr, "        **********************************************************\n");
                     MPI_Exit();
                  }
   
               } // i,j,k
//          } // if ( patch->ptr[0][lv][PID]->son == -1 )
            } // if ( true )
         } // for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
      } // OpenMP parallel region


//    compare the maximum CFL evaluated by different OMP threads
      for (int t=0; t<NT; t++)
      {
         if ( MaxCFL_OMP[t] > MaxCFL[lv] )   
         {
            MaxCFL[lv] = MaxCFL_OMP[t];

            MinDtVar_AllLv[lv][0] = MinDtVar_OMP[t][0];
            MinDtVar_AllLv[lv][1] = MinDtVar_OMP[t][1];
            MinDtVar_AllLv[lv][2] = MinDtVar_OMP[t][2];
            MinDtVar_AllLv[lv][3] = MinDtVar_OMP[t][3];
            MinDtVar_AllLv[lv][4] = MinDtVar_OMP[t][4];
         }
      }


#else // OOC

      OOC_Mis_GetMaxCFL( lv, MaxCFL, MinDtVar_AllLv[lv] );

#endif
   } // for (int lv=0; lv<NLEVEL; lv++)


   delete [] MaxCFL_OMP;
   delete [] MinDtVar_OMP;

} // FUNCTION : Hydro_GetMaxCFL



#endif // #if ( MODEL == HYDRO )
