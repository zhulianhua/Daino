
#include "DAINO.h"

static void Preparation_Step( const Solver_t TSolver, const int lv, const double PrepTime, const int NPG, 
                              const int *PID0_List, const int ArrayID );
static void Solver( const Solver_t TSolver, const int lv, const int NPG, const int ArrayID, const double dt,
                    const real Poi_Coeff );
static void Closing_Step( const Solver_t TSolver, const int lv, const int SaveSg, const int NPG,
                          const int *PID0_List, const int ArrayID );

extern Timer_t *Timer_Pre         [NLEVEL][4];
extern Timer_t *Timer_Sol         [NLEVEL][4];
extern Timer_t *Timer_Clo         [NLEVEL][4];
#ifdef GRAVITY
extern Timer_t *Timer_Poi_PreRho  [NLEVEL];
extern Timer_t *Timer_Poi_PreFlu  [NLEVEL];
extern Timer_t *Timer_Poi_PrePot_C[NLEVEL];
extern Timer_t *Timer_Poi_PrePot_F[NLEVEL];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  InvokeSolver
// Description :  Invoke the GPU (or CPU) solvers and enable the concurrent execution between CPU and GPU
//
// Note        :  a. Use the input parameter "TSolver" to control the targeted solver
//                b. Each solver involves three steps
//                   --> 1. preparation step : prepare the input data 
//                       2. execution   step : invoke the solvers --> advance solutions or evaluate potentials
//                       3. closing     step : store the updated data
//                c. Currently the fluid solver can only store the updated data in the different sandglass from 
//                   the input data
//                d. For LOAD_BALANCE, one can turn on the option "LB_INPUT__OVERLAP_MPI" to enable the 
//                   overlapping between MPI communication and CPU/GPU computation
//
// Parameter   :  TSolver        : Targeted solver
//                                 --> FLUID_SOLVER               : Fluid / ELBDM solver
//                                     POISSON_SOLVER             : Poisson solver
//                                     GRAVITY_SOLVER             : Gravity solver
//                                     POISSON_AND_GRAVITY_SOLVER : Poisson + Gravity solvers
//                lv             : Targeted refinement level 
//                PrepTime       : Targeted physical time to prepare the coarse-grid data
//                dt             : Time interval to advance solution (for the fluid and gravity solvers)
//                Poi_Coeff      : Coefficient in front of the RHS in the Poisson eq.
//                SaveSg         : Sandglass to store the updated data 
//                OverlapMPI     : true --> Overlap MPI time with CPU/GPU computation
//                Overlap_Sync   : true  --> Advance the patches which cannot be overlapped with MPI communication
//                                 false --> Advance the patches which can    be overlapped with MPI communication
//                                 (useful only if "OverlapMPI == true")
//-------------------------------------------------------------------------------------------------------
void InvokeSolver( const Solver_t TSolver, const int lv, const double PrepTime, const double dt,
                   const real Poi_Coeff, const int SaveSg, const bool OverlapMPI, const bool Overlap_Sync )
{

// currently the fluid solver can only store the updated data in the different sandglass from the input data
   if ( TSolver == FLUID_SOLVER  &&  SaveSg == patch->FluSg[lv] )
      Aux_Error( ERROR_INFO, "SaveSg (%d) == patch->FluSg (%d) in the fluid solver at level %d !!\n", 
                 SaveSg, patch->FluSg[lv], lv );

#  ifndef GRAVITY
   if ( TSolver != FLUID_SOLVER )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );
#  endif


// maximum number of patch groups to be updated at a time
#  ifdef GRAVITY
   const int NPG_Max = ( TSolver == FLUID_SOLVER ) ? FLU_GPU_NPGROUP : POT_GPU_NPGROUP;
#  else
   const int NPG_Max = FLU_GPU_NPGROUP;
#  endif

   int *PID0_List    = NULL;  // list recording the patch indicies with LocalID==0 to be udpated
   bool AllocateList = false; // whether to allocate PID0_List or not
   int  ArrayID      = 0;     // array index to load and store data ( 0 or 1 )
   int  NPG[2];               // number of patch groups to be updated at a time 
   int  NTotal;               // total number of patch groups to be updated 
   int  Disp;                 // index displacement in PID0_List
 
   if ( OverlapMPI )
   {
#     ifdef LOAD_BALANCE
      if ( TSolver == FLUID_SOLVER )  
      {
         if ( Overlap_Sync )
         {
            NTotal    = patch->LB->OverlapMPI_FluSyncN   [lv];
            PID0_List = patch->LB->OverlapMPI_FluSyncPID0[lv]; 
         }
         else
         {
            NTotal    = patch->LB->OverlapMPI_FluAsyncN   [lv];
            PID0_List = patch->LB->OverlapMPI_FluAsyncPID0[lv]; 
         }
      }

#     ifdef GRAVITY
      else              
      {
         if ( Overlap_Sync )
         {
            NTotal    = patch->LB->OverlapMPI_PotSyncN   [lv];
            PID0_List = patch->LB->OverlapMPI_PotSyncPID0[lv]; 
         }
         else
         {
            NTotal    = patch->LB->OverlapMPI_PotAsyncN   [lv];
            PID0_List = patch->LB->OverlapMPI_PotAsyncPID0[lv]; 
         }
      }
#     endif 

#     else // #ifdef LOAD_BALANCE ... else ...
      Aux_Error( ERROR_INFO, "MPI overlapping is NOT supported if LOAD_BALANCE is off !!\n" );
#     endif // #ifdef LOAD_BALANCE ... else ...

   } // if ( OverlapMPI )

   else
   {
      AllocateList = true;
      NTotal       = patch->NPatchComma[lv][1] / 8;
      PID0_List    = new int [NTotal];

      for (int t=0; t<NTotal; t++)  PID0_List[t] = 8*t;
   } // if ( OverlapMPI ) ... else ...

   NPG[ArrayID] = ( NPG_Max < NTotal ) ? NPG_Max : NTotal;


//-------------------------------------------------------------------------------------------------------------
   TIMING_SYNC(   Preparation_Step( TSolver, lv, PrepTime, NPG[ArrayID], PID0_List, ArrayID ),
                  Timer_Pre[lv][TSolver]  );
//-------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------
   TIMING_SYNC(   Solver( TSolver, lv, NPG[ArrayID], ArrayID, dt, Poi_Coeff ),
                  Timer_Sol[lv][TSolver]  );
//-------------------------------------------------------------------------------------------------------------


   for (Disp=NPG_Max; Disp<NTotal; Disp+=NPG_Max)
   {

      ArrayID      = 1 - ArrayID;
      NPG[ArrayID] = ( NPG_Max < NTotal-Disp ) ? NPG_Max : NTotal-Disp;


//-------------------------------------------------------------------------------------------------------------
      TIMING_SYNC(   Preparation_Step( TSolver, lv, PrepTime, NPG[ArrayID], PID0_List+Disp, ArrayID ),  
                     Timer_Pre[lv][TSolver]  );
//-------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------
#     ifdef GPU
      CUAPI_Synchronize();
#     endif
//-------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------
      TIMING_SYNC(   Solver( TSolver, lv, NPG[ArrayID], ArrayID, dt, Poi_Coeff ), 
                     Timer_Sol[lv][TSolver]  );
//-------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------
      TIMING_SYNC(   Closing_Step( TSolver, lv, SaveSg, NPG[1-ArrayID], PID0_List+Disp-NPG_Max, 1-ArrayID ),
                     Timer_Clo[lv][TSolver]  ); 
//-------------------------------------------------------------------------------------------------------------

   } // for (int Disp=NPG_Max; Disp<NTotal; Disp+=NPG_Max)


//-------------------------------------------------------------------------------------------------------------
#  ifdef GPU
   CUAPI_Synchronize();
#  endif
//-------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------
   TIMING_SYNC(   Closing_Step( TSolver, lv, SaveSg, NPG[ArrayID], PID0_List+Disp-NPG_Max, ArrayID ), 
                  Timer_Clo[lv][TSolver]  ); 
//-------------------------------------------------------------------------------------------------------------
     

   if ( AllocateList )  delete [] PID0_List;

} // FUNCTION : InvokeSolver



//-------------------------------------------------------------------------------------------------------
// Function    :  Preparation_Step
// Description :  Prepare the input data for CPU/GPU solvers 
//
// Note        :  Use the input parameter "TSolver" to control the targeted solver 
//
// Parameter   :  TSolver     : Targeted solver
//                              --> FLUID_SOLVER               : Fluid / ELBDM solver
//                                  POISSON_SOLVER             : Poisson solver
//                                  GRAVITY_SOLVER             : Gravity solver
//                                  POISSON_AND_GRAVITY_SOLVER : Poisson + Gravity solvers
//                lv          : Targeted refinement level 
//                PrepTime    : Targeted physical time to prepare the coarse-grid data
//                NPG         : Number of patch groups to be prepared at a time
//                PID0_List   : List recording the patch indicies with LocalID==0 to be udpated
//                ArrayID     : Array index to load and store data ( 0 or 1 )
//-------------------------------------------------------------------------------------------------------
void Preparation_Step( const Solver_t TSolver, const int lv, const double PrepTime, const int NPG, 
                       const int *PID0_List, const int ArrayID )
{

   switch ( TSolver )
   {
      case FLUID_SOLVER :  
         Flu_Prepare( lv, PrepTime, h_Flu_Array_F_In[ArrayID], NPG, PID0_List );
         break;

#     ifdef GRAVITY
      case POISSON_SOLVER :  
         TIMING_SYNC(   Poi_Prepare_Rho( lv, PrepTime, h_Rho_Array_P    [ArrayID], NPG, PID0_List ), 
                        Timer_Poi_PreRho[lv]   );

         TIMING_SYNC(   Poi_Prepare_Pot( lv, PrepTime, h_Pot_Array_P_In [ArrayID], NPG, PID0_List ),
                        Timer_Poi_PrePot_C[lv]   );
         break;

      case GRAVITY_SOLVER :  
         TIMING_SYNC(   Gra_Prepare_Flu( lv,           h_Flu_Array_G    [ArrayID], NPG, PID0_List ),
                        Timer_Poi_PreFlu[lv]   );

         TIMING_SYNC(   Gra_Prepare_Pot( lv, PrepTime, h_Pot_Array_P_Out[ArrayID], NPG, PID0_List ),
                        Timer_Poi_PrePot_F[lv]   );
         break;

      case POISSON_AND_GRAVITY_SOLVER :  
         TIMING_SYNC(   Poi_Prepare_Rho( lv, PrepTime, h_Rho_Array_P    [ArrayID], NPG, PID0_List ), 
                        Timer_Poi_PreRho[lv]   );

         TIMING_SYNC(   Poi_Prepare_Pot( lv, PrepTime, h_Pot_Array_P_In [ArrayID], NPG, PID0_List ),
                        Timer_Poi_PrePot_C[lv]   );

         TIMING_SYNC(   Gra_Prepare_Flu( lv,           h_Flu_Array_G    [ArrayID], NPG, PID0_List ),
                        Timer_Poi_PreFlu[lv]   );
         break;
#     endif

      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );

   } // switch ( TSolver )

} // FUNCTION : Preparation_Step



//-------------------------------------------------------------------------------------------------------
// Function    :  Solver
// Description :  Invoke the CPU/GPU solvers
//
// Note        :  Use the input parameter "TSolver" to control the targeted solver 
//
// Parameter   :  TSolver     : Targeted solver
//                              --> FLUID_SOLVER               : Fluid / ELBDM solver
//                                  POISSON_SOLVER             : Poisson solver
//                                  GRAVITY_SOLVER             : Gravity solver
//                                  POISSON_AND_GRAVITY_SOLVER : Poisson + Gravity solvers
//                lv          : Targeted refinement level 
//                NPG         : Number of patch groups to be upcated at a time
//                ArrayID     : Array index to load and store data ( 0 or 1 )
//                dt          : Time interval to advance solution (for the fluid and gravity solvers)
//                Poi_Coeff   : Coefficient in front of the RHS in the Poisson eq.
//-------------------------------------------------------------------------------------------------------
void Solver( const Solver_t TSolver, const int lv, const int NPG, const int ArrayID, const double dt,
             const real Poi_Coeff )
{

   const real dh          = patch->dh[lv];                                   

#  ifdef GRAVITY
   const bool POISSON_ON  = true;
   const bool GRAVITY_ON  = true;
   const bool POISSON_OFF = false;
   const bool GRAVITY_OFF = false;
#  endif // #ifdef GRAVITY


// define useless variables in different models
#  if ( MODEL != ELBDM )
   const real ETA = NULL_REAL;
#  endif

#  if ( MODEL != HYDRO )
   const bool  Flu_XYZ                  = true;
   const real  GAMMA                    = NULL_REAL;
   const real  MINMOD_COEFF             = NULL_REAL;
   const real  EP_COEFF                 = NULL_REAL;
   const LR_Limiter_t  OPT__LR_LIMITER  = LR_LIMITER_NONE;
   const WAF_Limiter_t OPT__WAF_LIMITER = WAF_LIMITER_NONE;
#  else
   const bool Flu_XYZ                   = 1 - ( AdvanceCounter[lv]%2 );    // forward/backward sweep
#  endif

#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#  error : ERROR : ADD THE MODEL-DEPENDENT USELESS VARIABLES FOR THE NEW MODELS HERE
#  endif


   switch ( TSolver )
   {
      case FLUID_SOLVER :  

#        ifdef GPU
         CUAPI_Asyn_FluidSolver( h_Flu_Array_F_In[ArrayID], h_Flu_Array_F_Out[ArrayID], h_Flux_Array[ArrayID], 
                                 h_MinDtInfo_Fluid_Array[ArrayID], NPG, dt, dh, GAMMA, OPT__FIXUP_FLUX, Flu_XYZ, 
                                 OPT__LR_LIMITER, MINMOD_COEFF, EP_COEFF, OPT__WAF_LIMITER, ETA, 
                                 OPT__ADAPTIVE_DT, GPU_NSTREAM );
#        else
         CPU_FluidSolver       ( h_Flu_Array_F_In[ArrayID], h_Flu_Array_F_Out[ArrayID], h_Flux_Array[ArrayID], 
                                 h_MinDtInfo_Fluid_Array[ArrayID], NPG, dt, dh, GAMMA, OPT__FIXUP_FLUX, Flu_XYZ, 
                                 OPT__LR_LIMITER, MINMOD_COEFF, EP_COEFF, OPT__WAF_LIMITER, ETA, 
                                 OPT__ADAPTIVE_DT );
#        endif
         break;


#     ifdef GRAVITY

      case POISSON_SOLVER :  

#        ifdef GPU     
         CUAPI_Asyn_PoissonGravitySolver( h_Rho_Array_P[ArrayID], h_Pot_Array_P_In[ArrayID], 
                                          h_Pot_Array_P_Out[ArrayID], NULL, 
                                          NPG, dt, dh, SOR_MIN_ITER, SOR_MAX_ITER, 
                                          SOR_OMEGA, MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, 
                                          MG_TOLERATED_ERROR, Poi_Coeff, OPT__POT_INT_SCHEME, 
                                          NULL_BOOL, ETA, POISSON_ON, GRAVITY_OFF, GPU_NSTREAM ); 
#        else
         CPU_PoissonGravitySolver       ( h_Rho_Array_P[ArrayID], h_Pot_Array_P_In[ArrayID], 
                                          h_Pot_Array_P_Out[ArrayID], NULL, 
                                          NPG, dt, dh, SOR_MIN_ITER, SOR_MAX_ITER, 
                                          SOR_OMEGA, MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, 
                                          MG_TOLERATED_ERROR, Poi_Coeff, OPT__POT_INT_SCHEME, 
                                          NULL_BOOL, ETA, POISSON_ON, GRAVITY_OFF ); 
#        endif
         break;


      case GRAVITY_SOLVER :  
              
#        ifdef GPU     
         CUAPI_Asyn_PoissonGravitySolver( NULL, NULL, 
                                          h_Pot_Array_P_Out[ArrayID], h_Flu_Array_G[ArrayID], 
                                          NPG, dt, dh, NULL_INT, NULL_INT, 
                                          NULL_REAL, NULL_INT, NULL_INT, NULL_INT, 
                                          NULL_REAL, NULL_REAL, (IntScheme_t)NULL_INT, 
                                          OPT__GRA_P5_GRADIENT, ETA, POISSON_OFF, GRAVITY_ON, GPU_NSTREAM ); 
#        else
         CPU_PoissonGravitySolver       ( NULL, NULL, 
                                          h_Pot_Array_P_Out[ArrayID], h_Flu_Array_G[ArrayID], 
                                          NPG, dt, dh, NULL_INT, NULL_INT, 
                                          NULL_REAL, NULL_INT, NULL_INT, NULL_INT, 
                                          NULL_REAL, NULL_REAL, (IntScheme_t)NULL_INT, 
                                          OPT__GRA_P5_GRADIENT, ETA, POISSON_OFF, GRAVITY_ON ); 
#        endif
         break;


      case POISSON_AND_GRAVITY_SOLVER :  
              
#        ifdef GPU     
         CUAPI_Asyn_PoissonGravitySolver( h_Rho_Array_P[ArrayID], h_Pot_Array_P_In[ArrayID], 
                                          h_Pot_Array_P_Out[ArrayID], h_Flu_Array_G[ArrayID], 
                                          NPG, dt, dh, SOR_MIN_ITER, SOR_MAX_ITER, 
                                          SOR_OMEGA, MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, 
                                          MG_TOLERATED_ERROR, Poi_Coeff, OPT__POT_INT_SCHEME, 
                                          OPT__GRA_P5_GRADIENT, ETA, POISSON_ON, GRAVITY_ON, GPU_NSTREAM ); 
#        else
         CPU_PoissonGravitySolver       ( h_Rho_Array_P[ArrayID], h_Pot_Array_P_In[ArrayID], 
                                          h_Pot_Array_P_Out[ArrayID], h_Flu_Array_G[ArrayID], 
                                          NPG, dt, dh, SOR_MIN_ITER, SOR_MAX_ITER, 
                                          SOR_OMEGA, MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, 
                                          MG_TOLERATED_ERROR, Poi_Coeff, OPT__POT_INT_SCHEME, 
                                          OPT__GRA_P5_GRADIENT, ETA, POISSON_ON, GRAVITY_ON ); 
#        endif
         break;

#     endif // #ifdef GRAVITY


      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );

   } // switch ( TSolver )

} // FUNCTION : Solver



//-------------------------------------------------------------------------------------------------------
// Function    :  Closing_Step
// Description :  Store the updated solutions back to the patch pointers 
//
// Note        :  Use the input parameter "TSolver" to control the targeted solver 
//
// Parameter   :  TSolver     : Targeted solver
//                              --> FLUID_SOLVER               : Fluid / ELBDM solver
//                                  POISSON_SOLVER             : Poisson solver
//                                  GRAVITY_SOLVER             : Gravity solver
//                                  POISSON_AND_GRAVITY_SOLVER : Poisson + Gravity solvers
//                lv          : Targeted refinement level 
//                SaveSg      : Sandglass to store the updated data 
//                NPG         : Number of patch groups to be evaluated at a time
//                PID0_List   : List recording the patch indicies with LocalID==0 to be udpated
//                ArrayID     : Array index to load and store data ( 0 or 1 )
//-------------------------------------------------------------------------------------------------------
void Closing_Step( const Solver_t TSolver, const int lv, const int SaveSg, const int NPG, const int *PID0_List,
                   const int ArrayID )
{

   switch ( TSolver )
   {
      case FLUID_SOLVER :   
         Flu_Close( lv, SaveSg, h_Flux_Array[ArrayID], h_Flu_Array_F_Out[ArrayID], 
                    h_MinDtInfo_Fluid_Array[ArrayID], NPG, PID0_List, OPT__ADAPTIVE_DT );
         break;

#     ifdef GRAVITY
      case POISSON_SOLVER : 
         Poi_Close( lv, SaveSg, h_Pot_Array_P_Out[ArrayID], NPG, PID0_List ); 
         break;

      case GRAVITY_SOLVER :  
         Gra_Close( lv, SaveSg, h_Flu_Array_G    [ArrayID], NPG, PID0_List );
         break;

      case POISSON_AND_GRAVITY_SOLVER :  
         Poi_Close( lv, SaveSg, h_Pot_Array_P_Out[ArrayID], NPG, PID0_List ); 
         Gra_Close( lv, SaveSg, h_Flu_Array_G    [ArrayID], NPG, PID0_List );
         break;
#     endif

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );

   } // switch ( TSolver )

} // FUNCTION : Closing_Step


