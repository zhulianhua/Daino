
#ifdef INDIVIDUAL_TIMESTEP



#include "DAINO.h"

extern Timer_t *Timer_Flu_Advance[NLEVEL];
extern Timer_t *Timer_FixUp      [NLEVEL];
extern Timer_t *Timer_Flag       [NLEVEL];
extern Timer_t *Timer_Refine     [NLEVEL];
extern Timer_t *Timer_Total      [NLEVEL];
extern Timer_t *Timer_GetBuf     [NLEVEL][6];
#ifdef GRAVITY
extern Timer_t *Timer_Gra_Advance[NLEVEL];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Integration_IndiviTimeStep
// Description :  Advance all physical attributes one global step by the individual time-step scheme 
// 
// Note        :  1. An alternative function "OOC_Integration_IndiviTimeStep" is called for the out-of-core 
//                   computation 
//                2. Each step contains TWO sub-steps, and each of which will evolve the solution 
//                   at lv by 0.5*dTime
//
// Parameter   :  lv    : Targeted refinement level
//                dTime : Time interval to update the physical time at this level
//-------------------------------------------------------------------------------------------------------
void Integration_IndiviTimeStep( const int lv, const double dTime )
{

   const double dTime_HalfStep = 0.5*dTime;
   double dt_HalfStep;

   int *FluSg = patch->FluSg;
#  ifdef GRAVITY
   int *PotSg = patch->PotSg;
#  endif


   for (int HalfStep=0; HalfStep<2; HalfStep++ )
   {

//    0. calculate the evolution time-step at this half-step
// ===============================================================================================
      dt_HalfStep = Mis_dTime2dt( Time[lv], dTime_HalfStep ); 


#     ifdef TIMING
      MPI_Barrier( MPI_COMM_WORLD );
      Timer_Total[lv]->Start();
#     endif

//    1. advance the current level by the hydrodynamic flux difference
// ===============================================================================================
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "   Lv %2d: Flu_AdvanceDt, counter = %8u ... ", lv, AdvanceCounter[lv] );

      if ( OPT__OVERLAP_MPI )
      {
//       enable OpenMP nested parallelism
#        ifdef OPENMP
         omp_set_nested( true );
#        endif

//       advance patches needed to be sent
         TIMING_FUNC(   Flu_AdvanceDt( lv, Time[lv], dt_HalfStep, 1-FluSg[lv], true, true ),
                        Timer_Flu_Advance[lv],   false   );

#        pragma omp parallel sections num_threads(2)
         {
#           pragma omp section
            {
//             transfer data simultaneously               
#              ifdef GRAVITY
               TIMING_FUNC(   Buf_GetBufferData( lv, 1-FluSg[lv], NULL_INT, DATA_GENERAL, _DENS, Rho_ParaBuf, 
                                                 USELB_YES ),
                              Timer_GetBuf[lv][0],   true   );
#              else
               TIMING_FUNC(   Buf_GetBufferData( lv, 1-FluSg[lv], NULL_INT, DATA_GENERAL, _FLU,  Flu_ParaBuf,
                                                 USELB_YES ),
                              Timer_GetBuf[lv][2],   true   );
#              endif
            }

#           pragma omp section
            {
//             advance patches not needed to be sent
               TIMING_FUNC(   Flu_AdvanceDt( lv, Time[lv], dt_HalfStep, 1-FluSg[lv], true, false ),
                              Timer_Flu_Advance[lv],   true   );
            }
         } // OpenMP parallel sections

//       disable OpenMP nested parallelism
#        ifdef OPENMP
         omp_set_nested( false );
#        endif
      } // if ( OPT__OVERLAP_MPI )

      else
      {
         TIMING_FUNC(   Flu_AdvanceDt( lv, Time[lv], dt_HalfStep, 1-FluSg[lv], false, false ),
                        Timer_Flu_Advance[lv],   true   );

#        ifdef GRAVITY
         TIMING_FUNC(   Buf_GetBufferData( lv, 1-FluSg[lv], NULL_INT, DATA_GENERAL, _DENS, Rho_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][0],   true   );
#        else
         TIMING_FUNC(   Buf_GetBufferData( lv, 1-FluSg[lv], NULL_INT, DATA_GENERAL, _FLU,  Flu_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][2],   true   );
#        endif
      } // if ( OPT__OVERLAP_MPI ) ... else ...

      FluSg[lv] = 1 - FluSg[lv];

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
// ===============================================================================================


//    2. advance the current level by the gravitational acceleration
// ===============================================================================================
#     ifdef GRAVITY
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
         Aux_Message( stdout, "   Lv %2d: Gra_AdvanceDt, counter = %8u ... ", lv, AdvanceCounter[lv] );

      if ( lv == 0 )
         Gra_AdvanceDt( lv, Time[lv]+dTime_HalfStep, dt_HalfStep, FluSg[lv], true, false, false );

      else // lv > 0
      {
         if ( OPT__OVERLAP_MPI )
         {
//          enable OpenMP nested parallelism
#           ifdef OPENMP
            omp_set_nested( true );
#           endif

//          advance patches needed to be sent
            TIMING_FUNC(   Gra_AdvanceDt( lv, Time[lv]+dTime_HalfStep, dt_HalfStep, FluSg[lv], 
                           true, true, true ),
                           Timer_Gra_Advance[lv],   false   );

#           pragma omp parallel sections num_threads(2)
            {
#              pragma omp section         
               {
//                transfer data simultaneously               
                  TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, FluSg[lv], POT_FOR_POISSON, _POTE,
                                                    Pot_ParaBuf, USELB_YES ),
                                 Timer_GetBuf[lv][1],   true   );

                  TIMING_FUNC(   Buf_GetBufferData( lv, FluSg[lv], NULL_INT, DATA_GENERAL,    _FLU,
                                                    Flu_ParaBuf, USELB_YES ),
                                 Timer_GetBuf[lv][2],   true   );
               }

#              pragma omp section
               {
//                advance patches not needed to be sent
                  TIMING_FUNC(   Gra_AdvanceDt( lv, Time[lv]+dTime_HalfStep, dt_HalfStep, FluSg[lv], 
                                 true, true, false),
                                 Timer_Gra_Advance[lv],   true   );
               }
            } // OpenMP parallel sections

//          disable OpenMP nested parallelism
#           ifdef OPENMP
            omp_set_nested( false );
#           endif
         } // if ( OPT__OVERLAP_MPI )

         else
         {
            TIMING_FUNC(   Gra_AdvanceDt( lv, Time[lv]+dTime_HalfStep, dt_HalfStep, FluSg[lv], 
                                          true, false, false ),
                           Timer_Gra_Advance[lv],   true   );

            TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, FluSg[lv], POT_FOR_POISSON, _POTE, Pot_ParaBuf,
                                              USELB_YES ),
                           Timer_GetBuf[lv][1],   true   );

            TIMING_FUNC(   Buf_GetBufferData( lv, FluSg[lv], NULL_INT, DATA_GENERAL,    _FLU,  Flu_ParaBuf,
                                              USELB_YES ),
                           Timer_GetBuf[lv][2],   true   );
         } // if ( OPT__OVERLAP_MPI ) ... else ...

         PotSg[lv] = FluSg[lv];
      } // if ( lv == 0 ) ... else ...

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
#     endif // #ifdef GRAVITY
// ===============================================================================================


      Time_Prev     [lv] = Time[lv];
      Time          [lv] += dTime_HalfStep;
      AdvanceCounter[lv] ++;

      if ( AdvanceCounter[lv] >= __UINT_MAX__ ) Aux_Message( stderr, "WARNING : AdvanceCounter overflow !!\n" );


      if ( lv != NLEVEL-1  &&  NPatchTotal[lv+1] != 0 )
      {

//       3. enter the next refinement level
// ===============================================================================================
#        ifdef TIMING
         MPI_Barrier( MPI_COMM_WORLD );
         Timer_Total[lv]->Stop( false );
#        endif

         Integration_IndiviTimeStep( lv+1, dTime_HalfStep );

#        ifdef TIMING
         MPI_Barrier( MPI_COMM_WORLD );
         Timer_Total[lv]->Start();
#        endif
// ===============================================================================================


//       4. correct the data at the current level by using the data at the next finer level
// ===============================================================================================
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d: Flu_FixUp %24s... ", lv, "" );

         if ( OPT__FIXUP_FLUX )
         TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, NULL_INT, COARSE_FINE_FLUX, _FLU, NULL_INT, USELB_YES ),
                        Timer_FixUp[lv],   false   );

         TIMING_FUNC(   Flu_FixUp( lv, dt_HalfStep ),   Timer_FixUp[lv],   true   );

         if ( OPT__FIXUP_FLUX  ||  OPT__FIXUP_RESTRICT )
         TIMING_FUNC(   Buf_GetBufferData( lv, FluSg[lv], NULL_INT, DATA_AFTER_FIXUP, _FLU, Flu_ParaBuf, USELB_YES  ),
                        Timer_GetBuf[lv][3],   true   );

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
// ===============================================================================================
//
      } // if ( lv != NLEVEL-1  &&  NPatchTotal[lv+1] != 0 )


//    5. flag the current level and create patches at the next finer level
// ===============================================================================================
      if ( lv != NLEVEL-1  &&  AdvanceCounter[lv] % REGRID_COUNT == 0 )
      {
//       flag
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d: Flag %29s... ", lv, "" );

#        ifdef LOAD_BALANCE
         TIMING_FUNC(   Flag_Real( lv, USELB_YES ),       Timer_Flag[lv],   true    );
#        else
         TIMING_FUNC(   Flag_Real( lv, USELB_NO ),        Timer_Flag[lv],   false   );

         TIMING_FUNC(   MPI_ExchangeBoundaryFlag( lv ),   Timer_Flag[lv],   false   );

         TIMING_FUNC(   Flag_Buffer( lv ),                Timer_Flag[lv],   true    );
#        endif

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


//       refine
         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d: Refine %27s... ", lv, "" );

         TIMING_FUNC(   Refine( lv ),   Timer_Refine[lv],   true   );

#        ifdef LOAD_BALANCE
         TIMING_FUNC(   Buf_GetBufferData( lv,   FluSg[lv  ], NULL_INT, DATA_AFTER_REFINE, _FLU,  Flu_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][4],    false   );
#        ifdef GRAVITY
         TIMING_FUNC(   Buf_GetBufferData( lv,   NULL_INT, PotSg[lv  ], POT_AFTER_REFINE,  _POTE, Pot_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][5],    false   );
#        endif
#        endif // #ifdef LOAD_BALANCE

         TIMING_FUNC(   Buf_GetBufferData( lv+1, FluSg[lv+1], NULL_INT, DATA_AFTER_REFINE, _FLU,  Flu_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][4],    true    );
#        ifdef GRAVITY
         TIMING_FUNC(   Buf_GetBufferData( lv+1, NULL_INT, PotSg[lv+1], POT_AFTER_REFINE,  _POTE, Pot_ParaBuf,
                                           USELB_YES ),
                        Timer_GetBuf[lv][5],    true    );
#        endif

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


         Time[lv+1] = Time[lv];

         if ( OPT__PATCH_COUNT > 2 )   Aux_PatchCount();

      } // if ( lv != NLEVEL-1  &&  AdvanceCounter[lv] % REGRID_COUNT == 0 )
// ===============================================================================================

#     ifdef TIMING
      else
      {
         Timer_Flag  [lv]   ->WorkingID ++;
         Timer_Refine[lv]   ->WorkingID ++;
         Timer_GetBuf[lv][4]->WorkingID ++;
         Timer_GetBuf[lv][5]->WorkingID ++;
      } // if ( lv != NLEVEL-1  &&  AdvanceCounter[lv] % REGRID_COUNT == 0 ) ... else ...
#     endif 

#     ifdef TIMING
      MPI_Barrier( MPI_COMM_WORLD );
      Timer_Total[lv]->Stop( true );
#     endif

   } // for (int HalfStep=0; HalfStep<2; HalfStep++ )


// synchronize the time array
   if ( lv > 0 )  Time[lv] = Time[lv-1];

} // FUNCTION : Integration_IndiviTimeStep



#endif // #ifdef INDIVIDUAL_TIMESTEP
