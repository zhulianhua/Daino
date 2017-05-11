
#include "DAINO.h"


#ifdef OOC 
#ifdef INDIVIDUAL_TIMESTEP
void Timing__IndividualTimestep_OOC( const char FileName[] );
#endif
#else // #ifdef OOC
#ifdef INDIVIDUAL_TIMESTEP
void Timing__IndividualTimestep( const char FileName[] );
#else 
void Timing__SharedTimestep( const char FileName[] );
#endif
#endif // #ifdef OOC ... else ...
#ifdef TIMING_SOLVER
void Timing__Solver( const char FileName[] );
#endif


// global timing variables
// ----------------------------------------------------------
extern Timer_t *Timer_Main[6];
extern Timer_t *Timer_Flu_Advance [NLEVEL];
extern Timer_t *Timer_Gra_Advance [NLEVEL];
extern Timer_t *Timer_FixUp       [NLEVEL];
extern Timer_t *Timer_Flag        [NLEVEL];
extern Timer_t *Timer_Refine      [NLEVEL];
extern Timer_t *Timer_GetBuf      [NLEVEL][6];

#ifdef INDIVIDUAL_TIMESTEP
extern Timer_t *Timer_Total       [NLEVEL];
#else
extern Timer_t *Timer_Flu_Total   [NLEVEL];
extern Timer_t *Timer_Gra_Restrict[NLEVEL];
#endif

#ifdef TIMING_SOLVER
extern Timer_t *Timer_Pre         [NLEVEL][4];
extern Timer_t *Timer_Sol         [NLEVEL][4];
extern Timer_t *Timer_Clo         [NLEVEL][4];
extern Timer_t *Timer_Poi_PreRho  [NLEVEL];
extern Timer_t *Timer_Poi_PreFlu  [NLEVEL];
extern Timer_t *Timer_Poi_PrePot_C[NLEVEL];
extern Timer_t *Timer_Poi_PrePot_F[NLEVEL];
#endif

#if ( defined OOC  &&  defined INDIVIDUAL_TIMESTEP )
extern Timer_t *Timer_OOC_Flu_Total [NLEVEL];
extern Timer_t *Timer_OOC_Flu_Load  [NLEVEL][2];
extern Timer_t *Timer_OOC_Flu_Dump  [NLEVEL][2];
extern Timer_t *Timer_OOC_Flu_GetBuf[NLEVEL][2];
extern Timer_t *Timer_OOC_Flu_MPIBuf[NLEVEL][2];
extern Timer_t *Timer_OOC_Flu_UpdBuf[NLEVEL];
extern 
extern Timer_t *Timer_OOC_Gra_Total  [NLEVEL];
extern Timer_t *Timer_OOC_Gra_Load   [NLEVEL][2];
extern Timer_t *Timer_OOC_Gra_Dump   [NLEVEL];
extern Timer_t *Timer_OOC_Gra_GetBuf [NLEVEL][2];
extern Timer_t *Timer_OOC_Gra_MPIBuf [NLEVEL][2];
extern Timer_t *Timer_OOC_Gra_UpdBuf [NLEVEL];
extern 
extern Timer_t *Timer_OOC_Ref_Total  [NLEVEL][2];
extern Timer_t *Timer_OOC_Ref_Load   [NLEVEL][5];
extern Timer_t *Timer_OOC_Ref_Dump   [NLEVEL][4];
extern Timer_t *Timer_OOC_Ref_GetBuf [NLEVEL][4];
extern Timer_t *Timer_OOC_Ref_MPIBuf [NLEVEL][3];
extern Timer_t *Timer_OOC_Ref_UpdBuf [NLEVEL][2];
extern Timer_t *Timer_OOC_Ref_GetFlag[NLEVEL];
extern Timer_t *Timer_OOC_Ref_MPIFlag[NLEVEL];
#endif // #if ( defined OOC  &&  defined INDIVIDUAL_TIMESTEP )




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_CreateTimer
// Description :  Create simulation timers by using the "Timer" structure
//-------------------------------------------------------------------------------------------------------
void Aux_CreateTimer()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_CreateTimer ... " );


   for (int t=0; t<6; t++)    Timer_Main[t] = new Timer_t( 1 );

   for (int lv=0; lv<NLEVEL; lv++)
   {
#     ifdef INDIVIDUAL_TIMESTEP

      const int N = 1<<(1+lv);

      Timer_Flu_Advance[lv] = new Timer_t( N );
      Timer_Gra_Advance[lv] = new Timer_t( N );
      Timer_FixUp      [lv] = new Timer_t( N );
      Timer_Flag       [lv] = new Timer_t( N );
      Timer_Refine     [lv] = new Timer_t( N );
      Timer_Total      [lv] = new Timer_t( N );

      for (int t=0; t<6; t++)    Timer_GetBuf[lv][t] = new Timer_t( N );

#     ifdef OOC
      Timer_OOC_Flu_Total  [lv] = new Timer_t( N );
      Timer_OOC_Flu_UpdBuf [lv] = new Timer_t( N );

      Timer_OOC_Gra_Total  [lv] = new Timer_t( N );
      Timer_OOC_Gra_UpdBuf [lv] = new Timer_t( N );
      Timer_OOC_Gra_Dump   [lv] = new Timer_t( N );

      Timer_OOC_Ref_GetFlag[lv] = new Timer_t( N );
      Timer_OOC_Ref_MPIFlag[lv] = new Timer_t( N );

      for (int t=0; t<2; t++)
      {
         Timer_OOC_Flu_Load  [lv][t] = new Timer_t( N );
         Timer_OOC_Flu_Dump  [lv][t] = new Timer_t( N );
         Timer_OOC_Flu_GetBuf[lv][t] = new Timer_t( N );
         Timer_OOC_Flu_MPIBuf[lv][t] = new Timer_t( N );

         Timer_OOC_Gra_Load  [lv][t] = new Timer_t( N );
         Timer_OOC_Gra_GetBuf[lv][t] = new Timer_t( N );
         Timer_OOC_Gra_MPIBuf[lv][t] = new Timer_t( N );
      }

      for (int t=0; t<2; t++)    Timer_OOC_Ref_Total [lv][t] = new Timer_t( N );
      for (int t=0; t<5; t++)    Timer_OOC_Ref_Load  [lv][t] = new Timer_t( N );
      for (int t=0; t<4; t++)    Timer_OOC_Ref_Dump  [lv][t] = new Timer_t( N );
      for (int t=0; t<4; t++)    Timer_OOC_Ref_GetBuf[lv][t] = new Timer_t( N );
      for (int t=0; t<3; t++)    Timer_OOC_Ref_MPIBuf[lv][t] = new Timer_t( N );
      for (int t=0; t<2; t++)    Timer_OOC_Ref_UpdBuf[lv][t] = new Timer_t( N );
#     endif // #ifdef OOC

#     else // #ifdef INDIVIDUAL_TIMESTEP

      Timer_Flu_Total   [lv] = new Timer_t( 2 );
      Timer_Flu_Advance [lv] = new Timer_t( 2 );
      Timer_Gra_Advance [lv] = new Timer_t( 1 );
      Timer_Gra_Restrict[lv] = new Timer_t( 1 );
      Timer_FixUp       [lv] = new Timer_t( 2 );
      Timer_Flag        [lv] = new Timer_t( 2 );
      Timer_Refine      [lv] = new Timer_t( 2 );

      for (int t=0; t<6; t++)    Timer_GetBuf[lv][t] = new Timer_t( 1 );

#     endif // #ifdef INDIVIDUAL_TIMESTEP ... else ...

#     ifdef TIMING_SOLVER
      for (int v=0; v<4; v++)
      {
         Timer_Pre[lv][v]    = new Timer_t( 1 );
         Timer_Sol[lv][v]    = new Timer_t( 1 );
         Timer_Clo[lv][v]    = new Timer_t( 1 );
      }
      Timer_Poi_PreRho  [lv] = new Timer_t( 1 );
      Timer_Poi_PreFlu  [lv] = new Timer_t( 1 );
      Timer_Poi_PrePot_C[lv] = new Timer_t( 1 );
      Timer_Poi_PrePot_F[lv] = new Timer_t( 1 );
#     endif // #ifdef TIMING_SOLVER
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : Aux_CreateTimer



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_DeleteTimer
// Description :  Delete simulation timers by using the "Timer" structure
//-------------------------------------------------------------------------------------------------------
void Aux_DeleteTimer()
{

   for (int t=0; t<6; t++)    delete Timer_Main[t];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      delete Timer_Flu_Advance [lv];
      delete Timer_Gra_Advance [lv];
      delete Timer_FixUp       [lv];
      delete Timer_Flag        [lv];
      delete Timer_Refine      [lv];

      for (int t=0; t<6; t++)    delete Timer_GetBuf[lv][t];

#     ifdef INDIVIDUAL_TIMESTEP

      delete Timer_Total       [lv];

#     ifdef OOC
      delete Timer_OOC_Flu_Total  [lv];
      delete Timer_OOC_Flu_UpdBuf [lv];

      delete Timer_OOC_Gra_Total  [lv];
      delete Timer_OOC_Gra_UpdBuf [lv];
      delete Timer_OOC_Gra_Dump   [lv];

      delete Timer_OOC_Ref_GetFlag[lv];
      delete Timer_OOC_Ref_MPIFlag[lv];

      for (int t=0; t<2; t++)
      {
         delete Timer_OOC_Flu_Load  [lv][t];
         delete Timer_OOC_Flu_Dump  [lv][t];
         delete Timer_OOC_Flu_GetBuf[lv][t];
         delete Timer_OOC_Flu_MPIBuf[lv][t];

         delete Timer_OOC_Gra_Load  [lv][t];
         delete Timer_OOC_Gra_GetBuf[lv][t];
         delete Timer_OOC_Gra_MPIBuf[lv][t];
      }

      for (int t=0; t<2; t++)    delete Timer_OOC_Ref_Total [lv][t];
      for (int t=0; t<5; t++)    delete Timer_OOC_Ref_Load  [lv][t];
      for (int t=0; t<4; t++)    delete Timer_OOC_Ref_Dump  [lv][t];
      for (int t=0; t<4; t++)    delete Timer_OOC_Ref_GetBuf[lv][t];
      for (int t=0; t<3; t++)    delete Timer_OOC_Ref_MPIBuf[lv][t];
      for (int t=0; t<2; t++)    delete Timer_OOC_Ref_UpdBuf[lv][t];
#     endif // #ifdef OOC

#     else // #ifdef INDIVIDUAL_TIMESTEP

      delete Timer_Flu_Total   [lv];
      delete Timer_Gra_Restrict[lv];

#     endif // #ifdef INDIVIDUAL_TIMESTEP ... else ...

#     ifdef TIMING_SOLVER
      for (int v=0; v<4; v++)
      {
         delete Timer_Pre      [lv][v];
         delete Timer_Sol      [lv][v];
         delete Timer_Clo      [lv][v];
      }
      delete Timer_Poi_PreRho  [lv];
      delete Timer_Poi_PreFlu  [lv];
      delete Timer_Poi_PrePot_C[lv];
      delete Timer_Poi_PrePot_F[lv];
#     endif // #ifdef TIMING_SOLVER
   }

} // FUNCTION : Aux_DeleteTimer



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_ResetTimer
// Description :  Reset simulation timers by using the "Timer" structure
//-------------------------------------------------------------------------------------------------------
void Aux_ResetTimer()
{

   for (int t=0; t<6; t++)    Timer_Main[t]->Reset();

   for (int lv=0; lv<NLEVEL; lv++)
   {
      Timer_Flu_Advance [lv]->Reset();
      Timer_Gra_Advance [lv]->Reset();
      Timer_FixUp       [lv]->Reset();
      Timer_Flag        [lv]->Reset();
      Timer_Refine      [lv]->Reset();

      for (int t=0; t<6; t++)    Timer_GetBuf[lv][t]->Reset();

#     ifdef INDIVIDUAL_TIMESTEP

      Timer_Total       [lv]->Reset();

#     ifdef OOC
      Timer_OOC_Flu_Total  [lv]->Reset();
      Timer_OOC_Flu_UpdBuf [lv]->Reset();

      Timer_OOC_Gra_Total  [lv]->Reset();
      Timer_OOC_Gra_UpdBuf [lv]->Reset();
      Timer_OOC_Gra_Dump   [lv]->Reset();

      Timer_OOC_Ref_GetFlag[lv]->Reset();
      Timer_OOC_Ref_MPIFlag[lv]->Reset();

      for (int t=0; t<2; t++)
      {
         Timer_OOC_Flu_Load  [lv][t]->Reset();
         Timer_OOC_Flu_Dump  [lv][t]->Reset();
         Timer_OOC_Flu_GetBuf[lv][t]->Reset();
         Timer_OOC_Flu_MPIBuf[lv][t]->Reset();

         Timer_OOC_Gra_Load  [lv][t]->Reset();
         Timer_OOC_Gra_GetBuf[lv][t]->Reset();
         Timer_OOC_Gra_MPIBuf[lv][t]->Reset();
      }

      for (int t=0; t<2; t++)    Timer_OOC_Ref_Total [lv][t]->Reset();
      for (int t=0; t<5; t++)    Timer_OOC_Ref_Load  [lv][t]->Reset();
      for (int t=0; t<4; t++)    Timer_OOC_Ref_Dump  [lv][t]->Reset();
      for (int t=0; t<4; t++)    Timer_OOC_Ref_GetBuf[lv][t]->Reset();
      for (int t=0; t<3; t++)    Timer_OOC_Ref_MPIBuf[lv][t]->Reset();
      for (int t=0; t<2; t++)    Timer_OOC_Ref_UpdBuf[lv][t]->Reset();
#     endif // #ifdef OOC

#     else // #ifdef INDIVIDUAL_TIMESTEP

      Timer_Flu_Total   [lv]->Reset();
      Timer_Gra_Restrict[lv]->Reset();

#     endif // #ifdef INDIVIDUAL_TIMESTEP ... else ...

#     ifdef TIMING_SOLVER
      for (int v=0; v<4; v++)
      {
         Timer_Pre      [lv][v]->Reset();
         Timer_Sol      [lv][v]->Reset();
         Timer_Clo      [lv][v]->Reset();
      }
      Timer_Poi_PreRho  [lv]->Reset();
      Timer_Poi_PreFlu  [lv]->Reset();
      Timer_Poi_PrePot_C[lv]->Reset();
      Timer_Poi_PrePot_F[lv]->Reset();
#     endif // #ifdef TIMING_SOLVER
   }

} // FUNCTION : Aux_ResetTimer



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_RecordTiming
// Description :  Record the timing results (in second)
//
// Note        :  The option "TIMING_SOLVER" record the MAXIMUM values of all ranks
//-------------------------------------------------------------------------------------------------------
void Aux_RecordTiming()
{

   const char FileName[] = "Record__Timing";


// only the root rank needs to output the timing results
   if ( MPI_Rank == 0 )
   {
//    check if file already exists
      static bool FirstTime = true;
      if ( MPI_Rank == 0  &&  FirstTime )
      {
         FILE *File_Check = fopen( FileName, "r" );
         if ( File_Check != NULL )  
         {
            Aux_Message( stderr, "WARNING : the file \"%s\" already exists !!\n", FileName );
            fclose( File_Check );
         }
         FirstTime = false;
      }


      FILE *File = fopen( FileName, "a" );

      fprintf( File, "Time : %13.7e -> %13.7e,     Step : %8ld -> %8ld\n\n", Time[0]-dTime_Base, Time[0],
                                                                             Step-1, Step );

//    1. main loop
      fprintf( File, "Main Loop\n" );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "--------------------------\n" );
      fprintf( File, "%9s%17s%15s%13s%13s%15s%15s\n", "Total", "TimeStep", "Integration", "Output", "Auxiliary",
               "LoadBalance", "Sum" );

      fprintf( File, "%9.4f%17.4f%15.4f%13.4f%13.4f%15.4f%15.4f\n",
               Timer_Main[0]->GetValue( 0 ),
               Timer_Main[1]->GetValue( 0 ),
               Timer_Main[2]->GetValue( 0 ),
               Timer_Main[3]->GetValue( 0 ),
               Timer_Main[4]->GetValue( 0 ),
               Timer_Main[5]->GetValue( 0 ),
               Timer_Main[1]->GetValue( 0 ) + Timer_Main[2]->GetValue( 0 ) + Timer_Main[3]->GetValue( 0 )
             + Timer_Main[4]->GetValue( 0 ) + Timer_Main[5]->GetValue( 0 ) );

      fprintf( File, "\n\n" );
      fclose( File );


#     ifdef OOC

#     ifdef INDIVIDUAL_TIMESTEP
//    2. for individual time-step + out-of-core
      Timing__IndividualTimestep_OOC( FileName );
#     endif

#     else // #ifdef OOC ... else ...

#     ifdef INDIVIDUAL_TIMESTEP
//    3. for individual time-step
      Timing__IndividualTimestep( FileName );
#     else
//    4. for shared time-step
      Timing__SharedTimestep( FileName );
#     endif

#     endif // #ifdef OOC ... else ...
   } // if ( MPI_Rank == 0 )


#  ifdef TIMING_SOLVER
// 5. for timing GPU/CPU solvers
   Timing__Solver( FileName );
#  endif


   if ( MPI_Rank == 0 )
   {
      FILE *File = fopen( FileName, "a" );
      fprintf( File, "========================================================================================" );
      fprintf( File, "=========================\n" );
      fprintf( File, "========================================================================================" );
      fprintf( File, "=========================\n\n\n\n" );
      fclose( File );
   }

} // FUNCTION : Aux_RecordTiming



#ifdef INDIVIDUAL_TIMESTEP
//-------------------------------------------------------------------------------------------------------
// Function    :  Timing__IndividualTimestep
// Description :  Record the timing results (in second) for the individual time-step integration
//-------------------------------------------------------------------------------------------------------
void Timing__IndividualTimestep( const char FileName[] )
{

   FILE *File = fopen( FileName, "a" );

   float Total[2][NLEVEL], Flu_Advance[2][NLEVEL], Gra_Advance[2][NLEVEL], FixUp[2][NLEVEL], Flag[2][NLEVEL];
   float Refine[2][NLEVEL], GetBuf[2][NLEVEL][6], Sum[2][NLEVEL];

// HS = Half-Step   
   for (int HS=0; HS<2; HS++)
   {
      fprintf( File, "Integration HalfStep %d\n", HS+1 );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "--------------------------\n" );
      fprintf( File, "%3s%9s%10s%9s%9s%9s%9s%9s%9s%9s%9s%11s%9s\n", 
               "Lv", "Total", "Flu_Adv", "Gra_Adv", "FixUp", "Flag", "Refine", 
               "Buf_Rho", "Buf_Pot", "Buf_Flu1", "Buf_Flu2", "Buf_Refine", "Sum" );

      for (int lv=0; lv<NLEVEL; lv++)
      {
         Total      [HS][lv]    = 0.0f;
         Flu_Advance[HS][lv]    = 0.0f;
         Gra_Advance[HS][lv]    = 0.0f;
         FixUp      [HS][lv]    = 0.0f;
         Flag       [HS][lv]    = 0.0f;
         Refine     [HS][lv]    = 0.0f;
         GetBuf     [HS][lv][0] = 0.0f;
         GetBuf     [HS][lv][1] = 0.0f;
         GetBuf     [HS][lv][2] = 0.0f;
         GetBuf     [HS][lv][3] = 0.0f;
         GetBuf     [HS][lv][4] = 0.0f;
         GetBuf     [HS][lv][5] = 0.0f;

         const int N = 1<<(1+lv);

         for (int t=HS; t<N; t+=2)
         {
            Total      [HS][lv]    += Timer_Total      [lv]   ->GetValue( t );
            Flu_Advance[HS][lv]    += Timer_Flu_Advance[lv]   ->GetValue( t );
            Gra_Advance[HS][lv]    += Timer_Gra_Advance[lv]   ->GetValue( t );
            FixUp      [HS][lv]    += Timer_FixUp      [lv]   ->GetValue( t );
            Flag       [HS][lv]    += Timer_Flag       [lv]   ->GetValue( t );
            Refine     [HS][lv]    += Timer_Refine     [lv]   ->GetValue( t );
            GetBuf     [HS][lv][0] += Timer_GetBuf     [lv][0]->GetValue( t );
            GetBuf     [HS][lv][1] += Timer_GetBuf     [lv][1]->GetValue( t );
            GetBuf     [HS][lv][2] += Timer_GetBuf     [lv][2]->GetValue( t );
            GetBuf     [HS][lv][3] += Timer_GetBuf     [lv][3]->GetValue( t );
            GetBuf     [HS][lv][4] += Timer_GetBuf     [lv][4]->GetValue( t );
            GetBuf     [HS][lv][5] += Timer_GetBuf     [lv][5]->GetValue( t );
         }

         Sum[HS][lv] =  Flu_Advance[HS][lv] + Gra_Advance[HS][lv] + FixUp[HS][lv] + Flag[HS][lv] + Refine[HS][lv] 
                      + GetBuf[HS][lv][0] + GetBuf[HS][lv][1] + GetBuf[HS][lv][2] + GetBuf[HS][lv][3] 
                      + GetBuf[HS][lv][4] + GetBuf[HS][lv][5]; 

         fprintf( File, "%3d%9.4f%10.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%11.4f%9.4f\n",
                  lv, Total[HS][lv], Flu_Advance[HS][lv], Gra_Advance[HS][lv], FixUp[HS][lv], Flag[HS][lv], 
                  Refine[HS][lv], GetBuf[HS][lv][0], GetBuf[HS][lv][1], GetBuf[HS][lv][2], GetBuf[HS][lv][3], 
                  GetBuf[HS][lv][4]+GetBuf[HS][lv][5], Sum[HS][lv] );

      } // for (int lv=0; lv<NLEVEL; lv++)

//    sum over all levels
      for (int lv=1; lv<NLEVEL; lv++)
      {
         Total      [HS][0]    += Total      [HS][lv];
         Flu_Advance[HS][0]    += Flu_Advance[HS][lv];
         Gra_Advance[HS][0]    += Gra_Advance[HS][lv];
         FixUp      [HS][0]    += FixUp      [HS][lv];
         Flag       [HS][0]    += Flag       [HS][lv];
         Refine     [HS][0]    += Refine     [HS][lv];
         GetBuf     [HS][0][0] += GetBuf     [HS][lv][0];
         GetBuf     [HS][0][1] += GetBuf     [HS][lv][1];
         GetBuf     [HS][0][2] += GetBuf     [HS][lv][2];
         GetBuf     [HS][0][3] += GetBuf     [HS][lv][3];
         GetBuf     [HS][0][4] += GetBuf     [HS][lv][4];
         GetBuf     [HS][0][5] += GetBuf     [HS][lv][5];
         Sum        [HS][0]    += Sum        [HS][lv];
      }

      fprintf( File, "%3s%9.4f%10.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%11.4f%9.4f\n",
               "Sum", Total[HS][0], Flu_Advance[HS][0], Gra_Advance[HS][0], FixUp[HS][0], Flag[HS][0], 
               Refine[HS][0], GetBuf[HS][0][0], GetBuf[HS][0][1], GetBuf[HS][0][2], GetBuf[HS][0][3], 
               GetBuf[HS][0][4]+GetBuf[HS][0][5], Sum[HS][0] );
      fprintf( File, "\n" );

   } // for (int HS=0; HS<2; HS++)


// summary
   float Everything, MPI, dt, Aux, Output, LB;
   float Flu_P, Gra_P, FixUp_P, Flag_P, Refine_P, Sum_P, MPI_P, dt_P, Aux_P, Output_P, LB_P;

   Everything = Timer_Main[0]->GetValue( 0 );
   dt         = Timer_Main[1]->GetValue( 0 );
   Output     = Timer_Main[3]->GetValue( 0 );
   Aux        = Timer_Main[4]->GetValue( 0 );
   LB         = Timer_Main[5]->GetValue( 0 );

// sum
   Flu_Advance[0][0]    += Flu_Advance[1][0];
   Gra_Advance[0][0]    += Gra_Advance[1][0];
   FixUp      [0][0]    += FixUp      [1][0];
   Flag       [0][0]    += Flag       [1][0];
   Refine     [0][0]    += Refine     [1][0];
   GetBuf     [0][0][0] += GetBuf     [1][0][0];
   GetBuf     [0][0][1] += GetBuf     [1][0][1];
   GetBuf     [0][0][2] += GetBuf     [1][0][2];
   GetBuf     [0][0][3] += GetBuf     [1][0][3];
   GetBuf     [0][0][4] += GetBuf     [1][0][4];
   GetBuf     [0][0][5] += GetBuf     [1][0][5];
   Sum        [0][0]    += Sum        [1][0];

   MPI = GetBuf[0][0][0]+GetBuf[0][0][1]+GetBuf[0][0][2]+GetBuf[0][0][3]+GetBuf[0][0][4]+GetBuf[0][0][5];
   Sum[0][0] += dt + Output + Aux + LB;

// percentage
   Flu_P    = 100.0*Flu_Advance[0][0]/Everything;
   Gra_P    = 100.0*Gra_Advance[0][0]/Everything;
   FixUp_P  = 100.0*FixUp      [0][0]/Everything;
   Flag_P   = 100.0*Flag       [0][0]/Everything;
   Refine_P = 100.0*Refine     [0][0]/Everything;
   MPI_P    = 100.0*MPI              /Everything;
   dt_P     = 100.0*dt               /Everything;
   Output_P = 100.0*Output           /Everything;
   Aux_P    = 100.0*Aux              /Everything;
   LB_P     = 100.0*LB               /Everything;
   Sum_P    = 100.0*Sum        [0][0]/Everything;


   fprintf( File, "\nSummary\n" );
   fprintf( File, "---------------------------------------------------------------------------------------" );
   fprintf( File, "--------------------------\n" );
   fprintf( File, "%4s%11s%9s%9s%9s%9s%9s%9s%9s%9s%9s%12s\n", 
            "", "Flu_Adv", "Gra_Adv", "FixUp", "Flag", "Refine", "MPI", "dt", "Output", "Aux", "LB", "Sum" );

   fprintf( File, "%4s%11.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%12.4f\n",
            "Time", Flu_Advance[0][0], Gra_Advance[0][0], FixUp[0][0], Flag[0][0], Refine[0][0], MPI,
            dt, Output, Aux, LB, Sum[0][0] );
   fprintf( File, "%4s%10.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%8.3f%%%11.3f%%\n",
            "Frac", Flu_P, Gra_P, FixUp_P, Flag_P, Refine_P, MPI_P, dt_P, Output_P, Aux_P, LB_P, Sum_P );
   fprintf( File, "\n" );

   fclose( File );

} // FUNCTION : Timing__IndividualTimestep
#endif // #ifdef INDIVIDUAL_TIMESTEP



#ifndef INDIVIDUAL_TIMESTEP
//-------------------------------------------------------------------------------------------------------
// Function    :  Timing__SharedTimestep
// Description :  Record the timing results (in second) for the shared time-step integration
//-------------------------------------------------------------------------------------------------------
void Timing__SharedTimestep( const char FileName[] )
{

   FILE *File = fopen( FileName, "a" );

   float Gra_Sum = 0.0f;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      Gra_Sum += Timer_Gra_Advance [lv]   ->GetValue( 0 );
      Gra_Sum += Timer_Gra_Restrict[lv]   ->GetValue( 0 );
      Gra_Sum += Timer_GetBuf      [lv][1]->GetValue( 0 );
      Gra_Sum += Timer_GetBuf      [lv][2]->GetValue( 0 );
   }

// a. overall
   fprintf( File, "Integration_SharedTimeStep\n" );
   fprintf( File, "-----------------------------------------------------------------------------------------" );
   fprintf( File, "------------------------\n" );
   fprintf( File, "%9s%15s%11s%11s%13s\n", "Total", "Fluid1", "Gravity", "Fluid2", "Sum" );
   fprintf( File, "%9.4f%15.4f%11.4f%11.4f%13.4f\n",
            Timer_Main[2]->GetValue( 0 ), Timer_Flu_Total[0]->GetValue( 0 ), Gra_Sum, 
            Timer_Flu_Total[0]->GetValue( 1 ), 
            Timer_Flu_Total[0]->GetValue( 0 ) + Timer_Flu_Total[0]->GetValue( 1 ) + Gra_Sum );
   fprintf( File, "\n" );


// b. hydro
   float Flu_Total;
   float Flu_Sum[NLEVEL];

   for (int t=0; t<2; t++)
   {
      fprintf( File, "Fluid %d\n", t+1 );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "--------------------------\n" );
#     ifdef GRAVITY
      if ( t == 0 )
      fprintf( File, "%5s%13s%13s%9s%9s%9s%9s%11s\n", 
               "Level", "Total", "Advance",  "FixUp", "Buf_Rho", "Flag", "Refine", "Sum" );
      else
      fprintf( File, "%5s%13s%13s%9s%9s%9s%9s%11s\n", 
               "Level", "Total", "Advance",  "FixUp", "Buf_Flu", "Flag", "Refine", "Sum" );
#     else
      fprintf( File, "%5s%13s%13s%9s%9s%9s%9s%11s\n", 
               "Level", "Total", "Advance",  "FixUp", "Buf_Flu", "Flag", "Refine", "Sum" );
#     endif

      for (int lv=0; lv<NLEVEL; lv++)
      {
         if ( lv != NLEVEL - 1 )
            Flu_Total = Timer_Flu_Total[lv]->GetValue( t ) - Timer_Flu_Total[lv+1]->GetValue( t );
         else
            Flu_Total = Timer_Flu_Total[lv]->GetValue( t );

         Flu_Sum[lv] =   Timer_Flu_Advance[lv]     ->GetValue( t ) + Timer_FixUp[lv]->GetValue( t )
                       + Timer_GetBuf     [lv][t*3]->GetValue( 0 ) + Timer_Flag [lv]->GetValue( t )
                       + Timer_Refine     [lv]     ->GetValue( t );

         fprintf( File, "%5d%13.4f%13.4f%9.4f%9.4f%9.4f%9.4f%11.4f\n",
                  lv, Flu_Total, 
                  Timer_Flu_Advance[lv]     ->GetValue( t ),
                  Timer_FixUp      [lv]     ->GetValue( t ),
                  Timer_GetBuf     [lv][t*3]->GetValue( 0 ),
                  Timer_Flag       [lv]     ->GetValue( t ),
                  Timer_Refine     [lv]     ->GetValue( t ),
                  Flu_Sum          [lv] );
      }

      for (int lv=1; lv<NLEVEL; lv++)  Flu_Sum[0] += Flu_Sum[lv];
      fprintf( File, "%78.4f\n", Flu_Sum[0] );

      fprintf( File, "\n" );

   } // for (int t=0; t<2; t++)


// c. gravity
#  ifdef GRAVITY

   float PoiGra_Sum[NLEVEL];

   fprintf( File, "Poisson + Gravity\n" );
   fprintf( File, "-----------------------------------------------------------------------------------------" );
   fprintf( File, "------------------------\n" );
   fprintf( File, "%5s%13s%11s%11s%11s%13s\n", 
            "Level", "Advance", "Restrict", "Buf_Pot", "Buf_Flu", "Sum" );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      PoiGra_Sum[lv] =   Timer_Gra_Advance [lv]   ->GetValue( 0 ) + Timer_Gra_Restrict[lv]   ->GetValue( 0 )
                       + Timer_GetBuf      [lv][1]->GetValue( 0 ) + Timer_GetBuf      [lv][2]->GetValue( 0 );

      fprintf( File, "%5d%13.4f%11.4f%11.4f%11.4f%13.4f\n", 
               lv, 
               Timer_Gra_Advance [lv]   ->GetValue( 0 ), Timer_Gra_Restrict[lv]   ->GetValue( 0 ),
               Timer_GetBuf      [lv][1]->GetValue( 0 ), Timer_GetBuf      [lv][2]->GetValue( 0 ),
               PoiGra_Sum[lv] );
   }

   for (int lv=1; lv<NLEVEL; lv++)  PoiGra_Sum[0] += PoiGra_Sum[lv];
   fprintf( File, "%64.4f\n", PoiGra_Sum[0] );

   fprintf( File, "\n" );

#  endif // #ifdef GRAVITY

   fclose( File );

} // FUNCTION : Timing__SharedTimestep
#endif // #ifndef INDIVIDUAL_TIMESTEP



#if ( defined OOC  &&  defined INDIVIDUAL_TIMESTEP )
//-------------------------------------------------------------------------------------------------------
// Function    :  Timing__IndividualTimestep_OOC 
// Description :  Record the timing results (in second) for the individual time-step integration with
//                out-of-core computing 
//-------------------------------------------------------------------------------------------------------
void Timing__IndividualTimestep_OOC( const char FileName[] )
{

   FILE *File = fopen( FileName, "a" );

// a. hydrodynamics
   float Flu_Tot[NLEVEL], Flu_Adv[NLEVEL], Flu_Load[NLEVEL][2], Flu_Dump[NLEVEL][2], Flu_GetBuf[NLEVEL][2];
   float Flu_MPIBuf[NLEVEL][2], Flu_UpdBuf[NLEVEL], Flu_Sum[NLEVEL];

   for (int HalfStep=0; HalfStep<2; HalfStep++)
   {
      fprintf( File, "Hydrodynamics : HalfStep %d\n", HalfStep+1 );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "--------------------------\n" );
      fprintf( File, "%3s%9s%10s%9s%9s%9s%9s%10s%9s%10s%9s%9s%9s\n", 
               "Lv", "Total", "Flu_Adv", "L_lv", "L_lv-1", "D_lv", "D_lv-1", "Buf_Flux", "Buf_Flu", 
               "MPI_Flux", "MPI_Flu", "UpdBuf", "Sum" );

      for (int lv=0; lv<NLEVEL; lv++)
      {
         Flu_Tot   [lv]    = 0.0f;
         Flu_Adv   [lv]    = 0.0f;
         Flu_Load  [lv][0] = 0.0f;
         Flu_Load  [lv][1] = 0.0f;
         Flu_Dump  [lv][0] = 0.0f;
         Flu_Dump  [lv][1] = 0.0f;
         Flu_GetBuf[lv][0] = 0.0f;
         Flu_GetBuf[lv][1] = 0.0f;
         Flu_MPIBuf[lv][0] = 0.0f;
         Flu_MPIBuf[lv][1] = 0.0f;
         Flu_UpdBuf[lv]    = 0.0f;

         const int N = 1<<(1+lv);

         for (int t=HalfStep; t<N; t+=2)
         {
            Flu_Tot   [lv]     += Timer_OOC_Flu_Total [lv]   ->GetValue( t );
            Flu_Adv   [lv]     += Timer_Flu_Advance   [lv]   ->GetValue( t );
            Flu_Load  [lv][0]  += Timer_OOC_Flu_Load  [lv][0]->GetValue( t );
            Flu_Load  [lv][1]  += Timer_OOC_Flu_Load  [lv][1]->GetValue( t );
            Flu_Dump  [lv][0]  += Timer_OOC_Flu_Dump  [lv][0]->GetValue( t );
            Flu_Dump  [lv][1]  += Timer_OOC_Flu_Dump  [lv][1]->GetValue( t );
            Flu_GetBuf[lv][0]  += Timer_OOC_Flu_GetBuf[lv][0]->GetValue( t );
            Flu_GetBuf[lv][1]  += Timer_OOC_Flu_GetBuf[lv][1]->GetValue( t );
            Flu_MPIBuf[lv][0]  += Timer_OOC_Flu_MPIBuf[lv][0]->GetValue( t );
            Flu_MPIBuf[lv][1]  += Timer_OOC_Flu_MPIBuf[lv][1]->GetValue( t );
            Flu_UpdBuf[lv]     += Timer_OOC_Flu_UpdBuf[lv]   ->GetValue( t );
         }
         
         Flu_Sum[lv] =   Flu_Adv[lv] + Flu_Load[lv][0]+ Flu_Load[lv][1] + Flu_Dump[lv][0] + Flu_Dump[lv][1]
                       + Flu_GetBuf[lv][0] + Flu_GetBuf[lv][1] + Flu_MPIBuf[lv][0] + Flu_MPIBuf[lv][1]
                       + Flu_UpdBuf[lv];
                   
         
         fprintf( File, "%3d%9.4f%10.4f%9.4f%9.4f%9.4f%9.4f%10.4f%9.4f%10.4f%9.4f%9.4f%9.4f\n",
                  lv, Flu_Tot[lv], Flu_Adv[lv], Flu_Load[lv][0], Flu_Load[lv][1], Flu_Dump[lv][0], 
                  Flu_Dump[lv][1], Flu_GetBuf[lv][0], Flu_GetBuf[lv][1], Flu_MPIBuf[lv][0], Flu_MPIBuf[lv][1],
                  Flu_UpdBuf[lv], Flu_Sum[lv] );

      } // for (int lv=0; lv<NLEVEL; lv++)

//    sum over all levels
      for (int lv=1; lv<NLEVEL; lv++)
      {
         Flu_Tot   [0]    += Flu_Tot   [lv];
         Flu_Adv   [0]    += Flu_Adv   [lv];
         Flu_Load  [0][0] += Flu_Load  [lv][0];
         Flu_Load  [0][1] += Flu_Load  [lv][1];
         Flu_Dump  [0][0] += Flu_Dump  [lv][0];
         Flu_Dump  [0][1] += Flu_Dump  [lv][1];
         Flu_GetBuf[0][0] += Flu_GetBuf[lv][0];
         Flu_GetBuf[0][1] += Flu_GetBuf[lv][1];
         Flu_MPIBuf[0][0] += Flu_MPIBuf[lv][0];
         Flu_MPIBuf[0][1] += Flu_MPIBuf[lv][1];
         Flu_UpdBuf[0]    += Flu_UpdBuf[lv];
         Flu_Sum   [0]    += Flu_Sum   [lv];
      }

      fprintf( File, "%3s%9.4f%10.4f%9.4f%9.4f%9.4f%9.4f%10.4f%9.4f%10.4f%9.4f%9.4f%9.4f\n",
               "Sum", Flu_Tot[0], Flu_Adv[0], Flu_Load[0][0], Flu_Load[0][1], Flu_Dump[0][0],
               Flu_Dump[0][1], Flu_GetBuf[0][0], Flu_GetBuf[0][1], Flu_MPIBuf[0][0], Flu_MPIBuf[0][1],
               Flu_UpdBuf[0], Flu_Sum[0] );

      fprintf( File, "\n" );

   } // for (int HalfStep=0; HalfStep<2; HalfStep++)


// b. self-gravity
   float Gra_Tot[NLEVEL], Gra_Adv[NLEVEL], Gra_Load[NLEVEL][2], Gra_Dump[NLEVEL], Gra_GetBuf[NLEVEL][2];
   float Gra_MPIBuf[NLEVEL][2], Gra_UpdBuf[NLEVEL], Gra_Sum[NLEVEL];

   for (int HalfStep=0; HalfStep<2; HalfStep++)
   {
      fprintf( File, "Self-Gravity : HalfStep %d\n", HalfStep+1 );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "--------------------------\n" );
      fprintf( File, "%3s%9s%10s%9s%9s%9s%9s%9s%9s%9s%9s%9s\n", 
               "Lv", "Total", "Gra_Adv", "L_lv", "L_lv-1", "D_lv", "Buf_Pot", "Buf_Flu", "MPI_Pot", "MPI_Flu",
               "UpdBuf", "Sum" );

      for (int lv=0; lv<NLEVEL; lv++)
      {
         Gra_Tot   [lv]    = 0.0f;
         Gra_Adv   [lv]    = 0.0f;
         Gra_Load  [lv][0] = 0.0f;
         Gra_Load  [lv][1] = 0.0f;
         Gra_Dump  [lv]    = 0.0f;
         Gra_GetBuf[lv][0] = 0.0f;
         Gra_GetBuf[lv][1] = 0.0f;
         Gra_MPIBuf[lv][0] = 0.0f;
         Gra_MPIBuf[lv][1] = 0.0f;
         Gra_UpdBuf[lv]    = 0.0f;

         const int N = 1<<(1+lv);

         for (int t=HalfStep; t<N; t+=2)
         {
            Gra_Tot   [lv]    += Timer_OOC_Gra_Total [lv]   ->GetValue( t );
            Gra_Adv   [lv]    += Timer_Gra_Advance   [lv]   ->GetValue( t );
            Gra_Load  [lv][0] += Timer_OOC_Gra_Load  [lv][0]->GetValue( t );
            Gra_Load  [lv][1] += Timer_OOC_Gra_Load  [lv][1]->GetValue( t );
            Gra_Dump  [lv]    += Timer_OOC_Gra_Dump  [lv]   ->GetValue( t );
            Gra_GetBuf[lv][0] += Timer_OOC_Gra_GetBuf[lv][0]->GetValue( t );
            Gra_GetBuf[lv][1] += Timer_OOC_Gra_GetBuf[lv][1]->GetValue( t );
            Gra_MPIBuf[lv][0] += Timer_OOC_Gra_MPIBuf[lv][0]->GetValue( t );
            Gra_MPIBuf[lv][1] += Timer_OOC_Gra_MPIBuf[lv][1]->GetValue( t );
            Gra_UpdBuf[lv]    += Timer_OOC_Gra_UpdBuf[lv]   ->GetValue( t );
         }
         
         Gra_Sum[lv] =   Gra_Adv[lv] + Gra_Load[lv][0]+ Gra_Load[lv][1] + Gra_Dump[lv] + Gra_GetBuf[lv][0]
                       + Gra_GetBuf[lv][1] + Gra_MPIBuf[lv][0] + Gra_MPIBuf[lv][1] + Gra_UpdBuf[lv];
                   
         
         fprintf( File, "%3d%9.4f%10.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
                  lv, Gra_Tot[lv], Gra_Adv[lv], Gra_Load[lv][0], Gra_Load[lv][1], Gra_Dump[lv], Gra_GetBuf[lv][0],
                  Gra_GetBuf[lv][1], Gra_MPIBuf[lv][0], Gra_MPIBuf[lv][1], Gra_UpdBuf[lv], Gra_Sum[lv] );

      } // for (int lv=0; lv<NLEVEL; lv++)

      for (int lv=1; lv<NLEVEL; lv++)
      {
         Gra_Tot   [0]    += Gra_Tot   [lv];
         Gra_Adv   [0]    += Gra_Adv   [lv];
         Gra_Load  [0][0] += Gra_Load  [lv][0];
         Gra_Load  [0][1] += Gra_Load  [lv][1];
         Gra_Dump  [0]    += Gra_Dump  [lv];
         Gra_GetBuf[0][0] += Gra_GetBuf[lv][0];
         Gra_GetBuf[0][1] += Gra_GetBuf[lv][1];
         Gra_MPIBuf[0][0] += Gra_MPIBuf[lv][0];
         Gra_MPIBuf[0][1] += Gra_MPIBuf[lv][1];
         Gra_UpdBuf[0]    += Gra_UpdBuf[lv];
         Gra_Sum   [0]    += Gra_Sum   [lv];
      }

      fprintf( File, "%3s%9.4f%10.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",
               "Sum", Gra_Tot[0], Gra_Adv[0], Gra_Load[0][0], Gra_Load[0][1], Gra_Dump[0], Gra_GetBuf[0][0],
               Gra_GetBuf[0][1], Gra_MPIBuf[0][0], Gra_MPIBuf[0][1], Gra_UpdBuf[0], Gra_Sum[0] );

      fprintf( File, "\n" );

   } // for (int HalfStep=0; HalfStep<2; HalfStep++)


// c. Fix-up + Flag
   float Ref_Tot[NLEVEL][2], Ref_FixUp[NLEVEL], Ref_Flag[NLEVEL], Ref_Refine[NLEVEL], Ref_Load[NLEVEL][5];
   float Ref_Dump[NLEVEL][4], Ref_GetBuf[NLEVEL][4], Ref_MPIBuf[NLEVEL][3], Ref_UpdBuf[NLEVEL][2];
   float Ref_GetFlag[NLEVEL], Ref_MPIFlag[NLEVEL], Ref_Sum[NLEVEL];

   for (int HalfStep=0; HalfStep<2; HalfStep++)
   {
      fprintf( File, "Fix-up + Flag : HalfStep %d\n", HalfStep+1 );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "--------------------------\n" );
      fprintf( File, "%3s%8s%8s%8s%8s%8s%8s%9s%8s%9s%8s%9s%8s%8s\n", 
               "Lv", "Total", "FixUp", "Flag", "L_lv", "L_lv+1", "D_lv", "Buf_Flux", "Buf_Flu", "Buf_Flag",
               "MPI_Flu", "MPI_Flag", "UpdBuf", "Sum" );

      for (int lv=0; lv<NLEVEL; lv++)
      {
         Ref_Tot    [lv][0] = 0.0f;
         Ref_FixUp  [lv]    = 0.0f;
         Ref_Flag   [lv]    = 0.0f;
         Ref_Load   [lv][0] = 0.0f;
         Ref_Load   [lv][1] = 0.0f;
         Ref_Dump   [lv][0] = 0.0f;
         Ref_GetBuf [lv][0] = 0.0f;
         Ref_GetBuf [lv][1] = 0.0f;
         Ref_GetFlag[lv]    = 0.0f;
         Ref_MPIBuf [lv][0] = 0.0f;
         Ref_MPIFlag[lv]    = 0.0f;
         Ref_UpdBuf [lv][0] = 0.0f;

         const int N = 1<<(1+lv);

         for (int t=HalfStep; t<N; t+=2)
         {
            Ref_Tot    [lv][0] += Timer_OOC_Ref_Total  [lv][0]->GetValue( t );
            Ref_FixUp  [lv]    += Timer_FixUp          [lv]   ->GetValue( t );
            Ref_Flag   [lv]    += Timer_Flag           [lv]   ->GetValue( t );
            Ref_Load   [lv][0] += Timer_OOC_Ref_Load   [lv][0]->GetValue( t );
            Ref_Load   [lv][1] += Timer_OOC_Ref_Load   [lv][1]->GetValue( t );
            Ref_Dump   [lv][0] += Timer_OOC_Ref_Dump   [lv][0]->GetValue( t );
            Ref_GetBuf [lv][0] += Timer_OOC_Ref_GetBuf [lv][0]->GetValue( t );
            Ref_GetBuf [lv][1] += Timer_OOC_Ref_GetBuf [lv][1]->GetValue( t );
            Ref_GetFlag[lv]    += Timer_OOC_Ref_GetFlag[lv]   ->GetValue( t );
            Ref_MPIBuf [lv][0] += Timer_OOC_Ref_MPIBuf [lv][0]->GetValue( t );
            Ref_MPIFlag[lv]    += Timer_OOC_Ref_MPIFlag[lv]   ->GetValue( t );
            Ref_UpdBuf [lv][0] += Timer_OOC_Ref_UpdBuf [lv][0]->GetValue( t );
         }
         
         Ref_Sum[lv] =   Ref_FixUp[lv] + Ref_Flag[lv] + Ref_Load[lv][0]+ Ref_Load[lv][1] + Ref_Dump[lv][0]
                       + Ref_GetBuf[lv][0] + Ref_GetBuf[lv][1] + Ref_GetFlag[lv] + Ref_MPIBuf[lv][0]
                       + Ref_MPIFlag[lv] + Ref_UpdBuf[lv][0];
                   
         
         fprintf( File, "%3d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%9.3f%8.3f%9.3f%8.3f%9.3f%8.3f%8.3f\n",
                  lv, Ref_Tot[lv][0], Ref_FixUp[lv], Ref_Flag[lv], Ref_Load[lv][0], Ref_Load[lv][1],
                  Ref_Dump[lv][0], Ref_GetBuf[lv][0], Ref_GetBuf[lv][1], Ref_GetFlag[lv], Ref_MPIBuf[lv][0],
                  Ref_MPIFlag[lv], Ref_UpdBuf[lv][0], Ref_Sum[lv] );

      } // for (int lv=0; lv<NLEVEL; lv++)

      for (int lv=1; lv<NLEVEL; lv++)
      {
         Ref_Tot    [0][0] += Ref_Tot    [lv][0];
         Ref_FixUp  [0]    += Ref_FixUp  [lv];
         Ref_Flag   [0]    += Ref_Flag   [lv];
         Ref_Load   [0][0] += Ref_Load   [lv][0];
         Ref_Load   [0][1] += Ref_Load   [lv][1];
         Ref_Dump   [0][0] += Ref_Dump   [lv][0];
         Ref_GetBuf [0][0] += Ref_GetBuf [lv][0];
         Ref_GetBuf [0][1] += Ref_GetBuf [lv][1];
         Ref_GetFlag[0]    += Ref_GetFlag[lv];
         Ref_MPIBuf [0][0] += Ref_MPIBuf [lv][0];
         Ref_MPIFlag[0]    += Ref_MPIFlag[lv];
         Ref_UpdBuf [0][0] += Ref_UpdBuf [lv][0];
         Ref_Sum    [0]    += Ref_Sum    [lv];
      }

      fprintf( File, "%3s%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%9.3f%8.3f%9.3f%8.3f%9.3f%8.3f%8.3f\n",
               "Sum", Ref_Tot[0][0], Ref_FixUp[0], Ref_Flag[0], Ref_Load[0][0], Ref_Load[0][1],
               Ref_Dump[0][0], Ref_GetBuf[0][0], Ref_GetBuf[0][1], Ref_GetFlag[0], Ref_MPIBuf[0][0],
               Ref_MPIFlag[0], Ref_UpdBuf[0][0], Ref_Sum[0] );

      fprintf( File, "\n" );

   } // for (int HalfStep=0; HalfStep<2; HalfStep++)


// d. Refine
   for (int HalfStep=0; HalfStep<2; HalfStep++)
   {
      fprintf( File, "Refine : HalfStep %d\n", HalfStep+1 );
      fprintf( File, "---------------------------------------------------------------------------------------" );
      fprintf( File, "--------------------------\n" );
      fprintf( File, "%3s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%7s\n", 
               "Lv", "Total", "Refine", "L_lv", "L_lv+1", "L_lv+2", "D_lv", "D_lv+1", "D_lv+2", "Buf_Flu",
               "Buf_Pot", "MPI_Flu", "MPI_Pot", "UpdBuf", "Sum" );

      for (int lv=0; lv<NLEVEL; lv++)
      {
         Ref_Tot   [lv][1] = 0.0f;
         Ref_Refine[lv]    = 0.0f;
         Ref_Load  [lv][2] = 0.0f;
         Ref_Load  [lv][3] = 0.0f;
         Ref_Load  [lv][4] = 0.0f;
         Ref_Dump  [lv][1] = 0.0f;
         Ref_Dump  [lv][2] = 0.0f;
         Ref_Dump  [lv][3] = 0.0f;
         Ref_GetBuf[lv][2] = 0.0f;
         Ref_GetBuf[lv][3] = 0.0f;
         Ref_MPIBuf[lv][1] = 0.0f;
         Ref_MPIBuf[lv][2] = 0.0f;
         Ref_UpdBuf[lv][1] = 0.0f;

         const int N = 1<<(1+lv);

         for (int t=HalfStep; t<N; t+=2)
         {
            Ref_Tot   [lv][1] += Timer_OOC_Ref_Total  [lv][1]->GetValue( t );
            Ref_Refine[lv]    += Timer_Refine         [lv]   ->GetValue( t );
            Ref_Load  [lv][2] += Timer_OOC_Ref_Load   [lv][2]->GetValue( t );
            Ref_Load  [lv][3] += Timer_OOC_Ref_Load   [lv][3]->GetValue( t );
            Ref_Load  [lv][4] += Timer_OOC_Ref_Load   [lv][4]->GetValue( t );
            Ref_Dump  [lv][1] += Timer_OOC_Ref_Dump   [lv][1]->GetValue( t );
            Ref_Dump  [lv][2] += Timer_OOC_Ref_Dump   [lv][2]->GetValue( t );
            Ref_Dump  [lv][3] += Timer_OOC_Ref_Dump   [lv][3]->GetValue( t );
            Ref_GetBuf[lv][2] += Timer_OOC_Ref_GetBuf [lv][2]->GetValue( t );
            Ref_GetBuf[lv][3] += Timer_OOC_Ref_GetBuf [lv][3]->GetValue( t );
            Ref_MPIBuf[lv][1] += Timer_OOC_Ref_MPIBuf [lv][1]->GetValue( t );
            Ref_MPIBuf[lv][2] += Timer_OOC_Ref_MPIBuf [lv][2]->GetValue( t );
            Ref_UpdBuf[lv][1] += Timer_OOC_Ref_UpdBuf [lv][1]->GetValue( t );
         }
         
         Ref_Sum[lv] =   Ref_Refine[lv] + Ref_Load[lv][2]+ Ref_Load[lv][3] + Ref_Load[lv][4] + Ref_Dump[lv][1]
                       + Ref_Dump[lv][2] + Ref_Dump[lv][3] + Ref_GetBuf[lv][2] + Ref_GetBuf[lv][3]
                       + Ref_MPIBuf[lv][1] + Ref_MPIBuf[lv][2] + Ref_UpdBuf[lv][1];
                   
         
         fprintf( File, "%3d%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%7.2f\n",
                  lv, Ref_Tot[lv][1], Ref_Refine[lv], Ref_Load[lv][2], Ref_Load[lv][3], Ref_Load[lv][4],
                  Ref_Dump[lv][1], Ref_Dump[lv][2], Ref_Dump[lv][3], Ref_GetBuf[lv][2], Ref_GetBuf[lv][3],
                  Ref_MPIBuf[lv][1], Ref_MPIBuf[lv][2], Ref_UpdBuf[lv][1], Ref_Sum[lv] );

      } // for (int lv=0; lv<NLEVEL; lv++)

      for (int lv=1; lv<NLEVEL; lv++)
      {
         Ref_Tot   [0][1] += Ref_Tot   [lv][1];
         Ref_Refine[0]    += Ref_Refine[lv];
         Ref_Load  [0][2] += Ref_Load  [lv][2];
         Ref_Load  [0][3] += Ref_Load  [lv][3];
         Ref_Load  [0][4] += Ref_Load  [lv][4];
         Ref_Dump  [0][1] += Ref_Dump  [lv][1];
         Ref_Dump  [0][2] += Ref_Dump  [lv][2];
         Ref_Dump  [0][3] += Ref_Dump  [lv][3];
         Ref_GetBuf[0][2] += Ref_GetBuf[lv][2];
         Ref_GetBuf[0][3] += Ref_GetBuf[lv][3];
         Ref_MPIBuf[0][1] += Ref_MPIBuf[lv][1];
         Ref_MPIBuf[0][2] += Ref_MPIBuf[lv][2];
         Ref_UpdBuf[0][1] += Ref_UpdBuf[lv][1];
         Ref_Sum   [0]    += Ref_Sum   [lv];
      }

      fprintf( File, "%3s%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%7.2f\n",
               "Sum", Ref_Tot[0][1], Ref_Refine[0], Ref_Load[0][2], Ref_Load[0][3], Ref_Load[0][4],
               Ref_Dump[0][1], Ref_Dump[0][2], Ref_Dump[0][3], Ref_GetBuf[0][2], Ref_GetBuf[0][3],
               Ref_MPIBuf[0][1], Ref_MPIBuf[0][2], Ref_UpdBuf[0][1], Ref_Sum[0] );

      fprintf( File, "\n" );

   } // for (int HalfStep=0; HalfStep<2; HalfStep++)

   fclose( File );

} // FUNCTION : Timing__IndividualTimestep_OOC
#endif // #if ( defined OOC  &&  defined INDIVIDUAL_TIMESTEP )



#ifdef TIMING_SOLVER
//-------------------------------------------------------------------------------------------------------
// Function    :  Timing__Solver
// Description :  Record the timing results (in second) for the option "TIMING_SOLVER"
//-------------------------------------------------------------------------------------------------------
void Timing__Solver( const char FileName[] )
{

// get the maximum values from all ranks
   int   ID;
   float Pre_loc[NLEVEL][4], Sol_loc[NLEVEL][4], Clo_loc[NLEVEL][4];
   float Pre_max[NLEVEL][4], Sol_max[NLEVEL][4], Clo_max[NLEVEL][4];
   float Poi_PreRho_loc[NLEVEL], Poi_PreFlu_loc[NLEVEL], Poi_PrePot_C_loc[NLEVEL], Poi_PrePot_F_loc[NLEVEL];
   float Poi_PreRho_max[NLEVEL], Poi_PreFlu_max[NLEVEL], Poi_PrePot_C_max[NLEVEL], Poi_PrePot_F_max[NLEVEL];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int v=0; v<4; v++)
      {
         Pre_loc[lv][v] = Timer_Pre[lv][v]->GetValue(0);
         Sol_loc[lv][v] = Timer_Sol[lv][v]->GetValue(0);
         Clo_loc[lv][v] = Timer_Clo[lv][v]->GetValue(0);
      }

      Poi_PreRho_loc  [lv] = Timer_Poi_PreRho  [lv]->GetValue(0);
      Poi_PreFlu_loc  [lv] = Timer_Poi_PreFlu  [lv]->GetValue(0);
      Poi_PrePot_C_loc[lv] = Timer_Poi_PrePot_C[lv]->GetValue(0);
      Poi_PrePot_F_loc[lv] = Timer_Poi_PrePot_F[lv]->GetValue(0);
   }

   MPI_Reduce( Pre_loc[0],       Pre_max[0],       NLEVEL*4, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Sol_loc[0],       Sol_max[0],       NLEVEL*4, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Clo_loc[0],       Clo_max[0],       NLEVEL*4, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD );

   MPI_Reduce( Poi_PreRho_loc,   Poi_PreRho_max,   NLEVEL,   MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Poi_PreFlu_loc,   Poi_PreFlu_max,   NLEVEL,   MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Poi_PrePot_C_loc, Poi_PrePot_C_max, NLEVEL,   MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Poi_PrePot_F_loc, Poi_PrePot_F_max, NLEVEL,   MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD );


   if ( MPI_Rank == 0 )
   {
      FILE *File = fopen( FileName, "a" );

      fprintf( File, "\nGPU/CPU solvers\n" );
      fprintf( File, "----------------------------------------------------------------------------------------" );
      fprintf( File, "-------------------------\n" );
      fprintf( File, "%3s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n", 
               "Lv", "Flu_Pre", "Flu_Sol", "Flu_Clo", "Poi_Pre", "PreRho", "PreFlu", "PrePot_C", "Pre_Pot_F", 
               "Poi_Sol", "Poi_Clo" );

      for (int lv=0; lv<NLEVEL; lv++)
      {
         if ( lv == 0 )    ID = 2;
         else              ID = 3;

         fprintf( File, "%3d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n", 
                  lv, 
                  Pre_max[lv][ 0],                                             Sol_max[lv][ 0], Clo_max[lv][ 0],
                  Pre_max[lv][ID], Poi_PreRho_max  [lv], Poi_PreFlu_max  [lv], 
                                   Poi_PrePot_C_max[lv], Poi_PrePot_F_max[lv], Sol_max[lv][ID], Clo_max[lv][ID] );
      } 

//    sum over all levels
      for (int lv=1; lv<NLEVEL; lv++)
      {
         Pre_max         [0][0] += Pre_max         [lv][0];
         Sol_max         [0][0] += Sol_max         [lv][0];
         Clo_max         [0][0] += Clo_max         [lv][0];
         Pre_max         [0][2] += Pre_max         [lv][3];
         Poi_PreRho_max  [0]    += Poi_PreRho_max  [lv];
         Poi_PreFlu_max  [0]    += Poi_PreFlu_max  [lv];
         Poi_PrePot_C_max[0]    += Poi_PrePot_C_max[lv];
         Poi_PrePot_F_max[0]    += Poi_PrePot_F_max[lv];
         Sol_max         [0][2] += Sol_max         [lv][3];
         Clo_max         [0][2] += Clo_max         [lv][3];
      }

      fprintf( File, "%3s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
               "Sum", 
               Pre_max[0][0],                                           Sol_max[0][0], Clo_max[0][0],
               Pre_max[0][2], Poi_PreRho_max  [0], Poi_PreFlu_max  [0], 
                              Poi_PrePot_C_max[0], Poi_PrePot_F_max[0], Sol_max[0][2], Clo_max[0][2] );
      fprintf( File, "\n" );

      fclose( File );
   } // if ( MPI_Rank == 0 )

} // FUNCTION : TimingSolver
#endif // TIMING_SOLVER



