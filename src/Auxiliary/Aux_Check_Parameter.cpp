
#include "DAINO.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_Parameter
// Description :  Check the initial parameter setting
//-------------------------------------------------------------------------------------------------------
void Aux_Check_Parameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_Check_Parameter ... \n" ); 


// general errors
// =======================================================================================
#  if ( PATCH_SIZE%2 != 0 )
#     error : ERROR : PATCH_SIZE must be an even number !!
#  endif

#  if ( defined TIMING_SOLVER  &&  !defined TIMING )
#     error : ERROR : option TIMING_SOLVER must work with the option TIMING !!
#  endif 

#  if ( defined OPENMP  &&  !defined _OPENMP )
#     error : ERROR : something is wrong in OpenMP, the macro "_OPENMP" is NOT defined !!
#  endif

#  if ( defined OVERLAP_MPI  &&  !defined LOAD_BALANCE )
#     error : ERROR : option OVERLAP_MPI must work with the option LOAD_BALANCE !!
#  endif 

#  ifndef INDIVIDUAL_TIMESTEP 
#     error : ERROR : currently the option "INDIVIDUAL_TIMESTEP" must be turned on !!
#  endif 

#  ifdef SERIAL
   int NRank = 1;
#  else
   int NRank, MPI_Thread_Status, MPI_Init_Status;

   MPI_Initialized( &MPI_Init_Status );
   if ( MPI_Init_Status == false )  Aux_Error( ERROR_INFO, "MPI_Init has not been called !!\n" );

   MPI_Query_thread( &MPI_Thread_Status );
   if ( OPT__OVERLAP_MPI  &&  MPI_Thread_Status == MPI_THREAD_SINGLE )
      Aux_Error( ERROR_INFO, "option \"%s\" is NOT supported if the level of MPI thread support == %s\n",
                 "OPT__OVERLAP_MPI", "MPI_THREAD_SINGLE" );

   MPI_Comm_size( MPI_COMM_WORLD, &NRank );
#  endif

   if ( MPI_NRank != NRank )  
      Aux_Error( ERROR_INFO, "MPI_NRank (%d) != MPI_Comm_size (%d) !!\n", MPI_NRank, NRank );

   if ( NX0_TOT[0] <= 0  ||  NX0_TOT[1] <= 0  ||  NX0_TOT[2] <= 0 )
      Aux_Error( ERROR_INFO, "incorrect number of base-level grids --> NX0_TOT[0/1/2] = [%d,%d,%d] !!\n", 
                 NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] );

   if ( NX0_TOT[0]%PS2 != 0  ||  NX0_TOT[1]%PS2 != 0  ||  NX0_TOT[2]%PS2 != 0 )
      Aux_Error( ERROR_INFO, "number of base-level patches in each direction must be \"%s\" !!\n",
                 "a multiple of TWO" );

#  ifdef LOAD_BALANCE
   if ( OPT__INIT != INIT_RESTART )
#  endif
   if ( NX0_TOT[0]%(PS2*DAINO_NRANK_X(0)) != 0  ||  NX0_TOT[1]%(PS2*DAINO_NRANK_X(1)) != 0  ||  
        NX0_TOT[2]%(PS2*DAINO_NRANK_X(2)) != 0 )
      Aux_Error( ERROR_INFO, "number of base-level patches in each direction in one rank must be \"%s\" !!\n",
                 "a multiple of TWO" );

   if ( OPT__INIT != INIT_STARTOVER  &&  OPT__INIT != INIT_RESTART  &&  OPT__INIT != INIT_UM )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__INIT = %d\" [1/2/3] !!\n", OPT__INIT );

   if ( MPI_NRank_X[0]*MPI_NRank_X[1]*MPI_NRank_X[2] != MPI_NRank )
      Aux_Error( ERROR_INFO, "MPI_NRank_X[0]*MPI_NRank_X[1]*MPI_NRank_X[2] != MPI_NRank !!\n" );

   if ( FLAG_BUFFER_SIZE > PATCH_SIZE )   Aux_Error( ERROR_INFO, "FLAG_BUFFER_SIZE > PATCH_SIZE !!\n" );

   if ( REGRID_COUNT <= 0 )   Aux_Error( ERROR_INFO, "REGRID_COUNT <= 0 !!\n" );

   if ( OPT__OUTPUT_MODE != OUTPUT_CONST_STEP  &&  OPT__OUTPUT_MODE != OUTPUT_CONST_DT  &&
        OPT__OUTPUT_MODE != OUTPUT_USE_TABLE )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__OUTPUT_MODE = %d\" [1/2/3] !!\n", OPT__OUTPUT_MODE );

   if ( OPT__OUTPUT_MODE == OUTPUT_CONST_STEP  &&  OUTPUT_STEP <= 0 )    
      Aux_Error( ERROR_INFO, "OUTPUT_STEP <= 0 !!\n" );

   if ( OPT__OUTPUT_MODE == OUTPUT_CONST_DT  &&  OUTPUT_DT <= 0 )      
      Aux_Error( ERROR_INFO, "OUTPUT_DT <= 0 !!\n" );

   if ( OPT__RESTART_HEADER != RESTART_HEADER_CHECK  &&  OPT__RESTART_HEADER != RESTART_HEADER_SKIP )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__RESTART_HEADER = %d\" [0/1] !!\n", 
                 OPT__RESTART_HEADER );

   if ( OPT__OUTPUT_TOTAL < 0  ||  OPT__OUTPUT_TOTAL > 2 ) 
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__OUTPUT_TOTAL = %d\" [0/1/2] !!\n", OPT__OUTPUT_TOTAL );

   if ( OPT__OUTPUT_PART != OUTPUT_NONE  &&  OPT__OUTPUT_PART != OUTPUT_DIAG  &&  
        OPT__OUTPUT_PART != OUTPUT_XY  &&  OPT__OUTPUT_PART != OUTPUT_YZ  &&  OPT__OUTPUT_PART != OUTPUT_XZ  &&  
        OPT__OUTPUT_PART != OUTPUT_X   &&  OPT__OUTPUT_PART != OUTPUT_Y   &&  OPT__OUTPUT_PART != OUTPUT_Z )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__OUTPUT_PART = %d\" [0 ~ 6] !!\n", OPT__OUTPUT_PART );

   if (  ( OPT__OUTPUT_PART == OUTPUT_YZ  ||  OPT__OUTPUT_PART == OUTPUT_Y  ||  OPT__OUTPUT_PART == OUTPUT_Z )  &&
         ( OUTPUT_PART_X < 0.0  ||  OUTPUT_PART_X >= patch->BoxSize[0] )  )
      Aux_Error( ERROR_INFO, "incorrect OUTPUT_PART_X (out of range [0<=X<%f]) !!\n", patch->BoxSize[0] );

   if (  ( OPT__OUTPUT_PART == OUTPUT_XZ  ||  OPT__OUTPUT_PART == OUTPUT_X  ||  OPT__OUTPUT_PART == OUTPUT_Z )  &&
         ( OUTPUT_PART_Y < 0.0  ||  OUTPUT_PART_Y >= patch->BoxSize[1] )  )
      Aux_Error( ERROR_INFO, "incorrect OUTPUT_PART_Y (out of range [0<=Y<%f]) !!\n", patch->BoxSize[1] );
      
   if (  ( OPT__OUTPUT_PART == OUTPUT_XY  ||  OPT__OUTPUT_PART == OUTPUT_X  ||  OPT__OUTPUT_PART == OUTPUT_Y )  &&
         ( OUTPUT_PART_Z < 0.0  ||  OUTPUT_PART_Z >= patch->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "incorrect OUTPUT_PART_Z (out of range [0<=Z<%f]) !!\n", patch->BoxSize[2] );

   if (  OPT__OUTPUT_PART == OUTPUT_DIAG  &&  ( NX0_TOT[0] != NX0_TOT[1] || NX0_TOT[0] != NX0_TOT[2] )  )
      Aux_Error( ERROR_INFO, "option \"%s\" only works with CUBIC domain !!\n", 
                 "OPT__OUTPUT_PART == 7 (OUTPUT_DIAG)" );

   if (  OPT__OUTPUT_BASEPS  &&  ( NX0_TOT[0] != NX0_TOT[1] || NX0_TOT[0] != NX0_TOT[2] )  )
      Aux_Error( ERROR_INFO, "option \"%s\" only works with CUBIC domain !!\n", "OPT__OUTPUT_BASEPS" );

   if (  OPT__INIT == INIT_UM  &&  ( OPT__UM_START_LEVEL < 0 || OPT__UM_START_LEVEL > NLEVEL-1 )  )
      Aux_Error( ERROR_INFO, "incorrect option \"OPT__UM_START_LEVEL = %d\" [0 ... NLEVEL-1] !!\n", 
                 OPT__UM_START_LEVEL );

   if (  OPT__INIT == INIT_UM  &&  ( OPT__UM_START_NVAR < 1 || OPT__UM_START_NVAR > NCOMP )  )
      Aux_Error( ERROR_INFO, "incorrect option \"OPT__UM_START_NVAR  = %d\" [1 ... NCOMP] !!\n",
                 OPT__UM_START_NVAR );

   if ( OPT__CK_CONSERVATION < 0  ||  OPT__CK_CONSERVATION > 2 )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__CK_CONSERVATION = %d\" [1/2/3] !!\n", 
                 OPT__CK_CONSERVATION );

   if ( OPT__CK_REFINE  &&  !OPT__FLAG_RHO )
      Aux_Error( ERROR_INFO, "currently the check \"%s\" must work with \"%s\" !!\n",
                 "OPT__CK_REFINE", "OPT__FLAG_RHO == 1" );

   if ( OPT__FLAG_LOHNER  &&  Flu_ParaBuf < 2 )
      Aux_Error( ERROR_INFO, "Lohner error estimator does NOT work when Flu_ParaBuf (%d) < 2 !!\n", Flu_ParaBuf );

   if ( OPT__FLAG_LOHNER < 0 )   Aux_Error( ERROR_INFO, "OPT__FLAG_LOHNER < 0 !!\n" );

#  if   ( MODEL == HYDRO )
   if ( OPT__FLAG_LOHNER > 2 )   Aux_Error( ERROR_INFO, "OPT__FLAG_LOHNER > 2 !!\n" );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   if ( OPT__FLAG_LOHNER > 1 )   Aux_Error( ERROR_INFO, "OPT__FLAG_LOHNER > 1 !!\n" );

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef GPU
   if ( OPT__GPUID_SELECT < 0  &&  OPT__GPUID_SELECT != -1  &&  OPT__GPUID_SELECT != -2 )
   {
#     ifdef LAOHU
      if ( OPT__GPUID_SELECT != -3 )
#     endif
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n", 
               "OPT__GPUID_SELECT", OPT__GPUID_SELECT );
   }
#  endif

#  ifdef SERIAL
   if ( MPI_NRank != 1 )   Aux_Error( ERROR_INFO, "\"MPI_NRank != 1\" in the serial code !!\n" );

   if ( MPI_NRank_X[0] != 1 )    Aux_Error( ERROR_INFO, "\"MPI_NRank_X[0] != 1\" in the serial code !!\n" );

   if ( MPI_NRank_X[1] != 1 )    Aux_Error( ERROR_INFO, "\"MPI_NRank_X[1] != 1\" in the serial code !!\n" );

   if ( MPI_NRank_X[2] != 1 )    Aux_Error( ERROR_INFO, "\"MPI_NRank_X[2] != 1\" in the serial code !!\n" );
#  endif // #ifdef SERIAL

#  ifndef OVERLAP_MPI
   if ( OPT__OVERLAP_MPI )
      Aux_Error( ERROR_INFO, "\"%s\" is NOT turned on in the Makefile for the option \"%s\" !!\n",
                 "OVERLAP_MPI", "OPT__OVERLAP_MPI" );
#  endif



// general warnings
// =======================================================================================
#  ifdef OVERLAP_MPI
#     warning : NOTE : make sure to link with multi-thread supported MPI and FFTW for the option "OVERLAP_MPI"
#  endif

   if ( MPI_Rank == 0 ) {

   if ( OPT__OUTPUT_TOTAL == 0  &&  OPT__OUTPUT_PART == OUTPUT_NONE  &&  !OPT__OUTPUT_ERROR  &&  !OPT__OUTPUT_BASEPS )
      Aux_Message( stderr, "WARNING : both %s, %s, %s, and %s are off --> no data will be output !!\n",
                   "OPT__OUTPUT_TOTAL", "OPT__OUTPUT_PART", "OPT__OUTPUT_ERROR", "OPT__OUTPUT_BASEPS" );

   if ( OPT__CK_REFINE )
      Aux_Message( stderr, "WARNING : \"%s\" check may fail due to the proper-nesting constraint !!\n",
                   "OPT__CK_REFINE" );

   if ( !OPT__INIT_RESTRICT )
      Aux_Message( stderr, "WARNING : option OPT__INIT_RESTRICT is NOT turned on !!\n" );

#  ifdef TIMING_SOLVER
   Aux_Message( stderr, "WARNING : option \"TIMING_SOLVER\" will disable the concurrent execution\n" );
   Aux_Message( stderr, "          between GPU and CPU and hence will decrease the overall performance !!\n" );
#  endif

   if ( OPT__RESTART_HEADER == RESTART_HEADER_SKIP )
   {
      Aux_Message( stderr, "WARNING : to skip the header check during restart, you must make sure that \n" );
      Aux_Message( stderr, "          everything is set correctly in both Input__Parameter and Makefile !!\n" );
   }

   bool Flag = ( OPT__FLAG_RHO_GRADIENT  ||  OPT__FLAG_USER );
#  if ( MODEL == HYDRO )
   Flag |= OPT__FLAG_PRES_GRADIENT;
#  endif
#  if ( MODEL == ELBDM )
   Flag |= OPT__FLAG_ENGY_DENSITY;
#  endif
   if ( OPT__CK_REFINE  &&  Flag )
      Aux_Message( stderr, "WARNING : currently the check \"%s\" only works with \"%s\" !!\n",
                   "OPT__CK_REFINE", "OPT__FLAG_RHO == 1" );

   Flag = ( OPT__FLAG_RHO  ||  OPT__FLAG_RHO_GRADIENT  ||  OPT__FLAG_LOHNER  ||  OPT__FLAG_USER );
#  if ( MODEL == HYDRO )
   Flag |= OPT__FLAG_PRES_GRADIENT;
#  endif
#  if ( MODEL == ELBDM )
   Flag |= OPT__FLAG_ENGY_DENSITY;
#  endif

   if ( !Flag )
   {
      Aux_Message( stderr, "WARNING : all flag criteria are turned off --> no refinement will be " );
      Aux_Message( stderr, "performed !!\n" );
   }

   if ( MAX_LEVEL < 0  ||  MAX_LEVEL > NLEVEL-1 )
      Aux_Message( stderr, "WARNING : MAX_LEVEL (%d) is not within the normal range [0 ... NLEVEL-1] !!\n", 
                   MAX_LEVEL );

   if ( OPT__OVERLAP_MPI )
   {
      Aux_Message( stderr, "WARNING : option \"%s\" is still experimental and is not fully optimized !!\n",
                   "OPT__OVERLAP_MPI" );

#     ifdef OPENMP
      omp_set_nested( true );

      if ( !omp_get_nested() )   
         Aux_Message( stderr, "WARNING : OpenMP nested parallelism is NOT supported for the option \"%s\" !!\n",
                      "OPT__OVERLAP_MPI" );

      omp_set_nested( false );
#     else
      Aux_Message( stderr, "WARNING : OpenMP is NOT turned on for the option \"%s\" !!\n", "OPT__OVERLAP_MPI" );
#     endif

      if ( OPT__TIMING_BARRIER )
         Aux_Message( stderr, "WARNING : option \"%s\" can seriously deteriorate the performance for \"%s\" !!\n",
                      "OPT__TIMING_BARRIER", "OPT__OVERLAP_MPI" );
   } // if ( OPT__OVERLAP_MPI )

   } // if ( MPI_Rank == 0 )



// load balance
// =======================================================================================
#ifdef LOAD_BALANCE

// errors
// ------------------------------
#  ifdef SERIAL
#     error : ERROR : options LOAD_BALANCE and SERIAL should NOT be turned on at the same time !!
#  endif

#  ifdef OOC
#     error : ERROR : options LOAD_BALANCE and OOC should NOT be turned on at the same time !!
#  endif

#  if ( LOAD_BALANCE != HILBERT )
#     error : ERROR : currently DAINO only supports "LOAD_BALANCE == HILBERT" !!
#  endif

// for sending fluid data fixed by coarse-fine fluxes correctly
   if ( OPT__FIXUP_FLUX  &&  Flu_ParaBuf >= PATCH_SIZE )
      Aux_Error( ERROR_INFO, "we must have \"%s\" for \"%s\" in LOAD_BALANCE !!\n",
                 "Flu_ParaBuf < PATCH_SIZE", "OPT__FIXUP_FLUX" );

// ensure that the variable "PaddedCr1D" will not overflow
   const int Padded              = 1<<NLEVEL;
   const int BoxNScale_Padded[3] = { patch->BoxScale[0]/PATCH_SIZE + 2*Padded,
                                     patch->BoxScale[1]/PATCH_SIZE + 2*Padded,
                                     patch->BoxScale[2]/PATCH_SIZE + 2*Padded };
   if (  log2( (double)BoxNScale_Padded[0]*(double)BoxNScale_Padded[1]*(double)BoxNScale_Padded[2] )
         > 8.0*sizeof(long)-1.0  )
   {
      Aux_Message( stderr, "ERROR : PaddedCr1D can overflow !!\n" );
      Aux_Message( stderr, "    --> Please either set NLEVEL to a smaller number or disable LOAD_BALANCE !!\n" );
      MPI_Exit();
   }


// warnings 
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( NX0_TOT[0] != NX0_TOT[1]  ||  NX0_TOT[0] != NX0_TOT[2] )
   {
      Aux_Message( stderr, "WARNING : LOAD_BALANCE has NOT been fully optimized for non-cubic simulation box\n" );
      Aux_Message( stderr, "          (NX0_TOT[0] != NX0_TOT[1] and/or NX0_TOT[0] != NX0_TOT[2]) !!\n" );
   }

   for (int d=0; d<3; d++)
   {
      if ( NX0_TOT[d] & (NX0_TOT[d]-1) )
      {
         Aux_Message( stderr, "WARNING : LOAD_BALANCE has NOT been fully optimized for non-power-of-two " );
         Aux_Message( stderr, "simulation box (NX0_TOT[%d] = %d) !!\n", d, NX0_TOT[d] );
      }
   }

   } // if ( MPI_Rank == 0 )

#else // #ifdef LOAD_BALANCE ... else ...

   if ( OPT__OVERLAP_MPI  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : currently the option \"%s\" is useful only in LOAD_BALANCE !!\n",
                   "OPT__OVERLAP_MPI" );

#endif // #ifdef LOAD_BALANCE ... else ...



// out-of-core
// =======================================================================================
#ifdef OOC

// errors
// ------------------------------
#  ifndef INDIVIDUAL_TIMESTEP 
#     error : ERROR : currently the option OOC does NOT work with the shared time-step scheme !!
#  endif

#  ifdef SERIAL
#     error : ERROR : currently the option OOC does NOT work with the option SERIAL !!
#  endif

#  ifdef LOAD_BALANCE
#     error : ERROR : currently the option OOC does NOT work with the option LOAD_BALANCE !!
#  endif

#  if ( MODEL != HYDRO )
#     error : ERROR : currently the option OOC only works with HYDRO !!
#  endif

   if ( ooc.NRank_X[0]*ooc.NRank_X[1]*ooc.NRank_X[2] != ooc.NRank )
      Aux_Error( ERROR_INFO, "ooc.NRank_X[0]*ooc.NRank_X[1]*ooc.NRank_X[2] != ooc.NRank !!\n" );

   if ( OPT__INIT == INIT_UM )
      Aux_Error( ERROR_INFO, "currently the out-of-core computing does NOT support \"%s\" !!\n",
                 "OPT__INIT == INIT_UM" );

// warnings 
// ------------------------------
#  ifndef OPENMP
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stderr, "WARNING : switching off the OPENMP option in the out-of-core computing will\n" );
      Aux_Message( stderr, "          disable the concurrent execution between data IO and computation !!\n" );
   }
#  endif

#endif // ifdef OOC



// comoving frame
// =======================================================================================
#ifdef COMOVING

// errors
// ------------------------------
   if ( OMEGA_M0 < 0.0  ||  OMEGA_M0 > 1.0 )
      Aux_Error( ERROR_INFO, "incorrect OMEGA_M0 (0.0 <= OMEGA_M0 <= 1.0) !!\n" );

   if ( A_INIT > 1.0  ||  A_INIT < 0.0 )
      Aux_Error( ERROR_INFO, "incorrect A_INIT (0.0 <= A_INIT <= 1.0) !!\n" );

#  if   ( MODEL == HYDRO )
   if ( fabs(GAMMA-5.0/3.0) > 1.e-4 )
      Aux_Error( ERROR_INFO, "GAMMA must be equal to 5.0/3.0 in cosmological simuluations !!\n" );
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif


// warnings 
// ------------------------------
#  ifndef GRAVITY
   if ( MPI_Rank == 0 )   
      Aux_Message( stderr, "WARNING : option \"%s\" is useless if the option \"%s\" is NOT turned on !!\n",
                   "COMOVING", "GRAVITY" );
#  endif

#endif // COMOVING



// fluid solver in all models
// =======================================================================================
  
// errors
// ------------------------------
   if ( Flu_ParaBuf > PATCH_SIZE )  
      Aux_Error( ERROR_INFO, "Flu_ParaBuf (%d) > PATCH_SIZE (%d) !!\n", Flu_ParaBuf, PATCH_SIZE );

   if ( GPU_NSTREAM < 1 )  Aux_Error( ERROR_INFO, "GPU_NSTREAM < 1 !!\n" );

   if ( FLU_GPU_NPGROUP % GPU_NSTREAM != 0 )    Aux_Error( ERROR_INFO, "FLU_GPU_NPGROUP%%GPU_NSTREAM != 0 !!\n" );

   if ( OPT__FLU_INT_SCHEME != INT_CENTRAL  &&  OPT__FLU_INT_SCHEME != INT_MINMOD  &&  
        OPT__FLU_INT_SCHEME != INT_VANLEER  &&  OPT__FLU_INT_SCHEME != INT_CQUAD   &&  
        OPT__FLU_INT_SCHEME != INT_QUAD     &&  OPT__FLU_INT_SCHEME != INT_CQUAR   &&
        OPT__FLU_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__FLU_INT_SCHEME", OPT__FLU_INT_SCHEME );

   if ( OPT__REF_FLU_INT_SCHEME != INT_CENTRAL  &&  OPT__REF_FLU_INT_SCHEME != INT_MINMOD  &&  
        OPT__REF_FLU_INT_SCHEME != INT_VANLEER  &&  OPT__REF_FLU_INT_SCHEME != INT_CQUAD   &&  
        OPT__REF_FLU_INT_SCHEME != INT_QUAD     &&  OPT__REF_FLU_INT_SCHEME != INT_CQUAR   &&
        OPT__REF_FLU_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__REF_FLU_INT_SCHEME", OPT__REF_FLU_INT_SCHEME );

   if ( OPT__FIXUP_FLUX  &&  !patch->WithFlux )
      Aux_Error( ERROR_INFO, "%s is turned on but patch->WithFlux is off !!\n", "OPT__FIXUP_FLUX" ); 

#  if ( NLEVEL > 1 )
   int Trash_RefFlu, NGhost_RefFlu;
   Int_Table( OPT__REF_FLU_INT_SCHEME, Trash_RefFlu, NGhost_RefFlu );
   if ( Flu_ParaBuf < NGhost_RefFlu )
      Aux_Error( ERROR_INFO, "Flu_ParaBuf (%d) < NGhost_RefFlu (%d) --> refinement will fail !!\n", 
                 Flu_ParaBuf, NGhost_RefFlu );
#  endif


// warnings 
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( OPT__CK_FLUX_ALLOCATE  &&  !patch->WithFlux )
      Aux_Message( stderr, "WARNING : option \"%s\" is useless since no flux is required !!\n", 
                   OPT__CK_FLUX_ALLOCATE );

   if ( DT__FLUID < 0.0  ||  DT__FLUID > 1.0 )
      Aux_Message( stderr, "WARNING : DT__FLUID (%14.7e) is not within the normal range [0...1] !!\n", 
                   DT__FLUID );

   } // if ( MPI_Rank == 0 )



// fluid solver in HYDRO
// =======================================================================================
#  if   ( MODEL == HYDRO )

// errors
// ------------------------------
#  if ( NCOMP != 5 )
#     error : ERROR : NCOMP != 5 in HYDRO !!
#  endif

   /*
#  if ( FLU_SCHEME != RTVD  &&  FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU  &&  \
        FLU_SCHEME != WAF )
#     error : ERROR : unsupported hydro scheme in the makefile !!
#  endif
   */

#  if ( defined LR_SCHEME  &&  LR_SCHEME != PLM  &&  LR_SCHEME != PPM )
#     error : ERROR : unsupported data reconstruction scheme (PLM/PPM) !!
#  endif

#  if ( defined RSOLVER  &&  RSOLVER != EXACT  &&  RSOLVER != ROE  &&  RSOLVER != HLLE  &&  RSOLVER != HLLC )
#     error : ERROR : unsupported Riemann solver (EXACT/ROE/HLLE/HLLC) !!
#  endif

#  if ( defined CHECK_INTERMEDIATE  &&  CHECK_INTERMEDIATE != EXACT  &&  CHECK_INTERMEDIATE != HLLE  &&  \
        CHECK_INTERMEDIATE != HLLC )
#     error : ERROR : unsupported option in CHECK_INTERMEDIATE (EXACT/HLLE/HLLC) !!
#  endif

   if ( OPT__CK_NEGATIVE < 0  ||  OPT__CK_NEGATIVE > 3 ) 
      Aux_Error( ERROR_INFO, "unsupported parameter \"%s = %d\" !!\n", "OPT__CK_NEGATIVE", OPT__CK_NEGATIVE );


// warnings 
// ------------------------------
   if ( MPI_Rank == 0 ) {

#  if ( !defined FERMI  &&  CHECK_INTERMEDIATE == EXACT  &&  \
        ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )
#     warning : WARNING : compilation for "CHECK_INTERMEDIATE == EXACT" will take long time in non-Fermi GPUs !!
#  endif

#  if ( FLU_SCHEME == MHM  &&  LR_SCHEME == PPM )
#     warning : WARNING : PPM is not recommended for the MHM scheme !!
      Aux_Message( stderr, "WARNING : PPM is not recommended for the MHM scheme !!\n" );
#  endif

#  if ( FLU_SCHEME == MHM_RP  &&  LR_SCHEME == PPM )
#     warning : WARNING : PPM is not recommended for MHM_RP scheme !!
      Aux_Message( stderr, "WARNING : PPM is not recommended for the MHM_RP scheme !!\n" );
#  endif

#  if ( defined RSOLVER  &&  RSOLVER == EXACT )
#     warning : WARNING : exact RSOLVER is not recommended since the vacuum solution has not been implemented
      Aux_Message( stderr, "WARNING : exact Riemann solver is not recommended since the vacuum solution " );
      Aux_Message( stderr,           "has not been implemented !!\n" );
#  endif

#  if ( defined CHAR_RECONSTRUCTION  &&  defined GRAVITY )
#     warning : WARNING : option "CHAR_RECONSTRUCTION" is less robust and can cause negative density/pressure !!
      Aux_Message( stderr, "WARNING : option \"CHAR_RECONSTRUCTION\" is less robust and can cause negative " );
      Aux_Message( stderr,           "density/pressure !!\n" );
#  endif

#  ifndef ENFORCE_POSITIVE
#     warning : WARNING : option "ENFORCE_POSITIVE" is turned off which can cause negative pressure 
      Aux_Message( stderr, "WARNING : option \"ENFORCE_POSITIVE\" is turned off and can cause negative " );
      Aux_Message( stderr,           "pressure !!\n" );
#  endif

   if ( OPT__FLU_INT_SCHEME == INT_CENTRAL  ||  OPT__FLU_INT_SCHEME == INT_VANLEER )
   {
      Aux_Message( stderr, "WARNING : interpolation scheme (%d) for fluid solver is ", OPT__FLU_INT_SCHEME );
      Aux_Message( stderr, "NOT monotonic and may introduce negative density !!\n" );
   }

   if ( OPT__REF_FLU_INT_SCHEME == INT_CENTRAL  ||  OPT__FLU_INT_SCHEME == INT_VANLEER )
   {
      Aux_Message( stderr, "WARNING : interpolation scheme (%d) for fluid refine is ", OPT__REF_FLU_INT_SCHEME );
      Aux_Message( stderr, "NOT monotonic and may introduce negative density !!\n" );
   }

   if ( !OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : option \"%s\" is disabled in HYDRO !!\n", "OPT__FIXUP_FLUX" );

   if ( !OPT__FIXUP_RESTRICT )
      Aux_Message( stderr, "WARNING : option \"%s\" is disabled in HYDRO !!\n", "OPT__FIXUP_RESTRICT" );

   if ( OPT__CK_FLUX_ALLOCATE  &&  !OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : option %s is useless since %s is off !!\n", 
                   OPT__CK_FLUX_ALLOCATE, OPT__FIXUP_FLUX );

   } // if ( MPI_Rank == 0 )


// check for MHM/MHM_RP/CTU
// ------------------------------
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

// errors
// ------------------------------
#  if ( LR_SCHEME == PPM )
   if ( OPT__LR_LIMITER == EXTPRE )
      Aux_Error( ERROR_INFO, "currently the PPM reconstruction does not support the \"%s\" limiter\n",
                 "extrema-preserving" );
#  endif

   if ( OPT__LR_LIMITER != VANLEER  &&  OPT__LR_LIMITER != GMINMOD  &&  OPT__LR_LIMITER != ALBADA  &&  
        OPT__LR_LIMITER != EXTPRE   &&  OPT__LR_LIMITER != VL_GMINMOD )
      Aux_Error( ERROR_INFO, "unsupported data reconstruction limiter (OPT__LR_IMITER = %d) !!\n", 
                 OPT__LR_LIMITER );


// warnings 
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if (  ( OPT__LR_LIMITER == GMINMOD || OPT__LR_LIMITER == EXTPRE || OPT__LR_LIMITER == VL_GMINMOD )  &&  
         ( MINMOD_COEFF < 1.0 || MINMOD_COEFF > 2.0 )  )
      Aux_Message( stderr, "WARNING : MinMod limiter coefficient (%14.7e) is not within the range [1...2] !!\n",
                   MINMOD_COEFF );

   if ( OPT__LR_LIMITER == EXTPRE  &&  EP_COEFF < 1.0 )
      Aux_Message( stderr, "WARNING : coefficient of the extrema-preserving limiter (%14.7e) < 1.0 !!\n",
                   EP_COEFF );

   } // if ( MPI_Rank == 0 )

#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )


// check for MHM/CTU
// ------------------------------
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == CTU )

#  if ( LR_SCHEME == PLM )
   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE < 3 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE > 3  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE < 2 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 2", "MHM/CTU scheme + PLM reconstruction + non-EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE > 2  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 2", "MHM/CTU scheme + PLM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PLM )

#  if ( LR_SCHEME == PPM )
   if ( FLU_GHOST_SIZE < 3 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PPM reconstruction + non-EXTPRE limiter" );

   if ( FLU_GHOST_SIZE > 3  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PPM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PPM )

#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == CTU )


// check for MHM_RP
// ------------------------------
#  if ( FLU_SCHEME == MHM_RP )

#  if ( LR_SCHEME == PLM )
   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE < 4 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE > 4  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE < 3 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 3", "MHM_RP scheme + PLM reconstruction + non-EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE > 3  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 3", "MHM_RP scheme + PLM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PLM )

#  if ( LR_SCHEME == PPM )
   if ( FLU_GHOST_SIZE < 4 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PPM reconstruction + non-EXTPRE limiter" );

   if ( FLU_GHOST_SIZE > 4  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PPM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PPM )

#  endif // #if ( FLU_SCHEME == MHM_RP )


// check for WAF
// ------------------------------
#  if ( FLU_SCHEME == WAF )

#  if ( RSOLVER == HLLE  ||  RSOLVER == HLLC )
#     error : ERROR : currently the WAF scheme does not support HLLE/HLLC Riemann solvers
#  endif

#  if ( FLU_GHOST_SIZE != 2 )
#     error : ERROR : please set FLU_GHOST_SIZE = 2 for the WAF scheme !!
#  endif

   if ( OPT__WAF_LIMITER != WAF_SUPERBEE  &&  OPT__WAF_LIMITER != WAF_VANLEER  &&  
        OPT__WAF_LIMITER != WAF_ALBADA    &&  OPT__WAF_LIMITER != WAF_MINBEE      )
      Aux_Error( ERROR_INFO, "unsupported WAF flux limiter (%d) !!\n", OPT__WAF_LIMITER );
#  endif // if ( FLU_SCHEME == WAF )


// check for RTVD 
// ------------------------------
#  if ( FLU_SCHEME == RTVD ) 
   
#  if ( FLU_GHOST_SIZE != 3 )
#     error : ERROR : please set FLU_GHOST_SIZE = 3 for the relaxing TVD scheme !!
#  endif

#  endif // if ( FLU_SCHEME == RTVD )



// fluid solver in MHD
// =======================================================================================
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

// errors
// ------------------------------

// warnings 
// ------------------------------
  


// fluid solver in ELBDM
// =======================================================================================
#  elif ( MODEL == ELBDM )

// errors
// ------------------------------
#  if ( NCOMP != 3 )
#     error : ERROR : NCOMP != 3 in ELBDM !!
#  endif

#  if ( FLU_NIN != 2 )
#     error : ERROR : FLU_NIN != 2 in ELBDM !!
#  endif

#  if ( FLU_NOUT != 3 )
#     error : ERROR : FLU_NOUT != 3 in ELBDM !!
#  endif

   if ( ELBDM_MASS   <= 0.0 )    Aux_Error( ERROR_INFO, "%s = %14.7e <= 0.0 !!\n", "ELBDM_MASS", ELBDM_MASS );
   if ( PLANCK_CONST <= 0.0 )    Aux_Error( ERROR_INFO, "%s = %14.7e <= 0.0 !!\n", "PLANCK_CONST", PLANCK_CONST );
   if ( ETA          <= 0.0 )    Aux_Error( ERROR_INFO, "%s = %14.7e <= 0.0 !!\n", "ETA", ETA );

   if ( OPT__INT_PHASE  &&
        OPT__FLU_INT_SCHEME != INT_CQUAD   &&  OPT__FLU_INT_SCHEME != INT_QUAD   &&
        OPT__FLU_INT_SCHEME != INT_CQUAR   &&  OPT__FLU_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported interpolation scheme \"%s = %d\" when OPT__INT_PHASE is on [4-7] !!\n",
                 "OPT__FLU_INT_SCHEME", OPT__FLU_INT_SCHEME );

   if ( OPT__INT_PHASE  &&
        OPT__REF_FLU_INT_SCHEME != INT_CQUAD   &&  OPT__REF_FLU_INT_SCHEME != INT_QUAD   &&
        OPT__REF_FLU_INT_SCHEME != INT_CQUAR   &&  OPT__REF_FLU_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported interpolation scheme \"%s = %d\" when OPT__INT_PHASE is on [4-7] !!\n",
                 "OPT__REF_FLU_INT_SCHEME", OPT__REF_FLU_INT_SCHEME );


// warnings 
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( DT__PHASE < 0.0  ||  DT__PHASE > 1.0 )
      Aux_Message( stderr, "WARNING : DT__PHASE (%14.7e) is not within the normal range [0...1] !!\n", 
                   DT__PHASE );

   if ( OPT__CK_FLUX_ALLOCATE )
      Aux_Message( stderr, "WARNING : option %s is useless in ELBDM !!\n", OPT__CK_FLUX_ALLOCATE );

   if ( OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : option %s is useless in ELBDM !!\n", OPT__FIXUP_FLUX );

   } // if ( MPI_Rank == 0 )

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL



// Poisson and Gravity solvers
// =======================================================================================
#ifdef GRAVITY

// errors
// ------------------------------
#  if ( POT_SCHEME != SOR  &&  POT_SCHEME != MG )
#     error : ERROR : unsupported Poisson solver in the makefile (SOR/MG) !!
#  endif

#  if ( POT_GHOST_SIZE < GRA_GHOST_SIZE )
      #error : ERROR : POT_GHOST_SIZE < GRA_GHOST_SIZE !!
#  endif

#  if ( POT_GHOST_SIZE < 1 )
#     error : ERROR : POT_GHOST_SIZE < 1 !!
#  endif

#  ifdef GPU
#     if ( PATCH_SIZE != 8 )
#        error : ERROR : PATCH_SIZE must == 8 for the GPU Poisson solver !!
#     endif

#     if ( POT_GHOST_SIZE > 5 )
#        error : ERROR : POT_GHOST_SIZE must <= 5 for the GPU Poisson solver !!
#     endif

#     if ( defined FLOAT8  &&  !defined FERMI )
#        error : ERROR : currently the double precision GPU Poisson solver only works with Fermi GPUs !!
#     endif
#  endif // GPU

#  ifndef LOAD_BALANCE
   if ( NX0_TOT[2]%MPI_NRank != 0 )
   {
      Aux_Message( stderr, "ERROR : NX0_TOT[2] %% MPI_NRank != 0 !!\n" );
      Aux_Message( stderr, "--> All MPI processes must have the same number of cells in the z direction for " );
      Aux_Message( stderr,     "the slab decomposition in FFTW 2.1.5 !!\n" );
      MPI_Exit();
   }
#  endif

   if ( Pot_ParaBuf > PATCH_SIZE )  
      Aux_Error( ERROR_INFO, "Pot_ParaBuf (%d) > PATCH_SIZE (%d) !!\n", Pot_ParaBuf, PATCH_SIZE );

   if ( Rho_ParaBuf > PATCH_SIZE )  
      Aux_Error( ERROR_INFO, "Rho_ParaBuf (%d) > PATCH_SIZE (%d) !!\n", Rho_ParaBuf, PATCH_SIZE );

#  if ( POT_SCHEME == SOR )
   if ( SOR_MIN_ITER < 3 )    Aux_Error( ERROR_INFO, "SOR_MIN_ITER < 3 !!\n" );
#  endif

   if ( POT_GPU_NPGROUP % GPU_NSTREAM != 0 )    
      Aux_Error( ERROR_INFO, "POT_GPU_NPGROUP %% GPU_NSTREAM != 0 !!\n" );

   if ( OPT__POT_INT_SCHEME != INT_CENTRAL  &&  OPT__POT_INT_SCHEME != INT_CQUAD  &&  
        OPT__POT_INT_SCHEME != INT_QUAD )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__POT_INT_SCHEME", OPT__POT_INT_SCHEME );

   if ( OPT__RHO_INT_SCHEME != INT_CENTRAL  &&  OPT__RHO_INT_SCHEME != INT_MINMOD  &&  
        OPT__RHO_INT_SCHEME != INT_VANLEER  &&  OPT__RHO_INT_SCHEME != INT_CQUAD   &&  
        OPT__RHO_INT_SCHEME != INT_QUAD     &&  OPT__RHO_INT_SCHEME != INT_CQUAR   &&
        OPT__RHO_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__RHO_INT_SCHEME", OPT__RHO_INT_SCHEME );

   if ( OPT__GRA_INT_SCHEME != INT_CENTRAL  &&  OPT__GRA_INT_SCHEME != INT_MINMOD  &&  
        OPT__GRA_INT_SCHEME != INT_VANLEER  &&  OPT__GRA_INT_SCHEME != INT_CQUAD   &&  
        OPT__GRA_INT_SCHEME != INT_QUAD     &&  OPT__GRA_INT_SCHEME != INT_CQUAR   &&
        OPT__GRA_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__GRA_INT_SCHEME", OPT__GRA_INT_SCHEME );

   if ( OPT__REF_POT_INT_SCHEME != INT_CENTRAL  &&  OPT__REF_POT_INT_SCHEME != INT_MINMOD  &&  
        OPT__REF_POT_INT_SCHEME != INT_VANLEER  &&  OPT__REF_POT_INT_SCHEME != INT_CQUAD   &&  
        OPT__REF_POT_INT_SCHEME != INT_QUAD     &&  OPT__REF_POT_INT_SCHEME != INT_CQUAR   &&
        OPT__REF_POT_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__REF_POT_INT_SCHEME", OPT__REF_POT_INT_SCHEME );

#  if ( NLEVEL > 1 )
   int Trash_RefPot, NGhost_RefPot;
   Int_Table( OPT__REF_POT_INT_SCHEME, Trash_RefPot, NGhost_RefPot );
   if ( Pot_ParaBuf < NGhost_RefPot )
      Aux_Error( ERROR_INFO, "Pot_ParaBuf (%d) < NGhost_RefPot (%d) --> refinement will fail !!\n", 
                 Pot_ParaBuf, NGhost_RefPot );
#  endif


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

#  if ( POT_SCHEME == MG  &&  PATCH_SIZE <= 8 )
   {
      Aux_Message( stderr, "WARNING : multigrid scheme gives lower performance than SOR for " );
      Aux_Message( stderr, "PATCH_SIZE <= 8 and hence is not recommended !!\n" );
   }
#  endif

   if ( DT__GRAVITY < 0.0  ||  DT__GRAVITY > 1.0 )
      Aux_Message( stderr, "WARNING : DT__GRAVITY (%14.7e) is not within the normal range [0...1] !!\n", 
                   DT__GRAVITY );

   } // if ( MPI_Rank == 0 )



// gravity solver in HYDRO
// =======================================================================================
#  if   ( MODEL == HYDRO )

// errors
// ------------------------------
#  if ( GRA_GHOST_SIZE < 1 )
#     error : ERROR : GRA_GHOST_SIZE must >= 1
#  endif

   if ( OPT__GRA_P5_GRADIENT  &&  GRA_GHOST_SIZE == 1 )
      Aux_Error( ERROR_INFO, "option \"%s\" requires \"%s\" !!\n",
                 "OPT__GRA_P5_GRADIENT", "GRA_GHOST_SIZE == 2" );


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( !OPT__GRA_P5_GRADIENT  &&  GRA_GHOST_SIZE == 2 )
   {
      Aux_Message( stderr, "WARNING : \"%s\" is useless if the option \"%s\" is NOT turned on !!\n",
                   "GRA_GHOST_SIZE == 2", "OPT__GRA_P5_GRADIENT" );
   }

   if ( GRA_GHOST_SIZE > 2 )  Aux_Message( stderr, "WARNING : \"GRA_GHOST_SIZE > 2\" !?\n" );

   } // if ( MPI_Rank == 0 )



// gravity solver in MHD
// =======================================================================================
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

// errors
// ------------------------------
  
// warnings 
// ------------------------------



// gravity solver in ELBDM
// =======================================================================================
#  elif ( MODEL == ELBDM )
  
// errors
// ------------------------------
  
// warnings 
// ------------------------------
#  if ( GRA_GHOST_SIZE != 0 )
#     warning : WARNING : GRA_GHOST_SIZE != 0 in ELBDM !!
#  endif


#  else
#  error : unsupported MODEL !!
#  endif // MODEL
  
#endif // GRAVITY


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_Check_Parameter ... done\n" ); 

} // FUNCTION : Aux_Check_Parameter
