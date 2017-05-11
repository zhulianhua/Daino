
#include "DAINO.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpData_Total
// Description :  Output all simulation data in the binary form, which can be used as a restart file
//
// Parameter   :  FileName : Name of the output file
//-------------------------------------------------------------------------------------------------------
void Output_DumpData_Total( const char *FileName )
{  

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_Check_Synchronization( Time[0], Time[lv], __FUNCTION__, true );


// check if the targeted file already exists
   if ( MPI_Rank == 0 )
   {
      FILE *File_Check = fopen( FileName, "r" );
      if ( File_Check != NULL )
      {
         Aux_Message( stderr, "WARNING : the file \"%s\" already exists and will be overwritten !!\n", FileName );
         fclose( File_Check );
      }
   }


// get the total number of patches that have no son
   int NDataPatch_Local[NLEVEL] = { 0 };
   int NDataPatch_Total[NLEVEL];

#  ifdef OOC
   for (int lv=0; lv<NLEVEL; lv++)    
   for (int r=0; r<ooc.NRank; r++)     NDataPatch_Local[lv] += ooc.NDataPatch[r][lv];
#  else
   for (int lv=0; lv<NLEVEL; lv++)    
   for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++) 
   {
      if ( patch->ptr[0][lv][PID]->son == -1 )  NDataPatch_Local[lv] ++;
   }
#  endif 

   MPI_Reduce( NDataPatch_Local, NDataPatch_Total, NLEVEL, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );


   FILE *File;

   if ( MPI_Rank == 0 )
   {
      File = fopen( FileName, "wb" );


//    check the size of different data types
//    =================================================================================================
      const int size_bool   = sizeof( bool   );
      const int size_int    = sizeof( int    );
      const int size_uint   = sizeof( uint   );
      const int size_long   = sizeof( long   );
      const int size_ulong  = sizeof( ulong  );
      const int size_real   = sizeof( real   );
      const int size_double = sizeof( double );

      if ( size_int != size_uint )
         Aux_Error( ERROR_INFO, "sizeof(int) = %d != sizeof(uint) = %d !!\n", size_int, size_uint );

      if ( size_long != size_ulong )
         Aux_Error( ERROR_INFO, "sizeof(long) = %d != sizeof(ulong) = %d !!\n", size_long, size_ulong );


//    set the size of the output buffers
  
//    **************************************************************************************************
//    ** if any of the following numbers are modified, please also correct the corresponding values   **
//    ** in the functions "Init_Reload" and "Load_Parameter_After_1200" in the file "Init_Reload.cpp" **
//    **************************************************************************************************

//    =================================================================================================
      const int NBuf_Format    =  256 -  0*size_bool -  5*size_int -  3*size_long -  0*size_real -  0*size_double;
      const int NBuf_Makefile  =  256 - 15*size_bool -  8*size_int -  0*size_long -  0*size_real -  0*size_double;
      const int NBuf_Constant  =  256 -  6*size_bool - 11*size_int -  0*size_long -  2*size_real -  0*size_double;
      const int NBuf_Parameter = 1024 - 18*size_bool - 35*size_int -  1*size_long - 12*size_real -  8*size_double;
      const int NBuf_Info      = 1024 -  0*size_bool - (1+3*NLEVEL)*size_int - 2*size_long 
                                      -  0*size_real - (1+NLEVEL)*size_double; // one size_long is for CheckCode

      if ( NBuf_Format    < 0 )  Aux_Error( ERROR_INFO, "%s = %d < 0 !!\n", "NBuf_Format",   NBuf_Format   );
      if ( NBuf_Makefile  < 0 )  Aux_Error( ERROR_INFO, "%s = %d < 0 !!\n", "NBuf_Makefile", NBuf_Makefile );
      if ( NBuf_Constant  < 0 )  Aux_Error( ERROR_INFO, "%s = %d < 0 !!\n", "NBuf_Constant", NBuf_Constant );
      if ( NBuf_Parameter < 0 )  Aux_Error( ERROR_INFO, "%s = %d < 0 !!\n", "NBuf_Format",   NBuf_Format   );
      if ( NBuf_Info      < 0 )  Aux_Error( ERROR_INFO, "%s = %d < 0 !!\n", "NBuf_Info",     NBuf_Info     );

      int MaxNBuf = -1;
      MaxNBuf = ( MaxNBuf >= NBuf_Format    ) ? MaxNBuf : NBuf_Format;
      MaxNBuf = ( MaxNBuf >= NBuf_Makefile  ) ? MaxNBuf : NBuf_Makefile;
      MaxNBuf = ( MaxNBuf >= NBuf_Constant  ) ? MaxNBuf : NBuf_Constant;
      MaxNBuf = ( MaxNBuf >= NBuf_Parameter ) ? MaxNBuf : NBuf_Parameter;
      MaxNBuf = ( MaxNBuf >= NBuf_Info      ) ? MaxNBuf : NBuf_Info;

      char *OutputBuf = new char [MaxNBuf];
      for (int t=0; t<MaxNBuf; t++)    OutputBuf[t] = 'b';


//    a. output the information of data format
//    =================================================================================================
      const long FormatVersion = 1201;
      const long HeaderSize    = 2048;          // it must be larger than output a+b+c+d
      const long CheckCode     = 123456789;

      fwrite( &FormatVersion,             sizeof(long),                    1,             File );
      fwrite( &HeaderSize,                sizeof(long),                    1,             File );
      fwrite( &CheckCode,                 sizeof(long),                    1,             File );
      fwrite( &size_bool,                 sizeof(int),                     1,             File );
      fwrite( &size_int,                  sizeof(int),                     1,             File );
      fwrite( &size_long,                 sizeof(int),                     1,             File );
      fwrite( &size_real,                 sizeof(int),                     1,             File );
      fwrite( &size_double,               sizeof(int),                     1,             File );

//    buffer space reserved for future usuage
      fwrite( OutputBuf,                  sizeof(char),          NBuf_Format,             File );


//    b. output the simulation options and parameters defined in the Makefile
//    =================================================================================================
#     ifdef MODEL
      const int  model               = MODEL;
#     else
      const int  model               = NULL_INT;
#     endif

#     ifdef GRAVITY
      const bool gravity             = true;
#     else
      const bool gravity             = false;
#     endif

#     ifdef POT_SCHEME
      const int  pot_scheme          = POT_SCHEME;
#     else
      const int  pot_scheme          = NULL_INT;
#     endif

#     ifdef INDIVIDUAL_TIMESTEP
      const bool individual_timestep = true;
#     else
      const bool individual_timestep = false;
#     endif

#     ifdef COMOVING
      const bool comoving            = true;
#     else
      const bool comoving            = false;
#     endif

#     ifdef FLU_SCHEME
      const int  flu_scheme          = FLU_SCHEME;
#     else
      const int  flu_scheme          = NULL_INT;
#     endif

#     ifdef LR_SCHEME
      const int  lr_scheme           = LR_SCHEME;
#     else
      const int  lr_scheme           = NULL_INT;
#     endif

#     ifdef RSOLVER
      const int  rsolver             = RSOLVER;
#     else
      const int  rsolver             = NULL_INT;
#     endif

#     ifdef GPU
      const bool gpu                 = true;
#     else
      const bool gpu                 = false;
#     endif

#     ifdef DAINO_OPTIMIZATION
      const bool DAINO_optimization  = true;
#     else
      const bool DAINO_optimization  = false;
#     endif

#     ifdef DAINO_DEBUG
      const bool DAINO_debug         = true;
#     else
      const bool DAINO_debug         = false;
#     endif

#     ifdef TIMING
      const bool timing              = true;
#     else
      const bool timing              = false;
#     endif

#     ifdef TIMING_SOLVER
      const bool timing_solver       = true;
#     else
      const bool timing_solver       = false;
#     endif

#     ifdef INTEL
      const bool intel               = true;
#     else
      const bool intel               = false;
#     endif

#     ifdef FLOAT8
      const bool float8              = true;
#     else
      const bool float8              = false;
#     endif

#     ifdef SERIAL
      const bool serial              = true;
#     else
      const bool serial              = false;
#     endif

#     ifdef OOC
      const bool ooc                 = true;
#     else
      const bool ooc                 = false;
#     endif

#     ifdef LOAD_BALANCE
      const int load_balance         = LOAD_BALANCE;
#     else
      const int load_balance         = NULL_INT;
#     endif

#     ifdef OVERLAP_MPI
      const bool overlap_mpi         = true;
#     else
      const bool overlap_mpi         = false;
#     endif

#     ifdef OPENMP
      const bool openmp              = true;
#     else
      const bool openmp              = false;
#     endif

#     ifdef FERMI
      const bool fermi               = true;
#     else
      const bool fermi               = false;
#     endif

      const int nlevel               = NLEVEL;
      const int max_patch            = MAX_PATCH;

      fwrite( &model,                     sizeof(int),                     1,             File );
      fwrite( &gravity,                   sizeof(bool),                    1,             File );
      fwrite( &pot_scheme,                sizeof(int),                     1,             File );
      fwrite( &individual_timestep,       sizeof(bool),                    1,             File );
      fwrite( &comoving,                  sizeof(bool),                    1,             File );
      fwrite( &flu_scheme,                sizeof(int),                     1,             File );
      fwrite( &lr_scheme,                 sizeof(int),                     1,             File );
      fwrite( &rsolver,                   sizeof(int),                     1,             File );
      fwrite( &gpu,                       sizeof(bool),                    1,             File );
      fwrite( &DAINO_optimization,        sizeof(bool),                    1,             File );
      fwrite( &DAINO_debug,               sizeof(bool),                    1,             File );
      fwrite( &timing,                    sizeof(bool),                    1,             File );
      fwrite( &timing_solver,             sizeof(bool),                    1,             File );
      fwrite( &intel,                     sizeof(bool),                    1,             File );
      fwrite( &float8,                    sizeof(bool),                    1,             File );
      fwrite( &serial,                    sizeof(bool),                    1,             File );
      fwrite( &ooc,                       sizeof(bool),                    1,             File );
      fwrite( &load_balance,              sizeof(int),                     1,             File );
      fwrite( &overlap_mpi,               sizeof(bool),                    1,             File );
      fwrite( &openmp,                    sizeof(bool),                    1,             File );
      fwrite( &fermi,                     sizeof(bool),                    1,             File );
      fwrite( &nlevel,                    sizeof(int),                     1,             File );
      fwrite( &max_patch,                 sizeof(int),                     1,             File );

//    buffer space reserved for future usuage
      fwrite( OutputBuf,                  sizeof(char),       NBuf_Makefile,              File );


//    c. output the symbolic constants defined in "Macro.h, CUPOT.h, and CUFLU.h"
//    =================================================================================================
      const int  ncomp                 = NCOMP;
      const int  patch_size            = PATCH_SIZE;
      const real min_value             = MIN_VALUE;
      const int  flu_ghost_size        = FLU_GHOST_SIZE;

#     ifdef GRAVITY
      const int  pot_ghost_size        = POT_GHOST_SIZE; 
      const int  gra_ghost_size        = GRA_GHOST_SIZE;
#     else
      const int  pot_ghost_size        = NULL_INT;
      const int  gra_ghost_size        = NULL_INT;
#     endif

#     ifdef ENFORCE_POSITIVE
      const bool enforce_positive      = true;
#     else
      const bool enforce_positive      = false;
#     endif

#     ifdef CHAR_RECONSTRUCTION
      const bool char_reconstruction   = true;
#     else
      const bool char_reconstruction   = false;
#     endif

#     ifdef CHECK_INTERMEDIATE
      const int  check_intermediate    = CHECK_INTERMEDIATE;;
#     else
      const int  check_intermediate    = NULL_INT;
#     endif

#     ifdef HLL_NO_REF_STATE
      const bool hll_no_ref_state      = true;
#     else
      const bool hll_no_ref_state      = false;
#     endif

#     ifdef HLL_INCLUDE_ALL_WAVES
      const bool hll_include_all_waves = true;
#     else
      const bool hll_include_all_waves = false;
#     endif

#     ifdef WAF_DISSIPATE
      const bool waf_dissipate         = true;
#     else
      const bool waf_dissipate         = false;
#     endif

#     ifdef MAX_ERROR
      const real max_error             = MAX_ERROR;
#     else
      const real max_error             = NULL_REAL;
#     endif

      const int  flu_block_size_x      = FLU_BLOCK_SIZE_X;
      const int  flu_block_size_y      = FLU_BLOCK_SIZE_Y;

#     ifdef USE_PSOLVER_10TO14
      const bool use_psolver_10to14    = true;
#     else
      const bool use_psolver_10to14    = false;
#     endif

#     ifdef POT_BLOCK_SIZE_X
      const int pot_block_size_x       = POT_BLOCK_SIZE_X;
#     else
      const int pot_block_size_x       = NULL_INT;
#     endif

#     ifdef POT_BLOCK_SIZE_Z
      const int pot_block_size_z       = POT_BLOCK_SIZE_Z;
#     else
      const int pot_block_size_z       = NULL_INT;
#     endif

#     ifdef GRA_BLOCK_SIZE_Z
      const int gra_block_size_z       = GRA_BLOCK_SIZE_Z;
#     else
      const int gra_block_size_z       = NULL_INT;
#     endif

      fwrite( &ncomp,                     sizeof(int),                     1,             File );
      fwrite( &patch_size,                sizeof(int),                     1,             File );
      fwrite( &min_value,                 sizeof(real),                    1,             File );
      fwrite( &flu_ghost_size,            sizeof(int),                     1,             File );
      fwrite( &pot_ghost_size,            sizeof(int),                     1,             File );
      fwrite( &gra_ghost_size,            sizeof(int),                     1,             File );
      fwrite( &enforce_positive,          sizeof(bool),                    1,             File );
      fwrite( &char_reconstruction,       sizeof(bool),                    1,             File );
      fwrite( &check_intermediate,        sizeof(int),                     1,             File );
      fwrite( &hll_no_ref_state,          sizeof(bool),                    1,             File );
      fwrite( &hll_include_all_waves,     sizeof(bool),                    1,             File );
      fwrite( &waf_dissipate,             sizeof(bool),                    1,             File );
      fwrite( &max_error,                 sizeof(real),                    1,             File );
      fwrite( &flu_block_size_x,          sizeof(int),                     1,             File );
      fwrite( &flu_block_size_y,          sizeof(int),                     1,             File );
      fwrite( &use_psolver_10to14,        sizeof(bool),                    1,             File );
      fwrite( &pot_block_size_x,          sizeof(int),                     1,             File );
      fwrite( &pot_block_size_z,          sizeof(int),                     1,             File );
      fwrite( &gra_block_size_z,          sizeof(int),                     1,             File );

//    buffer space reserved for future usuage
      fwrite( OutputBuf,                  sizeof(char),       NBuf_Constant,              File );


//    d. output the simulation parameters recorded in the file "Input__Parameter"
//    *** we always typecast "enum" to "int" when outputting data to ensure the data size ***
//    =================================================================================================
      const int    DAINO_nrank             = DAINO_NRANK;
      const int    DAINO_nrank_x[3]        = { DAINO_NRANK_X(0), DAINO_NRANK_X(1), DAINO_NRANK_X(2) };    
      const int    opt__output_part        = (int)OPT__OUTPUT_PART;
      const int    opt__output_mode        = (int)OPT__OUTPUT_MODE;
      const int    opt__flu_int_scheme     = (int)OPT__FLU_INT_SCHEME;
      const int    opt__ref_flu_int_scheme = (int)OPT__REF_FLU_INT_SCHEME;

#     ifdef OOC
      const int    ooc_nrank               = ooc.NRank;
      const int    ooc_nrank_x[3]          = { ooc.NRank_X[0], ooc.NRank_X[1], ooc.NRank_X[2] };
#     else
      const int    ooc_nrank               = NULL_INT;
      const int    ooc_nrank_x[3]          = { NULL_INT, NULL_INT, NULL_INT };
#     endif

#     ifndef COMOVING
      const double OMEGA_M0                = NULL_REAL;
      const double DT__MAX_DELTA_A         = NULL_REAL;
#     endif

#     ifdef GRAVITY
      const int    opt__pot_int_scheme     = (int)OPT__POT_INT_SCHEME;
      const int    opt__rho_int_scheme     = (int)OPT__RHO_INT_SCHEME;
      const int    opt__gra_int_scheme     = (int)OPT__GRA_INT_SCHEME;
      const int    opt__ref_pot_int_scheme = (int)OPT__REF_POT_INT_SCHEME;
#     else
      const double DT__GRAVITY             = NULL_REAL;
      const real   NEWTON_G                = NULL_REAL;
      const int    POT_GPU_NPGROUP         = NULL_INT;
      const bool   OPT__OUTPUT_POT         = NULL_BOOL;
      const bool   OPT__GRA_P5_GRADIENT    = NULL_BOOL;
      const real   SOR_OMEGA               = NULL_REAL;
      const int    SOR_MAX_ITER            = NULL_INT;
      const int    SOR_MIN_ITER            = NULL_INT;
      const real   MG_TOLERATED_ERROR      = NULL_REAL;
      const int    MG_MAX_ITER             = NULL_INT;
      const int    MG_NPRE_SMOOTH          = NULL_INT;
      const int    MG_NPOST_SMOOTH         = NULL_INT;
      const int    opt__pot_int_scheme     = NULL_INT;
      const int    opt__rho_int_scheme     = NULL_INT;
      const int    opt__gra_int_scheme     = NULL_INT;
      const int    opt__ref_pot_int_scheme = NULL_INT;
#     endif // #ifdef GRAVITY

#     ifdef LOAD_BALANCE
      const real   lb_wli_max              = patch->LB->WLI_Max;
#     else
      const real   lb_wli_max              = NULL_REAL;
#     endif

#     if ( MODEL == HYDRO )
      const int    opt__lr_limiter         = (int)OPT__LR_LIMITER;
      const int    opt__waf_limiter        = (int)OPT__WAF_LIMITER;
#     else
#     if ( MODEL == MHD )
#     warning : WAIT MHD !!!
#     endif
      const bool   OPT__FLAG_PRES_GRADIENT = NULL_BOOL;
      const real   GAMMA                   = NULL_REAL;
      const real   MINMOD_COEFF            = NULL_REAL;
      const real   EP_COEFF                = NULL_REAL;
      const int    opt__lr_limiter         = NULL_INT;
      const int    opt__waf_limiter        = NULL_INT;
#     endif

#     if ( MODEL != ELBDM )
      const double DT__PHASE               = NULL_REAL;
      const bool   OPT__FLAG_ENGY_DENSITY  = NULL_BOOL;
      const bool   OPT__INT_PHASE          = NULL_BOOL;
      const real   ELBDM_MASS              = NULL_REAL;
      const real   PLANCK_CONST            = NULL_REAL;
#     endif

      fwrite( &BOX_SIZE,                  sizeof(double),                  1,             File );
      fwrite(  NX0_TOT,                   sizeof(int),                     3,             File );
      fwrite( &DAINO_nrank,               sizeof(int),                     1,             File );
      fwrite(  DAINO_nrank_x,             sizeof(int),                     3,             File );
      fwrite( &OMP_NTHREAD,               sizeof(int),                     1,             File );
      fwrite( &END_T,                     sizeof(double),                  1,             File );
      fwrite( &END_STEP,                  sizeof(long),                    1,             File );
      fwrite( &ooc_nrank,                 sizeof(int),                     1,             File );
      fwrite(  ooc_nrank_x,               sizeof(int),                     3,             File );
      fwrite( &OMEGA_M0,                  sizeof(double),                  1,             File );
      fwrite( &DT__FLUID,                 sizeof(double),                  1,             File );
      fwrite( &DT__GRAVITY,               sizeof(double),                  1,             File );
      fwrite( &DT__PHASE,                 sizeof(double),                  1,             File );
      fwrite( &DT__MAX_DELTA_A,           sizeof(double),                  1,             File );
      fwrite( &OPT__ADAPTIVE_DT,          sizeof(bool),                    1,             File );
      fwrite( &OPT__DT_USER,              sizeof(bool),                    1,             File );
      fwrite( &REGRID_COUNT,              sizeof(int),                     1,             File );
      fwrite( &FLAG_BUFFER_SIZE,          sizeof(int),                     1,             File );
      fwrite( &MAX_LEVEL,                 sizeof(int),                     1,             File );
      fwrite( &OPT__FLAG_RHO,             sizeof(bool),                    1,             File );
      fwrite( &OPT__FLAG_RHO_GRADIENT,    sizeof(bool),                    1,             File );
      fwrite( &OPT__FLAG_PRES_GRADIENT,   sizeof(bool),                    1,             File );
      fwrite( &OPT__FLAG_ENGY_DENSITY,    sizeof(bool),                    1,             File );
      fwrite( &OPT__FLAG_USER,            sizeof(bool),                    1,             File );
      fwrite( &lb_wli_max,                sizeof(real),                    1,             File );
      fwrite( &GAMMA,                     sizeof(real),                    1,             File );
      fwrite( &MINMOD_COEFF,              sizeof(real),                    1,             File );
      fwrite( &EP_COEFF,                  sizeof(real),                    1,             File );
      fwrite( &opt__lr_limiter,           sizeof(int),                     1,             File );
      fwrite( &opt__waf_limiter,          sizeof(int),                     1,             File );
      fwrite( &ELBDM_MASS,                sizeof(real),                    1,             File );
      fwrite( &PLANCK_CONST,              sizeof(real),                    1,             File );
      fwrite( &FLU_GPU_NPGROUP,           sizeof(int),                     1,             File );
      fwrite( &GPU_NSTREAM,               sizeof(int),                     1,             File );
      fwrite( &OPT__FIXUP_FLUX,           sizeof(bool),                    1,             File );
      fwrite( &OPT__FIXUP_RESTRICT,       sizeof(bool),                    1,             File );
      fwrite( &OPT__OVERLAP_MPI,          sizeof(bool),                    1,             File );
      fwrite( &NEWTON_G,                  sizeof(real),                    1,             File );
      fwrite( &SOR_OMEGA,                 sizeof(real),                    1,             File );
      fwrite( &SOR_MAX_ITER,              sizeof(int),                     1,             File );
      fwrite( &SOR_MIN_ITER,              sizeof(int),                     1,             File );
      fwrite( &MG_MAX_ITER,               sizeof(int),                     1,             File );
      fwrite( &MG_NPRE_SMOOTH,            sizeof(int),                     1,             File );
      fwrite( &MG_NPOST_SMOOTH,           sizeof(int),                     1,             File );
      fwrite( &MG_TOLERATED_ERROR,        sizeof(real),                    1,             File );
      fwrite( &POT_GPU_NPGROUP,           sizeof(int),                     1,             File );
      fwrite( &OPT__GRA_P5_GRADIENT,      sizeof(bool),                    1,             File );
      fwrite( &OPT__INT_TIME,             sizeof(bool),                    1,             File );
      fwrite( &OPT__INT_PHASE,            sizeof(bool),                    1,             File );
      fwrite( &opt__flu_int_scheme,       sizeof(int),                     1,             File );
      fwrite( &opt__pot_int_scheme,       sizeof(int),                     1,             File );
      fwrite( &opt__rho_int_scheme,       sizeof(int),                     1,             File );
      fwrite( &opt__gra_int_scheme,       sizeof(int),                     1,             File );
      fwrite( &opt__ref_flu_int_scheme,   sizeof(int),                     1,             File );
      fwrite( &opt__ref_pot_int_scheme,   sizeof(int),                     1,             File );
      fwrite( &OPT__OUTPUT_TOTAL,         sizeof(int),                     1,             File );
      fwrite( &opt__output_part,          sizeof(int),                     1,             File );
      fwrite( &OPT__OUTPUT_ERROR,         sizeof(bool),                    1,             File );
      fwrite( &OPT__OUTPUT_BASE,          sizeof(bool),                    1,             File );
      fwrite( &OPT__OUTPUT_POT,           sizeof(bool),                    1,             File );
      fwrite( &opt__output_mode,          sizeof(int),                     1,             File );
      fwrite( &OUTPUT_STEP,               sizeof(int),                     1,             File );
      fwrite( &OUTPUT_DT,                 sizeof(double),                  1,             File );
      fwrite( &OUTPUT_PART_X,             sizeof(real),                    1,             File );
      fwrite( &OUTPUT_PART_Y,             sizeof(real),                    1,             File );
      fwrite( &OUTPUT_PART_Z,             sizeof(real),                    1,             File );
      fwrite( &OPT__TIMING_BARRIER,       sizeof(bool),                    1,             File );
      fwrite( &OPT__OUTPUT_BASEPS,        sizeof(bool),                    1,             File );

//    buffer space reserved for future usuage
      fwrite( OutputBuf,                  sizeof(char),       NBuf_Parameter,             File );


//    set the file position indicator to "HeaderSize" and record the check code
//    =================================================================================================
      fseek( File, HeaderSize, SEEK_SET );
      fwrite( &CheckCode,                 sizeof(long),                    1,             File );


//    e. output the simulation information
//    =================================================================================================
#     ifndef GRAVITY
      const double AveDensity = NULL_REAL;
#     endif

      fwrite( &DumpID,                    sizeof(int),                     1,             File );
      fwrite( Time,                       sizeof(double),             NLEVEL,             File );
      fwrite( &Step,                      sizeof(long),                    1,             File );
      fwrite( NPatchTotal,                sizeof(int),                NLEVEL,             File );
      fwrite( NDataPatch_Total,           sizeof(int),                NLEVEL,             File );
      fwrite( AdvanceCounter,             sizeof(uint),               NLEVEL,             File );
      fwrite( &AveDensity,                sizeof(double),                  1,             File );

//    buffer space reserved for future usuage
      fwrite( OutputBuf,                  sizeof(char),            NBuf_Info,             File );


      delete [] OutputBuf;

      fclose( File );

   } // if ( MPI_Rank == 0 )


// f. output the simulation data
// =================================================================================================
// array for re-ordering the fluid data from "xyzv" to "vxyz"
   real (*InvData_Flu)[PATCH_SIZE][PATCH_SIZE][NCOMP] = NULL;
   if ( OPT__OUTPUT_TOTAL == 1 )    InvData_Flu = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE][NCOMP];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
      {
         if ( MPI_Rank == TargetMPIRank )
         {
            File = fopen( FileName, "ab" );

#ifndef OOC

            for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
            {
//             f1. output the patch information 
//             (the father <-> son information will be re-constructed during the restart)
               fwrite(  patch->ptr[0][lv][PID]->corner, sizeof(int), 3, File );
               fwrite( &patch->ptr[0][lv][PID]->son,    sizeof(int), 1, File );


//             f2. output the patch data only if it has no son
               if ( patch->ptr[0][lv][PID]->son == -1 )
               {

//                f2-1. output the fluid variables
                  if ( OPT__OUTPUT_TOTAL == 1 )
                  {
                     for (int v=0; v<NCOMP; v++)
                     for (int k=0; k<PATCH_SIZE; k++)
                     for (int j=0; j<PATCH_SIZE; j++)
                     for (int i=0; i<PATCH_SIZE; i++)    
                        InvData_Flu[k][j][i][v] = patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

                     fwrite( InvData_Flu, sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NCOMP, File );
                  }
                  else
                     fwrite( patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid, sizeof(real), 
                             PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NCOMP, File );

#                 ifdef GRAVITY
//                f2-2. output the gravitational potential
                  if ( OPT__OUTPUT_POT )
                     fwrite( patch->ptr[ patch->PotSg[lv] ][lv][PID]->pot,   sizeof(real), 
                             PATCH_SIZE*PATCH_SIZE*PATCH_SIZE,       File );
#                 endif 

               } // if ( patch->ptr[0][lv][PID]->son == -1 )
            } // for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)

#else // OOC

            OOC_Output_DumpData_Total( lv, File, InvData_Flu );

#endif

            fclose( File );

         } // if ( MPI_Rank == TargetMPIRank )

         MPI_Barrier( MPI_COMM_WORLD );

      } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   } // for (int lv=0; lv<NLEVEL; lv++)


   if ( OPT__OUTPUT_TOTAL == 1 )    delete [] InvData_Flu;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_DumpData_Total


