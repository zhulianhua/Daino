
#include "DAINO.h"

static void ResetParameter();




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Load_Parameter
// Description :  Load the initial values of simulation parameters from the file "Input__Parameter"
//-------------------------------------------------------------------------------------------------------
void Init_Load_Parameter() 
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_Parameter ... \n" );


   const char FileName[] = "Input__Parameter";

   FILE *File = fopen( FileName, "r" );

   if ( File == NULL )
      Aux_Error( ERROR_INFO, "the file \"%s\" does not exist !!\n", FileName );


   int    temp_int;
   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;


// simulation scale
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &BOX_SIZE,                 string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &NX0_TOT[0],               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &NX0_TOT[1],               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &NX0_TOT[2],               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MPI_NRank,                string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MPI_NRank_X[0],           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MPI_NRank_X[1],           string );
   
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MPI_NRank_X[2],           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OMP_NTHREAD,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &END_T,                    string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%ld%s",  &END_STEP,                 string );

   getline( &input_line, &len, File );


// out-of-core computing
#  ifdef OOC
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &ooc.NRank,                string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &ooc.NRank_X[0],           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &ooc.NRank_X[1],           string );
   
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &ooc.NRank_X[2],           string );

#  else
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
#  endif

   getline( &input_line, &len, File );


// cosmology simulations (COMOVING frame)
#  ifdef COMOVING
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &A_INIT,                   string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OMEGA_M0,                 string );

#  else
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
#  endif

   getline( &input_line, &len, File );


// time-step determination
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &DT__FLUID,                string );

   getline( &input_line, &len, File );
#  ifdef GRAVITY
   sscanf( input_line, "%lf%s",  &DT__GRAVITY,              string );
#  endif
   
   getline( &input_line, &len, File );
#  if ( MODEL == ELBDM )
   sscanf( input_line, "%lf%s",  &DT__PHASE,                string );
#  endif

   getline( &input_line, &len, File );
#  ifdef COMOVING
   sscanf( input_line, "%lf%s",  &DT__MAX_DELTA_A,          string );
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__ADAPTIVE_DT = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__RECORD_DT = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__DT_USER = (bool)temp_int;

   getline( &input_line, &len, File );


// domain refinement
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &REGRID_COUNT,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &FLAG_BUFFER_SIZE,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MAX_LEVEL,                string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_RHO = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_RHO_GRADIENT = (bool)temp_int;

   getline( &input_line, &len, File );
#  if   ( MODEL == HYDRO ) 
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_PRES_GRADIENT = (bool)temp_int;
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif

   getline( &input_line, &len, File );
#  if ( MODEL == ELBDM ) 
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_ENGY_DENSITY = (bool)temp_int;
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OPT__FLAG_LOHNER,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLAG_USER = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OPT__PATCH_COUNT,         string );

   getline( &input_line, &len, File );


// load balance
#  ifdef LOAD_BALANCE
#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &LB_INPUT__WLI_MAX,        string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &LB_INPUT__WLI_MAX,        string );
#  endif

#  else // #ifdef LOAD_BALANCE ... else ...
   getline( &input_line, &len, File );
#  endif // #ifdef LOAD_BALANCE ... else ...

   getline( &input_line, &len, File );


// fluid solvers in HYDRO
#  if   ( MODEL == HYDRO )
#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &GAMMA,                    string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &MINMOD_COEFF,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &EP_COEFF,                 string );

#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &GAMMA,                    string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &MINMOD_COEFF,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &EP_COEFF,                 string );
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__LR_LIMITER = (LR_Limiter_t)temp_int;

   getline( &input_line, &len, File ); // skip two comment lines
   getline( &input_line, &len, File );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__WAF_LIMITER = (WAF_Limiter_t)temp_int;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  else
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
#  endif // MODEL

   getline( &input_line, &len, File );


// ELBDM solvers
#  if ( MODEL == ELBDM )
#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ELBDM_MASS,               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &PLANCK_CONST,             string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &ELBDM_MASS,               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &PLANCK_CONST,             string );
#  endif

#  else // #if ( MODEL == ELBDM ) ... else ...
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
#  endif // #if ( MODEL == ELBDM ) ... else ...

   getline( &input_line, &len, File );


// fluid solvers in both HYDRO/MHD/ELBDM
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &FLU_GPU_NPGROUP,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &GPU_NSTREAM,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FIXUP_FLUX = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FIXUP_RESTRICT = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OVERLAP_MPI = (bool)temp_int;

   getline( &input_line, &len, File );


// self-gravity
#  ifdef GRAVITY

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &NEWTON_G,                 string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &NEWTON_G,                 string );
#  endif

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SOR_OMEGA,                string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &SOR_OMEGA,                string );
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &SOR_MAX_ITER,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &SOR_MIN_ITER,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MG_MAX_ITER,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MG_NPRE_SMOOTH,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &MG_NPOST_SMOOTH,          string );

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &MG_TOLERATED_ERROR,       string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &MG_TOLERATED_ERROR,       string );
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &POT_GPU_NPGROUP,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__GRA_P5_GRADIENT = (bool)temp_int;

#  else // #ifdef GRAVITY ... else ...

   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );

#  endif // #ifdef GRAVITY ... else ...

   getline( &input_line, &len, File );


// initialization
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__INIT = (OptInit_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__RESTART_HEADER = (OptRestartH_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OPT__UM_START_LEVEL,      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OPT__UM_START_NVAR,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__INIT_RESTRICT = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OPT__GPUID_SELECT,        string );

   getline( &input_line, &len, File );



// interpolation schemes
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__INT_TIME = (bool)temp_int;

   getline( &input_line, &len, File );
#  if ( MODEL == ELBDM ) 
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__INT_PHASE = (bool)temp_int;
#  endif

   getline( &input_line, &len, File ); // skip one comment line

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__FLU_INT_SCHEME = (IntScheme_t)temp_int;

#  ifdef GRAVITY
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__POT_INT_SCHEME = (IntScheme_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__RHO_INT_SCHEME = (IntScheme_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__GRA_INT_SCHEME = (IntScheme_t)temp_int;

#  else
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__REF_FLU_INT_SCHEME = (IntScheme_t)temp_int;

   getline( &input_line, &len, File );
#  ifdef GRAVITY
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__REF_POT_INT_SCHEME = (IntScheme_t)temp_int;
#  endif

   getline( &input_line, &len, File );


// data dump
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OPT__OUTPUT_TOTAL,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_PART = (OptOutputPart_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_ERROR = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_BASEPS = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_BASE = (bool)temp_int;

   getline( &input_line, &len, File );
#  ifdef GRAVITY
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_POT = (bool)temp_int;
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__OUTPUT_MODE = (OptOutputMode_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OUTPUT_STEP,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OUTPUT_DT,                string );

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OUTPUT_PART_X,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OUTPUT_PART_Y,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OUTPUT_PART_Z,            string );

#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &OUTPUT_PART_X,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &OUTPUT_PART_Y,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &OUTPUT_PART_Z,            string );
#  endif

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &INIT_DUMPID,              string );

   getline( &input_line, &len, File );


// miscellaneous
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__VERBOSE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__TIMING_BARRIER = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__RECORD_MEMORY = (bool)temp_int;

   getline( &input_line, &len, File );


// simulation checks
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_REFINE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_PROPER_NESTING = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &OPT__CK_CONSERVATION,     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_RESTRICT = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_FINITE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_PATCH_ALLOCATE = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   OPT__CK_FLUX_ALLOCATE = (bool)temp_int;

   getline( &input_line, &len, File );
#  if   ( MODEL == HYDRO )
   sscanf( input_line, "%d%s",   &OPT__CK_NEGATIVE,         string );
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif // MODEL

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &OPT__CK_MEMFREE,          string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &OPT__CK_MEMFREE,          string );
#  endif

   fclose( File );

   if ( input_line != NULL )     free( input_line );


// reset parameters if necessary
   ResetParameter();


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_Parameter ... done\n" );

} // FUNCTION : Init_Load_Parameter



//-------------------------------------------------------------------------------------------------------
// Function    :  ResetParameter 
// Description :  1. Set parameters to their default values if they are not specified in the Input__Parameter 
//                2. Reset parameters which are either unsupported or useless
//-------------------------------------------------------------------------------------------------------
void ResetParameter() 
{

// set parameters to their default values
// ------------------------------------------------------------------------------------------------------
// (1) set the number of OpenMP threads and disable OpenMP nested parallelism by default
#  ifdef OPENMP
   const int OMP_Max_NThread = omp_get_max_threads();

   if ( OMP_NTHREAD <= 0 )  
   {
      OMP_NTHREAD = OMP_Max_NThread;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OMP_NTHREAD", OMP_NTHREAD );
   }

   else if ( OMP_NTHREAD > OMP_Max_NThread   &&  MPI_Rank == 0 )
   {
      Aux_Message( stderr, "WARNING : OMP_NTHREAD (%d) > omp_get_max_threads (%d) !!\n", 
                   OMP_NTHREAD, OMP_Max_NThread );
   }

   omp_set_num_threads( OMP_NTHREAD );
   omp_set_nested( false );

#  else 
   if ( OMP_NTHREAD != 1  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 since \"OPENMP\" is not turned on !!\n", 
                   "OMP_NTHREAD" );

   OMP_NTHREAD = 1;
#  endif


// (2) time-step factors
   if ( DT__FLUID < 0.0 )
   {
#     if   ( MODEL == HYDRO )
#     if   ( FLU_SCHEME == RTVD )
      DT__FLUID = 0.50;
#     elif ( FLU_SCHEME == WAF )
      DT__FLUID = 0.50;
#     elif ( FLU_SCHEME == MHM )
      DT__FLUID = 0.50;
#     elif ( FLU_SCHEME == MHM_RP )
      DT__FLUID = 0.50;
#     elif ( FLU_SCHEME == CTU )
      DT__FLUID = 0.50;
#     else
      DT__FLUID = 0.50;
#     endif

#     elif  ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif  ( MODEL == ELBDM )
#     ifdef GRAVITY
      DT__FLUID = 0.125;   // in order to minimize the numerical dissipation (amplitude error)
#     else
      DT__FLUID = 0.70;    // in order to make the scheme stable
#     endif

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

      if ( MPI_Rank == 0 ) Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n", 
                                        "DT__FLUID", DT__FLUID );
   } // if ( DT__FLUID < 0.0 )

#  ifdef GRAVITY
   if ( DT__GRAVITY < 0.0 )
   {
#     if   ( MODEL == HYDRO )
      DT__GRAVITY = 0.05;

#     elif  ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif  ( MODEL == ELBDM )
      DT__GRAVITY = 0.125;

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

      if ( MPI_Rank == 0 ) Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n", 
                                        "DT__GRAVITY", DT__GRAVITY );
   } // if ( DT__GRAVITY < 0.0 )
#  endif

#  if ( MODEL == ELBDM )
   if ( DT__PHASE < 0.0 )
   {
      DT__PHASE = 0.125;

      if ( MPI_Rank == 0 ) Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n", 
                                        "DT__PHASE", DT__PHASE );
   } // if ( DT__PHASE < 0.0 )
#  endif


// (3) SOR parameters
#  ifdef GRAVITY
#  if   ( POT_SCHEME == SOR )
   Init_Set_Default_SOR_Parameter( SOR_OMEGA, SOR_MAX_ITER, SOR_MIN_ITER );
#  elif ( POT_SCHEME == MG  )
   Init_Set_Default_MG_Parameter( MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, MG_TOLERATED_ERROR );
#  endif
#  endif // GRAVITY


// (4) set the GPU parameters to the default values when using CPU only (please set OMP_NTHREAD in advance)
#  ifndef GPU
   GPU_NSTREAM = 1;
   if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                      "GPU_NSTREAM", GPU_NSTREAM );

   if ( FLU_GPU_NPGROUP <= 0 )  
   {
#     ifdef OPENMP
      FLU_GPU_NPGROUP = OMP_NTHREAD*20;
#     else
      FLU_GPU_NPGROUP = 1;
#     endif

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "FLU_GPU_NPGROUP", FLU_GPU_NPGROUP );
   }

#  ifdef GRAVITY
   if ( POT_GPU_NPGROUP <= 0 )  
   {
#     ifdef OPENMP
      POT_GPU_NPGROUP = OMP_NTHREAD*20;
#     else
      POT_GPU_NPGROUP = 1;
#     endif

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "POT_GPU_NPGROUP", POT_GPU_NPGROUP );
   }
#  endif
#  endif // #ifndef GPU


// (5) grid size in different refinement levels and box size and scale in different directions
   int NX0_Max;
   NX0_Max = ( NX0_TOT[0] > NX0_TOT[1] ) ? NX0_TOT[0] : NX0_TOT[1];
   NX0_Max = ( NX0_TOT[2] > NX0_Max    ) ? NX0_TOT[2] : NX0_Max;

   for (int lv=0; lv<NLEVEL; lv++)     patch->dh[lv] = BOX_SIZE / (double)( NX0_Max*(1<<lv) );

   for (int d=0; d<3; d++)    
   {
      patch->BoxSize [d] = NX0_TOT[d]*patch->dh   [0];
      patch->BoxScale[d] = NX0_TOT[d]*patch->scale[0];
   }


// (6) whether of not to allocate fluxes at the coarse-fine boundaries
#  if   ( MODEL == HYDRO )
   if ( OPT__FIXUP_FLUX )  patch->WithFlux = true;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   patch->WithFlux = false;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// (7) ELBDM parameters
#  if ( MODEL == ELBDM )
#  ifdef COMOVING
   PLANCK_CONST = 1.9198e-26;    // overwrite the number of PLACNK_CONST in the cosmological simulations

   if ( MPI_Rank == 0 )
      Aux_Message( stdout, "NOTE : Planck constant is reset to %13.7e in the comological simulations\n",
                   PLANCK_CONST );
#  endif
   ETA = ELBDM_MASS / PLANCK_CONST;

   if ( MPI_Rank == 0 )    
   {
      Aux_Message( stdout, "NOTE : ETA is set to %13.7e\n", ETA );

#     ifdef COMOVING
      const real JeansK = pow( 6.0*A_INIT*ETA*ETA, 0.25 );
      Aux_Message( stdout, "Note : initial Jean's wavenumber = %13.7e (corresponding to %13.7e Mpc)\n",
                   JeansK, 2.0*M_PI/JeansK );
#     endif
   }
#  endif // #if ( MODEL == ELBDM )


// (8) interpolation scheme
// (8-1) Poisson/Gravity solvers and potential refinement
#  ifdef GRAVITY
   if ( OPT__POT_INT_SCHEME == INT_DEFAULT )
   {
      OPT__POT_INT_SCHEME = INT_QUAD;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__POT_INT_SCHEME", OPT__POT_INT_SCHEME );
   }

   if ( OPT__RHO_INT_SCHEME == INT_DEFAULT )
   {
//    OPT__RHO_INT_SCHEME = INT_MINMOD;
      OPT__RHO_INT_SCHEME = INT_CQUAD;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__RHO_INT_SCHEME", OPT__RHO_INT_SCHEME );
   }

   if ( OPT__GRA_INT_SCHEME == INT_DEFAULT )
   {
      OPT__GRA_INT_SCHEME = INT_QUAD;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__GRA_INT_SCHEME", OPT__GRA_INT_SCHEME );
   }

   if ( OPT__REF_POT_INT_SCHEME == INT_DEFAULT )
   {
      OPT__REF_POT_INT_SCHEME = INT_QUAD;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__REF_POT_INT_SCHEME", OPT__REF_POT_INT_SCHEME );
   }
#  endif // #ifdef GRAVITY

// (8-2) fluid solver and fluid refinement
#  if   ( MODEL == HYDRO )
   if ( OPT__FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__FLU_INT_SCHEME = INT_MINMOD;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__FLU_INT_SCHEME", OPT__FLU_INT_SCHEME );
   }

   if ( OPT__REF_FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__REF_FLU_INT_SCHEME = INT_MINMOD;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__REF_FLU_INT_SCHEME", OPT__REF_FLU_INT_SCHEME );
   }

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   if ( OPT__FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__FLU_INT_SCHEME = INT_CQUAR;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__FLU_INT_SCHEME", OPT__FLU_INT_SCHEME );
   }

   if ( OPT__REF_FLU_INT_SCHEME == INT_DEFAULT )
   {
      OPT__REF_FLU_INT_SCHEME = INT_CQUAR;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "OPT__REF_FLU_INT_SCHEME", OPT__REF_FLU_INT_SCHEME );
   }

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// (9) maximum refinement level
   if ( MAX_LEVEL < 0 )  
   {
      MAX_LEVEL = NLEVEL - 1;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "MAX_LEVEL", MAX_LEVEL );
   }


// (10) refinement frequency and the size of flag buffer
   if ( REGRID_COUNT < 0 )
   {
#     if   ( MODEL == HYDRO )
      REGRID_COUNT = 4;

#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      REGRID_COUNT = 4;

#     else
#     error : ERROR : PLEASE SET THE DEFAULT REGRID_COUNT FOR THE NEW MODEL !!
#     endif // MODEL

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "REGRID_COUNT", REGRID_COUNT );
   }

   if ( FLAG_BUFFER_SIZE < 0 )
   {
#     if   ( MODEL == HYDRO )
      FLAG_BUFFER_SIZE = PATCH_SIZE;

#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      FLAG_BUFFER_SIZE = PATCH_SIZE;

#     else
#     error : ERROR : PLEASE SET THE DEFAULT FLAG_BUFFER_SIZE FOR THE NEW MODEL !!
#     endif // MODEL

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "FLAG_BUFFER_SIZE", FLAG_BUFFER_SIZE );
   }


// (11) initial dump ID
   if ( INIT_DUMPID < 0 )  DumpID = 0;
   else                    DumpID = INIT_DUMPID;


// reset parameters and options which are either unsupported or useless
// ------------------------------------------------------------------------------------------------------
// (1) general
// (1-1) disable "OPT__ADAPTIVE_DT" (not supported yet)
   if ( OPT__ADAPTIVE_DT )
   {
      OPT__ADAPTIVE_DT = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled since it is not supported yet !!\n",
                      "OPT__ADAPTIVE_DT" );
   }

// (1-2) disable "OPT__OVERLAP_MPI" if "OVERLAP_MPI" is NOT turned on in the Makefile
#  ifndef OVERLAP_MPI
   if ( OPT__OVERLAP_MPI ) 
   {
      OPT__OVERLAP_MPI = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled since \"%s\" is off in the Makefile !!\n",
                      "OPT__OVERLAP_MPI", "OVERLAP_MPI" );
   }
#  endif

// (1-3) disable "OPT__CK_FLUX_ALLOCATE" if no flux arrays are going to be allocated
   if ( OPT__CK_FLUX_ALLOCATE  &&  !patch->WithFlux )
   {
      OPT__CK_FLUX_ALLOCATE = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled since no flux is required !!\n",
                      "OPT__CK_FLUX_ALLOCATE" );
   }


// (2) for shared time-step integration
#  ifndef INDIVIDUAL_TIMESTEP
// (2-1) disable "OPT__INT_TIME"
   if ( OPT__INT_TIME )   
   {
      OPT__INT_TIME = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled in the shared time-step scheme !!\n",
                      "OPT__INT_TIME" );
   }

// (2-2) turn off "OPT__OVERLAP_MPI"
   if ( OPT__OVERLAP_MPI ) 
   {
      OPT__OVERLAP_MPI = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled in the shared time-step scheme !!\n",
                      "OPT__OVERLAP_MPI" );
   }
#  endif


// (4) for serial mode
#  ifdef SERIAL
// (4-1) reset the MPI and OOC ranks
   if ( MPI_NRank != 1 )
   {
      MPI_NRank = 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "MPI_NRank" );
   }

   if ( MPI_NRank_X[0] != 1 )
   {
      MPI_NRank_X[0] = 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "MPI_NRank_X[0]");
   }

   if ( MPI_NRank_X[1] != 1 )
   {
      MPI_NRank_X[1] = 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "MPI_NRank_X[1]");
   }

   if ( MPI_NRank_X[2] != 1 )
   {
      MPI_NRank_X[2]= 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "MPI_NRank_X[2]");
   }

#  ifdef OOC
   if ( ooc.NRank != 1 )
   {
      ooc.NRank = 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "ooc.NRank" );
   }

   if ( ooc.NRank_X[0] != 1 )
   {
      ooc.NRank_X[0] = 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "ooc.NRank_X[0]");
   }

   if ( ooc.NRank_X[1] != 1 )
   {
      ooc.NRank_X[1] = 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "ooc.NRank_X[1]");
   }

   if ( ooc.NRank_X[2] != 1 )
   {
      ooc.NRank_X[2]= 1;
      Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 1 for the serial code !!\n", "ooc.NRank_X[2]");
   }
#  endif // ifdef OOC

// (4-2) turn off "OPT__OVERLAP_MPI"
   if ( OPT__OVERLAP_MPI ) 
   {
      OPT__OVERLAP_MPI = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled in the SERIAL mode !!\n",
                      "OPT__OVERLAP_MPI" );
   }
#  endif // ifdef SERIAL


// (5) for Load-balance simulation
#  ifdef LOAD_BALANCE
// (5-1) always turn on "OPT__PATCH_COUNT"
   if ( OPT__PATCH_COUNT <= 0 )
   {
      OPT__PATCH_COUNT = 2;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : parameter \"%s\" is reset to 2 for the LOAD_BALANCE simulation !!\n",
                      "OPT__PATCH_COUNT" );
   }

#  else
// (5-2) turn off "OPT__OVERLAP_MPI" if LOAD_BALANCE is not enabled
   if ( OPT__OVERLAP_MPI ) 
   {
      OPT__OVERLAP_MPI = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled since LOAD_BALANCE is NOT turned on !!\n",
                      "OPT__OVERLAP_MPI" );
   }
#  endif // #ifdef LOAD_BALANCE


// (6) always turn on "OPT__VERBOSE" in the debug mode
#  ifdef DAINO_DEBUG
   if ( !OPT__VERBOSE )
   {
      OPT__VERBOSE = true;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : parameter \"%s\" is turned on automatically in the debug mode !!\n",
                      "OPT__VERBOSE" );
   }
#  endif


// (7) for OpenMP
#  ifndef OPENMP
// (7-1) turn off "OPT__OVERLAP_MPI" if OPENMP is not enabled
   if ( OPT__OVERLAP_MPI ) 
   {
      OPT__OVERLAP_MPI = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled since OPENMP is NOT turned on !!\n",
                      "OPT__OVERLAP_MPI" );
   }
#  endif


// (8) for parallel mode
#  ifndef SERIAL
// (8-1) check the level of MPI thread support
   int MPI_Thread_Status;
   MPI_Query_thread( &MPI_Thread_Status );
   if ( OPT__OVERLAP_MPI  &&  MPI_Thread_Status == MPI_THREAD_SINGLE )
   {
      OPT__OVERLAP_MPI = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled if the level of MPI thread support == %s !!\n",
                      "OPT__OVERLAP_MPI", "MPI_THREAD_SINGLE" );
   }
#  endif


// (9) for different modes
#  if ( MODEL != HYDRO  &&  MODEL != MHD )
// (9-1) operations related to FLUX are useful in HYDRO and MHD only
   if ( OPT__FIXUP_FLUX ) 
   {
      OPT__FIXUP_FLUX = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is only supported in HYDRO/MHD and hence is disabled !!\n",
                      "OPT__FIXUP_FLUX" );
   }

   if ( OPT__CK_FLUX_ALLOCATE ) 
   {
      OPT__CK_FLUX_ALLOCATE = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is only supported in HYDRO/MHD and hence is disabled !!\n",
                      "OPT__CK_FLUX_ALLOCATE" );
   }
#  endif // #if ( MODEL == HYDRO  &&  MODEL != MHD )


// (9-2) turn off refinement criteria and checks related to density if "DENS" is not defined
#  ifndef DENS
   if ( OPT__FLAG_RHO )
   {
      OPT__FLAG_RHO = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled since the variable DENS is not defined !!\n", 
                      "OPT__FLAG_RHO" );
   }

   if ( OPT__FLAG_RHO_GRADIENT )
   {
      OPT__FLAG_RHO_GRADIENT = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled since the variable DENS is not defined !!\n", 
                      "OPT__FLAG_RHO_GRADIENT" );
   }

   if ( OPT__CK_REFINE )
   {
      OPT__CK_REFINE = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is disabled since the variable DENS is not defined !!\n", 
                      "OPT__CK_REFINE" );
   }
#  endif // #ifndef DENS


// (9-3) conservation check is supported only in the models HYDRO, MHD, and ELBDM
#  if ( MODEL != HYDRO  &&  MODEL != MHD  &&  MODEL != ELBDM )
   if ( OPT__CK_CONSERVATION )
   {
      OPT__CK_CONSERVATION = 0;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is not supported in this MODEL and hence is disabled !!\n",
                      "OPT__CK_CONSERVATION" );
   }
#  endif


// (9-4) set OPT__LR_LIMITER and OPT__WAF_LIMITER to NONE if they are useless (in HYDRO)
#  if ( MODEL == HYDRO )
#  if ( FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU )
   if ( OPT__LR_LIMITER != LR_LIMITER_NONE )    
   {
      OPT__LR_LIMITER = LR_LIMITER_NONE;

      if ( MPI_Rank == 0 )    
      {
         Aux_Message( stderr, "WARNING : \"%s\" is useless in the adopted hydro scheme ", "OPT__LR_LIMITER" );
         Aux_Message( stderr, "and has been set to NONE !!\n" );
      }
   }
#  endif

#  if ( FLU_SCHEME != WAF )
   if ( OPT__WAF_LIMITER != WAF_LIMITER_NONE )    
   {
      OPT__WAF_LIMITER = WAF_LIMITER_NONE;

      if ( MPI_Rank == 0 )    
      {
         Aux_Message( stderr, "WARNING : \"%s\" is useless in the adopted hydro scheme ", "OPT__WAF_LIMITER" );
         Aux_Message( stderr, "and has been set to NONE !!\n" );
      }
   }
#  endif
#  endif // #if ( MODEL == HYDRO )


// (10) currently OPT__OUTPUT_BASEPS is not supported if the self-gravity is disabled
#  ifndef GRAVITY 
   if ( OPT__OUTPUT_BASEPS )
   {
      OPT__OUTPUT_BASEPS = false;

      if ( MPI_Rank == 0 )    
         Aux_Message( stderr, "WARNING : option \"%s\" is not supported when GRAVITY is off and hence is disabled !!\n",
                      "OPT__OUTPUT_BASEPS" );
   }
#  endif

} // FUNCTION : ResetParameter
