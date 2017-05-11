
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init
// Description :  Initialize DAINO
//-------------------------------------------------------------------------------------------------------
void Init_DAINO( int *argc, char ***argv )
{

// initialize MPI
   Init_MPI( argc, argv );


// initialize the AMR structure
   patch = new AMR_t;


// load parameters from the file "Input__Parameter"
   Init_Load_Parameter();


// initialize parameters for the parallelization
   Init_Parallelization();


// initialize parameters for the out-of-core computing
#  ifdef OOC
   if ( ooc.NRank == 1 )   OOC_patch = patch;
   else                    OOC_patch = new AMR_t;

   Init_OOC();
#  endif


#  ifdef GRAVITY
// initialize FFTW
   Init_FFTW();
#  endif


// initialize test problem parameters
   Init_TestProb();


// set the GPU ID and several GPU parameters
#  ifdef GPU
#  ifndef GRAVITY
   int POT_GPU_NPGROUP = NULL_INT;
#  endif
   CUAPI_SetDevice( OPT__GPUID_SELECT );

   CUAPI_Set_Default_GPU_Parameter( GPU_NSTREAM, FLU_GPU_NPGROUP, POT_GPU_NPGROUP );
#  endif


// initialize the array recording the previous physical time as an arbitrary "negative" number
   for (int lv=0; lv<NLEVEL; lv++)     Time_Prev[lv] = -9999.9;


// verify the input parameters
   Aux_Check_Parameter();


// initialize the timer function
#  ifdef TIMING
   Aux_CreateTimer();
#  endif


#  ifdef COMOVING
// initialize the scale factor for the cosmological simulation (it will be overwritten during restart)
   for (int lv=0; lv<NLEVEL; lv++)     Time[lv] = A_INIT;


// reset the gravitational constant
#  ifdef GRAVITY
   NEWTON_G = real( 3.0*OMEGA_M0/8.0/M_PI );

   if ( MPI_Rank == 0 )   
      Aux_Message( stdout, "NOTE : gravitational constant is reset to %13.7e in the comological simulations\n",
                   NEWTON_G );
#  endif
#  endif // #ifdef COMOVING
   

// load the tables of the flag criteria from the input files "Input__FLAG_XXX"
   Init_Load_FlagCriteria();


// load the dump table from the input file "Input__DumpTable"
   if (  ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_ERROR || OPT__OUTPUT_BASEPS )  &&  
         OPT__OUTPUT_MODE == OUTPUT_USE_TABLE  )
   Init_Load_DumpTable();


// allocate memory for several global arrays
   Init_MemAllocate();


// create the initial condition
   switch ( OPT__INIT )
   {
      case INIT_STARTOVER:    Init_StartOver();    break;

      case INIT_RESTART :     Init_Reload();       break;

      case INIT_UM :          Init_UM();           break;

      default : Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "OPT__INIT", OPT__INIT );
   }


// get the total number of patches in all ranks
   for (int lv=0; lv<NLEVEL; lv++)     Mis_GetTotalPatchNumber( lv );


// improve load balance    
#  ifdef LOAD_BALANCE
   LB_Init_LoadBalance( OPT__INIT == INIT_RESTART );

// fill up the data for patches that are not leaf patches (for RESTART only)
   if ( OPT__INIT == INIT_RESTART )
   for (int lv=NLEVEL-2; lv>=0; lv--)
   {
      Flu_Restrict( lv, patch->FluSg[lv+1], patch->FluSg[lv], NULL_INT, NULL_INT, _FLU );

#     ifdef LOAD_BALANCE 
      LB_GetBufferData( lv, patch->FluSg[lv], NULL_INT, DATA_RESTRICT, _FLU, NULL_INT );
#     endif

      Buf_GetBufferData( lv, patch->FluSg[lv], NULL_INT, DATA_GENERAL, _FLU, Flu_ParaBuf, USELB_YES );
   }
#  endif // #ifdef LOAD_BALANCE


// sort all patches for the out-of-core computing
#  ifdef OOC
   OOC_Init_SortPatch();
#  endif


#  ifdef GRAVITY
// evaluate the average density if it is not set yet for the periodic Poisson solver
   if ( AveDensity <= 0.0 )   Poi_GetAverageDensity();


// evaluate the gravitational potential
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating gravitational potential" ); 

   for (int lv=0; lv<NLEVEL; lv++)     
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d ... ", lv );

#ifndef OOC

      Buf_GetBufferData( lv, patch->FluSg[lv], NULL_INT, DATA_GENERAL, _DENS, Rho_ParaBuf, USELB_YES );

      Gra_AdvanceDt( lv, Time[lv], 0.0, patch->PotSg[lv], false, false, false );

      if ( lv > 0 )  
      Buf_GetBufferData( lv, NULL_INT, patch->PotSg[lv], POT_FOR_POISSON, _POTE, Pot_ParaBuf, USELB_YES );

#else // OOC

      OOC_Init_DAINO_GetPot( lv );

#endif

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating gravitational potential" ); 
#  endif // #ifdef GARVITY


// initialize the array "MinDtInfo_Fluid" since the kernel "CUFLU_GetMaxCFL" will NOT work during initialization
   real MinDtVar_AllLv_Fluid[NLEVEL][NCOMP];
   if ( OPT__ADAPTIVE_DT )     
   {
#     if   ( MODEL == HYDRO )
      Hydro_GetMaxCFL( MinDtInfo_Fluid, MinDtVar_AllLv_Fluid );

#     elif ( MODEL == MHD )
#     error : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
//    nothing to do here

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL
   }

// initialize the array "MinDtInfo_Gravity" since the kernel "XXX" will NOT work during initialization
#  ifdef GRAVITY
   if ( OPT__ADAPTIVE_DT )     
   {
#     if   ( MODEL == HYDRO )
      Hydro_GetMaxAcc( MinDtInfo_Gravity );

#     elif ( MODEL == MHD )
#     error : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      ELBDM_GetMaxPot( MinDtInfo_Gravity );

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL
   }
#  endif

// initialize the array "MinDtInfo_Phase" since the kernel "XXX" will NOT work during initialization
#  if ( MODEL == ELBDM )
   real MinDtVar_AllLv_Phase[NLEVEL][3];

   if ( OPT__ADAPTIVE_DT )     
      ELBDM_GetMaxPhaseDerivative( MinDtInfo_Phase, MinDtVar_AllLv_Phase );
#  endif


} // FUNCTION : Init_DAINO
