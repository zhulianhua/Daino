#include "GetCube.h"

void Load_Parameter_Before_1200( FILE *File, const int FormatVersion, bool &DataOrder_xyzv, 
                                 bool &LoadPot, int *NX0_Tot, double &BoxSize, real &Gamma );
void Load_Parameter_After_1200( FILE *File, const int FormatVersion, bool &DataOrder_xyzv,
                                bool &LoadPot, int *NX0_Tot, double &BoxSize, real &Gamma );
void CompareVar( const char *VarName, const bool   RestartVar, const bool   RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const int    RestartVar, const int    RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const long   RestartVar, const long   RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const real   RestartVar, const real   RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const double RestartVar, const double RuntimeVar, const bool Fatal );




//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData
// Description :  Load an output data of DAINO
//
// Parameter   :  FileName       : The name of the input file
//                OldDataFormat  : Invoke the function "LoadData_Old" to load the old-format data
//                                 (no longer supported)
//-------------------------------------------------------------------------------------------------------
//void LoadData( const char *FileName, const bool OldDataFormat )
void LoadData( const char *FileName )
{

// invoke the function "LoadData_Old" to load file with the old data format
   /*
   if ( OldDataFormat )
   {
      LoadData_Old( FileName );

      return;
   }
   */


   if ( MyRank == 0 )   cout << "LoadData \"" << FileName << "\" ... "<< endl;


   FILE *File = fopen( FileName, "rb" );

   if ( File == NULL  &&  MyRank == 0 )
   {
      fprintf( stderr, "ERROR : the input file \"%s\" does not exist !!\n", FileName );
      MPI_Exit();
   }


// initialize the NPatchComma list as 0
   for (int lv=0; lv<NLEVEL; lv++)  for (int m=0; m<28; m++)      NPatchComma[lv][m] = 0;


// record the size of different data types
   const int size_bool   = sizeof( bool   );
   const int size_int    = sizeof( int    );
   const int size_uint   = sizeof( uint   );
   const int size_long   = sizeof( long   );
   const int size_ulong  = sizeof( ulong  );
   const int size_real   = sizeof( real   );
   const int size_double = sizeof( double );

   if ( size_int != size_uint )
   {
      fprintf( stderr, "ERROR : sizeof(int) = %d != sizeof(uint) = %d !!\n", size_int, size_uint );
      MPI_Exit();
   }

   if ( size_long != size_ulong )
   {
      fprintf( stderr, "ERROR : sizeof(long) = %d != sizeof(ulong) = %d !!\n", size_long, size_ulong );
      MPI_Exit();
   }


// a. load the information of data format
// =================================================================================================
   long FormatVersion, HeaderSize, CheckCode;

   fread( &FormatVersion, sizeof(long), 1, File );
   fread( &HeaderSize,    sizeof(long), 1, File );
   fread( &CheckCode,     sizeof(long), 1, File );


// verify the input data format version
   if ( MyRank == 0 )   fprintf( stdout, "   The version of the RESTART file's format = %ld\n", FormatVersion );

   if ( FormatVersion < 1100 )
   {
      fprintf( stderr, "ERROR : unsupported data format version (only support version >= 1100) !!\n" );
      MPI_Exit();
   }


// check if the size of different data types are consistent (only for version >= 1200)
   if ( FormatVersion >= 1200 )
   {
      int size_bool_restart, size_int_restart, size_long_restart, size_real_restart, size_double_restart;

      fread( &size_bool_restart,   sizeof(int), 1, File );
      fread( &size_int_restart,    sizeof(int), 1, File );
      fread( &size_long_restart,   sizeof(int), 1, File );
      fread( &size_real_restart,   sizeof(int), 1, File );
      fread( &size_double_restart, sizeof(int), 1, File );

      if ( size_bool_restart != size_bool )  
      {
         fprintf( stderr, "ERROR : sizeof(bool) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                  size_bool_restart, size_bool );
         MPI_Exit();
      }

      if ( size_int_restart != size_int )  
      {
         fprintf( stderr, "ERROR : sizeof(int) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                  size_int_restart, size_int );
         MPI_Exit();
      }

      if ( size_long_restart != size_long )  
      {
         fprintf( stderr, "ERROR : sizeof(long) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                  size_long_restart, size_long );
         MPI_Exit();
      }

      if ( size_real_restart != size_real )  
      {
         fprintf( stderr, "ERROR : sizeof(real) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                  size_real_restart, size_real );
         MPI_Exit();
      }

      if ( size_double_restart != size_double )  
      {
         fprintf( stderr, "ERROR : sizeof(double) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                  size_double_restart, size_double );
         MPI_Exit();
      }
   } // if ( FormatVersion >= 1200 )


// skip the buffer space
   const int NBuf_Format_1200 = 256 - 0*size_bool -  5*size_int -  3*size_long -  0*size_real -  0*size_double;
   const int NBuf_Format      = ( FormatVersion >= 1200 ) ? NBuf_Format_1200 : 40; 

   fseek( File, NBuf_Format, SEEK_CUR );



// b. load all simulation parameters
// =================================================================================================
   bool   DataOrder_xyzv;
   double BoxSize;

   if ( FormatVersion < 1200 )   
      Load_Parameter_Before_1200( File, FormatVersion, DataOrder_xyzv, OutputPot, NX0_TOT, BoxSize, GAMMA );
   else
      Load_Parameter_After_1200 ( File, FormatVersion, DataOrder_xyzv, OutputPot, NX0_TOT, BoxSize, GAMMA );


// set the file position indicator to "HeaderSize" and verify the check code
   long checkcode;
   fseek( File, HeaderSize, SEEK_SET );
   fread( &checkcode, sizeof(long), 1, File );

   if ( checkcode != CheckCode )
   {
      fprintf( stderr, "ERROR : incorrect check code in the RESTART file (input %ld <-> expeect %ld) !!\n",
               checkcode, CheckCode );
      MPI_Exit();
   }



// c. load the simulation information
// =================================================================================================
   uint AdvanceCounter[NLEVEL];
   int  NPatchTotal[NLEVEL], NDataPatch_Total[NLEVEL];

   fread( &DumpID,          sizeof(int),         1, File );
   fread( Time,             sizeof(double), NLEVEL, File );
   fread( &Step,            sizeof(long),        1, File );
   fread( NPatchTotal,      sizeof(int),    NLEVEL, File );
   fread( NDataPatch_Total, sizeof(int),    NLEVEL, File );
   fread( AdvanceCounter,   sizeof(uint),   NLEVEL, File );
   if ( FormatVersion >= 1200 )
   fseek( File, sizeof(double), SEEK_CUR );


// skip the buffer space
   const int NBuf_Info_1200 = 1024 - 0*size_bool - (1+3*NLEVEL)*size_int - 2*size_long 
                                   - 0*size_real - (1+NLEVEL)*size_double;
   const int NBuf_Info      = ( FormatVersion >= 1200 ) ? NBuf_Info_1200 : 80-size_double;

   fseek( File, NBuf_Info, SEEK_CUR );


// verify the size of the RESTART file
   long InfoSize, DataSize[NLEVEL], ExpectSize, InputSize, PatchDataSize;

   InfoSize =     sizeof(int   )*( 1 + 2*NLEVEL )
                + sizeof(long  )*( 2            )  // Step + checkcode
                + sizeof(uint  )*(       NLEVEL )
                + sizeof(double)*( 1 +   NLEVEL )
                + NBuf_Info;

   if ( OutputPot )  
   {
      NLoad ++;
      NOut  ++;
   }

   PatchDataSize = PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NLoad*sizeof(real);
   ExpectSize    = HeaderSize + InfoSize;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      DataSize[lv]  = 0;
      DataSize[lv] += NPatchTotal[lv]*4*sizeof(int);        // 4 = corner(3) + son(1)
      DataSize[lv] += NDataPatch_Total[lv]*PatchDataSize;

      ExpectSize   += DataSize[lv];
   }

   fseek( File, 0, SEEK_END );
   InputSize = ftell( File );

   if ( InputSize != ExpectSize )
   {
      fprintf( stderr, "ERROR : the size of the file <%s> is incorrect --> input = %ld <-> expect = %ld !!\n",
               FileName, InputSize, ExpectSize );
      MPI_Exit();
   }

   fclose( File );


// set up the simulation box   
   int NX0_Max;
   NX0_Max = ( NX0_TOT[0] > NX0_TOT[1] ) ? NX0_TOT[0] : NX0_TOT[1];
   NX0_Max = ( NX0_TOT[2] > NX0_Max    ) ? NX0_TOT[2] : NX0_Max;

   for (int lv=0; lv<NLEVEL; lv++)     patch.dh[lv] = BoxSize / (double)( NX0_Max*(1<<lv) );

   for (int d=0; d<3; d++)    
   {
      patch.BoxSize [d] = NX0_TOT[d]*patch.dh   [0];
      patch.BoxScale[d] = NX0_TOT[d]*patch.scale[0];
   }

   NX0[0] = NX0_TOT[0] / NGPU_X[0];
   NX0[1] = NX0_TOT[1] / NGPU_X[1];
   NX0[2] = NX0_TOT[2] / NGPU_X[2];



// d. set up and check parameters dedicated in this tool
// =================================================================================================
// initialize parameters related to the targeted domain
   Init_TargetDomain();

// allocate memory for the array "BaseP"
   Init_MemAllocate();

// verify the input parameters
   if ( MyRank == 0 )   CheckParameter();

// set the candidate box
   GetCandidateBox();

// set the range of the targeted sub-domain
   int TargetRange_Min[3], TargetRange_Max[3];

   for (int d=0; d<3; d++)
   {
      TargetRange_Min[d] = MyRank_X[d]*NX0[d]*patch.scale[0];
      TargetRange_Max[d] = TargetRange_Min[d] + NX0[d]*patch.scale[0];
   }



// e. load the simulation data
// =================================================================================================
   long int Offset = HeaderSize+InfoSize;
   int LoadCorner[3], LoadSon, PID;
   bool GotYou;

// array for re-ordering the fluid data from "xyzv" to "vxyz"
   real (*InvData_Flu)[PATCH_SIZE][PATCH_SIZE][NCOMP] = NULL;
   if ( DataOrder_xyzv )   InvData_Flu = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE][NCOMP];


   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int TargetRank=0; TargetRank<NGPU; TargetRank++)
      {
         if ( MyRank == 0 )
         {
            fprintf( stdout, "   Loading data: level %2d, MPI_Rank %3d ... ", lv, TargetRank ); 
            fflush( stdout );
         }

         if ( MyRank == TargetRank )
         {

            File = fopen( FileName, "rb" );
            fseek( File, Offset, SEEK_SET );

            for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)
            {

//             e1. load the patch information
               fread(  LoadCorner, sizeof(int), 3, File );
               fread( &LoadSon,    sizeof(int), 1, File );


//             verify that the loaded patch is within the targeted range
               if (  LoadCorner[0] >= TargetRange_Min[0]  &&  LoadCorner[0] < TargetRange_Max[0]  &&
                     LoadCorner[1] >= TargetRange_Min[1]  &&  LoadCorner[1] < TargetRange_Max[1]  &&
                     LoadCorner[2] >= TargetRange_Min[2]  &&  LoadCorner[2] < TargetRange_Max[2]     ) 
               {

//                verify that the loaded patch is within the candidate box
                  GotYou = WithinCandidateBox( LoadCorner, PATCH_SIZE*patch.scale[lv], CanBuf );

                  patch.pnew( lv, LoadCorner[0], LoadCorner[1], LoadCorner[2], -1, GotYou );

//                e2. load the physical data if it is a leaf patch
                  if ( LoadSon == -1 )
                  {
                     if ( GotYou )
                     {
                        PID = patch.num[lv] - 1;

//                      e2-1. load the fluid variables
                        if ( DataOrder_xyzv )
                        {
                           fread( InvData_Flu, sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NCOMP, File );

                           for (int v=0; v<NCOMP; v++)
                           for (int k=0; k<PATCH_SIZE; k++)
                           for (int j=0; j<PATCH_SIZE; j++)
                           for (int i=0; i<PATCH_SIZE; i++)    
                              patch.ptr[lv][PID]->fluid[v][k][j][i] = InvData_Flu[k][j][i][v];
                        }

                        else
                           fread( patch.ptr[lv][PID]->fluid, sizeof(real), 
                                  PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NCOMP, File );

//                      e2-2. load the gravitational potential
                        if ( OutputPot )
                           fread( patch.ptr[lv][PID]->pot, sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE, File );
                     }

                     else
                        fseek( File, PatchDataSize, SEEK_CUR );
                  }
               }

               else
               {
                  if ( LoadSon == -1 )    fseek( File, PatchDataSize, SEEK_CUR );
               }

            } // for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)

            fclose( File );

            Offset += DataSize[lv];

         } // if ( MyRank == TargetRank )

         MPI_Barrier( MPI_COMM_WORLD );

         if ( MyRank == 0 )
         {
            fprintf( stdout, "done\n" ); 
            fflush( stdout );
         }

      } // for (int TargetRank=0; TargetRank<NGPU; TargetRank++)
   } // for (int lv=0; lv<NLEVEL; lv++)


   if ( DataOrder_xyzv )  delete [] InvData_Flu;


// record the number of the real patches
   for (int lv=0; lv<NLEVEL; lv++)     
   for (int m=1; m<28; m++)               NPatchComma[lv][m] = patch.num[lv];


// complete all levels 
//=====================================================================================
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    construct the relation : father <-> son
      if ( lv > 0 )     FindFather( lv ); 

//    allocate the buffer patches 
      Buf_AllocateBufferPatch( lv );

//    set up the BaseP List
      if ( lv == 0 )    Init_RecordBasePatch();

//    set up the BounP_IDMap 
      Buf_RecordBoundaryPatch( lv );

//    construct the sibling relation
      SiblingSearch( lv );

//    get the IDs of patches for sending and receiving data between neighbor ranks
      Buf_RecordExchangeDataPatchID( lv );
   }


// fill up the data for patches that are not leaf patches
   for (int lv=NLEVEL-2; lv>=0; lv--)     Flu_Restrict( lv, OutputPot );


// fill up the data in the buffer patches
   for (int lv=0; lv<NLEVEL; lv++)        
   {
      Buf_GetBufferData( lv, 1, BUF_SIZE );

      if ( OutputPot )
      Buf_GetBufferData( lv, 2, BUF_SIZE );
   }


   if ( MyRank == 0 )   cout << "LoadData ... done" << endl;

} // FUNCTION : LoadData



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Parameter_Before_1200
// Description :  Load all simulation parameters from the RESTART file with format version < 1200
//
// Note        :  Only work in HYDRO 
//
// Parameter   :  File           : RESTART file pointer 
//                FormatVersion  : Format version of the RESTART file
//                DataOrder_xyzv : Order of data stored in the RESTART file (true/false --> xyzv/vxyz)
//                LoadPot        : Whether or not the RESTART file stores the potential data
//                NX0_Tot        : Total number of base-level cells along each direction
//                BoxSize        : Physical size of the simulation box
//                Gamma          : ratio of specific heat
//
// Return      :  DataOrder_xyzv, LoadPot, NX0_Tot, BoxSize, Gamma
//-------------------------------------------------------------------------------------------------------
void Load_Parameter_Before_1200( FILE *File, const int FormatVersion, bool &DataOrder_xyzv, 
                                 bool &LoadPot, int *NX0_Tot, double &BoxSize, real &Gamma )
{
   
// set the size of the output buffers
   const int NBuf_MakefileOption = 40  - 0*sizeof(int) - 0*sizeof(real) - 0*sizeof(double);
   const int NBuf_MakefileConst  = 80  - 0*sizeof(int) - 0*sizeof(real) - 0*sizeof(double);
   const int NBuf_Parameter      = 196 - 0*sizeof(int) - 0*sizeof(real) - 1*sizeof(double);


   if ( MyRank == 0 )   fprintf( stdout, "   Loading simulation parameters ...\n" );


// a. load the simulation options defined in the Makefile
// =================================================================================================
   bool gravity, comoving, float8;

   fread( &gravity,                    sizeof(bool),                    1,             File );
   fread( &comoving,                   sizeof(bool),                    1,             File );
   fread( &float8,                     sizeof(bool),                    1,             File );

// skip the buffer space
   fseek( File, NBuf_MakefileOption, SEEK_CUR );


// b. load the symbolic constants defined in the Makefile
// =================================================================================================
   int ncomp, patch_size, max_patch, nlevel, flu_ghost_size, pot_ghost_size, gra_ghost_size;

   fread( &ncomp,                      sizeof(int),                     1,             File );
   fread( &patch_size,                 sizeof(int),                     1,             File );
   fread( &max_patch,                  sizeof(int),                     1,             File );
   fread( &nlevel,                     sizeof(int),                     1,             File );
   fread( &flu_ghost_size,             sizeof(int),                     1,             File );
   fread( &pot_ghost_size,             sizeof(int),                     1,             File );
   fread( &gra_ghost_size,             sizeof(int),                     1,             File );

// skip the buffer space
   fseek( File, NBuf_MakefileConst, SEEK_CUR );


// c. load the simulation parameters recorded in the file "Input__Parameter"
// =================================================================================================
   int    nx0_tot[3], DAINO_nrank, DAINO_nrank_x[3], regrid_count, flag_buffer_size;
   int    opt__flu_int_scheme, opt__pot_int_scheme, opt__rho_int_scheme, opt__gra_int_scheme;
   int    opt__ref_flu_int_scheme, opt__ref_pot_int_scheme, opt__output_total;
   double omega_m0, dt__fluid, dt__gravity, dt__max_delta_a, box_size;
   real   gamma, newton_g;
   bool   opt__gra_p5_gradient, opt__output_pot;

   fread(  nx0_tot,                    sizeof(int),                     3,             File );
   fread( &DAINO_nrank,                sizeof(int),                     1,             File );
   fread(  DAINO_nrank_x,              sizeof(int),                     3,             File );
   fread( &omega_m0,                   sizeof(double),                  1,             File );
   fread( &dt__fluid,                  sizeof(double),                  1,             File );
   fread( &dt__gravity,                sizeof(double),                  1,             File );
   fread( &dt__max_delta_a,            sizeof(double),                  1,             File );
   fread( &regrid_count,               sizeof(int),                     1,             File );
   fread( &flag_buffer_size,           sizeof(int),                     1,             File );
   fread( &gamma,                      sizeof(real),                    1,             File );
   fread( &newton_g,                   sizeof(real),                    1,             File );
   fread( &opt__gra_p5_gradient,       sizeof(bool),                    1,             File );
   fread( &opt__flu_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__pot_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__rho_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__gra_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__ref_flu_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__ref_pot_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__output_pot,            sizeof(bool),                    1,             File );
   fread( &opt__output_total,          sizeof(int),                     1,             File );
   fread( &box_size,                   sizeof(double),                  1,             File );

// skip the buffer space
   fseek( File, NBuf_Parameter, SEEK_CUR );

// reset "box_size" for format version < 1102 (in which this parameter is not added yet)
   if ( FormatVersion < 1102 )
   {
      int nx0_max;
      nx0_max = ( nx0_tot[0] > nx0_tot[1] ) ? nx0_tot[0] : nx0_tot[1];
      nx0_max = ( nx0_tot[2] > nx0_max    ) ? nx0_tot[2] : nx0_max;

      box_size = nx0_max*( 1<<(nlevel-1) );

      if ( MyRank == 0 )  
         fprintf( stderr, "WARNING : loading data with format version < 1102 --> assuming BOX_SIZE = %f\n", 
                  box_size );
   }


   if ( MyRank == 0 )   fprintf( stdout, "   Loading simulation parameters ... done\n" );


// d. check parameters (before loading any size-dependent parameters)
// =================================================================================================
   if ( MyRank == 0 )   fprintf( stdout, "   Checking loaded parameters ...\n" );


   const bool Fatal    = true;
   const bool NonFatal = false;

#  ifdef FLOAT8
   if ( !float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "OFF", "ON" );
      MPI_Exit();
   }
#  else
   if (  float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "ON", "OFF" );
      MPI_Exit();
   }
#  endif

   CompareVar( "NLEVEL",                  nlevel,                 NLEVEL,                       Fatal );
   CompareVar( "NCOMP",                   ncomp,                  NCOMP,                        Fatal );
   CompareVar( "PATCH_SIZE",              patch_size,             PATCH_SIZE,                   Fatal );


   if ( MyRank == 0 )   fprintf( stdout, "   Checking loaded parameters ... done\n" );


// set the returned variables
   DataOrder_xyzv = ( FormatVersion >= 1101 ) ? ( (opt__output_total==1)?true:false ) : true;
   LoadPot        = opt__output_pot;
   BoxSize        = box_size;
   Gamma          = gamma;
   for (int d=0; d<3; d++)    
   NX0_Tot[d]     = nx0_tot[d];

} // FUNCTION : Load_Parameter_Before_1200



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Parameter_After_1200
// Description :  Load all simulation parameters from the RESTART file with format version >= 1200
//
// Note        :  "OPT__RESTART_HEADER == RESTART_HEADER_INHERIC" can be used in this function
//
// Parameter   :  File           : RESTART file pointer 
//                FormatVersion  : Format version of the RESTART file
//                DataOrder_xyzv : Order of data stored in the RESTART file (true/false --> xyzv/vxyz)
//                LoadPot        : Whether or not the RESTART file stores the potential data
//                NX0_Tot        : Total number of base-level cells along each direction
//                BoxSize        : Physical size of the simulation box
//                Gamma          : ratio of specific heat
//
// Return      :  DataOrder_xyzv, LoadPot, NX0_Tot, BoxSize, Gamma
//-------------------------------------------------------------------------------------------------------
void Load_Parameter_After_1200( FILE *File, const int FormatVersion, bool &DataOrder_xyzv,
                                bool &LoadPot, int *NX0_Tot, double &BoxSize, real &Gamma )
{

   const int size_bool      = sizeof( bool   );
   const int size_int       = sizeof( int    );
   const int size_long      = sizeof( long   );
   const int size_real      = sizeof( real   );
   const int size_double    = sizeof( double );

   const int NBuf_Makefile  =  256 - 15*size_bool -  8*size_int -  0*size_long -  0*size_real -  0*size_double;
   const int NBuf_Constant  =  256 -  6*size_bool - 11*size_int -  0*size_long -  2*size_real -  0*size_double;
   const int NBuf_Parameter = 1024 - 17*size_bool - 35*size_int -  1*size_long - 12*size_real -  8*size_double;


   if ( MyRank == 0 )   fprintf( stdout, "   Loading simulation parameters ...\n" );


// a. load the simulation options and parameters defined in the Makefile
// =================================================================================================
   bool gravity, individual_timestep, comoving, gpu, DAINO_optimization, DAINO_debug, timing, timing_solver;
   bool intel, float8, serial, ooc, overlap_mpi, openmp, fermi;
   int  model, pot_scheme, flu_scheme, lr_scheme, rsolver, load_balance, nlevel, max_patch;

   fread( &model,                      sizeof(int),                     1,             File );
   fread( &gravity,                    sizeof(bool),                    1,             File );
   fread( &pot_scheme,                 sizeof(int),                     1,             File );
   fread( &individual_timestep,        sizeof(bool),                    1,             File );
   fread( &comoving,                   sizeof(bool),                    1,             File );
   fread( &flu_scheme,                 sizeof(int),                     1,             File );
   fread( &lr_scheme,                  sizeof(int),                     1,             File );
   fread( &rsolver,                    sizeof(int),                     1,             File );
   fread( &gpu,                        sizeof(bool),                    1,             File );
   fread( &DAINO_optimization,         sizeof(bool),                    1,             File );
   fread( &DAINO_debug,                sizeof(bool),                    1,             File );
   fread( &timing,                     sizeof(bool),                    1,             File );
   fread( &timing_solver,              sizeof(bool),                    1,             File );
   fread( &intel,                      sizeof(bool),                    1,             File );
   fread( &float8,                     sizeof(bool),                    1,             File );
   fread( &serial,                     sizeof(bool),                    1,             File );
   fread( &ooc,                        sizeof(bool),                    1,             File );
   fread( &load_balance,               sizeof(int),                     1,             File );
   fread( &overlap_mpi,                sizeof(bool),                    1,             File );
   fread( &openmp,                     sizeof(bool),                    1,             File );
   fread( &fermi,                      sizeof(bool),                    1,             File );
   fread( &nlevel,                     sizeof(int),                     1,             File );
   fread( &max_patch,                  sizeof(int),                     1,             File );

// skip the buffer space
   fseek( File, NBuf_Makefile, SEEK_CUR );


// b. load the symbolic constants defined in "Macro.h, CUPOT.h, and CUFLU.h"
// =================================================================================================
   bool enforce_positive, char_reconstruction, hll_no_ref_state, hll_include_all_waves, waf_dissipate;
   bool use_psolver_10to14;
   int  ncomp, patch_size, flu_ghost_size, pot_ghost_size, gra_ghost_size, check_intermediate; 
   int  flu_block_size_x, flu_block_size_y, pot_block_size_x, pot_block_size_z, gra_block_size_z;
   real min_value, max_error;

   fread( &ncomp,                      sizeof(int),                     1,             File );
   fread( &patch_size,                 sizeof(int),                     1,             File );
   fread( &min_value,                  sizeof(real),                    1,             File );
   fread( &flu_ghost_size,             sizeof(int),                     1,             File );
   fread( &pot_ghost_size,             sizeof(int),                     1,             File );
   fread( &gra_ghost_size,             sizeof(int),                     1,             File );
   fread( &enforce_positive,           sizeof(bool),                    1,             File );
   fread( &char_reconstruction,        sizeof(bool),                    1,             File );
   fread( &check_intermediate,         sizeof(int),                     1,             File );
   fread( &hll_no_ref_state,           sizeof(bool),                    1,             File );
   fread( &hll_include_all_waves,      sizeof(bool),                    1,             File );
   fread( &waf_dissipate,              sizeof(bool),                    1,             File );
   fread( &max_error,                  sizeof(real),                    1,             File );
   fread( &flu_block_size_x,           sizeof(int),                     1,             File );
   fread( &flu_block_size_y,           sizeof(int),                     1,             File );
   fread( &use_psolver_10to14,         sizeof(bool),                    1,             File );
   fread( &pot_block_size_x,           sizeof(int),                     1,             File );
   fread( &pot_block_size_z,           sizeof(int),                     1,             File );
   fread( &gra_block_size_z,           sizeof(int),                     1,             File );

// skip the buffer space
   fseek( File, NBuf_Constant, SEEK_CUR );


// c. load the simulation parameters recorded in the file "Input__Parameter"
// =================================================================================================
   bool   opt__adaptive_dt, opt__dt_user, opt__flag_rho, opt__flag_rho_gradient, opt__flag_pres_gradient;
   bool   opt__flag_engy_density, opt__flag_user, opt__fixup_flux, opt__fixup_restrict, opt__overlap_mpi;
   bool   opt__gra_p5_gradient, opt__int_time, opt__output_error, opt__output_base, opt__output_pot; 
   bool   opt__timing_barrier, opt__int_phase;
   int    nx0_tot[3], DAINO_nrank, DAINO_nrank_x[3], omp_nthread, ooc_nrank, ooc_nrank_x[3], regrid_count;
   int    flag_buffer_size, max_level, opt__lr_limiter, opt__waf_limiter, flu_gpu_npgroup, gpu_nstream; 
   int    sor_max_iter, sor_min_iter, mg_max_iter, mg_npre_smooth, mg_npost_smooth, pot_gpu_npgroup;
   int    opt__flu_int_scheme, opt__pot_int_scheme, opt__rho_int_scheme;
   int    opt__gra_int_scheme, opt__ref_flu_int_scheme, opt__ref_pot_int_scheme;
   int    opt__output_total, opt__output_part, opt__output_mode, output_step;
   long   end_step; 
   real   lb_wli_max, gamma, minmod_coeff, ep_coeff, elbdm_mass, planck_const, newton_g, sor_omega;
   real   mg_tolerated_error, output_part_x, output_part_y, output_part_z;
   double box_size, end_t, omega_m0, dt__fluid, dt__gravity, dt__phase, dt__max_delta_a, output_dt;

   fread( &box_size,                   sizeof(double),                  1,             File );
   fread(  nx0_tot,                    sizeof(int),                     3,             File );
   fread( &DAINO_nrank,                sizeof(int),                     1,             File );
   fread(  DAINO_nrank_x,              sizeof(int),                     3,             File );
   fread( &omp_nthread,                sizeof(int),                     1,             File );
   fread( &end_t,                      sizeof(double),                  1,             File );
   fread( &end_step,                   sizeof(long),                    1,             File );
   fread( &ooc_nrank,                  sizeof(int),                     1,             File );
   fread(  ooc_nrank_x,                sizeof(int),                     3,             File );
   fread( &omega_m0,                   sizeof(double),                  1,             File );
   fread( &dt__fluid,                  sizeof(double),                  1,             File );
   fread( &dt__gravity,                sizeof(double),                  1,             File );
   fread( &dt__phase,                  sizeof(double),                  1,             File );
   fread( &dt__max_delta_a,            sizeof(double),                  1,             File );
   fread( &opt__adaptive_dt,           sizeof(bool),                    1,             File );
   fread( &opt__dt_user,               sizeof(bool),                    1,             File );
   fread( &regrid_count,               sizeof(int),                     1,             File );
   fread( &flag_buffer_size,           sizeof(int),                     1,             File );
   fread( &max_level,                  sizeof(int),                     1,             File );
   fread( &opt__flag_rho,              sizeof(bool),                    1,             File );
   fread( &opt__flag_rho_gradient,     sizeof(bool),                    1,             File );
   fread( &opt__flag_pres_gradient,    sizeof(bool),                    1,             File );
   fread( &opt__flag_engy_density,     sizeof(bool),                    1,             File );
   fread( &opt__flag_user,             sizeof(bool),                    1,             File );
   fread( &lb_wli_max,                 sizeof(real),                    1,             File );
   fread( &gamma,                      sizeof(real),                    1,             File );
   fread( &minmod_coeff,               sizeof(real),                    1,             File );
   fread( &ep_coeff,                   sizeof(real),                    1,             File );
   fread( &opt__lr_limiter,            sizeof(int),                     1,             File );
   fread( &opt__waf_limiter,           sizeof(int),                     1,             File );
   fread( &elbdm_mass,                 sizeof(real),                    1,             File );
   fread( &planck_const,               sizeof(real),                    1,             File );
   fread( &flu_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &gpu_nstream,                sizeof(int),                     1,             File );
   fread( &opt__fixup_flux,            sizeof(bool),                    1,             File );
   fread( &opt__fixup_restrict,        sizeof(bool),                    1,             File );
   fread( &opt__overlap_mpi,           sizeof(bool),                    1,             File );
   fread( &newton_g,                   sizeof(real),                    1,             File );
   fread( &sor_omega,                  sizeof(real),                    1,             File );
   fread( &sor_max_iter,               sizeof(int),                     1,             File );
   fread( &sor_min_iter,               sizeof(int),                     1,             File );
   fread( &mg_max_iter,                sizeof(int),                     1,             File );
   fread( &mg_npre_smooth,             sizeof(int),                     1,             File );
   fread( &mg_npost_smooth,            sizeof(int),                     1,             File );
   fread( &mg_tolerated_error,         sizeof(real),                    1,             File );
   fread( &pot_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &opt__gra_p5_gradient,       sizeof(bool),                    1,             File );
   fread( &opt__int_time,              sizeof(bool),                    1,             File );
   fread( &opt__int_phase,             sizeof(bool),                    1,             File );
   fread( &opt__flu_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__pot_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__rho_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__gra_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__ref_flu_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__ref_pot_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__output_total,          sizeof(int),                     1,             File );
   fread( &opt__output_part,           sizeof(int),                     1,             File );
   fread( &opt__output_error,          sizeof(bool),                    1,             File );
   fread( &opt__output_base,           sizeof(bool),                    1,             File );
   fread( &opt__output_pot,            sizeof(bool),                    1,             File );
   fread( &opt__output_mode,           sizeof(int),                     1,             File );
   fread( &output_step,                sizeof(int),                     1,             File );
   fread( &output_dt,                  sizeof(double),                  1,             File );
   fread( &output_part_x,              sizeof(real),                    1,             File );
   fread( &output_part_y,              sizeof(real),                    1,             File );
   fread( &output_part_z,              sizeof(real),                    1,             File );
   fread( &opt__timing_barrier,        sizeof(bool),                    1,             File );

// skip the buffer space
   fseek( File, NBuf_Parameter, SEEK_CUR );


   if ( MyRank == 0 )   fprintf( stdout, "   Loading simulation parameters ... done\n" );


// d. check parameters (before loading any size-dependent parameters)
// =================================================================================================
   if ( MyRank == 0 )   fprintf( stdout, "   Checking loaded parameters ...\n" );


   const bool Fatal    = true;
   const bool NonFatal = false;

#  ifdef FLOAT8
   if ( !float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "OFF", "ON" );
      MPI_Exit();
   }
#  else
   if (  float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "ON", "OFF" );
      MPI_Exit();
   }
#  endif

   CompareVar( "MODEL",                   model,                  MODEL,                        Fatal );
   CompareVar( "NLEVEL",                  nlevel,                 NLEVEL,                       Fatal );
   CompareVar( "NCOMP",                   ncomp,                  NCOMP,                        Fatal );
   CompareVar( "PATCH_SIZE",              patch_size,             PATCH_SIZE,                   Fatal );


   if ( MyRank == 0 )   fprintf( stdout, "   Checking loaded parameters ... done\n" );


// set the returned variables
   DataOrder_xyzv = ( opt__output_total == 1 ) ? true : false;
   LoadPot        = opt__output_pot;
   BoxSize        = box_size;
   Gamma          = gamma;
   for (int d=0; d<3; d++)    
   NX0_Tot[d]     = nx0_tot[d];

} // FUNCTION : Load_Parameter_After_1200



//-------------------------------------------------------------------------------------------------------
// Function    :  CompareVar 
// Description :  Compare the input variables 
//
// Note        :  This function is overloaded to work with different data types
//
// Parameter   :  VarName     : Name of the targeted variable 
//                RestartVar  : Variable loaded from the RESTART file
//                RuntimeVar  : Variable loaded from the Input__Parameter
//                Fatal       : Whether or not the difference between RestartVar and RuntimeVar is fatal
//                              --> true  : terminate the program if the input variables are different
//                                  false : display warning message if the input variables are different
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const bool RestartVar, const bool RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
      {
         fprintf( stderr, "ERROR : %s : RESTART file (%d) != runtime (%d) !!\n", 
                  VarName, RestartVar, RuntimeVar );
         MPI_Exit();
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%d) != runtime (%d) !!\n", 
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (bool)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "int"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const int RestartVar, const int RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
      {
         fprintf( stderr, "ERROR : %s : RESTART file (%d) != runtime (%d) !!\n", 
                  VarName, RestartVar, RuntimeVar );
         MPI_Exit();
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%d) != runtime (%d) !!\n", 
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (int)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "long"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const long RestartVar, const long RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
      {
         fprintf( stderr, "ERROR : %s : RESTART file (%ld) != runtime (%ld) !!\n", 
                  VarName, RestartVar, RuntimeVar );
         MPI_Exit();
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%ld) != runtime (%ld) !!\n", 
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (long)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "float"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const float RestartVar, const float RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
      {
         fprintf( stderr, "ERROR : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n", 
                  VarName, RestartVar, RuntimeVar );
         MPI_Exit();
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n", 
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (float)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "double"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const double RestartVar, const double RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
      {
         fprintf( stderr, "ERROR : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n", 
                  VarName, RestartVar, RuntimeVar );
         MPI_Exit();
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n", 
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (double)

