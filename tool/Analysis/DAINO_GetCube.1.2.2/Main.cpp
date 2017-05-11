#include "GetCube.h"

#define WRONG -9999999
#define PS    PATCH_SIZE

int       OutputXYZ         = WRONG;                     // output option (x,y,z,x-proj,y-proj,z-proj,3D)
int       InterScheme       = 0;                         // interpolation scheme (0/1:MinMod/Central)
char     *FileName_In       = NULL;                      // name of the input file
char     *Suffix            = NULL;                      // suffix attached to the output file name
//bool    OldDataFormat     = false;                     // true --> load the old-format data (no longer supported)
bool      OutputBinary      = false;                     // true --> output data in binary form
bool      OutputDivVel      = false;                     // true --> output divergence( velocity )
bool      OutputCurlVel     = false;                     // true --> output curl( velocity ) = vorticity
bool      OutputPot         = false;                     // true --> output gravitational potential
bool      InputScale        = false;                     // true --> input cell scales instead of coordinates
int       OMP_NThread       = -1;                        // number of OpenMP threads
int       TargetLevel       = 0;                         // targeted level
int       Scale_Start   [3] = { WRONG, WRONG, WRONG };   // targeted x, y, z starting scale
int       Scale_Size    [3] = { WRONG, WRONG, WRONG };   // targeted x, y, z scale size
int       Idx_Start     [3] = { WRONG, WRONG, WRONG };   // targeted x, y, z array indices
int       Idx_Size      [3] = { WRONG, WRONG, WRONG };   // targeted x, y, z array size
int       Idx_MyStart   [3] = { WRONG, WRONG, WRONG };   // starting x, y, z array indices of this process
int       Idx_MySize    [3] = { WRONG, WRONG, WRONG };   // array size of this process
real      PhyCoord_Start[3] = { WRONG, WRONG, WRONG };   // starting physical coordinates
real      PhyCoord_Size [3] = { WRONG, WRONG, WRONG };   // targeted size in physical coordinates
int       NGPU_X[3]         = { 1, 1, 1 };               // number of MPI ranks in each direction
int       CanBuf            = WRONG;                     // buffer size for the candidate box
int       NLoad             = NCOMP;                     // number of variables loaded from the input file
int       NOut              = NCOMP;                     // number of variables to be outputted
real     *OutputArray       = NULL;                      // array storing the output data
int       CanMin[3][3], CanMax[3][3];                    // range of the candidate box for loading patch data

DAINO_t   patch;
ParaVar_t ParaVar;
real      GAMMA;      
int       NX0_TOT[3], NX0[3], MyRank, MyRank_X[3], SibRank[26], NGPU, DumpID;
int       *BaseP, *BounP_IDMap[NLEVEL][26], SendP_NList[NLEVEL][26], *SendP_IDList[NLEVEL][26];  
int       RecvP_NList[NLEVEL][26], *RecvP_IDList[NLEVEL][26], NPatchComma[NLEVEL][28];
double    Time[NLEVEL];
long      Step;




//-------------------------------------------------------------------------------------------------------
// Function    :  ReadOption
// Description :  Read the command-line options
//-------------------------------------------------------------------------------------------------------
void ReadOption( int argc, char **argv )
{

   double Temp_Start[3] = { WRONG, WRONG, WRONG };
   double Temp_Size [3] = { WRONG, WRONG, WRONG };
   int c;

   while ( (c = getopt(argc, argv, "hbdvsi:o:l:n:x:y:z:X:Y:Z:p:q:r:I:t:")) != -1 )
   {
      switch ( c )
      {
         case 'i': FileName_In    = optarg;
                   break;
         case 'o': Suffix         = optarg;
                   break;
         case 'l': TargetLevel    = atoi(optarg);
                   break;
         case 'n': OutputXYZ      = atoi(optarg);
                   break;
         case 'x': Temp_Start[0]  = atof(optarg);
                   break;
         case 'y': Temp_Start[1]  = atof(optarg);
                   break;
         case 'z': Temp_Start[2]  = atof(optarg);
                   break;
         case 'X': Temp_Size[0]   = atof(optarg);
                   break;
         case 'Y': Temp_Size[1]   = atof(optarg);
                   break;
         case 'Z': Temp_Size[2]   = atof(optarg);
                   break;
         case 'p': NGPU_X[0]      = atoi(optarg);
                   break;
         case 'q': NGPU_X[1]      = atoi(optarg);
                   break;
         case 'r': NGPU_X[2]      = atoi(optarg);
                   break;
         case 'I': InterScheme    = atoi(optarg);
                   break;
         case 't': OMP_NThread    = atoi(optarg);
                   break;
//       case 'O': OldDataFormat  = true; 
//                 break;
         case 'b': OutputBinary   = true; 
                   break;
         case 'd': OutputDivVel   = true; 
                   break;
         case 'v': OutputCurlVel  = true; 
                   break;
         case 's': InputScale     = true; 
                   break;
         case 'h': 
         case '?': cerr << endl << "usage: " << argv[0] 
                        << " [-h (for help)] [-i input fileName] [-o suffix to the output file [none]]" 
                        << endl << "                      "
                        << " [-l targeted level [0]] [-I interpolation scheme (0/1: MinMod/central) [0]]"
                        << endl << "                      "
                        << " [-n output option (1~7 : X-slice, Y-slice, Z-slice, X-proj, Y-proj, Z-proj, 3D)" 
                        << endl << "                      "
                        << " [-x/y/z starting coordinate in x/y/z [0]] [-X/Y/Z target size in x/y/z [BoxSize]]"
                        << endl << "                      "
                        << " [-p/q/r # of ranks in the x/y/z directions [1,1,1]]"
                        << endl << "                      "
                        << " [-d (output div(vel)) [off]] [-v (output curl(vel)) [off]]"
                        << endl << "                      "
//                      << " [-O (load the old-format data) [off]] [-b (output binary file) [off]]"
                        << " [-b (output binary file) [off]]"
                        << endl << "                      "
                        << " [-s [input cell scales instead of physical coordinates to specify the range] [off]]"
                        << endl << "                      "
                        << " [-t number of OpenMP threads [omp_get_max_threads]]"
                        << endl << endl;
                   exit( 1 );

      } // switch ( c ) ...
   } // while ...


// reset the useless options
#  if ( MODEL != HYDRO  &&  MODEL != MHD )
   if ( OutputDivVel )
   {
      OutputDivVel = false;

      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -d (OutputDivVel) is useless in this MODEL (%d) !!\n", MODEL );
   }

   if ( OutputCurlVel )
   {
      OutputCurlVel = false;

      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -v (OutputCurlVel) is useless in this MODEL (%d) !!\n", MODEL );
   }
#  endif


// set up OpenMP
#  ifdef OPENMP
   const int OMP_Max_NThread = omp_get_max_threads();

   if ( OMP_NThread <= 0 )  
   {
      OMP_NThread = OMP_Max_NThread;

      if ( MyRank == 0 )  fprintf( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                   "OMP_NThread", OMP_NThread );
   }

   else if ( OMP_NThread > OMP_Max_NThread   &&  MyRank == 0 )
      fprintf( stderr, "WARNING : OMP_NThread (%d) > omp_get_max_threads (%d) !!\n", OMP_NThread, OMP_Max_NThread );

   omp_set_num_threads( OMP_NThread );
   omp_set_nested( false );

#  else 
   OMP_NThread = 1;
#  endif // #ifdef OPENMP ... else ...


// set target range
   if ( InputScale )
   {
      for (int d=0; d<3; d++)
      {
         Scale_Start[d] = (int)Temp_Start[d];
         Scale_Size [d] = (int)Temp_Size [d];
      }
   }
   else
   {
      for (int d=0; d<3; d++)
      {
         PhyCoord_Start[d] = (real)Temp_Start[d];
         PhyCoord_Size [d] = (real)Temp_Size [d];
      }
   }

   NGPU = NGPU_X[0]*NGPU_X[1]*NGPU_X[2];

   if ( OutputDivVel )     NOut ++;
   if ( OutputCurlVel )    NOut +=3;


// check the input command-line options
   if ( FileName_In == NULL )
   {
      cerr << "ERROR : please provide the name of the input file (-i FileName) !!" << endl;
      exit( 1 );
   }

   if ( OutputXYZ == WRONG )
   {
      cerr << "ERROR : please provide the targeted outupt data (-n output option) !!" << endl;
      exit( 1 );
   }

   if ( TargetLevel >= NLEVEL  ||  TargetLevel < 0 )
   {
      cerr << "ERROR : incorrect TargetLevel input (-l TargetLevel) !!" << endl;
      exit( 1 );
   }

   if ( OutputXYZ < 1  ||  OutputXYZ > 7 )
   {
      cerr << "ERROR : incorrect OutputXYZ input (-n output option) !!" << endl;
      exit( 1 );
   }

   if (  OutputBinary  &&  ( NGPU_X[0] != 1 || NGPU_X[1] != 1 )  )
   {
      cerr << "ERROR : only slab decomposition is allowed (-p 1 -q 1 -r NP) for binary output !!" << endl;
      exit( 1 );
   }

} // FUNCTION : ReadOption



//-------------------------------------------------------------------------------------------------------
// Function    :  TakeNote
// Description :  Output the simulation parameters 
//-------------------------------------------------------------------------------------------------------
void TakeNote()
{

   cout << "TakeNote ... " << endl;

   cout << "==========================================================" << endl;
   cout << "DumpID            = " << DumpID                             << endl;
   cout << "Time              = " << Time[0]                            << endl;
   cout << "Step              = " << Step                               << endl;
   cout << "NX0_TOT[0]        = " << NX0_TOT[0]                         << endl;
   cout << "NX0_TOT[1]        = " << NX0_TOT[1]                         << endl;
   cout << "NX0_TOT[2]        = " << NX0_TOT[2]                         << endl;
   cout << "BoxScale[0]       = " << patch.BoxScale[0]                  << endl;
   cout << "BoxScale[1]       = " << patch.BoxScale[1]                  << endl;
   cout << "BoxScale[2]       = " << patch.BoxScale[2]                  << endl;
   cout << "BoxSize[0]        = " << patch.BoxSize[0]                   << endl;
   cout << "BoxSize[1]        = " << patch.BoxSize[1]                   << endl;
   cout << "BoxSize[2]        = " << patch.BoxSize[2]                   << endl;
   cout << "----------------------------------------------------------" << endl;
#  if   ( MODEL == HYDRO )
   cout << "MODEL             = " << "HYDRO"                            << endl;
#  elif ( MODEL == MHD )
   cout << "MODEL             = " << "MHD"                              << endl;
#  elif ( MODEL == ELBDM )
   cout << "MODEL             = " << "ELBDM"                            << endl;
#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL
#  ifdef FLOAT8
   cout << "FLOAT8            = " << "ON"                               << endl;
#  else
   cout << "FLOAT8            = " << "OFF"                              << endl;
#  endif
#  ifdef OPENMP
   cout << "OPENMP            = " << "ON"                               << endl;
   cout << "Number of threads = " << OMP_NThread                        << endl;
#  else
   cout << "OPENMP            = " << "OFF"                              << endl;
#  endif
   cout << "BUF_SIZE          = " << BUF_SIZE                           << endl;
   cout << "NLoad             = " << NLoad                              << endl;
   cout << "NOut              = " << NOut                               << endl;
   cout << "OutputPot         = " << ( (OutputPot)     ? "YES" : "NO" ) << endl;
   cout << "OutputDivVel      = " << ( (OutputDivVel)  ? "YES" : "NO" ) << endl;
   cout << "OutputCurlVel     = " << ( (OutputCurlVel) ? "YES" : "NO" ) << endl;
   cout << "Output option     = " << OutputXYZ                          << endl;
   cout << "IntScheme         = " << InterScheme                        << endl;
   cout << "NGPU_X[0]         = " << NGPU_X[0]                          << endl;
   cout << "NGPU_X[1]         = " << NGPU_X[1]                          << endl;
   cout << "NGPU_X[2]         = " << NGPU_X[2]                          << endl;
   cout << "----------------------------------------------------------" << endl;
   cout << "TargetLv          = " << TargetLevel                        << endl;
   cout << "dh[TargetLv]      = " << patch.dh[TargetLevel]              << endl;
   cout << "dh[FinestLv]      = " << patch.dh[NLEVEL-1]                 << endl;
   cout << "scale[TargetLv]   = " << patch.scale[TargetLevel]           << endl;
   cout << "InputScale        = " << ( (InputScale) ? "YES" : "NO" )    << endl;
   cout << "Scale_Start[0]    = " << Scale_Start[0]                     << endl;
   cout << "Scale_Start[1]    = " << Scale_Start[1]                     << endl;
   cout << "Scale_Start[2]    = " << Scale_Start[2]                     << endl;
   cout << "PhyCoord_Start[0] = " << PhyCoord_Start[0]                  << endl;
   cout << "PhyCoord_Start[1] = " << PhyCoord_Start[1]                  << endl;
   cout << "PhyCoord_Start[2] = " << PhyCoord_Start[2]                  << endl;
   cout << "Idx_Start[0]      = " << Idx_Start[0]                       << endl;
   cout << "Idx_Start[1]      = " << Idx_Start[1]                       << endl;
   cout << "Idx_Start[2]      = " << Idx_Start[2]                       << endl;
   cout << "Scale_Size[0]     = " << Scale_Size[0]                      << endl;
   cout << "Scale_Size[1]     = " << Scale_Size[1]                      << endl;
   cout << "Scale_Size[2]     = " << Scale_Size[2]                      << endl;
   cout << "PhyCoord_Size[0]  = " << PhyCoord_Size[0]                   << endl;
   cout << "PhyCoord_Size[1]  = " << PhyCoord_Size[1]                   << endl;
   cout << "PhyCoord_Size[2]  = " << PhyCoord_Size[2]                   << endl;
   cout << "Idx_Size[0]       = " << Idx_Size[0]                        << endl;
   cout << "Idx_Size[1]       = " << Idx_Size[1]                        << endl;
   cout << "Idx_Size[2]       = " << Idx_Size[2]                        << endl;
   cout << "==========================================================" << endl;

} // FUNCTION : TakeNote



//-------------------------------------------------------------------------------------------------------
// Function    :  Output
// Description :  Output data 
//-------------------------------------------------------------------------------------------------------
void Output()
{

   if ( MyRank == 0 )   cout << "Output ... " << endl;


   const int  scale   = patch.scale[TargetLevel];
   const real scale_2 = 0.5*scale;
   const real dh_min  = patch.dh[NLEVEL-1];
   const long Size1v  = (long)Idx_MySize[0]*Idx_MySize[1]*Idx_MySize[2];   // total array size of one component

   char FileName_Out[200], FileName_Out_Binary[NOut][200];
   char Info[100];

   sprintf( Info, "_x%.3f-%.3f_y%.3f-%.3f_z%.3f-%.3f_lv%d", 
            Idx_Start[0]*scale*dh_min, (Idx_Start[0]+Idx_Size[0])*scale*dh_min,
            Idx_Start[1]*scale*dh_min, (Idx_Start[1]+Idx_Size[1])*scale*dh_min,
            Idx_Start[2]*scale*dh_min, (Idx_Start[2]+Idx_Size[2])*scale*dh_min,
            TargetLevel );

   switch ( OutputXYZ )
   {
      case 1:  sprintf( FileName_Out, "SliceX" );     break;
      case 2:  sprintf( FileName_Out, "SliceY" );     break;
      case 3:  sprintf( FileName_Out, "SliceZ" );     break;
      case 4:  sprintf( FileName_Out, "ProjX"  );     break;
      case 5:  sprintf( FileName_Out, "ProjY"  );     break;
      case 6:  sprintf( FileName_Out, "ProjZ"  );     break;
      case 7:  sprintf( FileName_Out, "Cube"   );     break;
   }

   strcat( FileName_Out, Info );

   if ( OutputBinary )
   {
      int  NextIdx = NCOMP;
      char comp[NOut][20];

#     if   ( MODEL == HYDRO )
      sprintf( comp[0], "_Binary_Dens" );
      sprintf( comp[1], "_Binary_MomX" );
      sprintf( comp[2], "_Binary_MomY" );
      sprintf( comp[3], "_Binary_MomZ" );
      sprintf( comp[4], "_Binary_Engy" );

#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      sprintf( comp[0], "_Binary_Dens" );
      sprintf( comp[1], "_Binary_Real"  );
      sprintf( comp[2], "_Binary_Imag"  );

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

      if ( OutputPot )     sprintf( comp[NextIdx++], "_Binary_Pot"   );
      if ( OutputDivVel )  sprintf( comp[NextIdx++], "_Binary_DivV"  );
      if ( OutputCurlVel )    
      {
         sprintf( comp[NextIdx++], "_Binary_CurlVx" );
         sprintf( comp[NextIdx++], "_Binary_CurlVy" );
         sprintf( comp[NextIdx++], "_Binary_CurlVz" );
      }

      for (int v=0; v<NOut; v++)
      {
         strcpy( FileName_Out_Binary[v], FileName_Out );
         strcat( FileName_Out_Binary[v], comp[v] );
      }
   }


   if ( Suffix != NULL )   
   {
      if ( OutputBinary )
         for (int v=0; v<NOut; v++)   strcat( FileName_Out_Binary[v], Suffix );
      else
         strcat( FileName_Out, Suffix );
   }


// clean the existing files
   if ( OutputBinary ) // binary files
   {
      if ( MyRank == 0  )
      {
         for (int v=0; v<NOut; v++)
         {
            if ( NULL != fopen(FileName_Out_Binary[v],"r") )  
            {
               fprintf( stderr, "Warning : the file \"%s\" already exists and will be overwritten !!\n", 
                        FileName_Out_Binary[v] );

               FILE *TempFile = fopen( FileName_Out_Binary[v], "w" );
               fclose( TempFile );
            }
         }
      }
   }

   else // text files
   {
      if ( MyRank == 0  &&  NULL != fopen(FileName_Out,"r") )  
      {
         fprintf( stderr, "Warning : the file \"%s\" already exists and will be overwritten !!\n", FileName_Out );

         FILE *TempFile = fopen( FileName_Out, "w" );
         fclose( TempFile );
      }
   }

   MPI_Barrier( MPI_COMM_WORLD );


// output data
   int  ii, jj, kk, NextIdx;
   long ID;
   real x, y, z, u[NOut];

   for (int TargetRank=0; TargetRank<NGPU; TargetRank++)
   {
//    for the output-slice operation, the useless rank will NOT output any data because one of the Idx_MySize[x]
//    will be equal to zero
      if (   MyRank == TargetRank  &&  (      OutputXYZ  < 4  
                                         ||   OutputXYZ == 7
                                         || ( OutputXYZ == 4 && MyRank_X[0] == 0 )
                                         || ( OutputXYZ == 5 && MyRank_X[1] == 0 ) 
                                         || ( OutputXYZ == 6 && MyRank_X[2] == 0 )  )   )
      {
//       output the binary file (different components will be outputted to different files)
         if ( OutputBinary )
         {
            for (int v=0; v<NOut; v++)    
            {
               FILE *File = fopen( FileName_Out_Binary[v], "ab" );

               fwrite( OutputArray+(long)v*Size1v, sizeof(real), Size1v, File );

               fclose( File );
            }
         }

//       output the text file (all components will be outputted to the same file)
         else
         {
            FILE *File = fopen( FileName_Out, "a" );

            if ( TargetRank == 0 )  
            {
#              if   ( MODEL == HYDRO )
               fprintf( File, "%6s %6s %6s %12s %12s %12s %13s %13s %13s %13s %13s %13s", 
                        "i", "j", "k", "x", "y", "z", 
                        "Density", "Momentum.x", "Momentum.y", "Momentum.z", "Energy", "Pressure" );  

               if ( OutputPot )        fprintf( File, " %13s", "Potential" );
               if ( OutputDivVel )     fprintf( File, " %13s", "Div(vel)" );
               if ( OutputCurlVel )    fprintf( File, " %13s %13s %13s", "Curl(vel).x", "Curl(vel).y", 
                                                                         "Curl(vel).z" );
#              elif ( MODEL == MHD )
#              warning : WAIT MHD !!!

#              elif ( MODEL == ELBDM )
               fprintf( File, "%6s %6s %6s %12s %12s %12s %13s %13s %13s", 
                        "i", "j", "k", "x", "y", "z", 
                        "Density", "Real", "Imag" );

               if ( OutputPot )  fprintf( File, " %13s", "Potential" );

#              else
#              error : ERROR : unsupported MODEL !!
#              endif // MODEL

               fprintf( File, "\n" );
            }


            for (int k=0; k<Idx_MySize[2]; k++)  {  kk = ( k + Idx_MyStart[2] )*scale;    z = (kk+scale_2)*dh_min;
            for (int j=0; j<Idx_MySize[1]; j++)  {  jj = ( j + Idx_MyStart[1] )*scale;    y = (jj+scale_2)*dh_min;
            for (int i=0; i<Idx_MySize[0]; i++)  {  ii = ( i + Idx_MyStart[0] )*scale;    x = (ii+scale_2)*dh_min;

               NextIdx = NCOMP; 
               ID      = ( (long)k*Idx_MySize[1] + j )*Idx_MySize[0] + i;

               for (int v=0; v<NOut; v++)    u[v] = OutputArray[ ID + (long)v*Size1v ];

#              if   ( MODEL == HYDRO )
               fprintf( File, "%6d %6d %6d %12.6e %12.6e %12.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e", 
                        ii, jj, kk, x, y, z,
                        u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], 
                        ( u[ENGY]-0.5*(u[MOMX]*u[MOMX]+u[MOMY]*u[MOMY]+u[MOMZ]*u[MOMZ])/u[DENS] )*(GAMMA-1.0) );

               if ( OutputPot )        fprintf( File, " %13.6e", u[NextIdx++] );
               if ( OutputDivVel )     fprintf( File, " %13.6e", u[NextIdx++] );
               if ( OutputCurlVel )    fprintf( File, " %13.6e %13.6e %13.6e", u[NextIdx++], u[NextIdx++], 
                                                                               u[NextIdx++] );

#              elif ( MODEL == MHD )
#              warning : WAIT MHD !!!

#              elif ( MODEL == ELBDM )
               fprintf( File, "%6d %6d %6d %12.6e %12.6e %12.6e %13.6e %13.6e %13.6e", 
                        ii, jj, kk, x, y, z,
                        u[DENS], u[REAL], u[IMAG] );

               if ( OutputPot )  fprintf( File, " %13.6e", u[NextIdx++] );

#              else
#              error : ERROR : unsupported MODEL !!
#              endif // MODEL

               fprintf( File, "\n" );

            }}} // i,j,k

            fclose( File );

         } // if ( OutputBinary ) ... else ...
      } // ( MyRank == TargetRank  &&  ... )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<NGPU; TargetRank++


   if ( MyRank == 0 )   cout << "Output ... done" << endl;

} // FUNCTION : Output



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckParameter
// Description :  Verify the input parameters 
//-------------------------------------------------------------------------------------------------------
void CheckParameter()
{

   if ( MyRank == 0 )   cout << "   CheckParameter ... " << endl;


#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#     error : ERROR : unsupported MODEL !!
#  endif

#  if ( BUF_SIZE != 2 )
#     error : ERROR : currently BUF_SIZE must == 2 !!
#  endif

#  if ( defined OPENMP  &&  !defined _OPENMP )
#     error : ERROR : something is wrong in OpenMP, the macro "_OPENMP" is NOT defined !!
#  endif

   int NRank;
   MPI_Comm_size( MPI_COMM_WORLD, &NRank );
   if ( NRank != NGPU ) 
   {
      fprintf( stderr, "ERROR : number of ranks (%d) != number of GPUs (%d) !!\n", NRank, NGPU );
      MPI_Exit();
   }

   for (int d=0; d<3; d++)
   {
      if ( PhyCoord_Size[d] < patch.dh[NLEVEL-1] )
      {
         fprintf( stderr, "ERROR : incorrect PhyCoord_Size[%d] (%14.7e) --> must be >= finest cell (%14.7e)!!\n", 
                  d, PhyCoord_Size[d], patch.dh[NLEVEL-1] );
         MPI_Exit();
      }

      if ( Scale_Size[d] <= 0 )
      {
         fprintf( stderr, "ERROR : incorrect Scale_Size[%d] (%d) --> must be positive !!\n", d, Scale_Size[d] );
         MPI_Exit();
      }

      if ( Idx_Size[d] <= 0 )
      {
         fprintf( stderr, "ERROR : incorrect Idx_Size[%d] (%d) --> must be positive !!\n", d, Idx_Size[d] );
         MPI_Exit();
      }

      if ( PhyCoord_Start[d] < 0.0  ||  PhyCoord_Start[d] >= patch.BoxSize[d] )
      {
         fprintf( stderr, "ERROR : incorrect PhyCoord_Start[%d] (%14.7e) --> out of range !!\n", 
                  d, PhyCoord_Start[d] );
         MPI_Exit();
      }

      if ( Scale_Start[d] < 0  ||  Scale_Start[d] >= patch.BoxScale[d] )
      {
         fprintf( stderr, "ERROR : incorrect Scale_Start[%d] (%d) --> out of range !!\n", d, Scale_Start[d] );
         MPI_Exit();
      }

      if ( Idx_Start[d] < 0  ||  Idx_Start[d] >= NX0_TOT[d]*(1<<TargetLevel) )
      {
         fprintf( stderr, "ERROR : incorrect Idx_Start[%d] (%d) --> out of range !!\n", d, Idx_Start[d] );
         MPI_Exit();
      }

      if ( PhyCoord_Start[d]+PhyCoord_Size[d] > patch.BoxSize[d] )
      {
         fprintf( stderr, "ERROR : incorrect PhyCoord_Start[%d]+PhyCoord_Size[%d] (%14.7e) --> out of range !!\n",
                  d, d, PhyCoord_Start[d]+PhyCoord_Size[d] );
         MPI_Exit();
      }

      if ( Scale_Start[d]+Scale_Size[d] > patch.BoxScale[d] )
      {
         fprintf( stderr, "ERROR : incorrect Scale_Start[%d]+Scale_Size[%d] (%d) --> out of range !!\n", 
                  d, d, Scale_Start[d]+Scale_Size[d] );
         MPI_Exit();
      }

      if ( Idx_Start[d]+Idx_Size[d] > NX0_TOT[d]*(1<<TargetLevel) )
      {
         fprintf( stderr, "ERROR : incorrect Idx_Start[%d]+Idx_Size[%d] (%d) --> out of range !!\n", 
                  d, d, Idx_Start[d]+Idx_Size[d] );
         MPI_Exit();
      }
   } // for (int d=0; d<3; d++)

   if ( NX0_TOT[0]%(PS2*NGPU_X[0]) != 0  ||  NX0_TOT[1]%(PS2*NGPU_X[1]) != 0  ||  NX0_TOT[2]%(PS2*NGPU_X[2]) != 0 )
   {
      fprintf( stderr, "ERROR : number of base-level patches in each direction in one rank must be \"%s\" !!\n",
                 "a multiple of TWO" );
      MPI_Exit();
   }


   if ( MyRank == 0 )   cout << "   CheckParameter ... done" << endl;

} // FUNCTION : CheckParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  WithinCandidateBox
// Description :  Check whether or not the input patch is within the candidate box 
//
// Parameter   :  Corner   : Pointer to the three corner coordinates  
//                Size     : Size of the targeted patch 
//                Buffer   : Buffer size attached to each side of the candidate box
//-------------------------------------------------------------------------------------------------------
bool WithinCandidateBox( const int *Corner, const int Size, const int Buf )
{

// if ( OldDataFormat )    return true;


   bool Inside[3];
   int cr1[3], cr2[3];

   for (int d=0; d<3; d++)
   {
      Inside[d] = false;
      cr1   [d] = Corner[d];
      cr2   [d] = Corner[d] + Size;
   }


   for (int d=0; d<3; d++)
   {
      for (int s=0; s<3; s++)
      {
         if (  cr1[d] < CanMax[d][s]+Buf  &&  cr2[d] > CanMin[d][s]-Buf  )
         {
            Inside[d] = true;

            break;
         }
      }
   }

   return ( Inside[0] && Inside[1] && Inside[2] );

} // FUNCTION : WithinCandidateBox



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCandidateBox
// Description :  Evaluate the range of the candidate box (considering the periodic boundary condition) 
//-------------------------------------------------------------------------------------------------------
void GetCandidateBox()
{

   if ( MyRank == 0 )   cout << "   GetCandidateBox ... " << flush;


   const int scale = patch.scale[TargetLevel];

   for (int d=0; d<3; d++)    
   for (int s=0; s<3; s++)
   {
      CanMin[d][s] = (Idx_Start[d]            )*scale + patch.BoxScale[d]*(s-1);
      CanMax[d][s] = (Idx_Start[d]+Idx_Size[d])*scale + patch.BoxScale[d]*(s-1);
   }


// set the buffer size to the size of two base-level patches
   CanBuf = 2*PATCH_SIZE*patch.scale[0];


   if ( MyRank == 0 )   cout << "done" << endl;

} // FUNCTION : GetCandidateBox



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TargetDomain
// Description :  Initialize the parameters related to parallelization and the targeted domain
//-------------------------------------------------------------------------------------------------------
void Init_TargetDomain()
{

   if ( MyRank == 0 )   cout << "   Init_TargetDomain ... " << endl;


// set parameters for parallelization
   MyRank_X[0] = MyRank % NGPU_X[0];
   MyRank_X[1] = (MyRank/NGPU_X[0]) % NGPU_X[1];
   MyRank_X[2] = (MyRank/NGPU_X[0]) / NGPU_X[1];

   SibRank[ 0] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]             )          *NGPU_X[0]           + 
                 ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

   SibRank[ 1] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]             )          *NGPU_X[0]           + 
                 ( MyRank_X[0]          +1 )%NGPU_X[0];

   SibRank[ 2] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]             );

   SibRank[ 3] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           +
                 ( MyRank_X[0]             );

   SibRank[ 4] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                 ( MyRank_X[1]             )          *NGPU_X[0]           +
                 ( MyRank_X[0]             );

   SibRank[ 5] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                 ( MyRank_X[1]             )          *NGPU_X[0]           +
                 ( MyRank_X[0]             );

   SibRank[ 6] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

   SibRank[ 7] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]          +1 )%NGPU_X[0];

   SibRank[ 8] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

   SibRank[ 9] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]          +1 )%NGPU_X[0];
   
   SibRank[10] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           +
                 ( MyRank_X[0]             );

   SibRank[11] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]             );

   SibRank[12] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           +
                 ( MyRank_X[0]             );             

   SibRank[13] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           +
                 ( MyRank_X[0]             );

   SibRank[14] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                 ( MyRank_X[1]             )          *NGPU_X[0]           + 
                 ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

   SibRank[15] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]             )*NGPU_X[0]                     + 
                 ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];
   
   SibRank[16] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]             )          *NGPU_X[0]           + 
                 ( MyRank_X[0]          +1 )%NGPU_X[0];
   
   SibRank[17] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                 ( MyRank_X[1]             )          *NGPU_X[0]           + 
                 ( MyRank_X[0]          +1 )%NGPU_X[0];
   
   SibRank[18] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];
   
   SibRank[19] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]          +1 )%NGPU_X[0];
   
   SibRank[20] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];
   
   SibRank[21] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]          +1 )%NGPU_X[0];
   
   SibRank[22] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];
   
   SibRank[23] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                 ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]          +1 )%NGPU_X[0];
   
   SibRank[24] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];
   
   SibRank[25] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] + 
                 ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           + 
                 ( MyRank_X[0]          +1 )%NGPU_X[0];


// set the default (start/size) of (cell scales/physical coordinates)
   if ( InputScale )
   {
      for (int d=0; d<3; d++)
      {
         if ( Scale_Start[d] == WRONG )   Scale_Start[d] = 0;
         if ( Scale_Size [d] == WRONG )   Scale_Size [d] = patch.BoxScale[d] - Scale_Start[d];

         PhyCoord_Start[d] = patch.BoxSize[d]*( (double)Scale_Start[d]/(double)patch.BoxScale[d] );
         PhyCoord_Size [d] = patch.BoxSize[d]*( (double)Scale_Size [d]/(double)patch.BoxScale[d] );
      }
   }

   else
   {
      for (int d=0; d<3; d++)
      {
         if ( PhyCoord_Start[d] == WRONG )   PhyCoord_Start[d] = 0.0;
         if ( PhyCoord_Size [d] == WRONG )   PhyCoord_Size [d] = patch.BoxSize[d] - PhyCoord_Start[d];

         Scale_Start[d] = int(patch.BoxScale[d]*PhyCoord_Start[d]/patch.BoxSize[d]);
         Scale_Size [d] = int(patch.BoxScale[d]*PhyCoord_Size [d]/patch.BoxSize[d]);
      }
   }

   switch ( OutputXYZ )
   {
      case 1 : Scale_Size[0] = 1;   PhyCoord_Size[0] = patch.dh[NLEVEL-1];    break;
      case 2 : Scale_Size[1] = 1;   PhyCoord_Size[1] = patch.dh[NLEVEL-1];    break;
      case 3 : Scale_Size[2] = 1;   PhyCoord_Size[2] = patch.dh[NLEVEL-1];    break;
   }


// set the targeted array indices and ranges
   for (int d=0; d<3; d++)    
   {
      Idx_Start[d] = Scale_Start[d] / patch.scale[TargetLevel];
      Idx_Size [d] = (Scale_Size[d]-1)/patch.scale[TargetLevel] + 1;
   }


// set the local array size and the starting indices
   int Min, Max, NGrid, X_End;

   for (int d=0; d<3; d++)
   {
      NGrid = NX0[d]*(1<<TargetLevel);
      Min   = (MyRank_X[d]  )*NGrid;
      Max   = (MyRank_X[d]+1)*NGrid-1;

//    check whether or not the sub-domain is within the targeted range
      if ( Min < Idx_Start[d]+Idx_Size[d]  &&  Max >= Idx_Start[d] )
      {
         Idx_MyStart[d] = ( Idx_Start[d]             > Min ) ? Idx_Start[d]-Min : 0;
         X_End          = ( Idx_Start[d]+Idx_Size[d] > Max ) ? NGrid-1 : NGrid-1-(Max-Idx_Start[d]-Idx_Size[d]+1);
         Idx_MySize [d] = X_End - Idx_MyStart[d] + 1;
         Idx_MyStart[d] += MyRank_X[d]*NGrid;

         if (   (  d == 0  &&  ( OutputXYZ == 1 || OutputXYZ == 4 )  )  ||
                (  d == 1  &&  ( OutputXYZ == 2 || OutputXYZ == 5 )  )  ||
                (  d == 2  &&  ( OutputXYZ == 3 || OutputXYZ == 6 )  )     )     Idx_MySize[d] = 1;


//       verify the variables "Idx_MyStart and Idx_MySize"
         if ( Idx_MyStart[d] < Idx_Start[d]  ||  Idx_MyStart[d] >= Idx_Start[d]+Idx_Size[d] )
         {
            fprintf( stderr, "ERROR : incorrect Idx_MyStart[%d] (%d) --> out of range !!\n", d, Idx_MyStart[d] );
            MPI_Exit();
         }

         if ( Idx_MySize[d] < 0  ||  Idx_MySize[d] > NX0[d]*(1<<TargetLevel) ) 
         {
            fprintf( stderr, "ERROR : incorrect Idx_MySize[%d] (%d) --> out of range !!\n", d, Idx_MySize[d] );
            MPI_Exit();
         }

         if ( Idx_MyStart[d]+Idx_MySize[d] > Idx_Start[d]+Idx_Size[d] )
         {
            fprintf( stderr, "ERROR : incorrect Idx_MyStart[%d]+Idx_MySize[%d] (%d) --> out of range (%d) !!\n", 
                     d, d, Idx_MyStart[d]+Idx_MySize[d], Idx_Start[d]+Idx_Size[d] );
            MPI_Exit();
         }

      } // if ( Min < Idx_Start[d]+Idx_Size[d]  &&  Max >= Idx_Start[d] )

      else
      {
         Idx_MyStart  [d] = 0;
         Idx_MySize[d] = ( OutputXYZ == d+4 ) ? 1 : 0;    // set Idx_MySize[X] = 1 for the function "SumOverRanks"
      } 
   } // for (int d=0; d<3; d++)


   if ( MyRank == 0 )   cout << "   Init_TargetDomain ... done" << endl;

} // FUNCTION : Init_TargetDomain



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MPI
// Description :  Initialize MPI and parameters for parallelization 
//-------------------------------------------------------------------------------------------------------
void Init_MPI( int *argc, char ***argv )
{

   if ( MPI_Init( argc, argv ) != MPI_SUCCESS )
   {
      cerr << "MPI_Init failed !!" << endl;
      exit( 1 );
   }

   if ( MPI_Comm_rank( MPI_COMM_WORLD, &MyRank ) != MPI_SUCCESS )
   {
      cerr << "MPI_Comm_rank failed !!" << endl;
      exit( 1 );
   }

   if ( MyRank == 0 )   cout << "Init_MPI ... done" << endl;;

} // FUNCTION : Init_MPI



//-------------------------------------------------------------------------------------------------------
// Function    :  PreparePatch
// Description :  Prepare the central and buffer data of a single patch 
//
// Note        :  If any of the sibling patch of PID does NOT exist, it will interpolate on the CData to fill up
//                the buffer data in FData
//
// Parameter   :  lv       : Refinement level of the targeted patch 
//                PID      : Patch index of the targeted patch
//                Buffer   : Size of buffer (ghost zone)
//                FData    : Array to store the prepared data ( width = PATCH_SIZE + 2*Buffer )
//                CData    : Array for interpolation ( width = PATCH_SIZE + 2*Buffer )
//-------------------------------------------------------------------------------------------------------
void PreparePatch( const int lv, const int PID, const int Buffer, real FData[], const real CData[] )
{

   if ( Buffer != 2 )
   {
      fprintf( stderr, "ERROR : currently Buffer must == 2 in the function \"%s\" !!\n", __FUNCTION__ );
      MPI_Exit();
   }

   const int  Size = PATCH_SIZE + 2*Buffer;
   const int  dii  = 1;
   const int  djj  = Size;
   const long dkk  = Size*Size;
   const long dvv  = Size*Size*Size;

   int  i0, j0, k0, ii0, jj0, kk0, i_loop, j_loop, k_loop;
   int  SibPID, LocalID;
   long ID, ID2;
   real Slope[3]={0,0,0}, LSlope[3], RSlope[3];

// prepare the central data
   if ( patch.ptr[lv][PID]->fluid == NULL )
   {
      fprintf( stderr, "ERROR : lv %d, PID %d, the fluid array is NOT allocated !!\n", lv, PID );
      MPI_Exit();
   }

   if ( OutputPot  &&  patch.ptr[lv][PID]->pot == NULL )
   {
      fprintf( stderr, "ERROR : lv %d, PID %d, the potential array is NOT allocated !!\n", lv, PID );
      MPI_Exit();
   }

   for (int v=0; v<NCOMP; v++)         
   for (int k=0, kk=Buffer; k<PATCH_SIZE; k++, kk++)
   for (int j=0, jj=Buffer; j<PATCH_SIZE; j++, jj++)
   for (int i=0, ii=Buffer; i<PATCH_SIZE; i++, ii++)
   {
      ID        = (long)v*dvv + kk*dkk + jj*djj + ii;
      FData[ID] = patch.ptr[lv][PID]->fluid[v][k][j][i];
   }

   if ( OutputPot )
   for (int k=0, kk=Buffer; k<PATCH_SIZE; k++, kk++)
   for (int j=0, jj=Buffer; j<PATCH_SIZE; j++, jj++)
   for (int i=0, ii=Buffer; i<PATCH_SIZE; i++, ii++)
   {
      ID        = (long)NCOMP*dvv + kk*dkk + jj*djj + ii; 
      FData[ID] = patch.ptr[lv][PID]->pot[k][j][i];
   }
   
   
// prepare the buffer data
   for (int sib=0; sib<26; sib++)
   {
      SibPID = patch.ptr[lv][PID]->sibling[sib];

      if ( SibPID != -1 )
      {
         if ( patch.ptr[lv][SibPID]->fluid == NULL )
         {
            fprintf( stderr, "ERROR : lv %d, SibPID %d, the fluid array is NOT allocated !!\n", lv, SibPID );
            MPI_Exit();
         }

         if ( OutputPot  &&  patch.ptr[lv][SibPID]->pot == NULL )
         {
            fprintf( stderr, "ERROR : lv %d, SibPID %d, the potential array is NOT allocated !!\n", lv, SibPID );
            MPI_Exit();
         }

         i0     = TABLE_01( sib, 'x', PATCH_SIZE-Buffer, 0, 0 );
         j0     = TABLE_01( sib, 'y', PATCH_SIZE-Buffer, 0, 0 );
         k0     = TABLE_01( sib, 'z', PATCH_SIZE-Buffer, 0, 0 );
         ii0    = TABLE_01( sib, 'x', 0, Buffer, PATCH_SIZE+Buffer );
         jj0    = TABLE_01( sib, 'y', 0, Buffer, PATCH_SIZE+Buffer );
         kk0    = TABLE_01( sib, 'z', 0, Buffer, PATCH_SIZE+Buffer );
         i_loop = TABLE_01( sib, 'x', Buffer, PATCH_SIZE, Buffer );
         j_loop = TABLE_01( sib, 'y', Buffer, PATCH_SIZE, Buffer );
         k_loop = TABLE_01( sib, 'z', Buffer, PATCH_SIZE, Buffer );
   
         for (int v=0; v<NCOMP; v++)
         for (int k=k0, kk=kk0; k<k0+k_loop; k++, kk++)
         for (int j=j0, jj=jj0; j<j0+j_loop; j++, jj++)
         for (int i=i0, ii=ii0; i<i0+i_loop; i++, ii++)
         {
            ID        = (long)v*dvv + kk*dkk + jj*djj + ii;
            FData[ID] = patch.ptr[lv][SibPID]->fluid[v][k][j][i];
         }

         if ( OutputPot )
         for (int k=k0, kk=kk0; k<k0+k_loop; k++, kk++)
         for (int j=j0, jj=jj0; j<j0+j_loop; j++, jj++)
         for (int i=i0, ii=ii0; i<i0+i_loop; i++, ii++)
         {
            ID        = (long)NCOMP*dvv + kk*dkk + jj*djj + ii;
            FData[ID] = patch.ptr[lv][SibPID]->pot[k][j][i];
         }
      }

      else
      {
//       the CData array must be properly prepared if any of the sibling patch does NOT exist
         if ( CData == NULL )
         {
            fprintf( stderr, "ERROR : CData == NULL, Rank = %d, lv = %d, PID = %d !!\n", MyRank, lv, PID );
            MPI_Exit();
         }


//###IMPROVE : Buffer != 2 will NOT work here !!
         LocalID = PID%8;
         i0      = TABLE_01( sib, 'x', 1, Buffer, Buffer+PS/2 ) + TABLE_02( LocalID, 'x', 0, PS/2 );
         j0      = TABLE_01( sib, 'y', 1, Buffer, Buffer+PS/2 ) + TABLE_02( LocalID, 'y', 0, PS/2 );
         k0      = TABLE_01( sib, 'z', 1, Buffer, Buffer+PS/2 ) + TABLE_02( LocalID, 'z', 0, PS/2 );
         ii0     = TABLE_01( sib, 'x', 0, Buffer, Buffer+PS );
         jj0     = TABLE_01( sib, 'y', 0, Buffer, Buffer+PS );
         kk0     = TABLE_01( sib, 'z', 0, Buffer, Buffer+PS );
         i_loop  = TABLE_01( sib, 'x', Buffer-1, PS/2, Buffer-1 );
         j_loop  = TABLE_01( sib, 'y', Buffer-1, PS/2, Buffer-1 );
         k_loop  = TABLE_01( sib, 'z', Buffer-1, PS/2, Buffer-1 );

         for (int v=0; v<NLoad; v++)
         for (int k=k0, kk=kk0; k<k0+k_loop; k++, kk+=2)
         for (int j=j0, jj=jj0; j<j0+j_loop; j++, jj+=2)
         for (int i=i0, ii=ii0; i<i0+i_loop; i++, ii+=2)
         {
            ID  = (long)v*dvv +  k*dkk +  j*djj +  i; 
            ID2 = (long)v*dvv + kk*dkk + jj*djj + ii; 

            switch ( InterScheme )
            {
               case 0 : // MinMod limiter
               {
                  LSlope[0] = CData[ID    ] - CData[ID-dii];
                  RSlope[0] = CData[ID+dii] - CData[ID    ];
                  if ( LSlope[0]*RSlope[0] <= 0.0 )  
                     Slope[0] = 0.0;
                  else                            
                     Slope[0] = 0.25*(  fabs( LSlope[0] ) < fabs( RSlope[0] ) ? LSlope[0] : RSlope[0]  );
                  
                  LSlope[1] = CData[ID    ] - CData[ID-djj];
                  RSlope[1] = CData[ID+djj] - CData[ID    ];
                  if ( LSlope[1]*RSlope[1] <= 0.0 )  
                     Slope[1] = 0.0;
                  else                            
                     Slope[1] = 0.25*(  fabs( LSlope[1] ) < fabs( RSlope[1] ) ? LSlope[1] : RSlope[1]  );


                  LSlope[2] = CData[ID    ] - CData[ID-dkk];
                  RSlope[2] = CData[ID+dkk] - CData[ID    ];
                  if ( LSlope[2]*RSlope[2] <= 0.0 )  
                     Slope[2] = 0.0;
                  else                            
                     Slope[2] = 0.25*(  fabs( LSlope[2] ) < fabs( RSlope[2] ) ? LSlope[2] : RSlope[2]  );
               }
               break;


               case 1 : // central interpolation
               {
                  Slope[0] = 0.125 * ( CData[ID+dii] - CData[ID-dii] );
                  Slope[1] = 0.125 * ( CData[ID+djj] - CData[ID-djj] );
                  Slope[2] = 0.125 * ( CData[ID+dkk] - CData[ID-dkk] );
               }
               break;


               default :
               {
                  cerr << "ERROR : unsupported interpolation scheme !!" << endl;
                  MPI_Exit();
               }

            } // switch ( InterScheme );


            FData[ID2            ] = CData[ID] - Slope[0] - Slope[1] - Slope[2];
            FData[ID2+dii        ] = CData[ID] + Slope[0] - Slope[1] - Slope[2];
            FData[ID2    +djj    ] = CData[ID] - Slope[0] + Slope[1] - Slope[2];
            FData[ID2        +dkk] = CData[ID] - Slope[0] - Slope[1] + Slope[2];
            FData[ID2+dii+djj    ] = CData[ID] + Slope[0] + Slope[1] - Slope[2];
            FData[ID2    +djj+dkk] = CData[ID] - Slope[0] + Slope[1] + Slope[2];
            FData[ID2+dii    +dkk] = CData[ID] + Slope[0] - Slope[1] + Slope[2];
            FData[ID2+dii+djj+dkk] = CData[ID] + Slope[0] + Slope[1] + Slope[2];

         } // i,j,k,v

      } // if ( SibPID != -1 ) ... else ...
   } // for (int sib=0; sib<26; sib++)

} // FUNCTION : PreparePatch



//-------------------------------------------------------------------------------------------------------
// Function    :  Interpolation
// Description :  Interpolate on the CData array and store the result in the FData array 
//
// Parameter   :  CSize : Width of the CData array 
//                CData : Array for interpolation 
//                FSize : Width of the FData array 
//                FData : Array to store the interpolated data 
//-------------------------------------------------------------------------------------------------------
void Interpolation( const int CSize, const real CData[], const int FSize, real FData[] )
{

   const int  Cdi = 1;
   const int  Cdj = CSize;
   const long Cdk = CSize*CSize;
   const long Cdv = CSize*CSize*CSize;
   const int  Fdi = 1;
   const int  Fdj = FSize;
   const long Fdk = FSize*FSize;
   const long Fdv = FSize*FSize*FSize;

   int  ii, jj, kk;
   long CID, FID;
   real Slope[3]={0,0,0}, LSlope[3], RSlope[3];


   for (int v=0; v<NLoad; v++)
   for (int k=1, kk=(k-1)*2; k<CSize-1; k++, kk=(k-1)*2)
   for (int j=1, jj=(j-1)*2; j<CSize-1; j++, jj=(j-1)*2)
   for (int i=1, ii=(i-1)*2; i<CSize-1; i++, ii=(i-1)*2)
   {
      CID = (long)v*Cdv +  k*Cdk +  j*Cdj +  i;
      FID = (long)v*Fdv + kk*Fdk + jj*Fdj + ii;

      switch ( InterScheme )
      {
         case 0 : // MinMod limiter
         {
            LSlope[0] = CData[CID    ] - CData[CID-Cdi];
            RSlope[0] = CData[CID+Cdi] - CData[CID    ];
            if ( LSlope[0]*RSlope[0] <= 0.0 )  
               Slope[0] = 0.0;
            else
               Slope[0] = 0.25*(  fabs( LSlope[0] ) < fabs( RSlope[0] ) ? LSlope[0] : RSlope[0]  );
            
            LSlope[1] = CData[CID    ] - CData[CID-Cdj];
            RSlope[1] = CData[CID+Cdj] - CData[CID    ];
            if ( LSlope[1]*RSlope[1] <= 0.0 )  
               Slope[1] = 0.0;
            else                                  
               Slope[1] = 0.25*(  fabs( LSlope[1] ) < fabs( RSlope[1] ) ? LSlope[1] : RSlope[1]  );


            LSlope[2] = CData[CID    ] - CData[CID-Cdk];
            RSlope[2] = CData[CID+Cdk] - CData[CID    ];
            if ( LSlope[2]*RSlope[2] <= 0.0 )  
               Slope[2] = 0.0;
            else                                  
               Slope[2] = 0.25*(  fabs( LSlope[2] ) < fabs( RSlope[2] ) ? LSlope[2] : RSlope[2]  );
         }
         break;


         case 1 : // central interpolation
         {
            Slope[0] = 0.125 * ( CData[CID+Cdi] - CData[CID-Cdi] );
            Slope[1] = 0.125 * ( CData[CID+Cdj] - CData[CID-Cdj] );
            Slope[2] = 0.125 * ( CData[CID+Cdk] - CData[CID-Cdk] );
         }
         break;


         default :
         {
            cerr << "ERROR : unsupported interpolation scheme !!" << endl;
            MPI_Exit();
         }

      } // switch ( InterScheme )


      FData[FID            ] = CData[CID] - Slope[0] - Slope[1] - Slope[2];
      FData[FID+Fdi        ] = CData[CID] + Slope[0] - Slope[1] - Slope[2];
      FData[FID    +Fdj    ] = CData[CID] - Slope[0] + Slope[1] - Slope[2];
      FData[FID        +Fdk] = CData[CID] - Slope[0] - Slope[1] + Slope[2];
      FData[FID+Fdi+Fdj    ] = CData[CID] + Slope[0] + Slope[1] - Slope[2];
      FData[FID    +Fdj+Fdk] = CData[CID] - Slope[0] + Slope[1] + Slope[2];
      FData[FID+Fdi    +Fdk] = CData[CID] + Slope[0] - Slope[1] + Slope[2];
      FData[FID+Fdi+Fdj+Fdk] = CData[CID] + Slope[0] + Slope[1] + Slope[2];

   } // i,j,k,v

} // FUNCTION : Interpolation



//-------------------------------------------------------------------------------------------------------
// Function    :  StoreData
// Description :  Store data in the array "Out"
//
// Note        :  Projection is also performed in this function
//
// Parameter   :  lv       : Refinement level of the input patch 
//                PID      : Pathc index of the input patch
//                FData    : Array storing the data to be outputted
//                Buffer   : Size of buffer
//                Out      : Array to store the results 
//-------------------------------------------------------------------------------------------------------
void StoreData( const int lv, const int PID, real FData[], const int Buffer, real *Out )
{

   const long Size1v    = Idx_MySize[0]*Idx_MySize[1]*Idx_MySize[2];  // total array size of one component
   const int  PatchSize = PATCH_SIZE*( 1<<(TargetLevel-lv) );
   const int  FSize     = PatchSize + 2*Buffer;
   const int  Corner[3] = { patch.ptr[lv][PID]->corner[0]/patch.scale[TargetLevel], 
                            patch.ptr[lv][PID]->corner[1]/patch.scale[TargetLevel], 
                            patch.ptr[lv][PID]->corner[2]/patch.scale[TargetLevel]  };

   int    ijk_min[3], ijk_max[3], NSum_local, ProjDir, ii, jj, kk, Stride[3]={0,0,0};
   long   Jump, ID1, ID2;
   double SumData, NAve;


// calculate the index range
   for (int d=0; d<3; d++)
   {
      if ( Corner[d] < Idx_Start[d]+Idx_Size[d]  &&  Corner[d]+PatchSize > Idx_Start[d] )
      {
         ijk_min[d] = ( Idx_Start[d] > Corner[d] ) ? Idx_Start[d]-Corner[d]+Buffer : Buffer;
         ijk_max[d] = ( Idx_Start[d]+Idx_Size[d] >= Corner[d]+PatchSize ) ? 
            Buffer+PatchSize-1 : Buffer+PatchSize-1-(Corner[d]+PatchSize-Idx_Start[d]-Idx_Size[d]);
      }

      else
         return;
   }


// projection
   if ( OutputXYZ == 4  ||  OutputXYZ == 5  ||  OutputXYZ == 6 )
   {
      ProjDir          = OutputXYZ - 4; 
      NSum_local       = ijk_max[ProjDir] - ijk_min[ProjDir] + 1; 
      ijk_max[ProjDir] = ijk_min[ProjDir];   
      NAve             = (double)Idx_Size[ProjDir];
      Jump             = 1;

      for (int t=0; t<ProjDir; t++)    Jump *= FSize;

      for (int v=0; v<NOut; v++)      
      for (int k=ijk_min[2]; k<=ijk_max[2]; k++)
      for (int j=ijk_min[1]; j<=ijk_max[1]; j++)
      for (int i=ijk_min[0]; i<=ijk_max[0]; i++)
      {
         SumData = 0.0;
         ID1     = (((long)v*FSize + k)*FSize + j)*FSize + i;

         for (int t=0; t<NSum_local; t++)    SumData += FData[ ID1 + (long)t*Jump ];

         FData[ID1] = SumData / NAve;
      }
   }


// store data in the OutputArray array
   switch ( OutputXYZ )
   {
      case 1 : case 4 :    Stride[0] = 0;    Stride[1] = 1;             Stride[2] = Idx_MySize[1];
                           break;

      case 2 : case 5 :    Stride[0] = 1;    Stride[1] = 0;             Stride[2] = Idx_MySize[0];
                           break;

      case 3 : case 6 :    Stride[0] = 1;    Stride[1] = Idx_MySize[0];  Stride[2] = 0;
                           break;

      case 7 :             Stride[0] = 1;    Stride[1] = Idx_MySize[0];  Stride[2] = Idx_MySize[1]*Idx_MySize[0];
                           break;
   }

   for (int v=0; v<NOut; v++)                   {  
   for (int k=ijk_min[2]; k<=ijk_max[2]; k++)   {  kk  = k - Buffer + Corner[2] - Idx_MyStart[2];
   for (int j=ijk_min[1]; j<=ijk_max[1]; j++)   {  jj  = j - Buffer + Corner[1] - Idx_MyStart[1];
   for (int i=ijk_min[0]; i<=ijk_max[0]; i++)   {  ii  = i - Buffer + Corner[0] - Idx_MyStart[0];
                                                   ID1 = (((long)v*FSize + k)*FSize + j)*FSize + i;
                                                   ID2 = (long)v*Size1v + kk*Stride[2] + jj*Stride[1] + ii*Stride[0];

      Out[ID2] += FData[ID1];

   }}}}

} // FUNCTION : StoreData



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCurlVel 
// Description :  Evaluate the curl( velocity ) == vorticity
//
// Parameter   :  FData    : Array to store the output data
//                FSize    : Width of the FData array 
//                Buffer   : Size of buffer
//                NextIdx  : Index to store the evaluated results
//                dh       : Cell size
//-------------------------------------------------------------------------------------------------------
void GetCurlVel( real FData[], const int FSize, const int Buffer, const int NextIdx, const real dh )
{

   const int   di  = 1;
   const int   dj  = FSize;
   const long  dk  = FSize*FSize;
   const long  dv  = FSize*FSize*FSize;
   const real _dh2 = 0.5/dh;

   const real *Rho  = FData +           0*dv;
   const real *MomX = FData +           1*dv;
   const real *MomY = FData +           2*dv;
   const real *MomZ = FData +           3*dv;
         real *Wx   = FData + (NextIdx+0)*dv;
         real *Wy   = FData + (NextIdx+1)*dv;
         real *Wz   = FData + (NextIdx+2)*dv;

   real *Vx = new real [dv];
   real *Vy = new real [dv];
   real *Vz = new real [dv];

   real _Rho, PxVy, PxVz, PyVx, PyVz, PzVx, PzVy;
   long Idx;


// get velocity
   for (int t=0; t<dv; t++)
   {
      _Rho  = 1.0/Rho[t];
      Vx[t] = _Rho*MomX[t];
      Vy[t] = _Rho*MomY[t];
      Vz[t] = _Rho*MomZ[t];
   }


//###OPTIMIZATION : perform calculations only for cells lying on the targeted slice
   for (int k=Buffer; k<FSize-Buffer; k++)
   for (int j=Buffer; j<FSize-Buffer; j++)
   for (int i=Buffer; i<FSize-Buffer; i++)
   {
      Idx  = (long)k*dk + j*dj + i*di;
      PxVy = Vy[Idx+di] - Vy[Idx-di];
      PxVz = Vz[Idx+di] - Vz[Idx-di];
      PyVx = Vx[Idx+dj] - Vx[Idx-dj];
      PyVz = Vz[Idx+dj] - Vz[Idx-dj];
      PzVx = Vx[Idx+dk] - Vx[Idx-dk];
      PzVy = Vy[Idx+dk] - Vy[Idx-dk];

      Wx[Idx] = _dh2*( PyVz - PzVy );
      Wy[Idx] = _dh2*( PzVx - PxVz );
      Wz[Idx] = _dh2*( PxVy - PyVx );
   }


   delete [] Vx;
   delete [] Vy;
   delete [] Vz;

} // FUNCTION : GetCurlVel



//-------------------------------------------------------------------------------------------------------
// Function    :  GetDivVel 
// Description :  Evaluate the divergence( velocity )
//
// Parameter   :  FData    : Array to store the output data
//                FSize    : Width of the FData array 
//                Buffer   : Size of buffer
//                NextIdx  : Index to store the evaluated results
//                dh       : Cell size
//-------------------------------------------------------------------------------------------------------
void GetDivVel( real FData[], const int FSize, const int Buffer, const int NextIdx, const real dh )
{

   const int   di  = 1;
   const int   dj  = FSize;
   const long  dk  = FSize*FSize;
   const long  dv  = FSize*FSize*FSize;
   const real _dh2 = 0.5/dh;

   const real *Rho  = FData +       0*dv;
   const real *MomX = FData +       1*dv;
   const real *MomY = FData +       2*dv;
   const real *MomZ = FData +       3*dv;
         real *Out  = FData + NextIdx*dv;

   real *Vx = new real [dv];
   real *Vy = new real [dv];
   real *Vz = new real [dv];

   real _Rho;
   long Idx;


// get velocity
   for (int t=0; t<dv; t++)
   {
      _Rho  = 1.0/Rho[t];
      Vx[t] = _Rho*MomX[t];
      Vy[t] = _Rho*MomY[t];
      Vz[t] = _Rho*MomZ[t];
   }


//###OPTIMIZATION : perform calculations only for cells lying on the targeted slice
   for (int k=Buffer; k<FSize-Buffer; k++)
   for (int j=Buffer; j<FSize-Buffer; j++)
   for (int i=Buffer; i<FSize-Buffer; i++)
   {
      Idx      = (long)k*dk + j*dj + i*di;
      Out[Idx] = _dh2*( ( Vz[Idx+dk] - Vz[Idx-dk] ) + ( Vy[Idx+dj] - Vy[Idx-dj] ) + ( Vx[Idx+di] - Vx[Idx-di] ) );
   }


   delete [] Vx;
   delete [] Vy;
   delete [] Vz;

} // FUNCTION : GetDivVel



//-------------------------------------------------------------------------------------------------------
// Function    :  Refine2TargetLevel
// Description :  Refine to the level "TargetLevel" patch by patch 
//-------------------------------------------------------------------------------------------------------
void Refine2TargetLevel()
{

   if ( MyRank == 0 )   cout << "Refine2TargetLevel ... " << endl;


   const int Buffer = BUF_SIZE;
   int   FaPID, FSize, CSize, NextIdx;
   real *CData = NULL, *FData = NULL;

#  ifdef OPENMP
   const long OutSize = (long)NOut*Idx_MySize[0]*Idx_MySize[1]*Idx_MySize[2];
   real **OutputArray_OMP = NULL;
   int TID;

   if ( OutputXYZ == 4  ||  OutputXYZ == 5  ||  OutputXYZ == 6 )
   {
      OutputArray_OMP = new real* [OMP_NThread];

      for (int t=0; t<OMP_NThread; t++)   
      {
         OutputArray_OMP[t] = new real [OutSize];

         for (long i=0; i<OutSize; i++)   OutputArray_OMP[t][i] = 0.0;
      }
   }
#  endif // #ifdef OPENMP


   for (int lv=0; lv<=TargetLevel; lv++)
   {
      if ( MyRank == 0 )   cout << "  Level = " << lv << " ... " << flush; 

#     pragma omp parallel private( FaPID, FSize, CSize, NextIdx, CData, FData, TID )
      {
#        ifdef OPENMP
         TID = omp_get_thread_num();
#        endif

#        pragma omp for
         for (int PID=0; PID<NPatchComma[lv][1]; PID++)
         {

            if (  ( patch.ptr[lv][PID]->son == -1 || lv == TargetLevel )  &&  
                  WithinCandidateBox( patch.ptr[lv][PID]->corner, PATCH_SIZE*patch.scale[lv], 0 )  )
            {

//             allocate memory
               CSize = PATCH_SIZE + 2*Buffer;
               FSize = PATCH_SIZE + 2*Buffer;
               CData = new real [ (long)NOut*CSize*CSize*CSize ];
               FData = new real [ (long)NOut*FSize*FSize*FSize ];


//###OPTIMIZATION : prepare the coarse-grid data "only once" for all son patches
//             prepare the coarse-grid data ( level = lv-1 )
               if ( lv != 0 )
               {
                  FaPID = patch.ptr[lv][PID]->father;

                  PreparePatch( lv-1, FaPID, Buffer, CData, NULL );
               } 


//             prepare the fine-grid data ( level = lv )
               PreparePatch( lv, PID, Buffer, FData, CData );


//             perform interpolation to reach the TagetLevel resolution
               for (int FineLv=lv+1; FineLv<=TargetLevel; FineLv++)         
               {
//                reset and reallocate CData and FData pointers
                  delete [] CData;
                  CSize = FSize;
                  CData = FData;
                  FSize = ( FSize - 2*Buffer )*2 + 2*Buffer;
                  FData = new real [ (long)NOut*FSize*FSize*FSize ];

//                interpolate from CData (level = FineLv-1) to FData (level = FineLv)
                  Interpolation( CSize, CData, FSize, FData );
               }

  
//             calculate the divergence(vel)
               NextIdx = NLoad;
               if ( OutputDivVel )     
               {
                  GetDivVel ( FData, FSize, Buffer, NextIdx, patch.dh[TargetLevel] );
                  NextIdx ++;
               }


//             calculate the curl(vel)
               if ( OutputCurlVel )    
               {
                  GetCurlVel( FData, FSize, Buffer, NextIdx, patch.dh[TargetLevel] );
                  NextIdx += 3;
               }


//             store the final result
#              ifdef OPENMP
               if ( OutputXYZ == 4  ||  OutputXYZ == 5  ||  OutputXYZ == 6 )  
               StoreData( lv, PID, FData, Buffer, OutputArray_OMP[TID] );
               else
               StoreData( lv, PID, FData, Buffer, OutputArray );

#              else
               StoreData( lv, PID, FData, Buffer, OutputArray );
#              endif


//             free memory
               delete [] CData;
               delete [] FData;
               CData = NULL;
               FData = NULL;

            } // if ( patch.ptr[lv][PID]->son == -1 ) ...
         } // for (int PID=0; PID<NPatchComma[lv][1]; PID++)
      } // OpenMP parallel region

      MPI_Barrier( MPI_COMM_WORLD );
      if ( MyRank == 0 )   cout << "done" << endl;

   } // for (int lv=0; lv<=TargetLevel; lv++)


#  ifdef OPENMP
   if ( OutputXYZ == 4  ||  OutputXYZ == 5  ||  OutputXYZ == 6 )
   {
      for (int t=0; t<OMP_NThread; t++)   
      {
         for (int i=0; i<OutSize; i++)    OutputArray[i] += OutputArray_OMP[t][i];

         delete [] OutputArray_OMP[t];
      }

      delete [] OutputArray_OMP;
   }
#  endif // #ifdef OPENMP


   if ( MyRank == 0 )   cout << "Refine2TargetLevel ... done" << endl;

} // FUNCTION : Refine2TargetLevel



//-------------------------------------------------------------------------------------------------------
// Function    :  AllocateOutputArray
// Description :  Allocate memory for the array "OutputArray"
//-------------------------------------------------------------------------------------------------------
void AllocateOutputArray() 
{

   const long Size = (long)NOut*Idx_MySize[0]*Idx_MySize[1]*Idx_MySize[2];

   OutputArray = new real [Size];

// initialize it as zero (necessary for the projection operation)
   for (long t=0; t<Size; t++)    OutputArray[t] = 0.0;

} // FUNCTION : AllocateOutputArray



//-------------------------------------------------------------------------------------------------------
// Function    :  SumOverRanks 
// Description :  When doing projection, this function will collect the projection results from different ranks 
//-------------------------------------------------------------------------------------------------------
void SumOverRanks()
{

   if ( MyRank == 0 )   cout << "SumOverRanks ... " << flush;


   const long Size = (long)NOut*Idx_MySize[0]*Idx_MySize[1]*Idx_MySize[2];

   int  RecvRank=0, NSend=0;
   int  *SendRank=NULL;
   real *RecvBuffer=NULL;

   switch ( OutputXYZ )
   {
      case 4:     NSend    = NGPU_X[0]-1;
                  RecvRank = MyRank_X[2]*NGPU_X[1]*NGPU_X[0] + MyRank_X[1]*NGPU_X[0];
                  SendRank = new int [NSend];

                  for (int ID=0; ID<NSend; ID++)   
                     SendRank[ID] = MyRank_X[2]*NGPU_X[1]*NGPU_X[0] + MyRank_X[1]*NGPU_X[0] + (ID+1);

                  break;

      case 5:     NSend    = NGPU_X[1]-1;
                  RecvRank = MyRank_X[2]*NGPU_X[1]*NGPU_X[0] + MyRank_X[0];
                  SendRank = new int [NSend];

                  for (int ID=0; ID<NSend; ID++)   
                     SendRank[ID] = MyRank_X[2]*NGPU_X[1]*NGPU_X[0] + (ID+1)*NGPU_X[0] + MyRank_X[0];

                  break;

      case 6:     NSend    = NGPU_X[2]-1;
                  RecvRank = MyRank_X[1]*NGPU_X[0] + MyRank_X[0];
                  SendRank = new int [NSend];

                  for (int ID=0; ID<NSend; ID++)   
                     SendRank[ID] = (ID+1)*NGPU_X[1]*NGPU_X[0] + MyRank_X[1]*NGPU_X[0] + MyRank_X[0];

                  break;

      default :
                  cerr << "ERROR : incorrect OutputXYZ in the function \"SumOverRanks\" !!" << endl;
                  MPI_Exit();
   }


   RecvBuffer = new real [Size];


   for (int SendID=0; SendID<NSend; SendID++)
   {
      if ( MyRank == RecvRank )
      {
#        ifdef FLOAT8
         MPI_Recv( RecvBuffer, Size, MPI_DOUBLE, SendRank[SendID], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
#        else
         MPI_Recv( RecvBuffer, Size, MPI_FLOAT,  SendRank[SendID], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
#        endif

         for (long t=0; t<Size; t++)  OutputArray[t] += RecvBuffer[t];
      }

      else if ( MyRank == SendRank[SendID] )
      {
#        ifdef FLOAT8
         MPI_Send( OutputArray, Size, MPI_DOUBLE, RecvRank, 0, MPI_COMM_WORLD );
#        else
         MPI_Send( OutputArray, Size, MPI_FLOAT,  RecvRank, 0, MPI_COMM_WORLD );
#        endif
      }

      MPI_Barrier( MPI_COMM_WORLD );
   }


   delete [] SendRank;
   delete [] RecvBuffer;


   if ( MyRank == 0 )   cout << "done" << endl;

} // FUNCTION : SumOverRanks



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :  
//-------------------------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{

   ReadOption( argc, argv );

   Init_MPI( &argc, &argv );

// LoadData( FileName_In, OldDataFormat );
   LoadData( FileName_In );

   if ( MyRank == 0 )   TakeNote();

   AllocateOutputArray();

   Refine2TargetLevel(); 

   if ( OutputXYZ == 4  ||  OutputXYZ == 5  ||  OutputXYZ == 6 )     SumOverRanks(); 

   End_MemFree();

   Output();

   delete [] OutputArray;

   MPI_Finalize();

   if ( MyRank == 0 )   cout << "Program terminated successfully" << endl << endl; 

   return 0;

} // FUNCTION : main


