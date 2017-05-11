#include "GetCube.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData_Old
// Description :  Load the old-format data
//
// Note        :  This function is no longer supported !!!
//
// Parameter   :  FileName : The name of the input file
//-------------------------------------------------------------------------------------------------------
void LoadData_Old( const char *FileName )
{

   if ( MyRank == 0 )   cout << "LoadData_Old \"" << FileName << "\" ... "<< endl;


   FILE *File = fopen( FileName, "rb" );

   if ( File == NULL  &&  MyRank == 0 )
   {
      fprintf( stderr, "Error : the file \"%s\" does not exist !!\n", FileName ); 
      MPI_Exit();
   }


// initialize the NPatchComma list as 0
   for (int lv=0; lv<NLEVEL; lv++)  for (int m=0; m<28; m++)      NPatchComma[lv][m] = 0;


// load the simulation parameters
   int patch_size, refine_ratio, flu_ghost_size, nlevel, block_size_y, max_patch, pot_ghost_size;
   int regrid_count, ngpu, ngpu_x[3];
   int NX0_AllRank[NGPU][3], NTruePatch_AllRank[NGPU][NLEVEL], NDataPatch_AllRank[NGPU];
   real g_newton_g;
   double dt__cfl, omega_m0;
   bool  gravity;

   fread( &gravity,              sizeof(bool),           1,             File );

   fread( &patch_size,           sizeof(int),            1,             File );
   fread( &refine_ratio,         sizeof(int),            1,             File );
   fread( &flu_ghost_size,       sizeof(int),            1,             File );
   fread( &nlevel,               sizeof(int),            1,             File );
   fread( &block_size_y,         sizeof(int),            1,             File );
   fread( &max_patch,            sizeof(int),            1,             File );
   fread( &pot_ghost_size,       sizeof(int),            1,             File );
   
   fread( &GAMMA,                sizeof(real),           1,             File );
   fread( NX0_TOT,               sizeof(int),            3,             File );
   fread( &dt__cfl,              sizeof(double),         1,             File );
   fread( &regrid_count,         sizeof(int),            1,             File );
   fread( &ngpu,                 sizeof(int),            1,             File );
   fread( ngpu_x,                sizeof(int),            3,             File );
   fread( &g_newton_g,           sizeof(real),           1,             File );
   fread( &omega_m0,             sizeof(double),         1,             File );
   fread( &OutputPot,            sizeof(bool),           1,             File );


// verify parameters before loading any size-dependent parameters
   if ( MyRank == 0 )
   {

      if ( patch_size != PATCH_SIZE )
      {
         fprintf( stderr, "Error : restart file PATCH_SIZE = %d, parameter file PATCH_SIZE = %d !!\n", 
                  patch_size, PATCH_SIZE );
         MPI_Exit();
      }

      if ( nlevel != NLEVEL )
      {
         fprintf( stderr, "Error : restart file NLEVEL = %d, parameter file NLEVEL = %d !!\n", nlevel, NLEVEL );
         MPI_Exit();
      }

      if ( ngpu != NGPU )
      {
         fprintf( stderr, "Error : restart file NGPU = %d, parameter file NGPU = %d !!\n", ngpu, NGPU );
         MPI_Exit();
      }

      if ( ngpu_x[0] != NGPU_X[0] )
      {
         fprintf( stderr, "Error : restart file NGPU_X[0] = %d, parameter file NGPU_X[0] = %d !!\n", 
                  ngpu_x[0], NGPU_X[0] );
         MPI_Exit();
      }

      if ( ngpu_x[1] != NGPU_X[1] )
      {
         fprintf( stderr, "Error : restart file NGPU_X[1] = %d, parameter file NGPU_X[1] = %d !!\n", 
                  ngpu_x[1], NGPU_X[1] );
         MPI_Exit();
      }

      if ( ngpu_x[2] != NGPU_X[2] )
      {
         fprintf( stderr, "Error : restart file NGPU_X[2] = %d, parameter file NGPU_X[2] = %d !!\n", 
                  ngpu_x[2], NGPU_X[2] );
         MPI_Exit();
      }
      
      if ( max_patch != MAX_PATCH )
      {
         fprintf( stderr, "Warning : restart file MAX_PATCH = %d, parameter file MAX_PATCH = %d !!\n", 
                  max_patch, MAX_PATCH );
      }

   } // if ( MyRank == 0 )


// load the remaining parameters
   uint evolve_counter[NLEVEL];

   fread( &DumpID,               sizeof(int),            1,             File );
   fread( Time,                  sizeof(double),         NLEVEL,        File );
   fread( &Step,                 sizeof(long int),       1,             File );
   fread( NTruePatch_AllRank,    sizeof(int),            NLEVEL*NGPU,   File );
   fread( NDataPatch_AllRank,    sizeof(int),            NGPU,          File );
   fread( NX0_AllRank,           sizeof(int),            3*NGPU,        File );
   fread( evolve_counter,        sizeof(uint),           NLEVEL,        File );

   fclose( File );


// ******************************************************************************
// set the NX0 for each rank
   for (int dim=0; dim<3; dim++)    NX0[dim] = NX0_AllRank[MyRank][dim];


// initialize parameters related to the targeted domain
   Init_TargetDomain();


// allocate memory for the array "BaseP"
   Init_MemAllocate();


// verify the input parameters
   if ( MyRank == 0 )   CheckParameter();


// set the candidate box
   GetCandidateBox();
// ******************************************************************************


// properly set the file position indicator for each process
   long int HeaderSize, DataSize[NGPU], offset;

   HeaderSize =   sizeof(int     )*( 16 + NLEVEL*NGPU + 4*NGPU )
                + sizeof(long int)*( 1                         )
                + sizeof(uint    )*( NLEVEL                    )
                + sizeof(real    )*( 2                         )
                + sizeof(double  )*( 2 + NLEVEL                )
                + sizeof(bool    )*( 2                         );

   offset = HeaderSize;

   if ( OutputPot )  
   {
      NLoad ++;
      NOut  ++;
   }

   for (int gpu=0; gpu<NGPU; gpu++)
   {
      DataSize[gpu] = 0;

      for (int lv=0; lv<NLEVEL; lv++)     DataSize[gpu] += NTruePatch_AllRank[gpu][lv]*5*sizeof(int);

      DataSize[gpu] += NDataPatch_AllRank[gpu]*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NLoad*sizeof(real);
   }

   for (int gpu=0; gpu<MyRank; gpu++)     offset += DataSize[gpu];


// begin to load data 
   int temp_corner[3], temp_father, temp_son;

// array for re-ordering the fluid data from "xyzv" to "vxyz"
   real (*InvData_Flu)[PATCH_SIZE][PATCH_SIZE][NCOMP] = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE][NCOMP];

   for (int TargetRank=0; TargetRank<NGPU; TargetRank++)
   {
      if ( MyRank == 0 )   
      {
         printf( "  Rank %2d loading ... ", TargetRank );   
         cout << flush << flush << flush;
      }

      if ( MyRank == TargetRank )
      {
         File = fopen( FileName, "rb" );
         fseek( File, offset, SEEK_SET );

         for (int lv=0; lv<NLEVEL; lv++)
         for (int n=0; n<NTruePatch_AllRank[MyRank][lv]; n++)
         {
            fread( temp_corner,  sizeof(int), 3, File );
            fread( &temp_father, sizeof(int), 1, File );
            fread( &temp_son,    sizeof(int), 1, File );

            patch.pnew( lv, temp_corner[0], temp_corner[1], temp_corner[2], temp_father, true );
            patch.ptr[lv][n]->son = temp_son;
         }

         for (int lv=0; lv<NLEVEL; lv++)
         for (int n=0; n<NTruePatch_AllRank[MyRank][lv]; n++)
         {
            if ( patch.ptr[lv][n]->son == -1 )
            {
//             load fluid variables
               fread( InvData_Flu, sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NCOMP, File );

               for (int v=0; v<NCOMP; v++)
               for (int k=0; k<PATCH_SIZE; k++)
               for (int j=0; j<PATCH_SIZE; j++)
               for (int i=0; i<PATCH_SIZE; i++)    
                  patch.ptr[lv][n]->fluid[v][k][j][i] = InvData_Flu[k][j][i][v];

//             load potential
               if ( OutputPot )
               fread( patch.ptr[lv][n]->pot, sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE, File );
            }
         }

         fclose( File );

      } // if ( MyRank == TargetRank )

      if ( MyRank == 0 )   cout << "done" << endl;

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<NGPU; TargetRank++)


   delete [] InvData_Flu;


// record the number of true patches
   for (int lv=0; lv<NLEVEL; lv++)     
   for (int m=1; m<28; m++)         NPatchComma[lv][m] = patch.num[lv];


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
   for (int lv=NLEVEL-2; lv>=0; lv--)     Flu_Restrict( lv );


// fill up the data in the buffer patches
   for (int lv=0; lv<NLEVEL; lv++)        
   {
      Buf_GetBufferData( lv, 1, BUF_SIZE );

      if ( OutputPot )
      Buf_GetBufferData( lv, 2, BUF_SIZE );
   }


   if ( MyRank == 0 )   cout << "LoadData_Old ... done" << endl;

}



