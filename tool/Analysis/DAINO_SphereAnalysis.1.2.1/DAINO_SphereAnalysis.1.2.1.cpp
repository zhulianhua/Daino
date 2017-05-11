#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "TypeDef.h"

using namespace std;

#define WRONG -999999

void ReadOption( int argc, char **argv );
void CheckParameter();
void LoadData();
//void LoadData_Old();
void SetMaxRhoPos( const int AveN );
void End();
void Init_ShellAve();
void Output_ShellAve();
void ShellAverage();
void GetRMS();
void GetMaxRho();
real GetMinShellWidth( const real TCen[], const real TCen_Map[] );
void Load_Parameter_Before_1200( FILE *File, const int FormatVersion, bool &DataOrder_xyzv, 
                                 bool &LoadPot, int *NX0_Tot, double &BoxSize, real &Gamma );
void Load_Parameter_After_1200( FILE *File, const int FormatVersion, bool &DataOrder_xyzv,
                                bool &LoadPot, int *NX0_Tot, double &BoxSize, real &Gamma );
void CompareVar( const char *VarName, const bool   RestartVar, const bool   RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const int    RestartVar, const int    RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const long   RestartVar, const long   RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const real   RestartVar, const real   RuntimeVar, const bool Fatal );
void CompareVar( const char *VarName, const double RestartVar, const double RuntimeVar, const bool Fatal );


// general-purpose variables
DAINO_t     patch;
int         DumpID, NX0_TOT[3];
long int    Step;
double      Time[NLEVEL];
real        GAMMA;      
real        Center_Map[3];             // the coordinates of the "mapped" sphere center for the periodic B.C.

char       *FileName_In    = NULL;
char       *Suffix         = NULL;
real        Center[3]      = { WRONG, WRONG, WRONG };    // the coordinates of the sphere center
real        MaxRadius      = WRONG;    // the maximum radius of the sphere
//bool      OldDataFormat  = false;    // (true/false) : (old,new)
bool        Periodic       = false;    // periodic B.C. (not the case by default)
bool        InputScale     = false;    // true --> input cell scales instead of coordinates
int         GetNShell      = 2;        // set the number of shells = "GetNShell" times the cell size around center
                                       // (<=0 --> disable)


// variables for the mode "shell average"
double      (*Average)[NCOMP] = NULL;  // the shell average in each shell ( rho, v_radial, v_transverse, 
                                       //                                   energy, pressure )
double      (*RMS)    [NCOMP] = NULL;  // the standard deviation (root-mean-squre) in each shell 
real        (*Max)    [NCOMP] = NULL;  // the maximum value in each shell
real        (*Min)    [NCOMP] = NULL;  // the minimum value in each shell
long int    *NCount           = NULL;  // the total number of grids in each shell
double      *Volume           = NULL;  // the total volume in each shell 
real        ShellWidth;                // the width of each shell

bool        Mode_ShellAve     = false; // true --> evaluate the shell average
int         NShell            = WRONG; // number of shells
int         UseMaxRhoPos      = 8;     // number of highest-density cells for determining the sphere center
                                       // (<= 0 --> disable)

// variables for the mode "maximum density"
bool        Mode_MaxRho       = false; // true --> get the cells exceeding the density threshold
real        RhoThres          = WRONG; // the density threshold (for the MaxRho mode)




//-------------------------------------------------------------------------------------------------------
// Function    :  GetMaxRho
// Description :  Output the mesh information if the density exceeds a given threshold "RhoThres" 
//-------------------------------------------------------------------------------------------------------
void GetMaxRho()
{

   cout << "GetMaxRho ..." << endl;


   const real dh_min = patch.dh[NLEVEL-1];
   real Radius, x, x1, x2, y, y1, y2, z, z1, z2, scale; // (x,y,z) : relative coordinates to the vector "Center"
   real xx, yy, zz;

#  if   ( MODEL == HYDRO )
   real Dens, VelX, VelY, VelZ, Engy, Pres;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   real Dens, Real, Imag;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// open the file
   char FileName[50] = "MaxRho";

   if ( Suffix != NULL )   strcat( FileName, Suffix );

   if ( NULL != fopen(FileName,"r") )  
         fprintf( stderr, "Warning : the file \"%s\" already exists and will be overwritten !!\n", FileName );

   FILE *File = fopen( FileName, "w" );

#  if   ( MODEL == HYDRO )
   fprintf( File, "%10s %10s %10s %2s %12s %12s %12s %12s %12s %12s\n", 
            "x", "y", "z", "Lv", "Dens", "VelX", "VelY", "VelZ", "Engy", "Pres" );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   fprintf( File, "%10s %10s %10s %2s %12s %12s %12s\n", 
            "x", "y", "z", "Lv", "Dens", "Real", "Imag" );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


   for (int lv=0; lv<NLEVEL; lv++)              
   { 
      scale = (real)patch.scale[lv];
      cout << "   Level " << lv << " ... ";

      for (int PID=0; PID<patch.num[lv]; PID++) 
      {
         for (int k=0; k<PATCH_SIZE; k++) {  zz = patch.ptr[lv][PID]->corner[2] + (k+0.5)*scale; 
                                             z1 = zz - Center    [2];
                                             z2 = zz - Center_Map[2];
                                             z  = ( fabs(z1) <= fabs(z2) ) ?  z1 : z2;
         for (int j=0; j<PATCH_SIZE; j++) {  yy = patch.ptr[lv][PID]->corner[1] + (j+0.5)*scale;  
                                             y1 = yy - Center    [1];
                                             y2 = yy - Center_Map[1];
                                             y  = ( fabs(y1) <= fabs(y2) ) ?  y1 : y2;
         for (int i=0; i<PATCH_SIZE; i++) {  xx = patch.ptr[lv][PID]->corner[0] + (i+0.5)*scale; 
                                             x1 = xx - Center    [0];
                                             x2 = xx - Center_Map[0];
                                             x  = ( fabs(x1) <= fabs(x2) ) ?  x1 : x2;

            Radius = sqrt( x*x + y*y + z*z );

//          output the mesh information if it's within the targeted sphere and the density exceeds "RhoThres"
            if ( Radius < MaxRadius  &&  patch.ptr[lv][PID]->fluid[DENS][k][j][i] >= RhoThres )
            {
#              if   ( MODEL == HYDRO )
               Dens = patch.ptr[lv][PID]->fluid[DENS][k][j][i];
               VelX = patch.ptr[lv][PID]->fluid[MOMX][k][j][i] / Dens;
               VelY = patch.ptr[lv][PID]->fluid[MOMY][k][j][i] / Dens;
               VelZ = patch.ptr[lv][PID]->fluid[MOMZ][k][j][i] / Dens;
               Engy = patch.ptr[lv][PID]->fluid[ENGY][k][j][i];
               Pres = (GAMMA-1.0) * ( Engy - 0.5*Dens*(VelX*VelX + VelY*VelY + VelZ*VelZ) ); 

               fprintf( File, "%10.4e %10.4e %10.4e %2d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                        xx*dh_min, yy*dh_min, zz*dh_min, lv, Dens, VelX, VelY, VelZ, Engy, Pres );

#              elif ( MODEL == MHD )
#              warning : WAIT MHD !!!

#              elif ( MODEL == ELBDM )
               Dens = patch.ptr[lv][PID]->fluid[DENS][k][j][i];
               Real = patch.ptr[lv][PID]->fluid[REAL][k][j][i];
               Imag = patch.ptr[lv][PID]->fluid[IMAG][k][j][i];

               fprintf( File, "%10.4e %10.4e %10.4e %2d %12.5e %12.5e %12.5e\n",
                        xx*dh_min, yy*dh_min, zz*dh_min, lv, Dens, Real, Imag );

#              else
#              error : ERROR : unsupported MODEL !!
#              endif // MODEL

            } // if ( Radius < MaxRadius )
         }}} // k, j, i
      } // for (int PID=0; PID<patch.num[lv]; PID++)

      cout << "done" << endl;

   } // for (int lv=0; lv<NLEVEL; lv++)


   fclose( File );

} // FUNCTION : GetMaxRho



//-------------------------------------------------------------------------------------------------------
// Function    :  GetRMS
// Description :  Evaluate the standard deviations from the average values 
//-------------------------------------------------------------------------------------------------------
void GetRMS()
{

   cout << "Evaluating the RMS ..." << endl;


   int  ShellID;
   real Radius, scale, dv;
   real x, x1, x2, y, y1, y2, z, z1, z2;     // (x,y,z) : relative coordinates to the vector "Center"

#  if   ( MODEL == HYDRO )
   real rho, vx, vy, vz, vr, vt, egy, pres;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   real Dens, Real, Imag;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


   for (int lv=0; lv<NLEVEL; lv++)              
   {  
      scale = (real)patch.scale[lv];
      dv    = scale*scale*scale;

      cout << "   Level " << lv << " ... ";

      for (int PID=0; PID<patch.num[lv]; PID++) 
      {
         for (int k=0; k<PATCH_SIZE; k++) {  z1 = patch.ptr[lv][PID]->corner[2] + (k+0.5)*scale - Center    [2];
                                             z2 = patch.ptr[lv][PID]->corner[2] + (k+0.5)*scale - Center_Map[2];
                                             z  = ( fabs(z1) <= fabs(z2) ) ? z1 : z2;
         for (int j=0; j<PATCH_SIZE; j++) {  y1 = patch.ptr[lv][PID]->corner[1] + (j+0.5)*scale - Center    [1];
                                             y2 = patch.ptr[lv][PID]->corner[1] + (j+0.5)*scale - Center_Map[1];
                                             y  = ( fabs(y1) <= fabs(y2) ) ? y1 : y2;
         for (int i=0; i<PATCH_SIZE; i++) {  x1 = patch.ptr[lv][PID]->corner[0] + (i+0.5)*scale - Center    [0];
                                             x2 = patch.ptr[lv][PID]->corner[0] + (i+0.5)*scale - Center_Map[0];
                                             x  = ( fabs(x1) <= fabs(x2) ) ? x1 : x2;

            Radius = sqrt( x*x + y*y + z*z );

            if ( Radius < MaxRadius )
            {
               ShellID = int( Radius / ShellWidth );

               if ( ShellID >= NShell )
               {
                  cerr << "ERROR : ShellID >= NShell !!" << endl;
                  exit( 1 );
               }

#              if   ( MODEL == HYDRO )
//             evaluate the values on the shell
               rho  = patch.ptr[lv][PID]->fluid[DENS][k][j][i];
               vx   = patch.ptr[lv][PID]->fluid[MOMX][k][j][i] / rho;
               vy   = patch.ptr[lv][PID]->fluid[MOMY][k][j][i] / rho;
               vz   = patch.ptr[lv][PID]->fluid[MOMZ][k][j][i] / rho;
               egy  = patch.ptr[lv][PID]->fluid[ENGY][k][j][i];

               pres = (GAMMA-1.0) * ( egy - 0.5*rho*(vx*vx + vy*vy + vz*vz) ); 
               vr   = ( x*vx + y*vy + z*vz ) / Radius;
               vt   = sqrt( fabs(vx*vx + vy*vy + vz*vz - vr*vr) );


//             evalute the square of deviation on the shell
               RMS[ShellID][0] += dv*pow( double(rho )-Average[ShellID][0], 2.0 );
               RMS[ShellID][1] += dv*pow( double(vr  )-Average[ShellID][1], 2.0 );
               RMS[ShellID][2] += dv*pow( double(vt  )-Average[ShellID][2], 2.0 );
               RMS[ShellID][3] += dv*pow( double(egy )-Average[ShellID][3], 2.0 );
               RMS[ShellID][4] += dv*pow( double(pres)-Average[ShellID][4], 2.0 );

#              elif ( MODEL == MHD )
#              warning : WAIT MHD !!!

#              elif ( MODEL == ELBDM )
               Dens = patch.ptr[lv][PID]->fluid[DENS][k][j][i];
               Real = patch.ptr[lv][PID]->fluid[REAL][k][j][i];
               Imag = patch.ptr[lv][PID]->fluid[IMAG][k][j][i];

//             evalute the square of deviation on the shell
               RMS[ShellID][0] += dv*pow( double(Dens)-Average[ShellID][0], 2.0 );
               RMS[ShellID][1] += dv*pow( double(Real)-Average[ShellID][1], 2.0 );
               RMS[ShellID][2] += dv*pow( double(Imag)-Average[ShellID][2], 2.0 );

#              else
#              error : ERROR : unsupported MODEL !!
#              endif // MODEL

            } // if ( Radius < MaxRadius )
         }}} // k, j, i
      } // for (int PID=0; PID<patch.num[lv]; PID++)

      cout << "done" << endl;

   } // for (int lv=0; lv<NLEVEL; lv++)


// get the root-mean-square at each level
   for (int n=0; n<NShell; n++)
   for (int v=0; v<NCOMP; v++)      RMS[n][v] = sqrt( RMS[n][v]/Volume[n] );

} // FUNCTION : GetRMS



//-------------------------------------------------------------------------------------------------------
// Function    :  ShellAverage
// Description :  Get the shell average of all variables 
//-------------------------------------------------------------------------------------------------------
void ShellAverage()
{

   cout << "Evaluating the shell average ..." << endl;


   int  ShellID;
   real Radius, scale, dv;
   real x, x1, x2, y, y1, y2, z, z1, z2;     // (x,y,z) : relative coordinates to the vector "Center"

#  if   ( MODEL == HYDRO )
   real px, py, pz, pr, pt, vr, vt, pres, rho, egy;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   real Dens, Real, Imag;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


   for (int lv=0; lv<NLEVEL; lv++)              
   { 
      scale = (real)patch.scale[lv];
      dv    = scale*scale*scale;

      cout << "   Level " << lv << " ... ";

      for (int PID=0; PID<patch.num[lv]; PID++) 
      {
         for (int k=0; k<PATCH_SIZE; k++) {  z1 = patch.ptr[lv][PID]->corner[2] + (k+0.5)*scale - Center    [2];
                                             z2 = patch.ptr[lv][PID]->corner[2] + (k+0.5)*scale - Center_Map[2];
                                             z  = ( fabs(z1) <= fabs(z2) ) ? z1 : z2;
         for (int j=0; j<PATCH_SIZE; j++) {  y1 = patch.ptr[lv][PID]->corner[1] + (j+0.5)*scale - Center    [1];
                                             y2 = patch.ptr[lv][PID]->corner[1] + (j+0.5)*scale - Center_Map[1];
                                             y  = ( fabs(y1) <= fabs(y2) ) ? y1 : y2;
         for (int i=0; i<PATCH_SIZE; i++) {  x1 = patch.ptr[lv][PID]->corner[0] + (i+0.5)*scale - Center    [0];
                                             x2 = patch.ptr[lv][PID]->corner[0] + (i+0.5)*scale - Center_Map[0];
                                             x  = ( fabs(x1) <= fabs(x2) ) ? x1 : x2;

            Radius = sqrt( x*x + y*y + z*z );

            if ( Radius < MaxRadius )
            {
               ShellID = int( Radius / ShellWidth );

               if ( ShellID >= NShell )
               {
                  cerr << "ERROR : ShellID >= NShell !!" << endl;
                  exit( 1 );
               }

#              if   ( MODEL == HYDRO )
//             evaluate the values on the shell
               rho  = patch.ptr[lv][PID]->fluid[DENS][k][j][i];
               px   = patch.ptr[lv][PID]->fluid[MOMX][k][j][i];
               py   = patch.ptr[lv][PID]->fluid[MOMY][k][j][i];
               pz   = patch.ptr[lv][PID]->fluid[MOMZ][k][j][i];
               egy  = patch.ptr[lv][PID]->fluid[ENGY][k][j][i];

               pr   = ( x*px + y*py + z*pz ) / Radius;
               pt   = sqrt( fabs(px*px + py*py + pz*pz - pr*pr) );
               pres = (GAMMA-1.0) * ( egy - 0.5*(px*px + py*py + pz*pz)/rho ); 
               vr   = pr/rho;
               vt   = pt/rho;


//             add up the values
               Average[ShellID][0] += (double)(dv*rho );
               Average[ShellID][1] += (double)(dv*pr  );
               Average[ShellID][2] += (double)(dv*pt  );
               Average[ShellID][3] += (double)(dv*egy );
               Average[ShellID][4] += (double)(dv*pres);


//             store the maximum and minimum values
               if ( rho  > Max[ShellID][0] )    Max[ShellID][0] = rho;
               if ( vr   > Max[ShellID][1] )    Max[ShellID][1] = vr;
               if ( vt   > Max[ShellID][2] )    Max[ShellID][2] = vt;
               if ( egy  > Max[ShellID][3] )    Max[ShellID][3] = egy;
               if ( pres > Max[ShellID][4] )    Max[ShellID][4] = pres;

               if ( rho  < Min[ShellID][0] )    Min[ShellID][0] = rho;
               if ( vr   < Min[ShellID][1] )    Min[ShellID][1] = vr;
               if ( vt   < Min[ShellID][2] )    Min[ShellID][2] = vt;
               if ( egy  < Min[ShellID][3] )    Min[ShellID][3] = egy;
               if ( pres < Min[ShellID][4] )    Min[ShellID][4] = pres;

#              elif ( MODEL == MHD )
#              warning : WAIT MHD !!!

#              elif ( MODEL == ELBDM )
               Dens = patch.ptr[lv][PID]->fluid[DENS][k][j][i];
               Real = patch.ptr[lv][PID]->fluid[REAL][k][j][i];
               Imag = patch.ptr[lv][PID]->fluid[IMAG][k][j][i];

//             add up the values
               Average[ShellID][0] += (double)(dv*Dens);
               Average[ShellID][1] += (double)(dv*Real);
               Average[ShellID][2] += (double)(dv*Imag);


//             store the maximum and minimum values
               if ( Dens > Max[ShellID][0] )    Max[ShellID][0] = Dens;
               if ( Real > Max[ShellID][1] )    Max[ShellID][1] = Real;
               if ( Imag > Max[ShellID][2] )    Max[ShellID][2] = Imag;

               if ( Dens < Min[ShellID][0] )    Min[ShellID][0] = Dens;
               if ( Real < Min[ShellID][1] )    Min[ShellID][1] = Real;
               if ( Imag < Min[ShellID][2] )    Min[ShellID][2] = Imag;

#              else
#              error : ERROR : unsupported MODEL !!
#              endif // MODEL

               Volume[ShellID] += dv;
               NCount[ShellID] ++;

            } // if ( Radius < MaxRadius )
         }}} // k, j, i
      } // for (int PID=0; PID<patch.num[lv]; PID++)

      cout << "done" << endl;

   } // for (int lv=0; lv<NLEVEL; lv++)


// get the average values
   for (int n=0; n<NShell; n++)
   for (int v=0; v<NCOMP; v++)      Average[n][v] /= Volume[n];


// get the average velocity ( define as Ave(momenum)/Ave(density) )
#  if   ( MODEL == HYDRO )
   for (int n=0; n<NShell; n++)
   {
      Average[n][1] /= Average[n][0];
      Average[n][2] /= Average[n][0];
   }

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif

} // FUNCTION : ShellAverage



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData
// Description :  Load data from the input file 
//-------------------------------------------------------------------------------------------------------
void LoadData()
{

// load data using the old data format
   /*
   if ( OldDataFormat )
   {
      LoadData_Old();

      return;
   }
   */


   cout << "Loading data \"" << FileName_In << "\" ..." << endl;


   FILE *File = fopen( FileName_In, "rb" );

   if ( File == NULL )
   {
      fprintf( stderr, "\nERROR : the file \"%s\" does not exist !!\n", FileName_In );
      exit( 1 );
   }


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
      exit( 1 );
   }

   if ( size_long != size_ulong )
   {
      fprintf( stderr, "ERROR : sizeof(long) = %d != sizeof(ulong) = %d !!\n", size_long, size_ulong );
      exit( 1 );
   }


// a. load the information of data format
// =================================================================================================
   long FormatVersion, HeaderSize, CheckCode;

   fread( &FormatVersion, sizeof(long), 1, File );
   fread( &HeaderSize,    sizeof(long), 1, File );
   fread( &CheckCode,     sizeof(long), 1, File );


// verify the input data format version
   fprintf( stdout, "   The version of the RESTART file's format = %ld\n", FormatVersion );

   if ( FormatVersion < 1100 )
   {
      fprintf( stderr, "ERROR : unsupported data format version (only support version >= 1100) !!\n" );
      exit( 1 );
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
         exit( 1 );
      }

      if ( size_int_restart != size_int )  
      {
         fprintf( stderr, "ERROR : sizeof(int) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                  size_int_restart, size_int );
         exit( 1 );
      }

      if ( size_long_restart != size_long )  
      {
         fprintf( stderr, "ERROR : sizeof(long) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                  size_long_restart, size_long );
         exit( 1 );
      }

      if ( size_real_restart != size_real )  
      {
         fprintf( stderr, "ERROR : sizeof(real) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                  size_real_restart, size_real );
         exit( 1 );
      }

      if ( size_double_restart != size_double )  
      {
         fprintf( stderr, "ERROR : sizeof(double) is inconsistent : RESTART file = %d, runtime = %d !!\n",
                  size_double_restart, size_double );
         exit( 1 );
      }
   } // if ( FormatVersion >= 1200 )


// skip the buffer space
   const int NBuf_Format_1200 = 256 - 0*size_bool -  5*size_int -  3*size_long -  0*size_real -  0*size_double;
   const int NBuf_Format      = ( FormatVersion >= 1200 ) ? NBuf_Format_1200 : 40; 

   fseek( File, NBuf_Format, SEEK_CUR );



// b. load all simulation parameters
// =================================================================================================
   bool   DataOrder_xyzv, LoadPot;
   double BoxSize;

   if ( FormatVersion < 1200 )   
      Load_Parameter_Before_1200( File, FormatVersion, DataOrder_xyzv, LoadPot, NX0_TOT, BoxSize, GAMMA );
   else
      Load_Parameter_After_1200 ( File, FormatVersion, DataOrder_xyzv, LoadPot, NX0_TOT, BoxSize, GAMMA );


// set the file position indicator to "HeaderSize" and verify the check code
   long checkcode;
   fseek( File, HeaderSize, SEEK_SET );
   fread( &checkcode, sizeof(long), 1, File );

   if ( checkcode != CheckCode )
   {
      fprintf( stderr, "ERROR : incorrect check code in the RESTART file (input %ld <-> expeect %ld) !!\n",
               checkcode, CheckCode );
      exit( 1 );
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


// verify the size of the input file
   long InfoSize, DataSize[NLEVEL], ExpectSize, InputSize, PatchDataSize;
   int  NVar;  // number of variables ( NCOMP or NCOMP+1 -> potential )

   InfoSize =     sizeof(int   )*( 1 + 2*NLEVEL )
                + sizeof(long  )*( 2            )  // Step + checkcode
                + sizeof(uint  )*(       NLEVEL )
                + sizeof(double)*( 1 +   NLEVEL )
                + NBuf_Info;

   NVar = ( LoadPot ) ? NCOMP+1 : NCOMP;

   PatchDataSize = PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NVar*sizeof(real);
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
               FileName_In, InputSize, ExpectSize );
      exit( 1 );
   }


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



// d. set up and check parameters dedicated in this tool
// =================================================================================================
   static bool FirstTime = true;

   if ( FirstTime )
   {
//    convert physical coordinates to cell scales
      if ( !InputScale )
      {
         const real _dh_min = 1.0 / patch.dh[NLEVEL-1];

         for (int d=0; d<3; d++)    
         {
            if ( Center[d] != WRONG )  Center[d] *= _dh_min;
         }

         if ( MaxRadius != WRONG )     MaxRadius *= _dh_min;
      }


//    initialize the sphere information if they are not provided by users
      if ( NShell == WRONG )   
      {
         NShell = ( NX0_TOT[0] < NX0_TOT[1] ) ? NX0_TOT[0] : NX0_TOT[1];
         NShell = ( NX0_TOT[2] < NShell     ) ? NX0_TOT[2] : NShell;
         NShell = NShell*patch.scale[0]/2;
      }

      if ( Center[0] == WRONG )  Center[0] = 0.5*patch.BoxScale[0]; 
      if ( Center[1] == WRONG )  Center[1] = 0.5*patch.BoxScale[1]; 
      if ( Center[2] == WRONG )  Center[2] = 0.5*patch.BoxScale[2]; 

      if ( MaxRadius == WRONG )  
      {
         MaxRadius = fmin( 0.5*patch.BoxScale[0], 0.5*patch.BoxScale[1] );
         MaxRadius = fmin( 0.5*patch.BoxScale[2], MaxRadius );
      }

      FirstTime = false;
   } // if ( FirstTime )


// set up the "mapped" sphere center for the periodic B.C.
   const real HalfBox[3] = { 0.5*patch.BoxScale[0], 0.5*patch.BoxScale[1], 0.5*patch.BoxScale[2] };

   if ( Periodic )
      for (int d=0; d<3; d++)
         Center_Map[d] = ( Center[d] > HalfBox[d] ) ? Center[d]-2.0*HalfBox[d] : Center[d]+2.0*HalfBox[d];

   else
      for (int d=0; d<3; d++)    
         Center_Map[d] = Center[d];


// set up the candidate box --> only patches with cells inside the candidate box will be allocated
   const real CanMax_x1 = Center    [0] + MaxRadius;
   const real CanMin_x1 = Center    [0] - MaxRadius;
   const real CanMax_y1 = Center    [1] + MaxRadius;
   const real CanMin_y1 = Center    [1] - MaxRadius;
   const real CanMax_z1 = Center    [2] + MaxRadius;
   const real CanMin_z1 = Center    [2] - MaxRadius;

   const real CanMax_x2 = Center_Map[0] + MaxRadius;
   const real CanMin_x2 = Center_Map[0] - MaxRadius;
   const real CanMax_y2 = Center_Map[1] + MaxRadius;
   const real CanMin_y2 = Center_Map[1] - MaxRadius;
   const real CanMax_z2 = Center_Map[2] + MaxRadius;
   const real CanMin_z2 = Center_Map[2] - MaxRadius;



// e. load the simulation data
// =================================================================================================
   int  LoadCorner[3], LoadSon, PID, cr1[3], cr2[3];
   bool GotYou;

// array for re-ordering the fluid data from "xyzv" to "vxyz"
   real (*InvData_Flu)[PATCH_SIZE][PATCH_SIZE][NCOMP] = NULL;
   if ( DataOrder_xyzv )   InvData_Flu = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE][NCOMP];

   fseek( File, HeaderSize+InfoSize, SEEK_SET );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      cout << "   Loading level " << lv << " ... ";

      for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)
      {

//       e1. load the patch information
         fread(  LoadCorner, sizeof(int), 3, File );
         fread( &LoadSon,    sizeof(int), 1, File );


//       e2. create the patch and load the physical data if it is a leaf patch
         if ( LoadSon == -1 )
         {
            GotYou = false;

            for (int dim=0; dim<3; dim++)
            {
               cr1[dim] = LoadCorner[dim];
               cr2[dim] = LoadCorner[dim] + PATCH_SIZE*patch.scale[lv];
            }

            if ( !Periodic )
            {
               if (  ( cr1[0]>=CanMin_x1 && cr1[0]<=CanMax_x1 ) || ( cr2[0]>=CanMin_x1 && cr2[0]<=CanMax_x1 )  )
               if (  ( cr1[1]>=CanMin_y1 && cr1[1]<=CanMax_y1 ) || ( cr2[1]>=CanMin_y1 && cr2[1]<=CanMax_y1 )  )
               if (  ( cr1[2]>=CanMin_z1 && cr1[2]<=CanMax_z1 ) || ( cr2[2]>=CanMin_z1 && cr2[2]<=CanMax_z1 )  )
                 GotYou = true;
            }

            else
            {
               if (  ( cr1[0]>=CanMin_x1 && cr1[0]<=CanMax_x1 ) || ( cr2[0]>=CanMin_x1 && cr2[0]<=CanMax_x1 ) ||
                     ( cr1[0]>=CanMin_x2 && cr1[0]<=CanMax_x2 ) || ( cr2[0]>=CanMin_x2 && cr2[0]<=CanMax_x2 )    )
               if (  ( cr1[1]>=CanMin_y1 && cr1[1]<=CanMax_y1 ) || ( cr2[1]>=CanMin_y1 && cr2[1]<=CanMax_y1 ) ||
                     ( cr1[1]>=CanMin_y2 && cr1[1]<=CanMax_y2 ) || ( cr2[1]>=CanMin_y2 && cr2[1]<=CanMax_y2 )    )
               if (  ( cr1[2]>=CanMin_z1 && cr1[2]<=CanMax_z1 ) || ( cr2[2]>=CanMin_z1 && cr2[2]<=CanMax_z1 ) ||
                     ( cr1[2]>=CanMin_z2 && cr1[2]<=CanMax_z2 ) || ( cr2[2]>=CanMin_z2 && cr2[2]<=CanMax_z2 )    )
                 GotYou = true;
            }


            if ( GotYou )
            {
               PID = patch.num[lv];

               patch.pnew( lv, LoadCorner[0], LoadCorner[1], LoadCorner[2] );

//             e2-1. load the fluid variables
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
                  fread( patch.ptr[lv][PID]->fluid, sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NCOMP, File );

//             e2-2. load the gravitational potential
               if ( LoadPot )
                  fread( patch.ptr[lv][PID]->pot,   sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE,       File );
            }

            else
            {
               fseek( File, PatchDataSize, SEEK_CUR );
            }

         } // if ( LoadSon == -1 )
      } // for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++) 

      cout << "done" << endl;

   } // for (int lv=0; lv<NLEVEL; lv++)


   if ( InvData_Flu != NULL )    delete [] InvData_Flu;

   fclose( File );



// f. set the shell width
// =================================================================================================
   if ( GetNShell > 0 )
   {
      ShellWidth = GetNShell * GetMinShellWidth( Center, Center_Map );
      NShell     = (int)ceil( MaxRadius / ShellWidth );
   }

   else
      ShellWidth = MaxRadius / (real)NShell;


   cout << "Loading data \"" << FileName_In << "\" ... done" << endl;

} // FUNCTION : LoadData



//-------------------------------------------------------------------------------------------------------
// Function    :  ReadOption
// Description :  Read the command line options 
//-------------------------------------------------------------------------------------------------------
void ReadOption( int argc, char **argv )
{

   int c;

   while ( (c = getopt(argc, argv, "hpsSMi:o:n:x:y:z:r:t:m:a:")) != -1 )
   {
      switch ( c )
      {
         case 'i': FileName_In   = optarg;
                   break;
         case 'o': Suffix        = optarg;
                   break;
         case 'n': NShell        = atoi(optarg);
                   break;
         case 'x': Center[0]     = atof(optarg);
                   break;
         case 'y': Center[1]     = atof(optarg);
                   break;
         case 'z': Center[2]     = atof(optarg);
                   break;
         case 'r': MaxRadius     = atof(optarg);
                   break;
         case 'p': Periodic      = true;
                   break;
         case 'm': UseMaxRhoPos  = atoi(optarg);
                   break;
         /*
         case 'O': OldDataFormat = true;
                   break;
         */
         case 'S': Mode_ShellAve = true;
                   break;
         case 'M': Mode_MaxRho   = true;
                   break;
         case 't': RhoThres      = atof(optarg);
                   break;
         case 's': InputScale    = true; 
                   break;
         case 'a': GetNShell     = atoi(optarg); 
                   break;
         case 'h': 
         case '?': cerr << endl << "usage: " << argv[0] 
                        << " [-h (for help)] [-i input fileName] [-o suffix to the output file [none]]" 
                        << endl << "                             " 
                        << " [-n # of shells [0.5*BoxSize]] [-(x,y,z) sphere center [box center]]" 
                        << endl << "                             " 
                        << " [-r sphere radius [0.5*BoxSize]] [-p periodic B.C. [off]]"
                        << endl << "                             " 
                        << " [-m # of highest-density cells for determining the sphere center (<=0: off) [8]"
                        << endl << "                             " 
                        << " [-a NCell (set shell width equal to NCell*cell_size around center, <=0: off) [2]]"
                        << endl << "                             " 
                        << " [-s [input cell scales instead of physical coordinates to specify the range] [off]]"
                        << endl << "                             " 
                        << " [-t density threshold [off]]" 
                        << endl << "                             " 
                        << " [-S (turn on the mode \"shell average\") [off]]"  
                        << endl << "                             " 
                        << " [-M (turn on the mode \"maximum density\") [off]]"
                        << endl << endl;
                   exit( 1 );

      } // switch ( c ) 
   } // while ..


// initial check
   if ( fopen( FileName_In, "rb" ) == NULL )
   {
      fprintf( stderr, "ERROR : the input file \"%s\" does not exist (-i input filename) !!\n", FileName_In );
      exit( 1 );
   }

   if ( !Mode_ShellAve  &&  !Mode_MaxRho )
   {
      fprintf( stderr, "ERROR : no operation mode is turned on (-S, -M) !!\n" );
      exit( 1 );
   }

   if ( Mode_MaxRho  &&  RhoThres == WRONG )
   {
      fprintf( stderr, "ERROR : please provide the density threshold (-t density threshold) !!\n" );
      exit( 1 );
   }

} // FUNCTION : ReadOption



//-------------------------------------------------------------------------------------------------------
// Function    :  GetMinShellWidth 
// Description :  Return the size (scale) of the cell closest to the input sphere center
//
// Parameter   :  TCen     : Center of the targeted sphere
//                TCen_Map : Center of the mapped targeted sphere (useful only if Periodic is enabled)
//-------------------------------------------------------------------------------------------------------
real GetMinShellWidth( const real TCen[], const real TCen_Map[] )
{

   cout << "   GetMinShellWidth ..." << endl;


   real x, x1, x2, xx, y, y1, y2, yy, z, z1, z2, zz; 
   real r, r_min, scale, ShellWidth_min;
   int  Lv_min;

   r_min  = __FLT_MAX__;
   Lv_min = 0;


   for (int lv=0; lv<NLEVEL; lv++)
   {
      scale = patch.scale[lv];

      for (int PID=0; PID<patch.num[lv]; PID++)
      {
         for (int k=0; k<PATCH_SIZE; k++) {  zz = patch.ptr[lv][PID]->corner[2] + (k+0.5)*scale; 
                                             z1 = zz - TCen    [2];
                                             z2 = zz - TCen_Map[2];
                                             z  = ( fabs(z1) <= fabs(z2) ) ?  z1 : z2;
         for (int j=0; j<PATCH_SIZE; j++) {  yy = patch.ptr[lv][PID]->corner[1] + (j+0.5)*scale;  
                                             y1 = yy - TCen    [1];
                                             y2 = yy - TCen_Map[1];
                                             y  = ( fabs(y1) <= fabs(y2) ) ?  y1 : y2;
         for (int i=0; i<PATCH_SIZE; i++) {  xx = patch.ptr[lv][PID]->corner[0] + (i+0.5)*scale; 
                                             x1 = xx - TCen    [0];
                                             x2 = xx - TCen_Map[0];
                                             x  = ( fabs(x1) <= fabs(x2) ) ?  x1 : x2;

            r = sqrt( x*x + y*y + z*z );

            if ( r < r_min )
            {
               r_min  = r;
               Lv_min = lv;
            }

         }}} // i,j,k
      } // for (int PID=0; PID<patch.num[lv]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)


   ShellWidth_min = patch.scale[Lv_min];

   printf( "      Minimum distance = %13.7e\n", r_min );
   printf( "      Level            = %d\n",     Lv_min );
   printf( "      Shell scale      = %13.7e\n", ShellWidth_min ); 


   cout << "   GetMinShellWidth ... done" << endl;

   return ShellWidth_min;

} // FUNCTION : GetMinShellWidth



//-------------------------------------------------------------------------------------------------------
// Function    :  Output_ShellAve 
// Description :  Output data for the mode "shell average"
//-------------------------------------------------------------------------------------------------------
void Output_ShellAve()
{

   cout << "Output_ShellAve ... " << flush;


// set output file names
#  if   ( MODEL == HYDRO )
   char FileName[NCOMP][50] = { "AveRho", "AveV_R", "AveV_T", "AveEgy", "AvePre" };

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   char FileName[NCOMP][50] = { "AveDens", "AveReal", "AveImag" };

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   if ( Suffix != NULL )       
      for (int v=0; v<NCOMP; v++)   strcat( FileName[v], Suffix );


// output data
   const real dh_min   = patch.dh[NLEVEL-1];
   const real VolConst = 4.0*M_PI/3.0;
   double AccMass      = 0.0;       // accumulated mass
   bool   OutputMass   = false;
   double v1, v2;

   for (int v=0; v<NCOMP; v++)
   {
//    determine whether or not to output the accumulated mass
#     ifdef DENS
      if ( v == DENS )  OutputMass = true;
      else              OutputMass = false;
#     endif 


//    header
      if ( NULL != fopen(FileName[v],"r") )  
         fprintf( stderr, "Warning : the file \"%s\" already exists and will be overwritten !!\n", FileName[v] );

      FILE *File = fopen( FileName[v], "w" );

      fprintf( File, "%12s  %10s  %13s  %13s  %13s  %13s", "Radius", "NCount", "Ave", "RMS", "Max", "Min" );
      if ( OutputMass )    fprintf( File, "  %13s  %13s", "Mass", "AveDens" );
      fprintf( File, "\n" );


//    data
      for (int n=0; n<NShell; n++)
      {
         fprintf( File, "%12.5e  %10ld  %13.6e  %13.6e  %13.6e  %13.6e", 
                  (n+0.5)*ShellWidth*dh_min, NCount[n], Average[n][v], RMS[n][v], Max[n][v], Min[n][v] );

         if ( OutputMass )
         {
            v1 = VolConst * pow( (double)(n+1)*ShellWidth*dh_min, 3.0 );
            v2 = VolConst * pow( (double)(n  )*ShellWidth*dh_min, 3.0 );

            if ( NCount[n] != 0 )   AccMass += Average[n][DENS]*(v1 - v2);

            fprintf( File, "  %13.6e  %13.6e", AccMass, AccMass/v1 );
         }

         fprintf( File, "\n" );
      }

      fclose( File );

   } // for (int v=0; v<NCOMP; v++)


   cout << "done" << endl;

} // FUNCTION : Output_ShellAve



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckParameter
// Description :  Verify parameters
//-------------------------------------------------------------------------------------------------------
void CheckParameter()
{

   cout << "CheckParameter ... " << flush;


#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#     error : ERROR : unsupported MODEL !!
#  endif

#  ifndef DENS
#     error : ERROR : symbolic constant DENS is not defined !!
#  endif

   if (  Center[0] > patch.BoxScale[0]  ||  Center[0] < 0.0  )
   {
      fprintf( stderr, "ERROR : the x coordinate of the sphere center lies outside the simulation box !!\n" );
      exit( 1 );
   }

   if (  Center[1] > patch.BoxScale[1]  ||  Center[1] < 0.0  )
   {
      fprintf( stderr, "ERROR : the y coordinate of the sphere center lies outside the simulation box !!\n" );
      exit( 1 );
   }

   if (  Center[2] > patch.BoxScale[2]  ||  Center[2] < 0.0  )
   {
      fprintf( stderr, "ERROR : the z coordinate of the sphere center lies outside the simulation box !!\n" );
      exit( 1 );
   }

   if ( MaxRadius > 0.5*patch.BoxScale[0] ) 
      fprintf( stderr, "Warning : the sphere radius exceeds the half box size in the x direction !!\n" );

   if ( MaxRadius > 0.5*patch.BoxScale[1] )
      fprintf( stderr, "Warning : the sphere radius exceeds the half box size in the y direction !!\n" );

   if ( MaxRadius > 0.5*patch.BoxScale[2] )
      fprintf( stderr, "Warning : the sphere radius exceeds the half box size in the z direction !!\n" );


// checks for the mode "shell average"
   if ( Mode_ShellAve )
   {
      if ( NShell <= 0 )
      {
         fprintf( stderr, "ERROR : please provide the correct number of shells (-n # of shells) !!\n" );
         exit( 1 );
      }
   }


   cout << "done" << endl;

} // FUNCTION : CheckParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  End
// Description :  Deallocate memory 
//-------------------------------------------------------------------------------------------------------
void End()
{

   int NPatch;

   for (int lv=NLEVEL-1; lv>=0; lv--)     
   {
      NPatch = patch.num[lv];

      for (int PID=0; PID<NPatch; PID++)  patch.pdelete( lv, PID );
   }

   if ( Mode_ShellAve )
   {
      if ( Average != NULL )  delete [] Average;
      if ( RMS     != NULL )  delete [] RMS;
      if ( Max     != NULL )  delete [] Max;
      if ( Min     != NULL )  delete [] Min;
      if ( NCount  != NULL )  delete [] NCount;
      if ( Volume  != NULL )  delete [] Volume;
   }

} // FUNCTION : End



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ShellAve
// Description :  Allocate memory and initialize all variables for the mode "shell average"
//-------------------------------------------------------------------------------------------------------
void Init_ShellAve()
{

   cout << "Init_ShellAve ... " << flush;


// allocate memory
   Average  = new double   [NShell][NCOMP];
   RMS      = new double   [NShell][NCOMP];
   Max      = new real     [NShell][NCOMP];
   Min      = new real     [NShell][NCOMP];
   NCount   = new long int [NShell];
   Volume   = new double   [NShell];


// initialize variables for storing average data
   for (int n=0; n<NShell; n++)
   {
      for (int v=0; v<NCOMP; v++)
      {
         Average[n][v] = 0.0;
         RMS    [n][v] = 0.0;
         Max    [n][v] = -__FLT_MAX__;
         Min    [n][v] = +__FLT_MAX__;
      }

      NCount[n] = 0;
      Volume[n] = 0.0;
   }


   cout << "done" << endl;

} // FUNCTION : Init_ShellAve



//-------------------------------------------------------------------------------------------------------
// Function    :  SetMaxRhoPos
// Description :  Set the sphere center to the center of mass of the AveN highest-density cells
//
// Parameter   :  AveN  : Number of highest-density cells for determining the sphere center
//-------------------------------------------------------------------------------------------------------
void SetMaxRhoPos( const int AveN )
{

   cout << "SetMaxRhoPos ..." << endl;


// check
   if ( AveN <= 0 )
   {
      fprintf( stderr, "ERROR : AveN (%d) <= 0 in %s !!\n", AveN, __FUNCTION__ );
      exit( 1 );
   }


   const real dh_min = patch.dh[NLEVEL-1];
   real  MinRho      = __FLT_MIN__;
   int   MinRho_ID   = 0;

   real MaxRho[AveN], MaxRho_Pos[AveN][3], Rho, scale, x, y, z;
   int  MaxRho_Lv[AveN];


// initialize the array recording the maximum density
   for (int t=0; t<AveN; t++)    MaxRho[t] = MinRho;


// begin to search the maximum density
   for (int lv=0; lv<NLEVEL; lv++)              
   {  
      scale = (real)patch.scale[lv];

      for (int PID=0; PID<patch.num[lv]; PID++) 
      {
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
         {
            Rho = patch.ptr[lv][PID]->fluid[DENS][k][j][i];

            if ( Rho > MinRho )
            {
//             record the density and coordinates 
               x = patch.ptr[lv][PID]->corner[0] + (i+0.5)*scale;
               y = patch.ptr[lv][PID]->corner[1] + (j+0.5)*scale;
               z = patch.ptr[lv][PID]->corner[2] + (k+0.5)*scale;

               MaxRho    [MinRho_ID]    = Rho;
               MaxRho_Lv [MinRho_ID]    = lv;
               MaxRho_Pos[MinRho_ID][0] = x;
               MaxRho_Pos[MinRho_ID][1] = y;
               MaxRho_Pos[MinRho_ID][2] = z;


//             find the density threshold 
               MinRho = __FLT_MAX__;

               for (int t=0; t<AveN; t++)
               {
                  if ( MaxRho[t] < MinRho )
                  {
                     MinRho    = MaxRho[t];
                     MinRho_ID = t;
                  } 
               }

            } // if ( Rho > MinRho )
         } // i, j, k
      } // for (int PID=0; PID<patch.num[lv]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)


// sorting
   real TempRho, TempPos[3];
   int  TempLv;
   for (int n=0; n<AveN-1; n++)
   for (int m=0; m<AveN-1-n; m++)
   {
      if ( MaxRho[m] < MaxRho[m+1] )
      {
         TempRho = MaxRho   [m];
         TempLv  = MaxRho_Lv[m];
         for (int dim=0; dim<3; dim++)    TempPos        [dim] = MaxRho_Pos[m  ][dim];

         MaxRho   [m] = MaxRho   [m+1];
         MaxRho_Lv[m] = MaxRho_Lv[m+1];
         for (int dim=0; dim<3; dim++)    MaxRho_Pos[m  ][dim] = MaxRho_Pos[m+1][dim];

         MaxRho   [m+1] = TempRho;
         MaxRho_Lv[m+1] = TempLv;
         for (int dim=0; dim<3; dim++)    MaxRho_Pos[m+1][dim] = TempPos        [dim];
      }
   }

   printf( "===================================================================================\n" );
   printf( "Maximum density\n" );
   printf( "%7s  %2s  (%13s, %13s, %13s)  %13s\n", "Ranking", "Lv", "x", "y", "z", "Density" );

   for (int t=0; t<AveN; t++)
      printf( "%7d  %2d  (%13.7e, %13.7e, %13.7e)  %13.7e\n", 
              t, MaxRho_Lv[t], MaxRho_Pos[t][0]*dh_min, MaxRho_Pos[t][1]*dh_min, MaxRho_Pos[t][2]*dh_min, 
              MaxRho[t] );
   printf( "===================================================================================\n" );


// set the new sphere center as the center of mass of AveN highest-density cells
   double PosSum[3] = { 0.0, 0.0, 0.0 };
   double TotalMass = 0.0;
   double dv, Mass;

   for (int t=0; t<AveN; t++)
   {
      dv   = pow( patch.dh[ MaxRho_Lv[t] ], 3.0 );
      Mass = MaxRho[t]*dv;

      for (int dim=0; dim<3; dim++)    PosSum[dim] += MaxRho_Pos[t][dim]*Mass;

      TotalMass += Mass;
   }

   for (int dim=0; dim<3; dim++)    Center[dim] = PosSum[dim] / TotalMass;


// remove the allocated patches (because the candidate box may be wrong due to the incorrect sphere center)
   int NPatch;

   for (int lv=NLEVEL-1; lv>=0; lv--)     
   {
      NPatch = patch.num[lv];

      for (int PID=0; PID<NPatch; PID++)  patch.pdelete( lv, PID );
   }


// reload the patches with the new sphere center
   LoadData();


   cout << "SetMaxRhoPos ... done" << endl;

} // FUNCTION : SetMaxRhoPos



//-------------------------------------------------------------------------------------------------------
// Function    :  TakeNote
// Description :  Output the simulation parameters 
//-------------------------------------------------------------------------------------------------------
void TakeNote()
{

   cout << "TakeNote ... " << endl;


   const real dh_min = patch.dh[NLEVEL-1];

   printf( "===============================================\n" );
#  if   ( MODEL == HYDRO )
   printf( "MODEL          = %14s\n",     "HYDRO"              );
#  elif ( MODEL == MHD )
   printf( "MODEL          = %14s\n",     "MHD"                );
#  elif ( MODEL == ELBDM )
   printf( "MODEL          = %14s\n",     "ELBDM"              );
#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL
   printf( "DumpID         = %14d\n",     DumpID               );
   printf( "Time           = %14.7e\n",   Time[0]              );
   printf( "Step           = %14ld\n",    Step                 );
   printf( "NX0_TOT[0]     = %14d\n",     NX0_TOT[0]           );
   printf( "NX0_TOT[1]     = %14d\n",     NX0_TOT[1]           );
   printf( "NX0_TOT[2]     = %14d\n",     NX0_TOT[2]           );
   printf( "Center_x       = %14.7e\n",   Center[0]*dh_min     );
   printf( "Center_y       = %14.7e\n",   Center[1]*dh_min     );
   printf( "Center_z       = %14.7e\n",   Center[2]*dh_min     );
   printf( "MappedCenter_x = %14.7e\n",   Center_Map[0]*dh_min );
   printf( "MappedCenter_y = %14.7e\n",   Center_Map[1]*dh_min );
   printf( "MappedCenter_z = %14.7e\n",   Center_Map[2]*dh_min );
   printf( "Radius         = %14.7e\n",   MaxRadius*dh_min     );
   printf( "UseMaxRhoPos   = %14d\n",     UseMaxRhoPos         );
   printf( "Periodic       = %14d\n",     Periodic             );
   printf( "InputScale     = %14d\n",     InputScale           );
   printf( "GetNShell      = %14d\n",     GetNShell            );
   printf( "NShell         = %14d\n",     NShell               );
   printf( "ShellWidth     = %14.7e\n",   ShellWidth*dh_min    );
   printf( "RhoThres       = %14.7e\n",   RhoThres             );
   printf( "===============================================\n" );


   cout << "TakeNote ... done" << endl;

} // FUNCTION : TakeNote



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


   fprintf( stdout, "   Loading simulation parameters ...\n" );


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

      fprintf( stderr, "WARNING : loading data with format version < 1102 --> assuming BOX_SIZE = %f\n", 
               box_size );
   }


   fprintf( stdout, "   Loading simulation parameters ... done\n" );


// d. check parameters (before loading any size-dependent parameters)
// =================================================================================================
   fprintf( stdout, "   Checking loaded parameters ...\n" );


   const bool Fatal    = true;
   const bool NonFatal = false;

#  ifdef FLOAT8
   if ( !float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "OFF", "ON" );
      exit( 1 );
   }
#  else
   if (  float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "ON", "OFF" );
      exit( 1 );
   }
#  endif

   CompareVar( "NLEVEL",                  nlevel,                 NLEVEL,                       Fatal );
   CompareVar( "NCOMP",                   ncomp,                  NCOMP,                        Fatal );
   CompareVar( "PATCH_SIZE",              patch_size,             PATCH_SIZE,                   Fatal );


   fprintf( stdout, "   Checking loaded parameters ... done\n" );


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


   fprintf( stdout, "   Loading simulation parameters ...\n" );


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


   fprintf( stdout, "   Loading simulation parameters ... done\n" );


// d. check parameters (before loading any size-dependent parameters)
// =================================================================================================
   fprintf( stdout, "   Checking loaded parameters ...\n" );


   const bool Fatal    = true;
   const bool NonFatal = false;

#  ifdef FLOAT8
   if ( !float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "OFF", "ON" );
      exit( 1 );
   }
#  else
   if (  float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "ON", "OFF" );
      exit( 1 );
   }
#  endif

   CompareVar( "MODEL",                   model,                  MODEL,                        Fatal );
   CompareVar( "NLEVEL",                  nlevel,                 NLEVEL,                       Fatal );
   CompareVar( "NCOMP",                   ncomp,                  NCOMP,                        Fatal );
   CompareVar( "PATCH_SIZE",              patch_size,             PATCH_SIZE,                   Fatal );


   fprintf( stdout, "   Checking loaded parameters ... done\n" );


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
         exit( 1 );
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
         exit( 1 );
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
         exit( 1 );
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
         exit( 1 );
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
         exit( 1 );
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n", 
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (double)



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :  
//-------------------------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{

   ReadOption( argc, argv );

   LoadData();

   if ( UseMaxRhoPos > 0 )    SetMaxRhoPos( UseMaxRhoPos );

   TakeNote();

   CheckParameter();

   if ( Mode_MaxRho )   GetMaxRho();

   if ( Mode_ShellAve )
   {
      Init_ShellAve();

      ShellAverage();

      GetRMS();

      Output_ShellAve();
   }

   End();


   cout << "Program terminated successfully" << endl;

   return 0;

} // FUNCTION : main
