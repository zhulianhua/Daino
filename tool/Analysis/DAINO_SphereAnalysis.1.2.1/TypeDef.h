#ifndef __TYPEDEF_H__
#define __TYPEDEF_H__



#ifdef FLOAT8
typedef double real;
#else
typedef float  real;
#endif


// models
#define HYDRO              1
#define MHD                2
#define ELBDM              3


#ifndef NULL
#define NULL               0
#endif


#ifndef __FLT_MIN__
#  define __FLT_MIN__      1.17549435e-38F
#endif


#ifndef __FLT_MAX__
#  define __FLT_MAX__      3.40282347e+38F
#endif


// patch size (number of cells of a single patch in the x/y/z directions)
#define PATCH_SIZE         8


// number of components in each cell and the variable indices in the array "fluid"
#if   ( MODEL == HYDRO )
#  define NCOMP            5
#  define DENS             0
#  define MOMX             1
#  define MOMY             2
#  define MOMZ             3
#  define ENGY             4

#elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  define NCOMP            8

#elif ( MODEL == ELBDM )
#  define NCOMP            3
#  define DENS             0
#  define REAL             1
#  define IMAG             2

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL




//-------------------------------------------------------------------------------------------------------
// Structure   :  patch_t 
// Description :  data structure of a single patch 
//
// Data Member :  fluid       : fluid variables (mass density, momentum density x, y ,z, energy density) 
//                pot         : potential
//                corner[3]   : physical coordinates of the patch corner
//
// Method      :  patch_t     : constructor 
//                ~patch_t    : destructor
//-------------------------------------------------------------------------------------------------------
struct patch_t
{

// data members
// ===================================================================================
   real (*fluid)[PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
   real (*pot  )[PATCH_SIZE][PATCH_SIZE];

   int  corner[3];



   //===================================================================================
   // Constructor :  patch_t 
   // Description :  constructor of the structure "patch_t"
   //
   // Note        :  initialize data members
   //
   // Parameter   :  x,y,z : physical coordinates of the patch corner
   //===================================================================================
   patch_t( const int x, const int y, const int z )
   {
      corner[0] = x; 
      corner[1] = y;
      corner[2] = z;
      
      fluid = new real [NCOMP][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
      pot   = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];

      fluid[0][0][0][0] = -1; 
   }



   //===================================================================================
   // Destructor  :  ~patch_t 
   // Description :  destructor of the structure "patch_t"
   //
   // Note        :  deallocate flux and data arrays
   //===================================================================================
   ~patch_t()
   {
      if ( fluid != NULL )    
      {
         delete [] fluid;
         fluid = NULL;
      }

      if ( pot != NULL )    
      {
         delete [] pot;
         pot   = NULL;
      }
   }

   
}; // struct patch_t





//-------------------------------------------------------------------------------------------------------
// Structure   :  DAINO_t
// Description :  data structure of DAINO
//
// Data Member :  ptr      : pointer of each patch
//                num      : number of patches (real patch + buffer patch) at each level
//                scale    : grid size at each level (normalize to the grid size at the finest level 
//                dh       : Grid size at each level
//                BoxSize  : Simulation box size
//                BoxScale : Simulation box scale
//
// Method      :  pnew     : allocate one patch
//                pdelete  : deallocate one patch
//-------------------------------------------------------------------------------------------------------
struct DAINO_t
{

// data members
// ===================================================================================
   patch_t *ptr[NLEVEL][MAX_PATCH];

   int    num  [NLEVEL];          
   int    scale[NLEVEL];        
   double dh   [NLEVEL];
   double BoxSize [3];
   int    BoxScale[3];
   


   //===================================================================================
   // Constructor :  DAINO_t
   // Description :  constructor of the structure "DAINO_t"
   //
   // Note        :  initialize data members
   //===================================================================================
   DAINO_t()
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
         num  [lv] = 0;
         scale[lv] = 1<<(NLEVEL-1-lv);
      }

      for (int lv=0; lv<NLEVEL; lv++)
      for (int PID=0; PID<MAX_PATCH; PID++)
         ptr[lv][PID] = NULL;
   }



   //===================================================================================
   // Method      :  pnew
   // Description :  allocate a single patch 
   //
   // Parameter   :  lv    : the targeted refinement level
   //                x,y,z : physical coordinates of the patch corner
   //===================================================================================
   void pnew( const int lv, const int x, const int y, const int z )
   {
      if ( ptr[lv][num[lv]] != NULL )
      {
         fprintf( stderr, "ERROR : \"allocate an existing patch !!\" : Lv %d, PID %d\n", lv, num[lv] );
         exit(-1);
      }

      ptr[lv][ num[lv] ] = new patch_t( x, y, z );

      num[lv] ++;

      if ( num[lv] > MAX_PATCH )
      {
         fprintf( stderr, "ERROR : exceed MAX_PATCH !!\n" );
         exit(-1);
      }
   }



   //===================================================================================
   // Method      :  pdelete
   // Description :  deallocate a single patch 
   //
   // Parameter   :  lv  : the targeted refinement level
   //                PID : the patch ID to be removed
   //===================================================================================
   void pdelete( const int lv, const int PID )
   {
      delete ptr[lv][PID];

      ptr[lv][PID] = NULL;

      num[lv] --;

      if ( num[lv] < 0 )
      {
         fprintf( stderr, "ERROR : num[%d] = %d < 0 !!\n", lv, num[lv] );
         exit(-1);
      }
   }


}; // struct DAINO_t



#endif // #ifndef __TYPEDEF_H__


