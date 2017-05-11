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


// patch size (number of cells of a single patch in the x/y/z directions)
#define PATCH_SIZE         8
#define PS1                (1*PATCH_SIZE)
#define PS2                (2*PATCH_SIZE)


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
//                sibling[26] : patch IDs of the 26 sibling patches
//                father      : patch ID of the father patch
//                son         : patch ID of the child patch
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
   int  sibling[26];
   int  father;
   int  son;



   //===================================================================================
   // Constructor :  patch_t 
   // Description :  constructor of the structure "patch_t"
   //
   // Note        :  initialize data members
   //
   // Parameter   :  x,y,z : physical coordinates of the patch corner
   //                FaPID : patch ID of the father patch    
   //                Data  : true --> allocate physical data (fluid + pot)
   //===================================================================================
   patch_t( const int x, const int y, const int z, const int FaPID, const bool Data )
   {
      corner[0] = x; 
      corner[1] = y;
      corner[2] = z;
      father    = FaPID;
      son       = -1;
      
      for (int s=0; s<26; s++ )     sibling[s] = -1;     // -1 <--> NO sibling

      fluid = NULL;
      pot   = NULL;

      if ( Data )
      {
         fluid = new real [NCOMP][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
         pot   = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];

         fluid[0][0][0][0] = -1; 
      }
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
// Data Member :  ptr      : Pointer of each patch
//                num      : Number of patches (real patch + buffer patch) at each level
//                scale    : Grid size at each level (normalize to the grid size at the finest level 
//                dh       : Grid size at each level
//                BoxSize  : Simulation box size
//                BoxScale : Simulation box scale
//
// Method      :  pnew    : Allocate one patch
//                pdelete : Deallocate one patch
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
   // Note        :  a. each patch contains two patch pointers --> SANDGLASS (Sg) = 0 / 1 
   //                b. Sg = 0 : store both data and relation (father,son.sibling,corner,flag,flux)
   //                   Sg = 1 : store only data 
   //
   // Parameter   :  lv    : the targeted refinement level
   //                x,y,z : physical coordinates of the patch corner
   //                FaPID : the patch ID of the parent patch at level "lv-1"
   //                Data  : true --> allocate physical data (fluid + pot)
   //===================================================================================
   void pnew( const int lv, const int x, const int y, const int z, const int FaPID, const bool Data )
   {
      if ( num[lv] == MAX_PATCH )
      {
         fprintf( stderr, "ERROR : \"exceeding the maximum number of patches [%d] !!\n", MAX_PATCH );
         exit(-1);
      }

      if ( ptr[lv][num[lv]] != NULL )
      {
         fprintf( stderr, "ERROR : \"allocate an existing patch !!\" : Lv %d, PID %d, FaPID %d\n", 
                  lv, num[lv], FaPID );
         exit(-1);
      }

      ptr[lv][ num[lv] ] = new patch_t( x, y, z, FaPID, Data );

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
   // Note        :  a. this function should NOT be applied to the base-level patches
   //                b. this function will also deallocate the flux arrays of the targeted patch
   //                c. delete a patch with son is forbidden
   //
   // Parameter   :  lv  : the targeted refinement level
   //                PID : the patch ID to be removed
   //===================================================================================
   void pdelete( const int lv, const int PID )
   {
      if ( lv == 0 )
      {
         fprintf( stderr, "ERROR : delete a base-level patch !!\n" );
         exit(-1);
      }

      if ( ptr[lv][PID]->son != -1 )
      {
         fprintf( stderr, "ERROR : \"delete a patch with son !!\" : Lv %d, PID %d, SonPID %d\n", 
                  lv, PID, ptr[lv][PID]->son );
         exit(-1);
      }

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





//-------------------------------------------------------------------------------------------------------
// Structure   :  ParaVar_t
// Description :  data structure collecting variables related to DAINO parallelization
//
// Data Member :  BounP_NList    : the number of boundary patches in 26 sibling directions
//                BounP_IDList   : the IDs of boundary patches in 26 sibling directions
//                BounP_PosList  : the positions of boundary patches recorded in "BounP_IDList"
//
// Method      :  
//-------------------------------------------------------------------------------------------------------
struct ParaVar_t
{

// data members
// ===================================================================================
   int  BounP_NList      [NLEVEL  ][26];
   int *BounP_IDList     [NLEVEL  ][26];
   int *BounP_PosList    [NLEVEL  ][26];

   int  SendP_NList      [NLEVEL  ][26];
   int *SendP_IDList     [NLEVEL  ][26];
   int  RecvP_NList      [NLEVEL  ][26];
   int *RecvP_IDList     [NLEVEL  ][26];



   //===================================================================================
   // Constructor :  ParaVar_t 
   // Description :  constructor of the structure "ParaVar_t"
   //
   // Note        :  initialize all pointers as NULL and all counters as 0
   //===================================================================================
   ParaVar_t()
   {
      for (int lv=0; lv<NLEVEL; lv++)
      for (int s=0; s<26; s++)
      {
         BounP_NList       [lv][s] = 0;
         BounP_IDList      [lv][s] = NULL;
         BounP_PosList     [lv][s] = NULL;

         SendP_NList       [lv][s] = 0;
         SendP_IDList      [lv][s] = NULL;
         RecvP_NList       [lv][s] = 0;
         RecvP_IDList      [lv][s] = NULL;
      }
   }


}; // struct ParaVar_t



#endif // #ifndef __TYPEDEF_H__
