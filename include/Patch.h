
#ifndef __PATCH_H__
#define __PATCH_H__



#include "Macro.h"

void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void Aux_Message( FILE *Type, const char *Format, ... );
#ifdef LOAD_BALANCE
long LB_Corner2Index( const int lv, const int Corner[], const Check_t Check );
long Mis_Idx3D2Idx1D( const int Size[], const int Idx3D[] );
#endif




//-------------------------------------------------------------------------------------------------------
// Structure   :  patch_t 
// Description :  Data structure of a single patch 
//
// Data Member :  fluid          : Fluid variables (mass density, momentum density x, y ,z, energy density) 
//                pot            : Potential
//                flux[6]        : Fluid flux (for the flux-correction operation)
//                flux_debug[6]  : Fluid flux for the debug mode (ensuring that the round-off errors are 
//                                 exactly the same in different parallelization parameters/strategies)
//                corner[3]      : Grid indices of the cell at patch corner
//                sibling[26]    : Patch IDs of the 26 sibling patches
//                father         : Patch ID of the father patch
//                son            : Patch ID of the child patch
//                flag           : Refinement flag
//                LB_Idx         : Space-filling-curve index for load balance
//                PaddedCr1D     : 1D corner coordiniate padded with two base-level patches on each side 
//                                 in each direction, normalized to the finest-level patch scale (PATCH_SIZE)
//                                 --> each PaddedCr1D defines a unique 3D position
//                                 --> patches at different levels with the same PaddedCr1D have the same 
//                                     3D corner coordinates
// Method      :  patch_t        : Constructor 
//               ~patch_t        : Destructor
//                fnew           : Allocate one flux array 
//                fdelete        : Deallocate one flux array
//                hnew           : Allocate hydrodynamic array
//                hdelete        : Deallocate hydrodynamic array
//                gnew           : Allocate potential array
//                gdelete        : Deallocate potential array
//-------------------------------------------------------------------------------------------------------
struct patch_t
{

// data members
// ===================================================================================
   real (*fluid)[PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
#  ifdef GRAVITY
   real (*pot  )[PATCH_SIZE][PATCH_SIZE];
#  endif

   real (*flux      [6])[PATCH_SIZE][PATCH_SIZE];
#  ifdef DAINO_DEBUG
   real (*flux_debug[6])[PATCH_SIZE][PATCH_SIZE];
#  endif

   int  corner[3];
   int  sibling[26];
   int  father;
   int  son;
   bool flag;
#  ifdef LOAD_BALANCE
   long LB_Idx;
   long PaddedCr1D;
#  endif



   //===================================================================================
   // Constructor :  patch_t 
   // Description :  Constructor of the structure "patch_t"
   //
   // Note        :  Initialize the data members
   //
   // Parameter   :  x,y,z    : Physical coordinates of the patch corner
   //                FaPID    : Patch ID of the father patch    
   //                FluData  : true --> Allocate hydrodynamic array "fluid"
   //                PotData  : true --> Allocate potential array "pot" (has no effect if "GRAVITY" is turned off)
   //                lv       : Refinement level of the newly created patch
   //                BoxScale : Simulation box scale
   //===================================================================================
   patch_t( const int x, const int y, const int z, const int FaPID, const bool FluData, const bool PotData,
            const int lv, const int BoxScale[] )
   {

      corner[0] = x; 
      corner[1] = y;
      corner[2] = z;
      father    = FaPID;
      son       = -1;

#     ifdef LOAD_BALANCE
      const int Padded              = 1<<NLEVEL;
      const int BoxNScale_Padded[3] = { BoxScale[0]/PATCH_SIZE + 2*Padded,
                                        BoxScale[1]/PATCH_SIZE + 2*Padded,
                                        BoxScale[2]/PATCH_SIZE + 2*Padded }; // normalized and padded box scale
      int Cr_Padded[3];
      for (int d=0; d<3; d++)    Cr_Padded[d] = corner[d]/PATCH_SIZE + Padded;

#     ifdef DAINO_DEBUG
      for (int d=0; d<3; d++)
      {
         if ( Cr_Padded[d] < 0  ||  Cr_Padded[d] >= BoxNScale_Padded[d] )
            Aux_Error( ERROR_INFO, "incorrect Cr_Padded[%d]=(%d) --> out of range [%d ... %d) !!\n",
                       d, Cr_Padded[d], 0, BoxNScale_Padded[d] );
      }
#     endif

      PaddedCr1D = Mis_Idx3D2Idx1D( BoxNScale_Padded, Cr_Padded );
      LB_Idx     = LB_Corner2Index( lv, corner, CHECK_OFF );   // this number can be wrong for external patches
#     endif // #ifdef LOAD_BALANCE
      
      for (int s=0; s<26; s++ )  sibling[s] = -1;     // -1 <--> NO sibling

      flag  = false;
      fluid = NULL;
#     ifdef GRAVITY
      pot   = NULL;
#     endif

      for (int s=0; s<6; s++)    
      {
         flux      [s] = NULL;
#        ifdef DAINO_DEBUG
         flux_debug[s] = NULL;
#        endif
      }

      if ( FluData )    hnew();
#     ifdef GRAVITY
      if ( PotData )    gnew();
#     endif
   } // METHOD : patch_t



   //===================================================================================
   // Destructor  :  ~patch_t 
   // Description :  Destructor of the structure "patch_t"
   //
   // Note        :  Deallocate flux and data arrays
   //===================================================================================
   ~patch_t()
   {
      fdelete();
      hdelete();
#     ifdef GRAVITY
      gdelete();
#     endif
   } // METHOD : ~patch_t



   //===================================================================================
   // Method      :  fnew 
   // Description :  Allocate flux array in the given direction
   //
   // Note        :  Flux array is initialized as zero
   //
   // Parameter   :  SibID : Targeted ID of the flux array (0,1,2,3,4,5) <--> (-x,+x,-y,+y,-z,+z) 
   //===================================================================================
   void fnew( const int SibID )
   {
#     ifdef DAINO_DEBUG
      if ( SibID < 0  ||  SibID > 5 )
         Aux_Error( ERROR_INFO, "incorrect input in the member function \"fnew\" !!\n" );

      if ( flux[SibID] != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing flux array (sibling = %d) !!\n", SibID );

      if ( flux_debug[SibID] != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing flux_debug array (sibling = %d) !!\n", SibID );
#     endif

      flux[SibID] = new real [NCOMP][PATCH_SIZE][PATCH_SIZE];
            
      for(int v=0; v<NCOMP; v++)
      for(int m=0; m<PATCH_SIZE; m++)
      for(int n=0; n<PATCH_SIZE; n++)
         flux[SibID][v][m][n] = 0.0;

#     ifdef DAINO_DEBUG
      flux_debug[SibID] = new real [NCOMP][PATCH_SIZE][PATCH_SIZE];
            
      for(int v=0; v<NCOMP; v++)
      for(int m=0; m<PATCH_SIZE; m++)
      for(int n=0; n<PATCH_SIZE; n++)
         flux_debug[SibID][v][m][n] = 0.0;
#     endif
   } // METHOD : fnew



   //===================================================================================
   // Method      :  fdelete
   // Description :  Deallocate all flux arrays allocated previously 
   //===================================================================================
   void fdelete()
   {
      for (int s=0; s<6; s++)
      {
         if ( flux[s] != NULL )
         {
            delete [] flux[s];
            flux[s] = NULL;

#           ifdef DAINO_DEBUG
            delete [] flux_debug[s];
            flux_debug[s] = NULL;
#           endif
         }
      }
   } // METHOD : fdelete



   //===================================================================================
   // Method      :  hnew 
   // Description :  Allocate hydrodynamic array
   //
   // Note        :  An error message will be displayed if array has already been allocated
   //===================================================================================
   void hnew()
   {
#     ifdef DAINO_DEBUG
      if ( fluid != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing fluid array !!\n" );
#     endif

      fluid = new real [NCOMP][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
      fluid[0][0][0][0] = -1; 
   } // METHOD : hnew



   //===================================================================================
   // Method      :  hdelete
   // Description :  Deallocate hydrodynamic array
   //===================================================================================
   void hdelete()
   {
      if ( fluid != NULL )    
      {
         delete [] fluid;
         fluid = NULL;
      }
   } // METHOD : hdelete



#  ifdef GRAVITY
   //===================================================================================
   // Method      :  gnew 
   // Description :  Allocate potential array
   //
   // Note        :  An error message will be displayed if array has already been allocated
   //===================================================================================
   void gnew()
   {
#     ifdef DAINO_DEBUG
      if ( pot != NULL )
         Aux_Error( ERROR_INFO, "allocate an existing pot array !!\n" );
#     endif

      pot = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
   } // METHOD : gnew



   //===================================================================================
   // Method      :  gdelete
   // Description :  Deallocate potential array
   //===================================================================================
   void gdelete()
   {
      if ( pot != NULL )    
      {
         delete [] pot;
         pot = NULL;
      }
   } // METHOD : gdelete
#  endif // #ifdef GRAVITY

   
}; // struct patch_t



#endif // #ifndef __PATCH_H__
