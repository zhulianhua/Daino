
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_BaseLevel
// Description :  Construct the base level during initialization 
//-------------------------------------------------------------------------------------------------------
void Init_BaseLevel()
{

   const int NPatch[3] = { NX0[0]/PATCH_SIZE, NX0[1]/PATCH_SIZE, NX0[2]/PATCH_SIZE };
   const int scale0    = patch->scale[0];
   const int Width     = PATCH_SIZE*scale0;

   int Cr[3];    


// allocate the real base-level patches
   for (int Pz=0; Pz<NPatch[2]; Pz+=2)    {  Cr[2] = DAINO_RANK_X(2)*NX0[2]*scale0 + Pz*PATCH_SIZE*scale0;
   for (int Py=0; Py<NPatch[1]; Py+=2)    {  Cr[1] = DAINO_RANK_X(1)*NX0[1]*scale0 + Py*PATCH_SIZE*scale0;
   for (int Px=0; Px<NPatch[0]; Px+=2)    {  Cr[0] = DAINO_RANK_X(0)*NX0[0]*scale0 + Px*PATCH_SIZE*scale0;

      patch->pnew( 0, Cr[0],       Cr[1],       Cr[2],       -1, true, true );
      patch->pnew( 0, Cr[0]+Width, Cr[1],       Cr[2],       -1, true, true );
      patch->pnew( 0, Cr[0],       Cr[1]+Width, Cr[2],       -1, true, true );
      patch->pnew( 0, Cr[0],       Cr[1],       Cr[2]+Width, -1, true, true );
      patch->pnew( 0, Cr[0]+Width, Cr[1]+Width, Cr[2],       -1, true, true );
      patch->pnew( 0, Cr[0],       Cr[1]+Width, Cr[2]+Width, -1, true, true );
      patch->pnew( 0, Cr[0]+Width, Cr[1],       Cr[2]+Width, -1, true, true );
      patch->pnew( 0, Cr[0]+Width, Cr[1]+Width, Cr[2]+Width, -1, true, true );

   }}}

   for (int m=1; m<28; m++)   patch->NPatchComma[0][m] = patch->num[0];


// record the number of real patches for the out-of-core computing
#  ifdef OOC
   ooc.NRealPatch[ooc.Rank][0] = patch->NPatchComma[0][1];
   ooc.NDataPatch[ooc.Rank][0] = 0;

   for (int PID=0; PID<patch->NPatchComma[0][1]; PID++)
      if ( patch->ptr[0][0][PID]->son == -1 )   ooc.NDataPatch[ooc.Rank][0]++;  
#  endif

// allocate the base-level buffer patches
#  ifdef OOC
   Buf_AllocateBufferPatch( patch, 0, 3, ooc.Rank );
#  else
   Buf_AllocateBufferPatch( patch, 0, 3, NULL_INT );
#  endif

// set up the BaseP List
   Init_RecordBasePatch();

// set up the BounP_IDMap[0]
   Buf_RecordBoundaryPatch( 0 );

// construct the sibling relation for the base level (including the buffer patches)
   SiblingSearch( 0 );

// get the IDs of patches for sending and receiving data between neighbor ranks
   Buf_RecordExchangeDataPatchID( 0 );

} // FUNCTION : Init_BaseLevel
