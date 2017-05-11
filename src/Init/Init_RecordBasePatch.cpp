
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_RecordBasePatch
// Description :  Record the IDs of base-level patches in the BaseP array 
//-------------------------------------------------------------------------------------------------------
void Init_RecordBasePatch()
{

   const int scale0      = patch->scale[0];
   const int NPatch1D[3] = { NX0[0]/PATCH_SIZE+4, NX0[1]/PATCH_SIZE+4, NX0[2]/PATCH_SIZE+4 };

   int order[3], P2[8];    // order[3] : (i,j,k)th patch in the (x,y,z) direction


   for (int P=0; P<patch->num[0]; P+=8)
   {
      for (int d=0; d<3; d++)
         order[d] = (patch->ptr[0][0][P]->corner[d] - DAINO_RANK_X(d)*NX0[d]*scale0) / (PATCH_SIZE*scale0) + 2;

      P2[0] = (order[2]+0)*NPatch1D[1]*NPatch1D[0] + (order[1]+0)*NPatch1D[0] + (order[0]+0);
      P2[1] = (order[2]+0)*NPatch1D[1]*NPatch1D[0] + (order[1]+0)*NPatch1D[0] + (order[0]+1);
      P2[2] = (order[2]+0)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]+0);
      P2[3] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]+0)*NPatch1D[0] + (order[0]+0);
      P2[4] = (order[2]+0)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]+1);
      P2[5] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]+0);
      P2[6] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]+0)*NPatch1D[0] + (order[0]+1);
      P2[7] = (order[2]+1)*NPatch1D[1]*NPatch1D[0] + (order[1]+1)*NPatch1D[0] + (order[0]+1);

//    record the patch ID in the BaseP array
      for (int m=0; m<8; m++)    BaseP[ P2[m] ] = P + m;
   }


// record the buffer patches by the periodic B.C. for the base-level sibling search in the serial code 
#  ifdef SERIAL
   int i2, j2, k2, Width[3], Disp1[3], Disp2[3], NP, ID1, ID2;

   for (int Sib=0; Sib<26; Sib++)
   {
      for (int d=0; d<3; d++)    
      {
         NP       = NX0[d]/PATCH_SIZE;

         Width[d] = TABLE_01( Sib, 'x'+d,  1, NP,    1 );
         Disp1[d] = TABLE_01( Sib, 'x'+d,  1,  2, 2+NP );
         Disp2[d] = TABLE_01( Sib, 'x'+d, NP,  0,  -NP );
      }

      for (int k1=Disp1[2]; k1<Disp1[2]+Width[2]; k1++)  {  k2 = k1 + Disp2[2]; 
      for (int j1=Disp1[1]; j1<Disp1[1]+Width[1]; j1++)  {  j2 = j1 + Disp2[1];
      for (int i1=Disp1[0]; i1<Disp1[0]+Width[0]; i1++)  {  i2 = i1 + Disp2[0];

         ID1 = ( k1*NPatch1D[1] + j1 )*NPatch1D[0] + i1;
         ID2 = ( k2*NPatch1D[1] + j2 )*NPatch1D[0] + i2;

         BaseP[ID1] = BaseP[ID2];

      }}}
   } // for (int Sib=0; Sib<26; Sib++)
#  endif

} // FUNCTION : Init_RecordBasePatch
