#include "GetCube.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_RecordBasePatch
// Description :  Record the IDs of base-level patches in the BasePID array 
//-------------------------------------------------------------------------------------------------------
void Init_RecordBasePatch()
{

// record the patch ID
   const int scale         = patch.scale[0];
   const int NPatch1D[3]   = { NX0[0]/PATCH_SIZE+4, NX0[1]/PATCH_SIZE+4, NX0[2]/PATCH_SIZE+4 };
   int ip, jp, kp, P2[8];              // (i,j,k)p : (i,j,k)th patch in (x,y,z) direction

   for (int P=0; P<patch.num[0]; P+=8)
   {
      ip = (patch.ptr[0][P]->corner[0] - MyRank_X[0]*NX0[0]*scale) / (PATCH_SIZE*scale) + 2;
      jp = (patch.ptr[0][P]->corner[1] - MyRank_X[1]*NX0[1]*scale) / (PATCH_SIZE*scale) + 2;
      kp = (patch.ptr[0][P]->corner[2] - MyRank_X[2]*NX0[2]*scale) / (PATCH_SIZE*scale) + 2;

      P2[0] = (kp+0)*NPatch1D[1]*NPatch1D[0] + (jp+0)*NPatch1D[0] + (ip+0);
      P2[1] = (kp+0)*NPatch1D[1]*NPatch1D[0] + (jp+0)*NPatch1D[0] + (ip+1);
      P2[2] = (kp+0)*NPatch1D[1]*NPatch1D[0] + (jp+1)*NPatch1D[0] + (ip+0);
      P2[3] = (kp+1)*NPatch1D[1]*NPatch1D[0] + (jp+0)*NPatch1D[0] + (ip+0);
      P2[4] = (kp+0)*NPatch1D[1]*NPatch1D[0] + (jp+1)*NPatch1D[0] + (ip+1);
      P2[5] = (kp+1)*NPatch1D[1]*NPatch1D[0] + (jp+1)*NPatch1D[0] + (ip+0);
      P2[6] = (kp+1)*NPatch1D[1]*NPatch1D[0] + (jp+0)*NPatch1D[0] + (ip+1);
      P2[7] = (kp+1)*NPatch1D[1]*NPatch1D[0] + (jp+1)*NPatch1D[0] + (ip+1);

//    record the patch ID in the BaseP array
      for (int m=0; m<8; m++)    BaseP[ P2[m] ] = P + m;
   }

}



