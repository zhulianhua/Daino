
#include "Macro.h"
#include "CUPOT.h"

#if ( defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )



#define GRA_NTHREAD  ( PATCH_SIZE*PATCH_SIZE*GRA_BLOCK_SIZE_Z )




//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_HydroGravitySolver
// Description :  GPU gravity solver which advances the momentum and energy density of a group of patches 
//                by gravitational acceleration
//
// Note        :  Prefix "g" for pointers pointing to the "Global" memory space
//                Prefix "s" for pointers pointing to the "Shared" memory space
//
// Parameter   :  g_Flu_Array : Global memory array to store the input and output fluid variables
//                g_Pot_Array : Global memory array storing the input potential for evaluating the 
//                              gravitational acceleration
//                Gra_Const   : 3-P stencil : dt / ( 2*dh) 
//                              5-P stencil : dt / (12*dh)
//                P5_Gradient : Use 5-points stecil to evaluate the potential gradient
//---------------------------------------------------------------------------------------------------
__global__ void CUPOT_HydroGravitySolver(       real g_Flu_Array[][NCOMP][ PATCH_SIZE*PATCH_SIZE*PATCH_SIZE ],
                                          const real g_Pot_Array[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                          const real Gra_Const, const bool P5_Gradient )
{

   const uint bx     = blockIdx.x;
   const uint tx     = threadIdx.x; 
   const uint ty     = threadIdx.y; 
   const uint tz     = threadIdx.z; 
   const uint ID     = __umul24( tz, PS1*PS1 ) + __umul24( ty, PS1 ) + tx;
   const uint NSlice = GRA_BLOCK_SIZE_Z;

   uint s_index =   __umul24( GRA_GHOST_SIZE+tz, GRA_NXT*GRA_NXT ) 
                  + __umul24( GRA_GHOST_SIZE+ty, GRA_NXT ) + (GRA_GHOST_SIZE+tx);
   uint g_index = ID; 

   uint ip1, jp1, kp1, im1, jm1, km1, t;
   real Acc[3], Eint, px, py, pz, rho, temp;

#  if ( GRA_GHOST_SIZE == 2 )
   uint ip2, jp2, kp2, im2, jm2, km2;
#  endif

   __shared__ real s_Pot[GRA_NXT*GRA_NXT*GRA_NXT];


// load the potential from the global memory to the shared memory 
   t = ID;
   do {  s_Pot[t] = g_Pot_Array[bx][t];   t += GRA_NTHREAD; }  while ( t < GRA_NXT*GRA_NXT*GRA_NXT );
   __syncthreads();

   
   for (uint z=tz; z<PS1; z+=NSlice)
   {
      ip1 = s_index + 1;
      jp1 = s_index + GRA_NXT;
      kp1 = s_index + GRA_NXT*GRA_NXT;
      im1 = s_index - 1;
      jm1 = s_index - GRA_NXT;
      km1 = s_index - GRA_NXT*GRA_NXT;

#     if ( GRA_GHOST_SIZE == 2 )
      ip2 = s_index + 2;
      jp2 = s_index + 2*GRA_NXT;
      kp2 = s_index + 2*GRA_NXT*GRA_NXT;
      im2 = s_index - 2;
      jm2 = s_index - 2*GRA_NXT;
      km2 = s_index - 2*GRA_NXT*GRA_NXT;
#     endif

      
//    evalute the gravitational acceleration
#     if ( GRA_GHOST_SIZE == 2 )
      if ( P5_Gradient )
      {
         Acc[0] = Gra_Const * ( - s_Pot[ip2] + (real)8.0*s_Pot[ip1] - (real)8.0*s_Pot[im1] + s_Pot[im2] );
         Acc[1] = Gra_Const * ( - s_Pot[jp2] + (real)8.0*s_Pot[jp1] - (real)8.0*s_Pot[jm1] + s_Pot[jm2] );
         Acc[2] = Gra_Const * ( - s_Pot[kp2] + (real)8.0*s_Pot[kp1] - (real)8.0*s_Pot[km1] + s_Pot[km2] );
      }
      else // P3_Gradient
#     endif
      {
         Acc[0] = Gra_Const * ( s_Pot[ip1] - s_Pot[im1] );
         Acc[1] = Gra_Const * ( s_Pot[jp1] - s_Pot[jm1] );
         Acc[2] = Gra_Const * ( s_Pot[kp1] - s_Pot[km1] );
      }


//    advance the fluid
      rho  = g_Flu_Array[bx][DENS][g_index];
      px   = g_Flu_Array[bx][MOMX][g_index];
      py   = g_Flu_Array[bx][MOMY][g_index];
      pz   = g_Flu_Array[bx][MOMZ][g_index];
      temp = (real)0.5/rho;
      Eint = g_Flu_Array[bx][ENGY][g_index] - temp*(px*px+py*py+pz*pz);
      
      px -= rho*Acc[0];
      py -= rho*Acc[1];
      pz -= rho*Acc[2];
      
      g_Flu_Array[bx][MOMX][g_index] = px;
      g_Flu_Array[bx][MOMY][g_index] = py;
      g_Flu_Array[bx][MOMZ][g_index] = pz;
      g_Flu_Array[bx][ENGY][g_index] = Eint + temp*(px*px+py*py+pz*pz);


      s_index += NSlice*GRA_NXT*GRA_NXT;
      g_index += NSlice*PS1*PS1;

   } // for (uint z=tz; z<PS1; z+=NSlice)

} // FUNCTION : CUPOT_HydroGravitySolver



#endif // #if ( defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )
