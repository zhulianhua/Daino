
#include "DAINO.h"

#if ( !defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )




//-----------------------------------------------------------------------------------------
// Function    :  CPU_HydroGravitySolver
// Description :  Use CPU to advance the momentum and energy density by gravitational acceleration
//
// Parameter   :  Flu_Array   : Array to store the input and output fluid variables
//                Pot_Array   : Array storing the input potential for evaluating the gravitational acceleration
//                NPatchGroup : Number of patch groups to be evaluated
//                dt          : Time interval to advance solution
//                dh          : Grid size
//                P5_Gradient : Use 5-points stencil to evaluate the potential gradient
//-----------------------------------------------------------------------------------------
void CPU_HydroGravitySolver(       real Flu_Array[][NCOMP][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE], 
                             const real Pot_Array[][GRA_NXT][GRA_NXT][GRA_NXT],
                             const int NPatchGroup, const real dt, const real dh, const bool P5_Gradient )
{ 

   const int NPatch     = NPatchGroup*8;
   const real Gra_Const = ( P5_Gradient ) ? dt/(12.0*dh ) : dt/(2.0*dh );
#  if ( GRA_GHOST_SIZE == 2 )       
   const real Const_8   = (real)8.0;
#  endif

   real Acc[3], Eint, px, py, pz, rho, temp;
   int ii, jj, kk;


// loop over all patches
#  pragma omp parallel for private( Acc, Eint, px, py, pz, rho, temp, ii, jj, kk )
   for (int P=0; P<NPatch; P++)
   {
      for (int k=GRA_GHOST_SIZE; k<GRA_NXT-GRA_GHOST_SIZE; k++)   {  kk = k - GRA_GHOST_SIZE;
      for (int j=GRA_GHOST_SIZE; j<GRA_NXT-GRA_GHOST_SIZE; j++)   {  jj = j - GRA_GHOST_SIZE;
      for (int i=GRA_GHOST_SIZE; i<GRA_NXT-GRA_GHOST_SIZE; i++)   {  ii = i - GRA_GHOST_SIZE;

//       evaluate the gravitational acceleration
#        if ( GRA_GHOST_SIZE == 2 )       
         if ( P5_Gradient )
         {
            Acc[0] = Gra_Const * ( -         Pot_Array[P][k  ][j  ][i+2] +         Pot_Array[P][k  ][j  ][i-2]
                                   + Const_8*Pot_Array[P][k  ][j  ][i+1] - Const_8*Pot_Array[P][k  ][j  ][i-1] );
            Acc[1] = Gra_Const * ( -         Pot_Array[P][k  ][j+2][i  ] +         Pot_Array[P][k  ][j-2][i  ]
                                   + Const_8*Pot_Array[P][k  ][j+1][i  ] - Const_8*Pot_Array[P][k  ][j-1][i  ] );
            Acc[2] = Gra_Const * ( -         Pot_Array[P][k+2][j  ][i  ] +         Pot_Array[P][k-2][j  ][i  ]
                                   + Const_8*Pot_Array[P][k+1][j  ][i  ] - Const_8*Pot_Array[P][k-1][j  ][i  ] );
         }

         else
#        endif
         {
            Acc[0] = Gra_Const * ( Pot_Array[P][k  ][j  ][i+1] - Pot_Array[P][k  ][j  ][i-1] );
            Acc[1] = Gra_Const * ( Pot_Array[P][k  ][j+1][i  ] - Pot_Array[P][k  ][j-1][i  ] );
            Acc[2] = Gra_Const * ( Pot_Array[P][k+1][j  ][i  ] - Pot_Array[P][k-1][j  ][i  ] );
         }


//       advance fluid
         rho  = Flu_Array[P][DENS][kk][jj][ii];
         px   = Flu_Array[P][MOMX][kk][jj][ii];
         py   = Flu_Array[P][MOMY][kk][jj][ii];
         pz   = Flu_Array[P][MOMZ][kk][jj][ii];
         temp = (real)0.5/rho;
         Eint = Flu_Array[P][ENGY][kk][jj][ii] - temp*(px*px+py*py+pz*pz);

         px -= rho*Acc[0];
         py -= rho*Acc[1];
         pz -= rho*Acc[2];

         Flu_Array[P][MOMX][kk][jj][ii] = px;
         Flu_Array[P][MOMY][kk][jj][ii] = py;
         Flu_Array[P][MOMZ][kk][jj][ii] = pz;
         Flu_Array[P][ENGY][kk][jj][ii] = Eint + temp*(px*px+py*py+pz*pz);

      }}} // i,j,k
   } // for (int P=0; P<NPatch; P++)

} // FUNCTION : CPU_HydroGravitySolver



#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )
