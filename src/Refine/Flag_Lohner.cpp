
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Lohner
// Description :  Check if the numerical error at the cell (i,j,k) estimated by Lohner's prescription
//                exceeds the given threshold
//
// Note        :  1. Invoked by the function "Flag_Check" 
//                2. Adopt the modified version in FLASH4 and MPI_AMRVAC
//
// Parameter   :  i,j,k       : Indices of the targeted cell in the array "Var1D"
//             :  Var1D       : Array storing the input variables for the Lohner error estimator
//                Slope1D     : Array storing the input slopes of Lohner_Var for the Lohner error estimator
//                NCell       : Size of the arrays Lohner_Var along one direction
//                NVar        : Number of variables stored in Lohner_Var and Lohner_Slope (HYDRO=1,ELBDM=2)
//                Threshold   : Refinement threshold
//                Filter      : Filter parameter for preventing refinement of small ripples
//                Soften      : Minimum number in the denominator --> error = sqrt( N/max(D,Soften) ), where
//                              N and D are numerator and denominator in the Lohner's formula, respectively
//
// Return      :  "true"  if the estimated error is larger           than the given threshold
//                "false" if the estimated error is equal or smaller than the given threshold
//-------------------------------------------------------------------------------------------------------
bool Flag_Lohner( const int i, const int j, const int k, const real *Var1D, const real *Slope1D, const int NCell, 
                  const int NVar, const double Threshold, const double Filter, const double Soften )
{

// check
#  ifdef DAINO_DEBUG
   if ( NCell != PS1 + 4 )    Aux_Error( ERROR_INFO, "NCell (%d) != %d !!\n", NCell, PS1+4 );

   if (  i < 2  ||  i >= PS1+2  ||  j < 2  ||  j >= PS1+2  ||  k < 2  ||  k >= PS1+2  )
      Aux_Error( ERROR_INFO, "incorrect index (i,j,k) = (%d,%d,%d) !!\n", i, j, k );

#  if   ( MODEL == HYDRO  ||  MODEL == MHD )
   if ( NVar != 1 )  Aux_Error( ERROR_INFO, "NVar (%d) != 1 !!\n", NVar );
#  elif ( MODEL == ELBDM )
   if ( NVar != 2 )  Aux_Error( ERROR_INFO, "NVar (%d) != 2 !!\n", NVar );
#  else
#  error : unsupported MODEL !!
#  endif // MODEL
#  endif // #ifdef DAINO_DEBUG


   const int NSlope = PS1 + 2;   // size of the slope array

// i,j,k: for the array "Var"
   const int ip1 = i + 1;   const int jp1 = j + 1;   const int kp1 = k + 1;
   const int im1 = i - 1;   const int jm1 = j - 1;   const int km1 = k - 1;
   const int ip2 = i + 2;   const int jp2 = j + 2;   const int kp2 = k + 2;
   const int im2 = i - 2;   const int jm2 = j - 2;   const int km2 = k - 2;

// ii,jj,kk: for the array "Slope"
   const int ii  =  i - 1;   const int jj  =  j - 1;   const int kk  =  k - 1;
   const int iip = ii + 1;   const int jjp = jj + 1;   const int kkp = kk + 1;
   const int iim = ii - 1;   const int jjm = jj - 1;   const int kkm = kk - 1;

   real Der2_xx, Der2_yy, Der2_zz, Der2_xy, Der2_yz, Der2_zx;              // grad X grad( Var) --> tensor
   real Der1_x, Der1_y, Der1_z;                                            // grad( Var )       --> vector
   real Filter_xx, Filter_yy, Filter_zz, Filter_xy, Filter_yz, Filter_zx;  // filters along different directions
   real Nume=0.0, Deno=0.0, Error;
   bool Flag;

// convert the 1D arrays
   real (*Var)     [NCell ][NCell ][NCell ] = ( real(*)   [NCell ][NCell ][NCell ] )  Var1D;
   real (*Slope)[3][NSlope][NSlope][NSlope] = ( real(*)[3][NSlope][NSlope][NSlope] )Slope1D;


   for (int v=0; v<NVar; v++)
   {
//    numerator
      Der2_xx = Slope[v][0][kk ][jj ][iip] - Slope[v][0][kk ][jj ][iim];
      Der2_yy = Slope[v][1][kk ][jjp][ii ] - Slope[v][1][kk ][jjm][ii ];
      Der2_zz = Slope[v][2][kkp][jj ][ii ] - Slope[v][2][kkm][jj ][ii ];
      Der2_xy = Slope[v][1][kk ][jj ][iip] - Slope[v][1][kk ][jj ][iim];
      Der2_yz = Slope[v][2][kk ][jjp][ii ] - Slope[v][2][kk ][jjm][ii ];
      Der2_zx = Slope[v][0][kkp][jj ][ii ] - Slope[v][0][kkm][jj ][ii ];

      Nume += SQR(Der2_xx) + SQR(Der2_yy) + SQR(Der2_zz) + (real)2.0*( SQR(Der2_xy) + SQR(Der2_yz) + SQR(Der2_zx) );


//    denominator
      Der1_x  = FABS( Slope[v][0][kk ][jj ][iip] ) + FABS( Slope[v][0][kk ][jj ][iim] );
      Der1_y  = FABS( Slope[v][1][kk ][jjp][ii ] ) + FABS( Slope[v][1][kk ][jjm][ii ] );
      Der1_z  = FABS( Slope[v][2][kkp][jj ][ii ] ) + FABS( Slope[v][2][kkm][jj ][ii ] );

      Filter_xx = Filter*(           FABS( Var[v][k  ][j  ][ip2] ) + 
                           (real)2.0*FABS( Var[v][k  ][j  ][i  ] ) + 
                                     FABS( Var[v][k  ][j  ][im2] )  );
      Filter_yy = Filter*(           FABS( Var[v][k  ][jp2][i  ] ) + 
                           (real)2.0*FABS( Var[v][k  ][j  ][i  ] ) + 
                                     FABS( Var[v][k  ][jm2][i  ] )  );
      Filter_zz = Filter*(           FABS( Var[v][kp2][j  ][i  ] ) + 
                           (real)2.0*FABS( Var[v][k  ][j  ][i  ] ) + 
                                     FABS( Var[v][km2][j  ][i  ] )  );

      Filter_xy = Filter*(  FABS( Var[v][k  ][jp1][ip1] ) + FABS( Var[v][k  ][jp1][im1] ) + 
                            FABS( Var[v][k  ][jm1][ip1] ) + FABS( Var[v][k  ][jm1][im1] )  );
      Filter_yz = Filter*(  FABS( Var[v][kp1][jp1][i  ] ) + FABS( Var[v][kp1][jm1][i  ] ) + 
                            FABS( Var[v][km1][jp1][i  ] ) + FABS( Var[v][km1][jm1][i  ] )  );
      Filter_zx = Filter*(  FABS( Var[v][kp1][j  ][ip1] ) + FABS( Var[v][kp1][j  ][im1] ) + 
                            FABS( Var[v][km1][j  ][ip1] ) + FABS( Var[v][km1][j  ][im1] )  );

      Deno += SQR( Der1_x + Filter_xx ) + SQR( Der1_x + Filter_xy ) + SQR( Der1_x + Filter_zx ) + 
              SQR( Der1_y + Filter_xy ) + SQR( Der1_y + Filter_yy ) + SQR( Der1_y + Filter_yz ) + 
              SQR( Der1_z + Filter_zx ) + SQR( Der1_z + Filter_yz ) + SQR( Der1_z + Filter_zz ); 

   } // for (int v=0; v<NVar; v++)

// check the flag
   Error = SQRT( Nume / MAX(Deno, Soften) );
   Flag  = Error > Threshold;

   return Flag;
 
} // FUNCTION : Flag_Lohner



//-------------------------------------------------------------------------------------------------------
// Function    :  GetSlope_for_Lohner 
// Description :  Evaluate slopes along x/y/z for the Lohner error estimator 
//
// Note        :  1. This function is called in "Flag_Real" before looping over all cells in the patch in order to 
//                   achieve higher performance
//                2. Evaluate slope by the discrete central difference: slope_x(i,j,k) = var(i+1,j,k) - var(i-1,j,k)
//                3. Do not take into account the physical size of each cell since the Lohner error estimator 
//                   is dimensionless
//
// Parameter   :  Var1D   : Array storing the input variables for the Lohner error estimator
//                Slope1D : Array to store the output slopes of Lohner_Var for the Lohner error estimator
//                NCell   : Size of the arrays Lohner_Var along one direction
//                NVar    : Number of variables stored in Lohner_Var and Lohner_Slope
//
// Return      :  None 
//-------------------------------------------------------------------------------------------------------
void GetSlope_for_Lohner( const real *Var1D, real *Slope1D, const int NCell, const int NVar )
{

// check
#  ifdef DAINO_DEBUG
   if ( NCell != PS1 + 4 )    Aux_Error( ERROR_INFO, "NCell (%d) != %d !!\n", NCell, PS1+4 );
#  endif


   const int NSlope = PS1 + 2;   // size of the slope array
   int ii, jj, kk, iim, jjm, kkm, iip, jjp, kkp;

// convert the 1D arrays
   real (*Var)     [NCell ][NCell ][NCell ] = ( real(*)   [NCell ][NCell ][NCell ] )  Var1D;
   real (*Slope)[3][NSlope][NSlope][NSlope] = ( real(*)[3][NSlope][NSlope][NSlope] )Slope1D;


   for (int v=0; v<NVar; v++)    {
   for (int k=0; k<NSlope; k++)  {  kk = k + 1;   kkp = kk + 1;   kkm = kk - 1;
   for (int j=0; j<NSlope; j++)  {  jj = j + 1;   jjp = jj + 1;   jjm = jj - 1;
   for (int i=0; i<NSlope; i++)  {  ii = i + 1;   iip = ii + 1;   iim = ii - 1;

      Slope[v][0][k][j][i] = Var[v][kk ][jj ][iip] - Var[v][kk ][jj ][iim];
      Slope[v][1][k][j][i] = Var[v][kk ][jjp][ii ] - Var[v][kk ][jjm][ii ];
      Slope[v][2][k][j][i] = Var[v][kkp][jj ][ii ] - Var[v][kkm][jj ][ii ];

   }}}}

} // FUNCTION : GetSlope_for_Lohner

