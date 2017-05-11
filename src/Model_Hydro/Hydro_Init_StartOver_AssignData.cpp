
#include "DAINO.h"

#if ( MODEL == HYDRO )

static void Init_Function_User( real fluid[], const real x, const real y, const real z, const double Time );
void (*Init_Function_Ptr)( real fluid[], const real x, const real y, const real z, const double Time ) = Init_Function_User;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Function_User 
// Description :  Function to initialize the fluid field 
//
// Note        :  Invoked by "Hydro_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Init_Function_User( real fluid[], const real x, const real y, const real z, const double Time )
{

   const real Gamma2  = real( 1.0/GAMMA/(GAMMA-1.0) );
   const real C1[3]   = { 0.5*patch->BoxSize[0], 
                          0.5*patch->BoxSize[1], 
                          0.5*patch->BoxSize[2] };

   const real Cs      =  1.0;
   const real Height1 = 10.0;
   const real Width1  =  0.20;


   fluid[DENS] = 1.0 + Height1*exp(  -( SQR(x-C1[0])+ SQR(y-C1[1]) + SQR(z-C1[2]) ) / SQR(Width1)  );
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = Cs*Cs*fluid[DENS]*Gamma2 + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

} // FUNCTION : Init_Function_User



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Init_StartOver_AssignData
// Description :  Construct the initial condition in HYDRO
//
// Note        :  1. Work for the option "OPT__INIT == INIT_STARTOVER"
//                2. The initialization function should be specified in "Init_Function_Ptr", which is a function
//                   pointer pointing to either "Init_Function_User" or the test problem specified function
//                   (e.g., Hydro_TestProbSol_Riemann) set in "Init_TestProb"
//
// Parameter   :  lv : Targeted refinement level
//-------------------------------------------------------------------------------------------------------
void Hydro_Init_StartOver_AssignData( const int lv )
{

   const real scale  = (real)patch->scale[lv];
   const real dh     = patch->dh[lv];

   real fluid[NCOMP], fluid_sub[NCOMP];
   real x, y, z, x0, y0, z0;

   for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
   for (int k=0; k<PS1; k++)  {  z = ( patch->ptr[0][lv][PID]->corner[2]/scale + k + 0.5 )*dh;
   for (int j=0; j<PS1; j++)  {  y = ( patch->ptr[0][lv][PID]->corner[1]/scale + j + 0.5 )*dh;
   for (int i=0; i<PS1; i++)  {  x = ( patch->ptr[0][lv][PID]->corner[0]/scale + i + 0.5 )*dh;

      Init_Function_Ptr( fluid, x, y, z, Time[lv] );

      for (int v=0; v<NCOMP; v++)   patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[v][k][j][i] = fluid[v];

   }}}

} // FUNCTION : Hydro_Init_StartOver_AssignData



#endif // #if ( MODEL == HYDRO )
