
#ifndef __CUFLU_FLUUTILITY_CU__
#define __CUFLU_FLUUTILITY_CU__



#include "Macro.h"
#include "CUFLU.h"

static __device__ FluVar CUFLU_Pri2Con( const FluVar Pri, const real _Gamma_m1 );
static __device__ FluVar CUFLU_Con2Pri( const FluVar Con, const real Gamma_m1 );
static __device__ FluVar CUFLU_Con2Flux( const FluVar Input, const real Gamma_m1, const int XYZ );
static __device__ FluVar CUFLU_Rotate3D( const FluVar In, const int XYZ, const bool Forward );
static __device__ void CUFLU_Con2Pri_AllGrids( const real g_Fluid_In[][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                               real g_PriVar[][5][ FLU_NXT*FLU_NXT*FLU_NXT ], const real Gamma );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Rotate3D
// Description :  Rotate the input 5-element fluid variables properly to simplify the 3D calculation
//
// Note        :  x : (0,1,2,3,4) <--> (0,1,2,3,4)   
//                y : (0,1,2,3,4) <--> (0,2,3,1,4)
//                z : (0,1,2,3,4) <--> (0,3,1,2,4)
//
// Parameter   :  In       : Input variables to be rotated
//                XYZ      : Targeted spatial direction : (0/1/2) --> (x/y/z)
//                Forward  : (true/false) <--> (forward/backward)
//-------------------------------------------------------------------------------------------------------
__device__ FluVar CUFLU_Rotate3D( const FluVar In, const int XYZ, const bool Forward )
{

   if ( XYZ == 0 )   return In;


   FluVar Out;

   if ( Forward )
   {
      switch ( XYZ )
      {
         case 1 : Out.Rho=In.Rho;  Out.Px=In.Py;  Out.Py=In.Pz;  Out.Pz=In.Px;  Out.Egy=In.Egy;    break;
         case 2 : Out.Rho=In.Rho;  Out.Px=In.Pz;  Out.Py=In.Px;  Out.Pz=In.Py;  Out.Egy=In.Egy;    break;
      }
   }

   else
   {
      switch ( XYZ )
      {
         case 1 : Out.Rho=In.Rho;  Out.Px=In.Pz;  Out.Py=In.Px;  Out.Pz=In.Py;  Out.Egy=In.Egy;    break;
         case 2 : Out.Rho=In.Rho;  Out.Px=In.Py;  Out.Py=In.Pz;  Out.Pz=In.Px;  Out.Egy=In.Egy;    break;
      }
   }

   return Out;

} // FUNCTION : CUFLU_Rotate3D



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Con2Flux
// Description :  Conserved variables --> fluxes
//
// Parameter   :  Input    : Input conserved variables
//                Gamma_m1 : Gamma - 1
//                XYZ      : Targeted spatial direction : (0/1/2) --> (x/y/z)
//-------------------------------------------------------------------------------------------------------
__device__ FluVar CUFLU_Con2Flux( const FluVar Input, const real Gamma_m1, const int XYZ )
{

   FluVar Var, TempFlux, Output; 
   real Pres, _Rho, Vx;

   Var = CUFLU_Rotate3D( Input, XYZ, true );

   _Rho = (real)1.0 / Var.Rho;
   Pres = Gamma_m1 * (  Var.Egy - (real)0.5*( Var.Px*Var.Px + Var.Py*Var.Py + Var.Pz*Var.Pz )*_Rho  );
   Vx   = _Rho*Var.Px;

   TempFlux.Rho = Var.Px;
   TempFlux.Px  = Vx*Var.Px + Pres;
   TempFlux.Py  = Vx*Var.Py;
   TempFlux.Pz  = Vx*Var.Pz;
   TempFlux.Egy = Vx*( Var.Egy + Pres );

   Output = CUFLU_Rotate3D( TempFlux, XYZ, false );

   return Output;

} // FUNCTION : CUFLU_Con2Flux



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Pri2Con
// Description :  Primitive variables --> conserved variables
//
// Parameter   :  Pri         : Input primitive variables
//                _Gamma_m1   : 1.0 / (Gamma-1.0) 
//-------------------------------------------------------------------------------------------------------
__device__ FluVar CUFLU_Pri2Con( const FluVar Pri, const real _Gamma_m1 )
{

   FluVar Con;

   Con.Rho = Pri.Rho;
   Con.Egy = _Gamma_m1*Pri.Egy + (real)0.5*Pri.Rho*( Pri.Px*Pri.Px + Pri.Py*Pri.Py + Pri.Pz*Pri.Pz );
   Con.Px  = Pri.Rho*Pri.Px;
   Con.Py  = Pri.Rho*Pri.Py;
   Con.Pz  = Pri.Rho*Pri.Pz;

   return Con;

} // FUNCTION : CUFLU_Pri2Con



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Con2Pri
// Description :  Conserved variables --> primitive variables 
//
// Parameter   :  Con      : Input conserved variables
//                Gamma_m1 : Gamma - 1.0
//-------------------------------------------------------------------------------------------------------
__device__ FluVar CUFLU_Con2Pri( const FluVar Con, const real Gamma_m1 )
{

   const real _Rho = (real)1.0/Con.Rho;
   FluVar Pri;

   Pri.Rho = Con.Rho;
   Pri.Px  = Con.Px*_Rho;
   Pri.Py  = Con.Py*_Rho;
   Pri.Pz  = Con.Pz*_Rho;
   Pri.Egy = Gamma_m1*(  Con.Egy - (real)0.5*( Con.Px*Con.Px + Con.Py*Con.Py + Con.Pz*Con.Pz )*_Rho  );

#  ifdef ENFORCE_POSITIVE
   Pri.Egy = FMAX( Pri.Egy, MIN_VALUE );
#  endif

   return Pri;

} // FUNCTION : CUFLU_Con2Pri



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Con2Pri_AllGrids
// Description :  Conserved variables --> primitive variables for all grids
//
// Parameter   :  g_Fluid_In  : Global memory array storing the input fluid variables
//                g_PriVar    : Global memory array to store the output primitive variables
//                Gamma       : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_Con2Pri_AllGrids( const real g_Fluid_In[][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                        real g_PriVar[][5][ FLU_NXT*FLU_NXT*FLU_NXT ], const real Gamma )
{

   const uint bx       = blockIdx.x;
   const uint tx       = threadIdx.x;
   const uint dID      = blockDim.x;
   const real Gamma_m1 = Gamma - (real)1.0;

   uint   ID = tx;
   FluVar Var;

// loop over all cells
   while ( ID < FLU_NXT*FLU_NXT*FLU_NXT )
   {
      Var.Rho = g_Fluid_In[bx][0][ID];
      Var.Px  = g_Fluid_In[bx][1][ID];
      Var.Py  = g_Fluid_In[bx][2][ID];
      Var.Pz  = g_Fluid_In[bx][3][ID];
      Var.Egy = g_Fluid_In[bx][4][ID];

      Var = CUFLU_Con2Pri( Var, Gamma_m1 );

      g_PriVar[bx][0][ID] = Var.Rho;
      g_PriVar[bx][1][ID] = Var.Px;
      g_PriVar[bx][2][ID] = Var.Py;
      g_PriVar[bx][3][ID] = Var.Pz;
      g_PriVar[bx][4][ID] = Var.Egy;

      ID += dID;
   }

} // FUNCTION : CUFLU_Con2Pri_AllGrids



#endif // #ifndef __CUFLU_FLUUTILITY_CU__
