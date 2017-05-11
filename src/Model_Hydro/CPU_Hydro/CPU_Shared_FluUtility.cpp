
#include "DAINO.h"
#include "CUFLU.h"

#if (  !defined GPU  &&  MODEL == HYDRO  &&  \
       ( FLU_SCHEME == WAF || FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Rotate3D
// Description :  Rotate the input 5-element fluid variables properly to simplify the 3D calculation
//
// Note        :  x : (0,1,2,3,4) <--> (0,1,2,3,4)   
//                y : (0,1,2,3,4) <--> (0,2,3,1,4)
//                z : (0,1,2,3,4) <--> (0,3,1,2,4)
//
// Parameter   :  InOut    : Array storing both the input and output data
//                XYZ      : Targeted spatial direction : (0/1/2) --> (x/y/z)
//                Forward  : (true/false) <--> (forward/backward)
//-------------------------------------------------------------------------------------------------------
void CPU_Rotate3D( real InOut[], const int XYZ, const bool Forward )
{
   
   if ( XYZ == 0 )   return;


   real Temp[5];
   for (int v=0; v<5; v++)    Temp[v] = InOut[v];

   if ( Forward )
   {
      switch ( XYZ )
      {
         case 1 : InOut[1] = Temp[2];  InOut[2] = Temp[3];  InOut[3] = Temp[1];     break;
         case 2 : InOut[1] = Temp[3];  InOut[2] = Temp[1];  InOut[3] = Temp[2];     break;
      }
   }

   else // backward
   {
      switch ( XYZ )
      {
         case 1 : InOut[1] = Temp[3];  InOut[2] = Temp[1];  InOut[3] = Temp[2];     break;
         case 2 : InOut[1] = Temp[2];  InOut[2] = Temp[3];  InOut[3] = Temp[1];     break;
      }
   }

} // FUNCTION : CPU_Rotate3D



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Con2Pri 
// Description :  Convert the conserved variables to the primitive variables 
//
// Parameter   :  In       : Array storing the input conserved variables
//                Out      : Array to store the output primitive variables
//                Gamma_m1 : Gamma - 1
//-------------------------------------------------------------------------------------------------------
void CPU_Con2Pri( const real In[], real Out[], const real Gamma_m1 )
{

   const real _Rho = (real)1.0 / In[0];
   
   Out[0] = In[0];
   Out[1] = In[1]*_Rho;
   Out[2] = In[2]*_Rho;
   Out[3] = In[3]*_Rho;
   Out[4] = (  In[4] - (real)0.5*Out[0]*( Out[1]*Out[1] + Out[2]*Out[2] + Out[3]*Out[3] )  )*Gamma_m1;

#  ifdef ENFORCE_POSITIVE
   Out[4] = FMAX( Out[4], MIN_VALUE );
#  endif

} // FUNCTION : CPU_Con2Pri



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Pri2Con 
// Description :  Convert the primitive variables to the conserved variables 
//
// Parameter   :  In       : Array storing the input primitive variables
//                Out      : Array to store the output conserved variables
//                Gamma_m1 : 1 / (Gamma - 1)
//-------------------------------------------------------------------------------------------------------
void CPU_Pri2Con( const real In[], real Out[], const real _Gamma_m1 )
{

   Out[0] = In[0];
   Out[1] = In[0]*In[1];
   Out[2] = In[0]*In[2];
   Out[3] = In[0]*In[3];
   Out[4] = In[4]*_Gamma_m1 + (real)0.5*In[0]*( In[1]*In[1] + In[2]*In[2] + In[3]*In[3] );

} // FUNCTION : CPU_Pri2Con



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Con2Flux
// Description :  Evaluate the hydrodynamic fluxes by the input conserved variables
//
// Parameter   :  XYZ   : Targeted spatial direction : (0/1/2) --> (x/y/z) 
//                Flux  : Array to store the output fluxes
//                Input : Array storing the input conserved variables
//                Gamma : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
void CPU_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma )
{

   real Var[5];
   real Pres, _Rho, Vx;

   for (int v=0; v<5; v++)    Var[v] = Input[v];

   CPU_Rotate3D( Var, XYZ, true );

   _Rho = (real)1.0 / Var[0];
   Pres = (Gamma-(real)1.0) * (  Var[4] - (real)0.5*( Var[1]*Var[1] + Var[2]*Var[2] + Var[3]*Var[3] )*_Rho  );
   Vx   = _Rho*Var[1];

   Flux[0] = Var[1];
   Flux[1] = Vx*Var[1] + Pres;
   Flux[2] = Vx*Var[2];
   Flux[3] = Vx*Var[3];
   Flux[4] = Vx*( Var[4] + Pres );

   CPU_Rotate3D( Flux, XYZ, false );

} // FUNCTION : CPU_Con2Flux



#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  (FLU_SCHEME == WAF || MHM || MHM_RP || CTU) )
