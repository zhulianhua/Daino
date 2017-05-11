
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_Restrict
// Description :  Replace the data at level "FaLv" by the average data at level "FaLv+1" 
//
// Note        :  Use the input parameter "TVar" to determine the targeted variables, which can be any
//                subset of (_FLU | _POTE)
//
// Parameter   :  FaLv     : Targeted refinement level at which the data are going to be replaced
//                SonFluSg : Fluid sandglass at level "FaLv+1"
//                FaFluSg  : Fluid sandglass at level "FaLv"
//                SonPotSg : Potential sandglass at level "FaLv+1"
//                FaPotSg  : Potential sandglass at level "FaLv"
//                TVar     : Targeted variables
//                           --> Supported variables in different models:
//                               HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY, _FLU [, _POTE]
//                               MHD   : 
//                               ELBDM : _DENS, _REAL, _IMAG, _FLU [, _POTE]
//-------------------------------------------------------------------------------------------------------
void Flu_Restrict( const int FaLv, const int SonFluSg, const int FaFluSg, const int SonPotSg, const int FaPotSg,
                   const int TVar )
{

   const int SonLv = FaLv + 1;

// check   
   if ( FaLv < 0  ||  FaLv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaLv", FaLv );

   if ( FaLv == NLEVEL-1 )
   {
      Aux_Message( stderr, "WARNING : applying \"%s\" to the maximum level is meaningless !! \n", __FUNCTION__ ); 
      return;
   }

   if (  ( TVar & _FLU )  &&  ( SonFluSg != 0 && SonFluSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SonFluSg", SonFluSg );

   if (  ( TVar & _FLU )  &&  ( FaFluSg != 0 && FaFluSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaFluSg", FaFluSg );

#  ifdef GRAVITY
   if (  ( TVar & _POTE )  &&  ( SonPotSg != 0 && SonPotSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SonPotSg", SonPotSg );

   if (  ( TVar & _POTE )  &&  ( FaPotSg != 0 && FaPotSg != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaPotSg", FaPotSg );

   if (  !( TVar & (_FLU|_POTE) )  )
      Aux_Error( ERROR_INFO, "no suitable targeted variable is found --> missing (_FLU | _POTE) !!\n" );
#  else
   if (  !( TVar & _FLU )  )
      Aux_Error( ERROR_INFO, "no suitable targeted variable is found --> missing _FLU !!\n" );
#  endif

// nothing to do if there are no real patches at lv+1
   if ( patch->NPatchComma[SonLv][1] == 0 )  return;

// check the synchronization
   Mis_Check_Synchronization( Time[FaLv], Time[SonLv], __FUNCTION__, true );


   const bool ResFlu = TVar & _FLU;
#  ifdef GRAVITY
   const bool ResPot = TVar & _POTE;
#  endif
   int SonPID, FaPID, Disp_i, Disp_j, Disp_k, ii, jj, kk, I, J, K, Ip, Jp, Kp; 
   int NVar_Flu, NVar_Tot, TFluVarIdx, TFluVarIdxList[NCOMP];


// determine the components to be restricted (TFluVarIdx : targeted fluid variable indices ( = [0 ... NCOMP-1] )
   NVar_Flu = 0;

   for (int v=0; v<NCOMP; v++)
      if ( TVar & (1<<v) )    TFluVarIdxList[ NVar_Flu++ ] = v;

// check again
   NVar_Tot = NVar_Flu;
#  ifdef GRAVITY
   if ( ResPot )   NVar_Tot ++; 
#  endif
   if ( NVar_Tot == 0 )
   {
      Aux_Message( stderr, "WARNING : no targeted variable is found !!\n" );
      return;
   }


// restrict
#  pragma omp parallel for private( SonPID, FaPID, Disp_i, Disp_j, Disp_k, ii, jj, kk, I, J, K, Ip, Jp, Kp, \
                                    TFluVarIdx )
   for (int SonPID0=0; SonPID0<patch->NPatchComma[SonLv][1]; SonPID0+=8)
   {
      FaPID = patch->ptr[0][SonLv][SonPID0]->father;

//    check
#     ifdef DAINO_DEBUG
      if ( FaPID < 0 )
         Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 %d has no father patch (FaPID = %d) !!\n", 
                    SonLv, SonPID0, FaPID );

      if ( ResFlu  &&  patch->ptr[FaFluSg][FaLv][FaPID]->fluid == NULL )
         Aux_Error( ERROR_INFO, "FaFluSg %d, FaLv %d, FaPID %d has no fluid array allocated !!\n", 
                    FaFluSg, FaLv, FaPID );

#     ifdef GRAVITY
      if ( ResPot  &&  patch->ptr[FaPotSg][FaLv][FaPID]->pot == NULL )
         Aux_Error( ERROR_INFO, "FaPotSg %d, FaLv %d, FaPID %d has no potential array allocated !!\n", 
                    FaPotSg, FaLv, FaPID );
#     endif
#     endif // #ifdef DAINO_DEBUG


//    loop over eight sons
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         SonPID = SonPID0 + LocalID;
         Disp_i = TABLE_02( LocalID, 'x', 0, PATCH_SIZE/2 ); 
         Disp_j = TABLE_02( LocalID, 'y', 0, PATCH_SIZE/2 ); 
         Disp_k = TABLE_02( LocalID, 'z', 0, PATCH_SIZE/2 ); 

//       check         
#        ifdef DAINO_DEBUG
         if ( ResFlu  &&  patch->ptr[SonFluSg][SonLv][SonPID]->fluid == NULL )
            Aux_Error( ERROR_INFO, "SonFluSg %d, SonLv %d, SonPID %d has no fluid array allocated !!\n", 
                       SonFluSg, SonLv, SonPID );

#        ifdef GRAVITY
         if ( ResPot  &&  patch->ptr[SonPotSg][SonLv][SonPID]->pot == NULL )
            Aux_Error( ERROR_INFO, "SonPotSg %d, SonLv %d, SonPID %d has no potential array allocated !!\n", 
                       SonPotSg, SonLv, SonPID );
#        endif
#        endif // #ifdef DAINO_DEBUG


//       restrict the fluid data
         if ( ResFlu )
         for (int v=0; v<NVar_Flu; v++)
         {  
            TFluVarIdx = TFluVarIdxList[v]; 

            for (int k=0; k<PATCH_SIZE/2; k++)  {  K = k*2;    Kp = K+1;   kk = k + Disp_k;
            for (int j=0; j<PATCH_SIZE/2; j++)  {  J = j*2;    Jp = J+1;   jj = j + Disp_j;
            for (int i=0; i<PATCH_SIZE/2; i++)  {  I = i*2;    Ip = I+1;   ii = i + Disp_i;

               patch->ptr[FaFluSg][FaLv][FaPID]->fluid[TFluVarIdx][kk][jj][ii] 
                  = 0.125 * ( patch->ptr[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][K ][J ][I ] + 
                              patch->ptr[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][K ][J ][Ip] +
                              patch->ptr[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][K ][Jp][I ] +
                              patch->ptr[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][Kp][J ][I ] +
                              patch->ptr[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][K ][Jp][Ip] +
                              patch->ptr[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][Kp][Jp][I ] +
                              patch->ptr[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][Kp][J ][Ip] +
                              patch->ptr[SonFluSg][SonLv][SonPID]->fluid[TFluVarIdx][Kp][Jp][Ip]   );
            }}}
         }


//       restrict the potential data
#        ifdef GRAVITY
         if ( ResPot )
         {  
            for (int k=0; k<PATCH_SIZE/2; k++)  {  K = k*2;    Kp = K+1;   kk = k + Disp_k;
            for (int j=0; j<PATCH_SIZE/2; j++)  {  J = j*2;    Jp = J+1;   jj = j + Disp_j;
            for (int i=0; i<PATCH_SIZE/2; i++)  {  I = i*2;    Ip = I+1;   ii = i + Disp_i;

               patch->ptr[FaPotSg][FaLv][FaPID]->pot[kk][jj][ii] 
                  = 0.125 * ( patch->ptr[SonPotSg][SonLv][SonPID]->pot[K ][J ][I ] + 
                              patch->ptr[SonPotSg][SonLv][SonPID]->pot[K ][J ][Ip] +
                              patch->ptr[SonPotSg][SonLv][SonPID]->pot[K ][Jp][I ] +
                              patch->ptr[SonPotSg][SonLv][SonPID]->pot[Kp][J ][I ] +
                              patch->ptr[SonPotSg][SonLv][SonPID]->pot[K ][Jp][Ip] +
                              patch->ptr[SonPotSg][SonLv][SonPID]->pot[Kp][Jp][I ] +
                              patch->ptr[SonPotSg][SonLv][SonPID]->pot[Kp][J ][Ip] +
                              patch->ptr[SonPotSg][SonLv][SonPID]->pot[Kp][Jp][Ip]   );
            }}}
         }
#        endif
      } // for (int LocalID=0; LocalID<8; LocalID++)


//    rescale real and imaginary parts to get the correct density in ELBDM
#     if ( MODEL == ELBDM )
      real Real, Imag, Rho_Wrong, Rho_Corr, Rescale;

      if (  ( TVar & _DENS )  &&  ( TVar & _REAL )  &&  (TVar & _IMAG )  )
      for (int k=0; k<PATCH_SIZE; k++)
      for (int j=0; j<PATCH_SIZE; j++)
      for (int i=0; i<PATCH_SIZE; i++)
      {
         Real      = patch->ptr[FaFluSg][FaLv][FaPID]->fluid[REAL][k][j][i];
         Imag      = patch->ptr[FaFluSg][FaLv][FaPID]->fluid[IMAG][k][j][i];
         Rho_Wrong = Real*Real + Imag*Imag;
         Rho_Corr  = patch->ptr[FaFluSg][FaLv][FaPID]->fluid[DENS][k][j][i];
         Rescale   = SQRT( Rho_Corr/Rho_Wrong );

         patch->ptr[FaFluSg][FaLv][FaPID]->fluid[REAL][k][j][i] *= Rescale;
         patch->ptr[FaFluSg][FaLv][FaPID]->fluid[IMAG][k][j][i] *= Rescale;
      }
#     endif

   } // for (int SonPID0=0; SonPID0<patch->NPatchComma[SonLv][1]; SonPID0+=8)

} // FUNCTION : Flu_Restrict
