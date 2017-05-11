
#include "DAINO.h"

void Flag_Grandson( const int lv, const int PID, const int LocalID );
void GetSlope_for_Lohner( const real *Var1D, real *Slope1D, const int NCell, const int NVar );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Real
// Description :  Flag the real patches at level "lv" according to the given refinement criteria
//
// Note        :  1. Flag operation of the buffer patches is performed by the function "Flag_Buffer"
//                2. In this function, the buffer patches may still be flagged due to the FLAG_BUFFER_SIZE
//                   extension and the grandson check
//                3. We must assume that the refinement criteria do NOT depend on data in the sibling patches 
//                   (even if we adopt the gradient(XXX) as the criteria)
//                   --> Otherwise, we must rewrite all procedures to ensure that the buffer data are already
//                       assigned or received BEFORE the flag operation
//                   --> Related function : Init_UM, Integration_IndiviTimeStep, Integration_SharedTimeStep,
//                                          Init_Start
//                4. To add new refinement criteria, please edit the function "Flag_Check"
//                5. Definition of the function "GetSlope_for_Lohner" is put in the file "Flag_Lohner"
//
// Parameter   :  lv          : Targeted refinement level to be flagged
//                UseLBFunc   : Use the load-balance alternative functions for the grandson check and exchanging
//                              the buffer flags (useless if LOAD_BALANCE is off)
//                              --> USELB_YES : use the load-balance alternative functions
//                                  USELB_NO  : do not use the load-balance alternative functions
//-------------------------------------------------------------------------------------------------------
void Flag_Real( const int lv, const UseLBFunc_t UseLBFunc )
{

// check
   if ( lv == NLEVEL-1 )
      Aux_Error( ERROR_INFO, "function <%s> should NOT be applied to the finest level\" !!\n", __FUNCTION__ );


// initialize all flags as false
#  pragma omp parallel for
   for (int PID=0; PID<patch->num[lv]; PID++)   patch->ptr[0][lv][PID]->flag = false;


// set sibling IDs
   const int SibID_Array[3][3][3]     = {  { {18, 10, 19}, {14,   4, 16}, {20, 11, 21} }, 
                                           { { 6,  2,  7}, { 0, 999,  1}, { 8,  3,  9} }, 
                                           { {22, 12, 23}, {15,   5, 17}, {24, 13, 25} }  };
   const bool IntPhase_No             = false;                  // for invoking "Prepare_PatchGroupData"
   const int  NPG                     = 1;                      // for invoking "Prepare_PatchGroupData"
   const int  Lohner_NGhost           = 2;                      // number of ghost cells for the Lohner error estimator
   const int  Lohner_NCell            = PS1 + 2*Lohner_NGhost;  // size of the variable array for Lohner
   const int  Lohner_NSlope           = Lohner_NCell - 2;       // size of the slope array for Lohner
   const IntScheme_t Lohner_IntScheme = INT_MINMOD;             // interpolation scheme for Lohner

#  if   ( MODEL == HYDRO )
   const int  Lohner_NVar             = 1;                      // number of variables for Lohner (HYDRO)
   const int  Lohner_TVar             = ( OPT__FLAG_LOHNER == 1 ) ? _DENS : _ENGY; // target variable for Lohner (HYDRO)

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   const int  Lohner_NVar             = 2;                      // number of variables for Lohner (ELBDM)
   const int  Lohner_TVar             = _REAL | _IMAG;          // target variable for Lohner (ELBDM)

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

   const int  Lohner_Stride           = Lohner_NVar*Lohner_NCell*Lohner_NCell*Lohner_NCell;  // stride of array for one patch
      

#  pragma omp parallel
   {
      const real (*Fluid)[PS1][PS1][PS1] = NULL;
      real (*Pres)[PS1][PS1]             = NULL;
      real (*Pot )[PS1][PS1]             = NULL;
      real (*Lohner_Var)                 = NULL;   // array storing the variables for Lohner
      real (*Lohner_Slope)               = NULL;   // array storing the slopes of Lohner_Var for Lohner
      int  i_start, i_end, j_start, j_end, k_start, k_end, SibID, PID;
      bool ProperNesting, NextPatch;

#     if   ( MODEL == HYDRO )
      if ( OPT__FLAG_PRES_GRADIENT )   Pres = new real [PS1][PS1][PS1];
#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!
#     endif // MODEL

      if ( OPT__FLAG_LOHNER )    
      { 
         Lohner_Var   = new real [ 8*Lohner_NVar*Lohner_NCell *Lohner_NCell *Lohner_NCell  ]; // 8:# of local patches
         Lohner_Slope = new real [ 3*Lohner_NVar*Lohner_NSlope*Lohner_NSlope*Lohner_NSlope ]; // 3: X/Y/Z of 1 patch
      }


//    loop over all REAL patches (the buffer patches will be flagged only due to the FLAG_BUFFER_SIZE
//    extension or the grandson check )
#     pragma omp for
      for (int PID0=0; PID0<patch->NPatchComma[lv][1]; PID0+=8)
      {
//       prepare the ghost-zone data for Lohner
         if ( OPT__FLAG_LOHNER )    
            Prepare_PatchGroupData( lv, Time[lv], Lohner_Var, Lohner_NGhost, NPG, &PID0, Lohner_TVar, 
                                    Lohner_IntScheme, UNIT_PATCH, NSIDE_26, IntPhase_No );


//       loop over all local patches within the same patch group
         for (int LocalID=0; LocalID<8; LocalID++)
         {
            PID = PID0 + LocalID;

            ProperNesting = true;

            for (int sib=0; sib<26; sib++)
            {
               if ( patch->ptr[0][lv][PID]->sibling[sib] == -1 )
               {
                  ProperNesting = false;
                  break;
               }
            }


//          do flag check only if 26 siblings all exist (proper-nesting constraint)
            if ( ProperNesting )
            {
               NextPatch = false;
               Fluid     = patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid;
#              ifdef GRAVITY
               Pot       = patch->ptr[ patch->PotSg[lv] ][lv][PID]->pot;
#              endif


//             evaluate pressure
#              if   ( MODEL == HYDRO )
               if ( OPT__FLAG_PRES_GRADIENT )
               {
                  const real Gamma_m1 = GAMMA - (real)1.0;
                  real (*FluData)[PS1][PS1][PS1] = patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid;
                  real Ek; 

                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
                     Ek = (real)0.5*( FluData[MOMX][k][j][i]*FluData[MOMX][k][j][i] + 
                                      FluData[MOMY][k][j][i]*FluData[MOMY][k][j][i] +
                                      FluData[MOMZ][k][j][i]*FluData[MOMZ][k][j][i] ) / FluData[DENS][k][j][i];

                     Pres[k][j][i] = Gamma_m1 * ( FluData[ENGY][k][j][i] - Ek );
                  }
               }

#              elif ( MODEL == MHD )
#              warning : WAIT MHD !!!
#              endif // MODEL


//             evaluate the slopes along x/y/z for Lohner
               if ( OPT__FLAG_LOHNER )    GetSlope_for_Lohner( Lohner_Var+LocalID*Lohner_Stride, Lohner_Slope, Lohner_NCell, 
                                                               Lohner_NVar );


//             loop over all cells within the target patch
               for (int k=0; k<PS1; k++)  {  if ( NextPatch )  break;
                                             k_start = ( k - FLAG_BUFFER_SIZE < 0    ) ? 0 : 1;
                                             k_end   = ( k + FLAG_BUFFER_SIZE >= PS1 ) ? 2 : 1;

               for (int j=0; j<PS1; j++)  {  if ( NextPatch )  break; 
                                             j_start = ( j - FLAG_BUFFER_SIZE < 0    ) ? 0 : 1;
                                             j_end   = ( j + FLAG_BUFFER_SIZE >= PS1 ) ? 2 : 1;

               for (int i=0; i<PS1; i++)  {  if ( NextPatch )  break;
                                             i_start = ( i - FLAG_BUFFER_SIZE < 0    ) ? 0 : 1;
                                             i_end   = ( i + FLAG_BUFFER_SIZE >= PS1 ) ? 2 : 1;

//                check if the targeted cell satisfies the refinement criteria (useless pointers are always == NULL)
                  if (  lv < MAX_LEVEL  &&  Flag_Check( lv, PID, i, j, k, Fluid, Pot, Pres, Lohner_Var+LocalID*Lohner_Stride, 
                                                        Lohner_Slope, Lohner_NCell, Lohner_NVar )  )
                  {
//                   flag itself
                     patch->ptr[0][lv][PID]->flag = true;

//                   flag sibling patches according to the size of FLAG_BUFFER_SIZE
                     for (int kk=k_start; kk<=k_end; kk++)
                     for (int jj=j_start; jj<=j_end; jj++)
                     for (int ii=i_start; ii<=i_end; ii++)
                     {
                        SibID = SibID_Array[kk][jj][ii];

                        if ( SibID != 999 )  
                           patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[SibID] ]->flag = true;
                     }

//                   for FLAG_BUFFER_SIZE == PATCH_SIZE, once an cell is flagged, all 26 siblings will be flagged
                     if ( FLAG_BUFFER_SIZE == PS1 )   NextPatch = true;
                  }

               }}} // k, j, i
            } // if ( ProperNesting )
         } // for (int LocalID=0; LocalID<8; LocalID++)
      } // for (int PID0=0; PID0<patch->NPatchComma[lv][1]; PID0+=8)


#     if   ( MODEL == HYDRO )
      if ( OPT__FLAG_PRES_GRADIENT )   delete [] Pres;
#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!
#     endif // MODEL

      if ( OPT__FLAG_LOHNER )    
      {
         delete [] Lohner_Var;
         delete [] Lohner_Slope;
      }

   } // OpenMP parallel region


// apply the proper-nesting constraint (should also apply to the buffer patches)
#  pragma omp parallel for
   for (int PID=0; PID<patch->num[lv]; PID++)
   {
      for (int sib=0; sib<26; sib++)
      {
         if ( patch->ptr[0][lv][PID]->sibling[sib] == -1 )
         {
            patch->ptr[0][lv][PID]->flag = false;
            break;
         }
      }
   }


// invoke the load-balance functions 
#  ifdef LOAD_BALANCE
   if ( UseLBFunc == USELB_YES )
   {
//    grandson check
      if ( lv < NLEVEL-2 )    LB_GrandsonCheck( lv );

//    exchange the flagged buffer patches 
      LB_ExchangeFlaggedBuffer( lv );    
    
      return;
   }
#  endif


   int SonPID;

// grandson check
   if ( lv < NLEVEL-2 )
   {
//###ISSUE: use atomic ??      
//#     pragma omp parallel for private( SonPID )
      for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
      {
         SonPID = patch->ptr[0][lv][PID]->son;

         if ( SonPID != -1 )
         {
            for (int LocalID=0; LocalID<8; LocalID++)
            {
               if ( patch->ptr[0][lv+1][SonPID+LocalID]->son != -1 )    // if grandson exists
               {
//                flag the corresponding siblings of patch PID
                  Flag_Grandson( lv, PID, LocalID );
   
//                flag the patch PID
                  patch->ptr[0][lv][PID]->flag = true;
               }
            }
         }
      }
   }  // if ( lv < NLEVEL-2 )


// set up the BounFlag_NList and BounFlag_PosList 
   Buf_RecordBoundaryFlag( lv );

} // FUNCTION : Flag_Real



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Grandson
// Description :  Properly flag the siblings of patch "PID" at level "lv" to ensure that the proper-nesting
//                condition is satisfied at level "lv+2"
//
// Parameter   :  lv      : Targeted level to be flagged
//                PID     : Targeted patch ID at level "lv"
//                LocalID : Index of son (0~7) which has its own son (grandson of patch "PID" at level "lv")
//-------------------------------------------------------------------------------------------------------
void Flag_Grandson( const int lv, const int PID, const int LocalID )
{

   switch ( LocalID )
   {
      case 0:
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 0] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 2] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 4] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 6] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[14] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[18] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[10] ]->flag = true;
         break;
      
      case 1:
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 1] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 7] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[19] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[16] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 2] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 4] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[10] ]->flag = true;
         break;
      
      case 2:
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 3] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 8] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[11] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[20] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 0] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[14] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 4] ]->flag = true;
         break;
      
      case 3:
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 5] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[15] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[12] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[22] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 0] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 6] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 2] ]->flag = true;
         break;
      
      case 4:
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 1] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[21] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 9] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[16] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[11] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 3] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 4] ]->flag = true;
         break;
      
      case 5:
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 0] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[24] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 8] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[15] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 3] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[13] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 5] ]->flag = true;
         break;
      
      case 6:
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 1] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[23] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 7] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[17] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 2] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[12] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 5] ]->flag = true;
         break;
       
      case 7:
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 1] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 3] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 5] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[25] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[ 9] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[17] ]->flag = true;
         patch->ptr[0][lv][ patch->ptr[0][lv][PID]->sibling[13] ]->flag = true;
         break;
      
      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "LocalID", LocalID );
   }

} // FUNCTION : Flag_Grandson
