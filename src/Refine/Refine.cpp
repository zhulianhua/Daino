
#include "DAINO.h"

#if ( MODEL == ELBDM  &&  defined DAINO_DEBUG )
void ELBDM_GetPhase_DebugOnly( real *CData, const int CSize );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Refine
// Description :  Construct patches at level "lv+1" according to the flagging result at level "lv" 
//
// Note        :  1. This function will also construct buffer patches at level "lv+1" by calling the function
//                   "Refine_Buffer"
//                2. Data of all sibling-buffer patches must be prepared in advance for creating new 
//                   fine-grid patches by spatial interpolation
//                3. If LOAD_BALANCE is turned on, this function will invoke "LB_Refine" and then return
//
// Parameter   :  lv : Targeted refinement level to be refined
//-------------------------------------------------------------------------------------------------------
void Refine( const int lv )
{

// invoke the load-balance refine function
#  ifdef LOAD_BALANCE
   LB_Refine( lv );
   return;
#  endif


// check
   if ( lv == NLEVEL-1 )
   {
      Aux_Message( stderr, "WARNING : function <%s> should NOT be applied to the finest level !!\n", 
                   __FUNCTION__ );
      return;
   }

// check the synchronization
   if ( NPatchTotal[lv+1] != 0 )    Mis_Check_Synchronization( Time[lv], Time[lv+1], __FUNCTION__, true );


// for the share time-step scheme, we always have PotSg = !FluSg after one step (except during the initialization)
#  ifndef INDIVIDUAL_TIMESTEP
   patch->FluSg[lv+1]  =  patch->FluSg[lv];
   patch->PotSg[lv+1]  = !patch->FluSg[lv];
#  endif

   const int CFluSg = patch->FluSg[lv  ];                // sandglass of fluid variables at level "lv"
   const int FFluSg = patch->FluSg[lv+1];                // sandglass of fluid variables at level "lv+1"
#  ifdef GRAVITY
   const int CPotSg = patch->PotSg[lv  ];                // sandglass of potential       at level "lv"
   const int FPotSg = patch->PotSg[lv+1];                // sandglass of potential       at level "lv+1"
#  endif
   const int Width  = PATCH_SIZE * patch->scale[lv+1];   // scale of a single patch at level "lv+1"

   int *Cr            = NULL;    // corner coordinates
   int *BufGrandTable = NULL;    // table recording the patch IDs of grandson buffer patches
   int *BufSonTable   = NULL;    // table recording the linking index of each buffer father patch to BufGrandTable
   patch_t *Pedigree  = NULL;    // pointer storing the relation of the targeted patch at level "lv"


// parameters for spatial interpolation
   const int CRange[3]     = { PATCH_SIZE, PATCH_SIZE, PATCH_SIZE };
   const int FSize         = 2*PATCH_SIZE;
   const int FStart[3]     = { 0, 0, 0 };

   int NSide_Flu, CGhost_Flu;
   Int_Table( OPT__REF_FLU_INT_SCHEME, NSide_Flu, CGhost_Flu );

   const int CSize_Flu     = PATCH_SIZE + 2*CGhost_Flu;
   const int CStart_Flu[3] = { CGhost_Flu, CGhost_Flu, CGhost_Flu }; 

   real Flu_CData[NCOMP][CSize_Flu][CSize_Flu][CSize_Flu];  // coarse-grid fluid array for interpolation
   real Flu_FData[NCOMP][FSize][FSize][FSize];  // fine-grid fluid array storing the interpolation result

#  ifdef GRAVITY
   int NSide_Pot, CGhost_Pot;
   Int_Table( OPT__REF_POT_INT_SCHEME, NSide_Pot, CGhost_Pot );

   const int CSize_Pot     = PATCH_SIZE + 2*CGhost_Pot;
   const int CStart_Pot[3] = { CGhost_Pot, CGhost_Pot, CGhost_Pot }; 

   real Pot_CData[CSize_Pot][CSize_Pot][CSize_Pot];         // coarse-grid potential array for interpolation
   real Pot_FData[FSize][FSize][FSize];         // fine-grid potential array storing the interpolation result
#  endif



// a. record the tables "BufGrandTable" and "BufFathTable"
// *** the grandson patch indices must be stored in advance before we remove all buffer patches at level "lv+1"
// ------------------------------------------------------------------------------------------------
   if ( lv < NLEVEL-2 )
   {
      const int NBufSon = patch->NPatchComma[lv+1][27] - patch->NPatchComma[lv+1][1];
      const int NBufFa  = patch->NPatchComma[lv  ][27] - patch->NPatchComma[lv  ][1];

      BufGrandTable = new int [NBufSon];
      BufSonTable   = new int [NBufFa ];

//    initialize the table BufSonTable as -1
//#     pragma omp parallel for
      for (int t=0; t<NBufFa; t++)  BufSonTable[t] = -1;

//#     pragma omp parallel for
      for (int m=0; m<NBufSon; m+=8)
      {
//       record the grandson patch ID
         for (int n=m; n<m+8; n++)
         {
            const int BufSonPID = patch->NPatchComma[lv+1][1] + n;

            BufGrandTable[n] = patch->ptr[0][lv+1][BufSonPID]->son; 
         }

//       record the index of BufGrandTable array for father patches with son
         const int BufFaID = patch->ptr[0][lv+1][ patch->NPatchComma[lv+1][1] + m ]->father 
                             - patch->NPatchComma[lv][1];

         BufSonTable[BufFaID] = m;
      }

   } // if ( lv < NLEVEL-2 )



// b. deallocate all buffer patches at level "lv+1"
// ------------------------------------------------------------------------------------------------
//#  pragma omp parallel for
   for (int PID=patch->NPatchComma[lv+1][1]; PID<patch->NPatchComma[lv+1][27]; PID++)
   {
      patch->ptr[0][lv+1][PID]->son = -1;
      patch->pdelete( lv+1, PID );
   }

//#  pragma omp parallel for
   for (int PID=patch->NPatchComma[lv][1]; PID<patch->NPatchComma[lv][27]; PID++)
      patch->ptr[0][lv][PID]->son = -1; 



// c. check the refinement flags for all real patches at level "lv"
// ------------------------------------------------------------------------------------------------
//#  pragma omp parallel for private( ??? )
   for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)
   {
      Pedigree = patch->ptr[0][lv][PID];

//    (c1) construct new child patches if they are newly born (one patch group is allocated at a time)
// ================================================================================================
      if ( Pedigree->flag  &&  Pedigree->son == -1 )
      {

//       (c1.1) construct relation : father -> child
         Pedigree->son = patch->num[lv+1];


//       (c1.2) allocate child patches and construct relation : child -> father
         Cr = Pedigree->corner;

         patch->pnew( lv+1, Cr[0],       Cr[1],       Cr[2],       PID, true, true );
         patch->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2],       PID, true, true );
         patch->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2],       PID, true, true );
         patch->pnew( lv+1, Cr[0],       Cr[1],       Cr[2]+Width, PID, true, true );
         patch->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2],       PID, true, true );
         patch->pnew( lv+1, Cr[0],       Cr[1]+Width, Cr[2]+Width, PID, true, true );
         patch->pnew( lv+1, Cr[0]+Width, Cr[1],       Cr[2]+Width, PID, true, true );
         patch->pnew( lv+1, Cr[0]+Width, Cr[1]+Width, Cr[2]+Width, PID, true, true );


//       (c1.3) assign data to child patches by spatial interpolation
//       (c1.3.1) fill up the central region of CData
         int I, J, K;

         for (int v=0; v<NCOMP; v++)         {
         for (int k=0; k<PATCH_SIZE; k++)    {  K = k + CGhost_Flu;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = j + CGhost_Flu;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = i + CGhost_Flu;

            Flu_CData[v][K][J][I] = patch->ptr[CFluSg][lv][PID]->fluid[v][k][j][i];

         }}}}

#        ifdef GRAVITY
         for (int k=0; k<PATCH_SIZE; k++)    {  K = k + CGhost_Pot;
         for (int j=0; j<PATCH_SIZE; j++)    {  J = j + CGhost_Pot;
         for (int i=0; i<PATCH_SIZE; i++)    {  I = i + CGhost_Pot;

            Pot_CData[K][J][I] = patch->ptr[CPotSg][lv][PID]->pot[k][j][i];

         }}}
#        endif


//       (c1.3.2) fill up the ghost zone of CData (no interpolation is required)
         int Loop_i, Loop_j, Loop_k, Disp_i, Disp_j, Disp_k, Disp_i2, Disp_j2, Disp_k2, I2, J2, K2;
         int SibPID;

         for (int sib=0; sib<NSide_Flu; sib++)
         {
            SibPID = Pedigree->sibling[sib];

//          it will violate the proper-nesting condition if the flagged patch is NOT surrounded by sibling patches
#           ifdef DAINO_DEBUG
            if ( SibPID == -1 )
               Aux_Error( ERROR_INFO, "no sibling patch is found for FaLv %d, FaPID %d, sib %d !!\n",
                          lv, PID, sib );
#           endif

            Loop_i  = TABLE_01( sib, 'x', CGhost_Flu, PATCH_SIZE, CGhost_Flu );
            Loop_j  = TABLE_01( sib, 'y', CGhost_Flu, PATCH_SIZE, CGhost_Flu );
            Loop_k  = TABLE_01( sib, 'z', CGhost_Flu, PATCH_SIZE, CGhost_Flu );
            Disp_i  = TABLE_01( sib, 'x', 0, CGhost_Flu, CGhost_Flu+PATCH_SIZE );
            Disp_j  = TABLE_01( sib, 'y', 0, CGhost_Flu, CGhost_Flu+PATCH_SIZE );
            Disp_k  = TABLE_01( sib, 'z', 0, CGhost_Flu, CGhost_Flu+PATCH_SIZE );
            Disp_i2 = TABLE_01( sib, 'x', PATCH_SIZE-CGhost_Flu, 0, 0 );
            Disp_j2 = TABLE_01( sib, 'y', PATCH_SIZE-CGhost_Flu, 0, 0 );
            Disp_k2 = TABLE_01( sib, 'z', PATCH_SIZE-CGhost_Flu, 0, 0 );
           
            for (int v=0; v<NCOMP; v++)   {
            for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k;   K2 = k + Disp_k2;
            for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j;   J2 = j + Disp_j2;
            for (int i=0; i<Loop_i; i++)  {  I = i + Disp_i;   I2 = i + Disp_i2;

               Flu_CData[v][K][J][I] = patch->ptr[CFluSg][lv][SibPID]->fluid[v][K2][J2][I2];

            }}}}

         } // for (int sib=0; sib<NSide_Flu; sib++)


#        ifdef GRAVITY
         for (int sib=0; sib<NSide_Pot; sib++)
         {
            SibPID = Pedigree->sibling[sib];

//          it will violate the proper-nesting condition if the flagged patch is NOT surrounded by sibling patches
#           ifdef DAINO_DEBUG
            if ( SibPID == -1 )
               Aux_Error( ERROR_INFO, "no sibling patch is found for FaLv %d, FaPID %d, sib %d !!\n",
                          lv, PID, sib );
#           endif

            Loop_i  = TABLE_01( sib, 'x', CGhost_Pot, PATCH_SIZE, CGhost_Pot );
            Loop_j  = TABLE_01( sib, 'y', CGhost_Pot, PATCH_SIZE, CGhost_Pot );
            Loop_k  = TABLE_01( sib, 'z', CGhost_Pot, PATCH_SIZE, CGhost_Pot );
            Disp_i  = TABLE_01( sib, 'x', 0, CGhost_Pot, CGhost_Pot+PATCH_SIZE );
            Disp_j  = TABLE_01( sib, 'y', 0, CGhost_Pot, CGhost_Pot+PATCH_SIZE );
            Disp_k  = TABLE_01( sib, 'z', 0, CGhost_Pot, CGhost_Pot+PATCH_SIZE );
            Disp_i2 = TABLE_01( sib, 'x', PATCH_SIZE-CGhost_Pot, 0, 0 );
            Disp_j2 = TABLE_01( sib, 'y', PATCH_SIZE-CGhost_Pot, 0, 0 );
            Disp_k2 = TABLE_01( sib, 'z', PATCH_SIZE-CGhost_Pot, 0, 0 );
           
            for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k;   K2 = k + Disp_k2;
            for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j;   J2 = j + Disp_j2;
            for (int i=0; i<Loop_i; i++)  {  I = i + Disp_i;   I2 = i + Disp_i2;

               Pot_CData[K][J][I] = patch->ptr[CPotSg][lv][SibPID]->pot[K2][J2][I2];

            }}}

         } // for (int sib=0; sib<NSide_Pot; sib++)
#        endif // #ifdef GRAVITY


//       (c1.3.3) perform spatial interpolation
         const int  CSize_Flu_Temp[3]    = { CSize_Flu, CSize_Flu, CSize_Flu };
         const int  FSize_Temp    [3]    = { FSize, FSize, FSize };
         const bool PhaseUnwrapping_Yes  = true;
         const bool PhaseUnwrapping_No   = false;
         const bool EnsurePositivity_Yes = true;
         const bool EnsurePositivity_No  = false;

//       (c1.3.3.1) determine the variables which must be positive
         bool Positivity[NCOMP];

         for (int v=0; v<NCOMP; v++)
         {
#           if ( MODEL == HYDRO )
            if ( v == DENS  ||  v == ENGY )  Positivity[v] = EnsurePositivity_Yes;
            else                             Positivity[v] = EnsurePositivity_No;

#           elif ( MODEL == MHD )
#           warning : WAIT MHD !!!

#           elif ( MODEL == ELBDM )
            if ( v == DENS )                 Positivity[v] = EnsurePositivity_Yes;
            else                             Positivity[v] = EnsurePositivity_No;

#           else
#           warning : WARNING : DO YOU WANT TO ENSURE THE POSITIVITY OF INTERPOLATION ??          
#           endif // MODEL
         }

//       (c1.3.2.2) interpolation
#        if ( MODEL == ELBDM )
         if ( OPT__INT_PHASE )
         {
//          get the wrapped phase (store in the REAL component)
#           ifdef DAINO_DEBUG
            ELBDM_GetPhase_DebugOnly( &Flu_CData[0][0][0][0], CSize_Flu );
#           else
            for (int k=0; k<CSize_Flu; k++)
            for (int j=0; j<CSize_Flu; j++)
            for (int i=0; i<CSize_Flu; i++)
               Flu_CData[REAL][k][j][i] = ATAN2( Flu_CData[IMAG][k][j][i], Flu_CData[REAL][k][j][i] );
#           endif

//          interpolate density 
            Interpolate( &Flu_CData[DENS][0][0][0], CSize_Flu_Temp, CStart_Flu, CRange, &Flu_FData[DENS][0][0][0],
                         FSize_Temp, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No, 
                         EnsurePositivity_Yes );

//          interpolate phase
            Interpolate( &Flu_CData[REAL][0][0][0], CSize_Flu_Temp, CStart_Flu, CRange, &Flu_FData[REAL][0][0][0],
                         FSize_Temp, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_Yes,
                         EnsurePositivity_No );
         }

         else // if ( OPT__INT_PHASE )
         {
            for (int v=0; v<NCOMP; v++)
            Interpolate( &Flu_CData[v][0][0][0], CSize_Flu_Temp, CStart_Flu, CRange, &Flu_FData[v][0][0][0], 
                         FSize_Temp, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No, 
                         Positivity[v] );
         }

         if ( OPT__INT_PHASE )
         {
//          retrieve real and imaginary parts
            real Amp, Phase, Rho;

            for (int k=0; k<FSize; k++)
            for (int j=0; j<FSize; j++)
            for (int i=0; i<FSize; i++)
            {
               Phase                    = Flu_FData[REAL][k][j][i];
               Rho                      = Flu_FData[DENS][k][j][i];

               if ( Rho < 0.0 )    
                  Aux_Error( ERROR_INFO, "negative density (%14.7e) is obtained in %s !!\n", Rho, __FUNCTION__ );

               Amp                      = SQRT( Rho );
               Flu_FData[REAL][k][j][i] = Amp*COS( Phase );
               Flu_FData[IMAG][k][j][i] = Amp*SIN( Phase );
            }
         }

#        else // #if ( MODEL == ELBDM )

         for (int v=0; v<NCOMP; v++)
         Interpolate( &Flu_CData[v][0][0][0], CSize_Flu_Temp, CStart_Flu, CRange, &Flu_FData[v][0][0][0], 
                      FSize_Temp, FStart, 1, OPT__REF_FLU_INT_SCHEME, PhaseUnwrapping_No, 
                      Positivity[v] );

#        endif // #if ( MODEL == ELBDM ) ... else 


#        ifdef GRAVITY
         const int CSize_Pot_Temp[3] = { CSize_Pot, CSize_Pot, CSize_Pot };

         Interpolate( &Pot_CData[0][0][0], CSize_Pot_Temp, CStart_Pot, CRange, &Pot_FData[0][0][0],
                      FSize_Temp, FStart, 1, OPT__REF_POT_INT_SCHEME, PhaseUnwrapping_No,
                      EnsurePositivity_No );
#        endif


//       (c1.3.4) copy data from IntData to patch pointers
         int SonPID;

         for (int LocalID=0; LocalID<8; LocalID++)
         {
            SonPID = patch->num[lv+1] - 8 + LocalID;
            Disp_i = TABLE_02( LocalID, 'x', 0, PATCH_SIZE ); 
            Disp_j = TABLE_02( LocalID, 'y', 0, PATCH_SIZE ); 
            Disp_k = TABLE_02( LocalID, 'z', 0, PATCH_SIZE ); 
               
//          fluid data
            for (int v=0; v<NCOMP; v++)         {
            for (int k=0; k<PATCH_SIZE; k++)    {  K = k + Disp_k;
            for (int j=0; j<PATCH_SIZE; j++)    {  J = j + Disp_j;
            for (int i=0; i<PATCH_SIZE; i++)    {  I = i + Disp_i;

               patch->ptr[FFluSg][lv+1][SonPID]->fluid[v][k][j][i] = Flu_FData[v][K][J][I];

            }}}}

//          potential data
#           ifdef GRAVITY
            for (int k=0; k<PATCH_SIZE; k++)    {  K = k + Disp_k;
            for (int j=0; j<PATCH_SIZE; j++)    {  J = j + Disp_j;
            for (int i=0; i<PATCH_SIZE; i++)    {  I = i + Disp_i;

               patch->ptr[FPotSg][lv+1][SonPID]->pot[k][j][i] = Pot_FData[K][J][I];

            }}}
#           endif

//          rescale real and imaginary parts to get the correct density in ELBDM if OPT__INT_PHASE is off
#           if ( MODEL == ELBDM )
            real Real, Imag, Rho_Wrong, Rho_Corr, Rescale;

            if ( !OPT__INT_PHASE )
            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++)
            {
               Real = patch->ptr[FFluSg][lv+1][SonPID]->fluid[REAL][k][j][i];
               Imag = patch->ptr[FFluSg][lv+1][SonPID]->fluid[IMAG][k][j][i];

//             patch->ptr[FFluSg][lv+1][SonPID]->fluid[DENS][k][j][i] = Real*Real + Imag*Imag;

               Rho_Wrong = Real*Real + Imag*Imag;
               Rho_Corr  = patch->ptr[FFluSg][lv+1][SonPID]->fluid[DENS][k][j][i];
               Rescale   = SQRT( Rho_Corr/Rho_Wrong );

               if ( Rho_Corr < 0.0 )    
                  Aux_Error( ERROR_INFO, "negative density (%14.7e) is obtained in %s !!\n", 
                             Rho_Corr, __FUNCTION__ );

               patch->ptr[FFluSg][lv+1][SonPID]->fluid[REAL][k][j][i] *= Rescale;
               patch->ptr[FFluSg][lv+1][SonPID]->fluid[IMAG][k][j][i] *= Rescale;
            }
#           endif
         }

      } // if ( Pedigree->flag  &&  Pedigree->son == -1 )


//    (c2) remove unflagged child patches if they originally existed (one patch group is removed at a time)
// ================================================================================================
      else if ( !Pedigree->flag  &&  Pedigree->son != -1 )
      {

//       (c2.1) deallocate the unflagged child patches
         const int SonPID0 = Pedigree->son;
         const int NewPID0 = SonPID0;
         const int OldPID0 = patch->num[lv+1] - 8;

         for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)     patch->pdelete( lv+1, SonPID );


//       (c2.2) construct relation : father -> son
         Pedigree->son = -1;


//       (c2.3) relink the child patch pointers so that no patch indices are skipped
         if ( NewPID0 != OldPID0 )  
         {
            int NewPID, OldPID, GrandPID0, FaPID;

            for (int t=0; t<8; t++)
            {
               NewPID = NewPID0 + t;
               OldPID = OldPID0 + t;

//             relink pointers
               patch->ptr[0][lv+1][NewPID] = patch->ptr[0][lv+1][OldPID];
               patch->ptr[1][lv+1][NewPID] = patch->ptr[1][lv+1][OldPID];

///            set redundant patch pointers as NULL
               patch->ptr[0][lv+1][OldPID] = NULL; 
               patch->ptr[1][lv+1][OldPID] = NULL; 

//             re-construct relation : grandson -> son
               GrandPID0 = patch->ptr[0][lv+1][NewPID]->son;
               if ( GrandPID0 != -1 )
               {
                  for (int GrandPID=GrandPID0; GrandPID<GrandPID0+8; GrandPID++)
                     patch->ptr[0][lv+2][GrandPID]->father = NewPID;
               }
            }

//          re-construct relation : father -> son
            FaPID = patch->ptr[0][lv+1][NewPID0]->father;
            patch->ptr[0][lv][FaPID]->son = NewPID0;

         } // if ( NewPID0 != OldPID0 )

      } // else if ( !Pedigree->flag  &&  Pedigree->son != -1 )
   } // for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)


// initialize the patch->NPatchComma list for the buffer patches
   for (int m=1; m<28; m++)   patch->NPatchComma[lv+1][m] = patch->num[lv+1];


// record the numbers of real and data patches for the out-of-core computing
#  ifdef OOC
   ooc.NDataPatch[ooc.Rank][lv  ] = 0;
   ooc.NDataPatch[ooc.Rank][lv+1] = 0;
   ooc.NRealPatch[ooc.Rank][lv+1] = patch->NPatchComma[lv+1][1];

   for (int PID=0; PID<patch->NPatchComma[lv  ][1]; PID++)
      if ( patch->ptr[0][lv  ][PID]->son == -1 )      ooc.NDataPatch[ooc.Rank][lv  ]++; 

   for (int PID=0; PID<patch->NPatchComma[lv+1][1]; PID++)
      if ( patch->ptr[0][lv+1][PID]->son == -1 )      ooc.NDataPatch[ooc.Rank][lv+1]++;
#  endif


// d. refine buffer patches
// ------------------------------------------------------------------------------------------------
   Refine_Buffer( lv, BufSonTable, BufGrandTable );


// deallocate tables 
   if ( lv < NLEVEL-2 )
   {
      delete [] BufGrandTable;
      delete [] BufSonTable;
   }



// e. re-construct tables and sibling relations
// ------------------------------------------------------------------------------------------------

// set up the BounP_IDMap for the level just created
   Buf_RecordBoundaryPatch( lv+1 );


// construct relation : siblings
   SiblingSearch( lv+1 );


// allocate flux arrays for level "lv"
   if ( patch->WithFlux )
   Flu_AllocateFluxArray( lv );


// allocate flux arrays for level "lv+1"
   if ( lv < NLEVEL-2  &&  patch->WithFlux )
      Flu_AllocateFluxArray( lv+1 );


// get the IDs of patches for sending and receiving data between neighbor ranks
   Buf_RecordExchangeDataPatchID( lv+1 );


// get the total number of patches at lv+1   
#  ifndef OOC
   Mis_GetTotalPatchNumber( lv+1 );
#  endif

} // FUNCTION : Refine



#if ( MODEL == ELBDM  &&  defined DAINO_DEBUG )
//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetPhase_DebugOnly
// Description :  Alternative function to calculate phase in the debug mode so that the functions "Refine" and 
//                "LB_Refine_AllocateNewPatch" will give EXACTLY THE SAME RESULTS (even the round-off errors
//                are the same)
//
// Note        :  It is found that it is necessary to use this alternative function to calculate phase 
//                in order to make runs with "SERIAL", "NO LOAD_BALANCE", and "LOAD_BALANCE" have exactly
//                the same results
//
// Parameter   :  CData : Coarse-grid array
//                CSize : Size of CData in each direction 
//-------------------------------------------------------------------------------------------------------
void ELBDM_GetPhase_DebugOnly( real *CData, const int CSize )
{

   const int CSize_1v = CSize*CSize*CSize;

   real *const CData_Dens = CData + DENS*CSize_1v;
   real *const CData_Real = CData + REAL*CSize_1v;
   real *const CData_Imag = CData + IMAG*CSize_1v;

   for (int t=0; t<CSize_1v; t++)   CData_Real[t] = ATAN2( CData_Imag[t], CData_Real[t] );

} // FUNCTION : 
#endif // #if ( MODEL == ELBDM  &&  defined DAINO_DEBUG )
