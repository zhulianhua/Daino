
#include "DAINO.h"

void InterpolateGhostZone( const int lv, const int PID, real IntData[], const int SibID, const bool IntTime, 
                           const int GhostSize, const int FluSg, const int PotSg, 
                           const IntScheme_t IntScheme, const int NTSib[], int *TSib[],
                           const int NVar_Flu, const int TFluVarIdxList[], const bool PrepPot,
                           const bool IntPhase );
void SetTargetSibling( int NTSib[], int* TSib[] );
static int Table_01( const int SibID, const char dim, const int Count, const int GhostSize );
static int Table_02( const int lv, const int PID, const int Side );




//-------------------------------------------------------------------------------------------------------
// Function    :  Prepare_PatchGroupData
// Description :  Prepare uniform data including ghost zone for the targeted patch groups
//
// Note        :  1. Use the input parameter "TVar" to control the targeted variables
//                   --> TVar can be any combination of the symbolic constants defined in "Macro.h" 
//                       (e.g., "TVar = _DENS", "TVar = _MOMX|ENGY", or "TVar = _FLU|_POTE")
//                2. If "GhostSize != 0" --> the function "InterpolateGhostZone" will be used to fill up the
//                   ghost-zone values by spatial interpolation if the corresponding sibling patches do
//                   NOT exist
//                3. The parameter "PrepTime" is used to determine whether or not the "temporal interpolation"
//                   is necessary
//
// Parameter   :  lv             : Targeted refinement level
//                PrepTime       : Targeted physical time to prepare data
//                h_Input_Array  : Host array to store the prepared data
//                GhostSize      : Number of ghost zones to be prepared
//                NPG            : Number of patch groups prepared at a time
//                PID0_List      : List recording the patch indicies with LocalID==0 to be prepared
//                TVar           : Targeted variables to be prepared
//                                 --> Supported variables in different models:
//                                     HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY, _FLU [, _POTE]
//                                     MHD   : 
//                                     ELBDM : _DENS, _REAL, _IMAG [, _POTE]
//                IntScheme      : Interpolation scheme
//                                 --> currently supported schemes include
//                                     INT_CENTRAL : central 
//                                     INT_MINMOD  : MinMod 
//                                     INT_VANLEER : vanLeer
//                                     INT_CQUAD   : conservative quadratic
//                                     INT_QUAD    : quadratic
//                PrepUnit       : Whether or not to separate the prepared data into individual patches 
//                                 --> UNIT_PATCH      : prepare data "patch by patch"
//                                     UNIT_PATCHGROUP : prepare data "patch group by patch group"
//                NSide          : Number of sibling directions to prepare data
//                                 --> NSIDE_06 (=  6) : prepare only sibling directions 0~5
//                                     NSIDE_26 (= 26) : prepare all sibling directions 0~25
//                IntPhase       : true --> Perform interpolation on rho/phase instead of real/imag parts in ELBDM
//                                      --> TVar must contain _REAL and _IMAG
//-------------------------------------------------------------------------------------------------------
void Prepare_PatchGroupData( const int lv, const double PrepTime, real *h_Input_Array, const int GhostSize, 
                             const int NPG, const int *PID0_List, const int TVar, const IntScheme_t IntScheme,
                             const PrepUnit_t PrepUnit, const NSide_t NSide, const bool IntPhase )
{

// check
#  ifdef GRAVITY
   if (  TVar & ~( _FLU | _POTE )  )  
      Aux_Error( ERROR_INFO, "unsupported parameter %s = %d !!\n", "TVar", TVar );
#  else
   if ( TVar & ~_FLU )   
      Aux_Error( ERROR_INFO, "unsupported parameter %s = %d !!\n", "TVar", TVar );
#  endif

   if ( IntPhase )
   {
#     if ( MODEL == ELBDM )
      if (  !(TVar & _REAL)  ||  !(TVar & _IMAG)  )
      Aux_Error( ERROR_INFO, "real and/or imag parts are not found for phase interpolation in ELBDM !!\n" );
#     else
      Aux_Error( ERROR_INFO, "\"interpolation on phase\" is useful only in ELBDM !!\n" );
#     endif
   }


   const int PGSize1D = 2*( PATCH_SIZE + GhostSize );    // size of a single patch group including the ghost zone
   const int PGSize3D = PGSize1D*PGSize1D*PGSize1D;
   const int FluSg    = patch->FluSg[lv];

#  ifdef GRAVITY
   const int PotSg    = patch->PotSg[lv];
   const bool PrepPot = ( TVar & _POTE ) ? true : false;
#  endif

// TFluVarIdxList : List recording the targeted fluid variable indices ( = [0 ... NCOMP-1] )
   int  NTSib[26], *TSib[26], NVar_Flu, NVar_Tot, TFluVarIdxList[NCOMP];


// set up the targeted sibling indices for the function "InterpolateGhostZone"
   SetTargetSibling( NTSib, TSib );


// determine the components to be prepared
   NVar_Flu = 0;

   for (int v=0; v<NCOMP; v++)
      if ( TVar & (1<<v) )    TFluVarIdxList[ NVar_Flu++ ] = v;

   NVar_Tot = NVar_Flu;
#  ifdef GRAVITY
   if ( PrepPot )    NVar_Tot ++; 
#  endif

   if ( NVar_Tot == 0 )    
   {
      Aux_Message( stderr, "WARNING : no targeted variable is found !!\n" );
      return;
   }


// start to prepare data
#  pragma omp parallel
   {
      int J, K, I2, J2, K2, Idx1, Idx2, PID0, TFluVarIdx;

//    Array : array to store the prepared data of one patch group (including the ghost-zone data) 
      real *Array_Ptr = NULL;
      real *Array     = new real [ NVar_Tot*PGSize3D ];

      
//    prepare eight nearby patches (one patch group) at a time 
#     pragma omp for
      for (int TID=0; TID<NPG; TID++)
      {
         PID0 = PID0_List[TID];
         
//       a. fill up the central region of Array (ghost zone is not filled up yet)
// ------------------------------------------------------------------------------------------------------------
         for (int LocalID=0; LocalID<8; LocalID++ )
         {
            const int PID    = PID0 + LocalID;
            const int Disp_i = TABLE_02( LocalID, 'x', GhostSize, GhostSize+PATCH_SIZE );
            const int Disp_j = TABLE_02( LocalID, 'y', GhostSize, GhostSize+PATCH_SIZE );
            const int Disp_k = TABLE_02( LocalID, 'z', GhostSize, GhostSize+PATCH_SIZE );

            Array_Ptr = Array;
            
//          fluid data            
            for (int v=0; v<NVar_Flu; v++)         
            {  
               TFluVarIdx = TFluVarIdxList[v]; 

               for (int k=0; k<PATCH_SIZE; k++)    {  K    = k + Disp_k;
               for (int j=0; j<PATCH_SIZE; j++)    {  J    = j + Disp_j;   
                                                      Idx1 = IDX321( Disp_i, J, K, PGSize1D, PGSize1D );
               for (int i=0; i<PATCH_SIZE; i++)    {
               
                  Array_Ptr[ Idx1 ++ ] = patch->ptr[FluSg][lv][PID]->fluid[TFluVarIdx][k][j][i];
               
               }}}

               Array_Ptr += PGSize3D;
            }


//          potential data            
#           ifdef GRAVITY
            if ( PrepPot )
            {
               for (int k=0; k<PATCH_SIZE; k++)    {  K    = k + Disp_k;
               for (int j=0; j<PATCH_SIZE; j++)    {  J    = j + Disp_j;
                                                      Idx1 = IDX321( Disp_i, J, K, PGSize1D, PGSize1D );
               for (int i=0; i<PATCH_SIZE; i++)    {
               
                  Array_Ptr[ Idx1 ++ ] = patch->ptr[PotSg][lv][PID]->pot[k][j][i];
               
               }}}
            } // if ( PrepPot )
#           endif
         } // for (int LocalID=0; LocalID<8; LocalID++ )


//       b. fill up the ghost zone of Array
// ------------------------------------------------------------------------------------------------------------
         for (int Side=0; Side<NSide; Side++)
         {
//          nothing to do if no ghost zone is required
            if ( GhostSize == 0 )   break;


            const int SibPID0 = Table_02( lv, PID0, Side );    // the 0th patch of the sibling patch group

//          (b1) if the targeted sibling patch exists --> just copy data from the nearby patch at the same level
            if ( SibPID0 != -1 )
            {
               const int Loop_i  = TABLE_01( Side, 'x', GhostSize, PATCH_SIZE, GhostSize );
               const int Loop_j  = TABLE_01( Side, 'y', GhostSize, PATCH_SIZE, GhostSize );
               const int Loop_k  = TABLE_01( Side, 'z', GhostSize, PATCH_SIZE, GhostSize );
               const int Disp_i2 = TABLE_01( Side, 'x', PATCH_SIZE-GhostSize, 0, 0 );
               const int Disp_j2 = TABLE_01( Side, 'y', PATCH_SIZE-GhostSize, 0, 0 );
               const int Disp_k2 = TABLE_01( Side, 'z', PATCH_SIZE-GhostSize, 0, 0 );
                  
//###OPTIMIZATION: simplify TABLE_03 and TABLE_04               
               for (int Count=0; Count<TABLE_04( Side ); Count++)
               {
                  const int SibPID = TABLE_03( Side, Count ) + SibPID0;
                  const int Disp_i = Table_01( Side, 'x', Count, GhostSize );
                  const int Disp_j = Table_01( Side, 'y', Count, GhostSize );
                  const int Disp_k = Table_01( Side, 'z', Count, GhostSize );
                  
                  Array_Ptr = Array;

                  for (int v=0; v<NVar_Flu; v++)  
                  {
                     TFluVarIdx = TFluVarIdxList[v]; 

                     for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k;   K2 = k + Disp_k2;
                     for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j;   J2 = j + Disp_j2;
                                                      Idx1 = IDX321( Disp_i, J, K, PGSize1D, PGSize1D );
                     for (I2=Disp_i2; I2<Disp_i2+Loop_i; I2++) {

                        Array_Ptr[ Idx1 ++ ] = patch->ptr[FluSg][lv][SibPID]->fluid[TFluVarIdx][K2][J2][I2];

                     }}}

                     Array_Ptr += PGSize3D;
                  }


#                 ifdef GRAVITY
                  if ( PrepPot )
                  {
                     for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k;   K2 = k + Disp_k2;
                     for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j;   J2 = j + Disp_j2;
                                                      Idx1 = IDX321( Disp_i, J, K, PGSize1D, PGSize1D );
                     for (I2=Disp_i2; I2<Disp_i2+Loop_i; I2++) {

                        Array_Ptr[ Idx1 ++ ] = patch->ptr[PotSg][lv][SibPID]->pot[K2][J2][I2];

                     }}}
                  }
#                 endif
               } // for (int Count=0; Count<TABLE_04( Side ); Count++)
            } // if ( SibPID0 != -1 )  


//          (b2) if the targeted sibling patch does not exist --> interpolate from patches at level lv-1
            else
            {
//             interpolation should never be applied to the base level
               if ( lv == 0 )    Aux_Error( ERROR_INFO, "performing interpolation at the base level !!\n" );


//             allocate the array to store the interpolation result
               const int GhostSize_Padded = GhostSize + (GhostSize&1);

               int FSize[3];
               for (int d=0; d<3; d++)    
                  FSize[d] = TABLE_01( Side, 'x'+d, GhostSize_Padded, 2*PATCH_SIZE, GhostSize_Padded );

               real *IntData_Ptr = NULL;
               real *IntData     = new real [ NVar_Tot*FSize[0]*FSize[1]*FSize[2] ];


//             determine the parameters for the spatial and temporal interpolations
               const int FaPID    = patch->ptr[0][lv][PID0]->father;
               const int FaSibPID = patch->ptr[0][lv-1][FaPID]->sibling[Side];

               bool IntTime; 
               int  IntFluSg; 
#              ifdef GRAVITY
               int  IntPotSg;
#              endif

#              ifdef INDIVIDUAL_TIMESTEP
               if      (  Mis_Check_Synchronization( PrepTime,      Time[lv-1], NULL, false )  )
               {
                  IntTime  = false;
                  IntFluSg = patch->FluSg[lv-1];
#                 ifdef GRAVITY
                  IntPotSg = patch->PotSg[lv-1];
#                 endif
               }
               else if (  Mis_Check_Synchronization( PrepTime, Time_Prev[lv-1], NULL, false )  )
               {
                  IntTime  = false;
                  IntFluSg = 1 - patch->FluSg[lv-1];
#                 ifdef GRAVITY
                  IntPotSg = 1 - patch->PotSg[lv-1];
#                 endif
               }
               else
               {
                  IntTime  = ( OPT__INT_TIME ) ? true : false;
                  IntFluSg = 1 - patch->FluSg[lv-1];
#                 ifdef GRAVITY
                  IntPotSg = 1 - patch->PotSg[lv-1];
#                 endif
               }
#              else // #ifdef INDIVIDUAL_TIMESTEP ... else ...

               IntTime  = false;
               IntFluSg = patch->FluSg[lv];
#              ifdef GRAVITY
               IntPotSg = patch->PotSg[lv];
#              endif
#              endif // #ifdef INDIVIDUAL_TIMESTEP ... else ...


//             verify if the targeted time for temporal interpolation is within the permitted range
               if ( IntTime )
               {
                  if ( PrepTime <= Time_Prev[lv-1]  ||  PrepTime-Time[lv-1] >= 1.e-12*fabs(Time[lv-1]) )
                  {
                     Aux_Message( stderr, "ERROR : targeted time for temporal interpolation is incorrect !!\n" );
                     Aux_Message( stderr, "        (lv %d, T_Prep %20.14e, T_Min %20.14e, T_Max %20.14e)\n", 
                                  PrepTime, Time_Prev[lv-1], Time[lv-1] );
                     MPI_Exit();
                  }
               }


//             perform interpolation and store the results in IntData
#              ifdef GRAVITY
               InterpolateGhostZone( lv-1, FaSibPID, IntData, Side, IntTime, GhostSize, IntFluSg, IntPotSg, 
                                     IntScheme, NTSib, TSib, NVar_Flu, TFluVarIdxList, PrepPot, IntPhase );
#              else
               InterpolateGhostZone( lv-1, FaSibPID, IntData, Side, IntTime, GhostSize, IntFluSg, NULL_INT, 
                                     IntScheme, NTSib, TSib, NVar_Flu, TFluVarIdxList, false, IntPhase );
#              endif


//             properly copy data from IntData array to Array
               const int NUseless = GhostSize & 1;
               const int Loop_i   = TABLE_01( Side, 'x', GhostSize, 2*PATCH_SIZE, GhostSize );
               const int Loop_j   = TABLE_01( Side, 'y', GhostSize, 2*PATCH_SIZE, GhostSize );
               const int Loop_k   = TABLE_01( Side, 'z', GhostSize, 2*PATCH_SIZE, GhostSize );
               const int Disp_i1  = TABLE_01( Side, 'x', 0, GhostSize, GhostSize+2*PATCH_SIZE );
               const int Disp_j1  = TABLE_01( Side, 'y', 0, GhostSize, GhostSize+2*PATCH_SIZE );
               const int Disp_k1  = TABLE_01( Side, 'z', 0, GhostSize, GhostSize+2*PATCH_SIZE );
               const int Disp_i2  = TABLE_01( Side, 'x', NUseless, 0, 0 );
               const int Disp_j2  = TABLE_01( Side, 'y', NUseless, 0, 0 );
               const int Disp_k2  = TABLE_01( Side, 'z', NUseless, 0, 0 );

               Array_Ptr   = Array;
               IntData_Ptr = IntData;

               for (int v=0; v<NVar_Tot; v++)  
               {
                  for (int k=0; k<Loop_k; k++)  {  K = k + Disp_k1;  K2 = k + Disp_k2;
                  for (int j=0; j<Loop_j; j++)  {  J = j + Disp_j1;  J2 = j + Disp_j2;
                                                   Idx1 = IDX321( Disp_i1, J,  K,  PGSize1D, PGSize1D );
                                                   Idx2 = IDX321( Disp_i2, J2, K2, FSize[0], FSize[1] );
                  for (int i=0; i<Loop_i; i++)  {

                     Array_Ptr[ Idx1 ++ ] = IntData_Ptr[ Idx2 ++ ];

                  }}}

                  Array_Ptr   += PGSize3D;
                  IntData_Ptr += FSize[0]*FSize[1]*FSize[2];
               }

               delete [] IntData;

            } // if ( SibPID0 != -1 ) ... else ...
         } // for (int Side=0; Side<NSide; Side++)
            

//       c. copy data from Array to h_Input_Array
// ------------------------------------------------------------------------------------------------------------
         if ( PrepUnit == UNIT_PATCH ) // separate the prepared patch group data into individual patches
         {
            const int PSize1D = PATCH_SIZE + 2*GhostSize;  // size of a single patch including the ghost zone
            const int PSize3D = PSize1D*PSize1D*PSize1D;
            real *InArray_Ptr = NULL;

            for (int LocalID=0; LocalID<8; LocalID++)
            {
               const int N      = 8*TID + LocalID;
               const int Disp_i = TABLE_02( LocalID, 'x', 0, PATCH_SIZE ); 
               const int Disp_j = TABLE_02( LocalID, 'y', 0, PATCH_SIZE ); 
               const int Disp_k = TABLE_02( LocalID, 'z', 0, PATCH_SIZE ); 

               Array_Ptr   = Array;
               InArray_Ptr = h_Input_Array + N*NVar_Tot*PSize3D;
               Idx2        = 0;
               
               for (int v=0; v<NVar_Tot; v++)
               {
                  for (int k=Disp_k; k<Disp_k+PSize1D; k++)
                  for (int j=Disp_j; j<Disp_j+PSize1D; j++)
                  {
                     Idx1 = IDX321( Disp_i, j, k, PGSize1D, PGSize1D );

                     for (int i=0; i<PSize1D; i++)    InArray_Ptr[ Idx2 ++ ] = Array_Ptr[ Idx1 ++ ];
                  }

                  Array_Ptr += PGSize3D;
               }
            }
         } // if ( PatchByPatch )

         else if ( PrepUnit == UNIT_PATCHGROUP )
            memcpy( h_Input_Array + TID*NVar_Tot*PGSize3D, Array, NVar_Tot*PGSize3D*sizeof(real) );

         else
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "PrepUnit", PrepUnit );

      } // for (int TID=0; TID<NPG; TID++)

      delete [] Array;

   } // OpenMP parallel region


// free memroy
   for (int s=0; s<26; s++)   delete [] TSib[s];

} // FUNCTION : Prepare_PatchGroupData



// ============
// |  Tables  | 
// ============

//-------------------------------------------------------------------------------------------------------
// Function    :  Table_01 
// Description :  Return the displacement for the function "Prepare_PatchGroupData"
//
// Parameter   :  SibID     : Sibling index (0~25) 
//                dim       : Targeted spatial direction (x/y/z)
//                Count     : Patch counter (0~3)
//                GhostSize : Number of ghost zones
//-------------------------------------------------------------------------------------------------------
int Table_01( const int SibID, const char dim, const int Count, const int GhostSize )
{

   switch ( dim )
   {
      case 'x':
      {
         switch ( SibID )
         {
            case 0:case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
               return 0;
               
            case 2: case 3:
            {
               switch ( Count )
               {
                  case 0: case 2:   return GhostSize;
                  case 1: case 3:   return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
               
            case 4: case 5: 
            {
               switch ( Count )
               {
                  case 0: case 1:   return GhostSize;
                  case 2: case 3:   return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
               
            case 10: case 11: case 12: case 13:
            {
               switch ( Count )
               {
                  case 0:  return GhostSize;
                  case 1:  return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
               
            case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
               return GhostSize + 2*PATCH_SIZE;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SibID", SibID );

         } // switch ( SibID )
      } // case 'x':


      case 'y':
      {
         switch ( SibID )
         {
            case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
               return 0;
               
            case 0: case 1:
            {
               switch ( Count )
               {
                  case 0: case 1:   return GhostSize;
                  case 2: case 3:   return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }

            case 4: case 5:
            {
               switch ( Count )
               {
                  case 0: case 2:   return GhostSize;
                  case 1: case 3:   return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }

            case 14: case 15: case 16: case 17:
            {
               switch ( Count ) 
               {
                  case 0:  return GhostSize;
                  case 1:  return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
                  
            case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
               return GhostSize + 2*PATCH_SIZE; 

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SibID", SibID );

         } // switch ( SibID )
      } // case 'y':


      case 'z':
      {
         switch ( SibID )
         {
            case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
               return 0;
               
            case 0: case 1:
            {
               switch ( Count )
               {
                  case 0: case 2:   return GhostSize;
                  case 1: case 3:   return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
               
            case 2: case 3:
            {
               switch ( Count )
               {
                  case 0: case 1:   return GhostSize;
                  case 2: case 3:   return GhostSize + PATCH_SIZE;                     
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
               
            case 6: case 7: case 8: case 9:
            {
               switch ( Count )
               {
                  case 0:  return GhostSize;
                  case 1:  return GhostSize + PATCH_SIZE;
                  default:
                     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d !!\n", 
                                "SibID", SibID, "Count", Count );
               }
            }
                  
            case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
               return GhostSize + 2*PATCH_SIZE;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SibID", SibID );

         } // switch ( SibID )
      } // case 'z':


      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %c !!\n", "dim", dim );
         exit(1);
   } // switch ( dim )

   return NULL_INT;

} // FUNCTION : Table_01



//-------------------------------------------------------------------------------------------------------
// Function    :  Table_02 
// Description :  Return the patch ID of the 0th patch (local ID = 0) of the sibling patch group 
//
// Note        :  Work for the function "Prepare_PatchGroupData" 
//
// Parameter   :  lv    : Targeted refinement level 
//                PID   : Targeted patch ID to find its sibling patches
//                Side  : Sibling index (0~25) 
//-------------------------------------------------------------------------------------------------------
int Table_02( const int lv, const int PID, const int Side )
{

   int Sib;

   switch ( Side )
   {
      case 0: 
         Sib = patch->ptr[0][lv][PID  ]->sibling[0];
         if ( Sib != -1 )  return Sib-1;
         else              return -1;

      case 1:
         Sib = patch->ptr[0][lv][PID+1]->sibling[1];
         if ( Sib != -1 )  return Sib;
         else              return -1;

      case 2:
         Sib = patch->ptr[0][lv][PID  ]->sibling[2];
         if ( Sib != -1 )  return Sib-2;
         else              return -1;

      case 3:
         Sib = patch->ptr[0][lv][PID+2]->sibling[3];
         if ( Sib != -1 )  return Sib;
         else              return -1;

      case 4:
         Sib = patch->ptr[0][lv][PID  ]->sibling[4];
         if ( Sib != -1 )  return Sib-3;
         else              return -1;

      case 5: 
         Sib = patch->ptr[0][lv][PID+3]->sibling[5];
         if ( Sib != -1 )  return Sib;
         else              return -1;

      case 6:
         Sib = patch->ptr[0][lv][PID  ]->sibling[6];
         if ( Sib != -1 )  return Sib-4;
         else              return -1;

      case 7:
         Sib = patch->ptr[0][lv][PID+1]->sibling[7];
         if ( Sib != -1 )  return Sib-2;
         else              return -1;

      case 8:
         Sib = patch->ptr[0][lv][PID+2]->sibling[8];
         if ( Sib != -1 )  return Sib-1;
         else              return -1;

      case 9:
         Sib = patch->ptr[0][lv][PID+4]->sibling[9];
         if ( Sib != -1 )  return Sib;
         else              return -1;

      case 10: 
         Sib = patch->ptr[0][lv][PID  ]->sibling[10];
         if ( Sib != -1 )  return Sib-5;
         else              return -1;

      case 11:
         Sib = patch->ptr[0][lv][PID+2]->sibling[11];
         if ( Sib != -1 )  return Sib-3;
         else              return -1;

      case 12:
         Sib = patch->ptr[0][lv][PID+3]->sibling[12];
         if ( Sib != -1 )  return Sib-2;
         else              return -1;

      case 13:
         Sib = patch->ptr[0][lv][PID+5]->sibling[13];
         if ( Sib != -1 )  return Sib;
         else              return -1;

      case 14:
         Sib = patch->ptr[0][lv][PID  ]->sibling[14];
         if ( Sib != -1 )  return Sib-6;
         else              return -1;

      case 15: 
         Sib = patch->ptr[0][lv][PID+3]->sibling[15];
         if ( Sib != -1 )  return Sib-1;
         else              return -1;

      case 16:
         Sib = patch->ptr[0][lv][PID+1]->sibling[16];
         if ( Sib != -1 )  return Sib-3;
         else              return -1;

      case 17:
         Sib = patch->ptr[0][lv][PID+6]->sibling[17];
         if ( Sib != -1 )  return Sib;
         else              return -1;

      case 18:
         Sib = patch->ptr[0][lv][PID  ]->sibling[18];
         if ( Sib != -1 )  return Sib-7;
         else              return -1;

      case 19:
         Sib = patch->ptr[0][lv][PID+1]->sibling[19];
         if ( Sib != -1 )  return Sib-5;
         else              return -1;

      case 20: 
         Sib = patch->ptr[0][lv][PID+2]->sibling[20];
         if ( Sib != -1 )  return Sib-6;
         else              return -1;

      case 21:
         Sib = patch->ptr[0][lv][PID+4]->sibling[21];
         if ( Sib != -1 )  return Sib-3;
         else              return -1;

      case 22:
         Sib = patch->ptr[0][lv][PID+3]->sibling[22];
         if ( Sib != -1 )  return Sib-4;
         else              return -1;

      case 23:
         Sib = patch->ptr[0][lv][PID+6]->sibling[23];
         if ( Sib != -1 )  return Sib-2;
         else              return -1;

      case 24:
         Sib = patch->ptr[0][lv][PID+5]->sibling[24];
         if ( Sib != -1 )  return Sib-1;
         else              return -1;

      case 25: 
         Sib = patch->ptr[0][lv][PID+7]->sibling[25];
         if ( Sib != -1 )  return Sib;
         else              return -1;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Side", Side );
         exit(1);

   } // switch ( Side )

   return NULL_INT;

} // FUNCTION : Table_02



//-------------------------------------------------------------------------------------------------------
// Function    :  SetTargetSibling 
// Description :  Set the targeted sibling directions for preparing ghost-zone data at coarse-grid level
//
// Note        :  1. Works for the functions "Prepare_PatchGroupData" and "LB_RecordExchangeDataPatchID.cpp"
//                2. TSib needs to be deallocated manually
// 
// Parameter   :  NTSib : Number of targeted sibling patches along different sibling directions 
//                TSib  : Targeted sibling indices along different sibling directions
//-------------------------------------------------------------------------------------------------------
void SetTargetSibling( int NTSib[], int* TSib[] )
{

   for (int t= 0; t< 6; t++)  NTSib[t] = 17;
   for (int t= 6; t<18; t++)  NTSib[t] = 11;
   for (int t=18; t<26; t++)  NTSib[t] =  7;

   for (int s=0; s<26; s++)   TSib[s] = new int [ NTSib[s] ];
   
   TSib[ 0][ 0] = 10;
   TSib[ 0][ 1] = 19;
   TSib[ 0][ 2] =  4;
   TSib[ 0][ 3] = 16;
   TSib[ 0][ 4] = 11;
   TSib[ 0][ 5] = 21;
   TSib[ 0][ 6] =  2;
   TSib[ 0][ 7] =  7;
   TSib[ 0][ 8] =  1;
   TSib[ 0][ 9] =  3;
   TSib[ 0][10] =  9;
   TSib[ 0][11] = 12;
   TSib[ 0][12] = 23;
   TSib[ 0][13] =  5;
   TSib[ 0][14] = 17;
   TSib[ 0][15] = 13;
   TSib[ 0][16] = 25;
   
   TSib[ 1][ 0] = 18;
   TSib[ 1][ 1] = 10;
   TSib[ 1][ 2] = 14;
   TSib[ 1][ 3] =  4;
   TSib[ 1][ 4] = 20;
   TSib[ 1][ 5] = 11;
   TSib[ 1][ 6] =  6;
   TSib[ 1][ 7] =  2;
   TSib[ 1][ 8] =  0;
   TSib[ 1][ 9] =  8;
   TSib[ 1][10] =  3;
   TSib[ 1][11] = 22;
   TSib[ 1][12] = 12;
   TSib[ 1][13] = 15;
   TSib[ 1][14] =  5;
   TSib[ 1][15] = 24;
   TSib[ 1][16] = 13;
   
   TSib[ 2][ 0] = 14;
   TSib[ 2][ 1] =  4;
   TSib[ 2][ 2] = 16;
   TSib[ 2][ 3] = 20;
   TSib[ 2][ 4] = 11;
   TSib[ 2][ 5] = 21;
   TSib[ 2][ 6] =  0;
   TSib[ 2][ 7] =  1;
   TSib[ 2][ 8] =  8;
   TSib[ 2][ 9] =  3;
   TSib[ 2][10] =  9;
   TSib[ 2][11] = 15;
   TSib[ 2][12] =  5;
   TSib[ 2][13] = 17;
   TSib[ 2][14] = 24;
   TSib[ 2][15] = 13;
   TSib[ 2][16] = 25;
   
   TSib[ 3][ 0] = 18;
   TSib[ 3][ 1] = 10;
   TSib[ 3][ 2] = 19;
   TSib[ 3][ 3] = 14;
   TSib[ 3][ 4] =  4;
   TSib[ 3][ 5] = 16;
   TSib[ 3][ 6] =  6;
   TSib[ 3][ 7] =  2;
   TSib[ 3][ 8] =  7;
   TSib[ 3][ 9] =  0;
   TSib[ 3][10] =  1;
   TSib[ 3][11] = 22;
   TSib[ 3][12] = 12;
   TSib[ 3][13] = 23;
   TSib[ 3][14] = 15;
   TSib[ 3][15] =  5;
   TSib[ 3][16] = 17;
   
   TSib[ 4][ 0] =  6;
   TSib[ 4][ 1] =  2;
   TSib[ 4][ 2] =  7;
   TSib[ 4][ 3] =  0;
   TSib[ 4][ 4] =  1;
   TSib[ 4][ 5] =  8;
   TSib[ 4][ 6] =  3;
   TSib[ 4][ 7] =  9;
   TSib[ 4][ 8] = 22;
   TSib[ 4][ 9] = 12;
   TSib[ 4][10] = 23;
   TSib[ 4][11] = 15;
   TSib[ 4][12] =  5;
   TSib[ 4][13] = 17;
   TSib[ 4][14] = 24;
   TSib[ 4][15] = 13;
   TSib[ 4][16] = 25;
   
   TSib[ 5][ 0] = 18;
   TSib[ 5][ 1] = 10;
   TSib[ 5][ 2] = 19;
   TSib[ 5][ 3] = 14;
   TSib[ 5][ 4] =  4;
   TSib[ 5][ 5] = 16;
   TSib[ 5][ 6] = 20;
   TSib[ 5][ 7] = 11;
   TSib[ 5][ 8] = 21;
   TSib[ 5][ 9] =  6;
   TSib[ 5][10] =  2;
   TSib[ 5][11] =  7;
   TSib[ 5][12] =  0;
   TSib[ 5][13] =  1;
   TSib[ 5][14] =  8;
   TSib[ 5][15] =  3;
   TSib[ 5][16] =  9;
   
   TSib[ 6][ 0] =  4;
   TSib[ 6][ 1] = 16;
   TSib[ 6][ 2] = 11;
   TSib[ 6][ 3] = 21;
   TSib[ 6][ 4] =  1;
   TSib[ 6][ 5] =  3;
   TSib[ 6][ 6] =  9;
   TSib[ 6][ 7] =  5;
   TSib[ 6][ 8] = 17;
   TSib[ 6][ 9] = 13;
   TSib[ 6][10] = 25;
   
   TSib[ 7][ 0] = 14;
   TSib[ 7][ 1] =  4;
   TSib[ 7][ 2] = 20;
   TSib[ 7][ 3] = 11;
   TSib[ 7][ 4] =  0;
   TSib[ 7][ 5] =  8;
   TSib[ 7][ 6] =  3;
   TSib[ 7][ 7] = 15;
   TSib[ 7][ 8] =  5;
   TSib[ 7][ 9] = 24;
   TSib[ 7][10] = 13;
   
   TSib[ 8][ 0] = 10;
   TSib[ 8][ 1] = 19;
   TSib[ 8][ 2] =  4;
   TSib[ 8][ 3] = 16;
   TSib[ 8][ 4] =  2;
   TSib[ 8][ 5] =  7;
   TSib[ 8][ 6] =  1;
   TSib[ 8][ 7] = 12;
   TSib[ 8][ 8] = 23;
   TSib[ 8][ 9] =  5;
   TSib[ 8][10] = 17;
   
   TSib[ 9][ 0] = 18;
   TSib[ 9][ 1] = 10;
   TSib[ 9][ 2] = 14;
   TSib[ 9][ 3] =  4;
   TSib[ 9][ 4] =  6;
   TSib[ 9][ 5] =  2;
   TSib[ 9][ 6] =  0;
   TSib[ 9][ 7] = 22;
   TSib[ 9][ 8] = 12;
   TSib[ 9][ 9] = 15;
   TSib[ 9][10] =  5;
   
   TSib[10][ 0] =  0;
   TSib[10][ 1] =  1;
   TSib[10][ 2] =  8;
   TSib[10][ 3] =  3;
   TSib[10][ 4] =  9;
   TSib[10][ 5] = 15;
   TSib[10][ 6] =  5;
   TSib[10][ 7] = 17;
   TSib[10][ 8] = 24;
   TSib[10][ 9] = 13;
   TSib[10][10] = 25;
   
   TSib[11][ 0] =  6;
   TSib[11][ 1] =  2;
   TSib[11][ 2] =  7;
   TSib[11][ 3] =  0;
   TSib[11][ 4] =  1;
   TSib[11][ 5] = 22;
   TSib[11][ 6] = 12;
   TSib[11][ 7] = 23;
   TSib[11][ 8] = 15;
   TSib[11][ 9] =  5;
   TSib[11][10] = 17;
   
   TSib[12][ 0] = 14;
   TSib[12][ 1] =  4;
   TSib[12][ 2] = 16;
   TSib[12][ 3] = 20;
   TSib[12][ 4] = 11;
   TSib[12][ 5] = 21;
   TSib[12][ 6] =  0;
   TSib[12][ 7] =  1;
   TSib[12][ 8] =  8;
   TSib[12][ 9] =  3;
   TSib[12][10] =  9;
   
   TSib[13][ 0] = 18;
   TSib[13][ 1] = 10;
   TSib[13][ 2] = 19;
   TSib[13][ 3] = 14;
   TSib[13][ 4] =  4;
   TSib[13][ 5] = 16;
   TSib[13][ 6] =  6;
   TSib[13][ 7] =  2;
   TSib[13][ 8] =  7;
   TSib[13][ 9] =  0;
   TSib[13][10] =  1;
   
   TSib[14][ 0] =  2;
   TSib[14][ 1] =  7;
   TSib[14][ 2] =  1;
   TSib[14][ 3] =  3;
   TSib[14][ 4] =  9;
   TSib[14][ 5] = 12;
   TSib[14][ 6] = 23;
   TSib[14][ 7] =  5;
   TSib[14][ 8] = 17;
   TSib[14][ 9] = 13;
   TSib[14][10] = 25;
   
   TSib[15][ 0] = 10;
   TSib[15][ 1] = 19;
   TSib[15][ 2] =  4;
   TSib[15][ 3] = 16;
   TSib[15][ 4] = 11;
   TSib[15][ 5] = 21;
   TSib[15][ 6] =  2;
   TSib[15][ 7] =  7;
   TSib[15][ 8] =  1;
   TSib[15][ 9] =  3;
   TSib[15][10] =  9;
   
   TSib[16][ 0] =  6;
   TSib[16][ 1] =  2;
   TSib[16][ 2] =  0;
   TSib[16][ 3] =  8;
   TSib[16][ 4] =  3;
   TSib[16][ 5] = 22;
   TSib[16][ 6] = 12;
   TSib[16][ 7] = 15;
   TSib[16][ 8] =  5;
   TSib[16][ 9] = 24;
   TSib[16][10] = 13;
   
   TSib[17][ 0] = 18;
   TSib[17][ 1] = 10;
   TSib[17][ 2] = 14;
   TSib[17][ 3] =  4;
   TSib[17][ 4] = 20;
   TSib[17][ 5] = 11;
   TSib[17][ 6] =  6;
   TSib[17][ 7] =  2;
   TSib[17][ 8] =  0;
   TSib[17][ 9] =  8;
   TSib[17][10] =  3;
   
   TSib[18][ 0] =  1;
   TSib[18][ 1] =  3;
   TSib[18][ 2] =  9;
   TSib[18][ 3] =  5;
   TSib[18][ 4] = 17;
   TSib[18][ 5] = 13;
   TSib[18][ 6] = 25;
   
   TSib[19][ 0] =  0;
   TSib[19][ 1] =  8;
   TSib[19][ 2] =  3;
   TSib[19][ 3] = 15;
   TSib[19][ 4] =  5;
   TSib[19][ 5] = 24;
   TSib[19][ 6] = 13;
   
   TSib[20][ 0] =  2;
   TSib[20][ 1] =  7;
   TSib[20][ 2] =  1;
   TSib[20][ 3] = 12;
   TSib[20][ 4] = 23;
   TSib[20][ 5] =  5;
   TSib[20][ 6] = 17;
   
   TSib[21][ 0] =  6;
   TSib[21][ 1] =  2;
   TSib[21][ 2] =  0;
   TSib[21][ 3] = 22;
   TSib[21][ 4] = 12;
   TSib[21][ 5] = 15;
   TSib[21][ 6] =  5;
   
   TSib[22][ 0] =  4;
   TSib[22][ 1] = 16;
   TSib[22][ 2] = 11;
   TSib[22][ 3] = 21;
   TSib[22][ 4] =  1;
   TSib[22][ 5] =  3;
   TSib[22][ 6] =  9;
   
   TSib[23][ 0] = 14;
   TSib[23][ 1] =  4;
   TSib[23][ 2] = 20;
   TSib[23][ 3] = 11;
   TSib[23][ 4] =  0;
   TSib[23][ 5] =  8;
   TSib[23][ 6] =  3;
   
   TSib[24][ 0] = 10;
   TSib[24][ 1] = 19;
   TSib[24][ 2] =  4;
   TSib[24][ 3] = 16;
   TSib[24][ 4] =  2;
   TSib[24][ 5] =  7;
   TSib[24][ 6] =  1;
   
   TSib[25][ 0] = 18;
   TSib[25][ 1] = 10;
   TSib[25][ 2] = 14;
   TSib[25][ 3] =  4;
   TSib[25][ 4] =  6;
   TSib[25][ 5] =  2;
   TSib[25][ 6] =  0;

} // FUNCTION : SetTargetSibling
