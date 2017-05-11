
#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__



#include "Macro.h"
#include "Typedef.h"
#include "AMR.h"

#ifdef OOC
#include "OOC.h"
#endif


// Auxiliary
void Aux_Check_MemFree( const real MinMemFree_Total, const char *comment );
void Aux_Check_Conservation( const bool Output2File, const char *comment );
void Aux_Check();
void Aux_Check_Finite( const int lv, const char *comment );
void Aux_Check_FluxAllocate( const int lv, const char *comment );
void Aux_Check_Parameter();
void Aux_Check_PatchAllocate( const int lv, const char *comment );
void Aux_Check_ProperNesting( const int lv, const char *comment );
void Aux_Check_Refinement( const int lv, const char *comment );
void Aux_Check_Restrict( const int lv, const char *comment );
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void Aux_GetCPUInfo( const char *FileName );
void Aux_GetMemInfo();
void Aux_Message( FILE *Type, const char *Format, ... );
void Aux_PatchCount();
void Aux_TakeNote();
void Aux_RecordTiming();
void Aux_CreateTimer();
void Aux_DeleteTimer();
void Aux_ResetTimer();
#ifndef SERIAL
void Aux_RecordBoundaryPatch( const int lv, int *NList, int **IDList, int **PosList );
#endif


// Buffer
#ifndef SERIAL
void Buf_AllocateBufferPatch_Base( AMR_t *Tpatch, const int OOC_MyRank );
void Buf_AllocateBufferPatch( AMR_t *Tpatch, const int lv, const int Mode, const int OOC_MyRank );
void Buf_GetBufferData( const int lv, const int FluSg, const int PotSg, const GetBufMode_t GetBufMode, 
                        const int TVar, const int ParaBuffer, const UseLBFunc_t UseLBFunc );
void Buf_RecordBoundaryFlag( const int lv );
void Buf_RecordBoundaryPatch_Base();
void Buf_RecordBoundaryPatch( const int lv );
void Buf_RecordExchangeDataPatchID( const int lv );
void Buf_RecordExchangeFluxPatchID( const int lv );
void Buf_ResetBufferFlux( const int lv );
void Buf_SortBoundaryPatch( const int NPatch, int *IDList, int *PosList );
#endif // #ifndef SERIAL


// Hydrodynamics
void CPU_FluidSolver( real h_Flu_Array_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                      real h_Flu_Array_Out[][FLU_NOUT][ PS2*PS2*PS2 ], 
                      real h_Flux_Array[][9][NCOMP   ][ PS2*PS2 ], 
                      real h_MinDtInfo_Array[],
                      const int NPatchGroup, const real dt, const real dh, const real Gamma, const bool StoreFlux,
                      const bool XYZ, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const real EP_Coeff,
                      const WAF_Limiter_t WAF_Limiter, const real Eta, const bool GetMinDtInfo );
void Flu_AdvanceDt( const int lv, const double PrepTime, const double dt, const int SaveSg,
                    const bool OverlapMPI, const bool Overlap_Sync );
void Flu_AllocateFluxArray( const int lv );
void Flu_Close( const int lv, const int SaveSg, const real h_Flux_Array[][9][NCOMP][4*PATCH_SIZE*PATCH_SIZE],
                const real h_Flu_Array_F_Out[][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE], 
                const real h_MinDtInfo_Array[], const int NPG, const int *PID0_List, const bool GetMinDtInfo );
void Flu_FixUp( const int lv, const double dt );
void Flu_Prepare( const int lv, const double PrepTime, real h_Flu_Array_F_In[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT], 
                  const int NPG, const int *PID0_List );
void Flu_Restrict( const int FaLv, const int SonFluSg, const int FaFluSg, const int SonPotSg, const int FaPotSg,
                   const int TVar );
#ifndef SERIAL
void Flu_AllocateFluxArray_Buffer( const int lv );
#endif


// DAINO
void Integration_IndiviTimeStep( const int lv, const double dTime );
void InvokeSolver( const Solver_t TSolver, const int lv, const double PrepTime, const double dt,
                   const real Poi_Coeff, const int SaveSg, const bool OverlapMPI, const bool Overlap_Sync );
void Prepare_PatchGroupData( const int lv, const double PrepTime, real *h_Input_Array, const int GhostSize, 
                             const int NPG, const int *PID0_List, const int TVar, const IntScheme_t IntScheme,
                             const PrepUnit_t PrepUnit, const NSide_t NSide, const bool IntPhase );


// Init
void End_DAINO();
void End_MemFree();
void End_MemFree_Fluid();
void End_StopManually( int &Terminate_global );
void Init_BaseLevel();
void Init_DAINO( int *argc, char ***argv );
void Init_Load_DumpTable();
void Init_Load_FlagCriteria();
void Init_Load_Parameter();
void Init_MemAllocate();
void Init_MemAllocate_Fluid( const int Flu_NPatchGroup );
void Init_Parallelization();
void Init_RecordBasePatch();
void Init_Refine( const int lv );
void Init_Reload();
void Init_Reload_OldFormat();
void Init_StartOver();
void Init_TestProb();
void Init_UM();


// Interpolation
void Int_Table( const IntScheme_t IntScheme, int &NSide, int &NGhost );
void Interpolate( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
                  real FData [], const int FSize[3], const int FStart[3], 
                  const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase, 
                  const bool EnsurePositivity );


// Miscellaneous
void   Mis_Idx1D2Idx3D( const int Size[], const int  Idx1D, int Idx3D[] );
void   Mis_Idx1D2Idx3D( const int Size[], const long Idx1D, int Idx3D[] );
long   Mis_Idx3D2Idx1D( const int Size[], const int Idx3D[] );
int    Mis_BinarySearch( const int  Array[], int Min, int Max, const  int Key );
int    Mis_BinarySearch( const long Array[], int Min, int Max, const long Key );
bool   Mis_Check_Synchronization( const double Time1, const double Time2, const char *comment, 
                                  const bool Verbose );
void   Mis_GetTimeStep();
void   Mis_GetTimeStep_UserCriteria( double &dt, double &dTime, const double dt_dTime );
double Mis_dTime2dt( const double Time_In, const double dTime_In );
void   Mis_GetTotalPatchNumber( const int lv );
void   Mis_Heapsort( const int N, int  Array[], int IdxTable[] );
void   Mis_Heapsort( const int N, long Array[], int IdxTable[] );
int    Mis_Matching( const int N, const int  Array[], const int M, const int  Key[], char Match[] );
int    Mis_Matching( const int N, const long Array[], const int M, const long Key[], char Match[] );
int    Mis_Matching( const int N, const int  Array[], const int M, const int  Key[], int  Match[] );
int    Mis_Matching( const int N, const long Array[], const int M, const long Key[], int  Match[] );


// MPI
#ifndef SERIAL
void Init_MPI( int *argc, char ***argv );
void MPI_ExchangeBoundaryFlag( const int lv );
void MPI_ExchangeBufferPosition( int NSend[26], int NRecv[26], int *Send_PosList[26], int *Recv_PosList[26] );
void MPI_ExchangeData( const int TargetRank[2], const int SendSize[2], const int RecvSize[2],
                       real *SendBuffer[2], real *RecvBuffer[2] );
void MPI_Exit();
#endif // #ifndef SERIAL


// Output
void Output_DumpData( const int Stage );
void Output_DumpData_Part( const OptOutputPart_t Part, const bool BaseOnly, const real x, const real y, 
                           const real z, const char *FileName );
void Output_DumpData_Total( const char *FileName );
void Output_DumpManually( int &Dump_global );
void Output_FlagMap( const int lv, const int xyz, const char *comment );
void Output_Flux( const int lv, const int PID, const int Sib, const char *comment );
void Output_PatchCorner( const int lv, const char *comment );
void Output_Patch( const int lv, const int PID, const int FluSg, const int PotSg, const char *comment );
void Output_PatchMap( const int lv, const int PID, const int TSg, const int Comp, const char *comment );
void Output_PreparedPatch_Fluid( const int TLv, const int TPID, 
                                 const real h_Flu_Array[][FLU_NIN][FLU_NXT*FLU_NXT*FLU_NXT], 
                                 const int NPG, const int *PID0_List, const int CLv, const char *comment );
void Output_TestProbErr( const bool BaseOnly );
void Output_BasePowerSpectrum( const char *FileName );
#ifndef SERIAL
void Output_ExchangePatchMap( const int lv, const int xyz, const char *comment );
void Output_ExchangeFluxPatchList( const int option, const int lv, const char *comment );
void Output_ExchangeDataPatchList( const int option, const int lv, const char *comment );
void Output_BoundaryFlagList( const int option, const int lv, const char *comment );
#endif


// Refine
void FindFather( const int lv, const int Mode );
void Flag_Real( const int lv, const UseLBFunc_t UseLBFunc );
bool Flag_Check( const int lv, const int PID, const int i, const int j, const int k, 
                 const real Fluid[][PS1][PS1][PS1], const real Pot[][PS1][PS1], const real Pres[][PS1][PS1],
                 const real *Lohner_Var, const real *Lohner_Slope, const int Lohner_NCell, const int Lohner_NVar );
bool Flag_UserCriteria( const int i, const int j, const int k, const int lv, const int PID, const real Threshold);
bool Flag_Lohner( const int i, const int j, const int k, const real *Var1D, const real *Slope1D, const int NCell, 
                  const int NVar, const double Threshold, const double Filter, const double Soften );
void Refine( const int lv );
void SiblingSearch( const int lv );
void SiblingSearch_Base();
#ifndef SERIAL
void Flag_Buffer( const int lv );
void Refine_Buffer( const int lv, const int *SonTable, const int *GrandTable );
#endif


// SelfGravity
#ifdef GRAVITY
void CPU_PoissonGravitySolver( const real h_Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT], 
                                     real h_Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                     real h_Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                                     real h_Flu_Array    [][GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE], 
                               const int NPatchGroup, const real dt, const real dh, const int SOR_Min_Iter, 
                               const int SOR_Max_Iter, const real SOR_Omega, const int MG_Max_Iter,
                               const int MG_NPre_Smooth, const int MG_NPost_Smooth, const real MG_Tolerated_Error,
                               const real Poi_Coeff, const IntScheme_t IntScheme, const bool P5_Gradient,
                               const real Eta, const bool Poisson, const bool GraAcc );
void CPU_PoissonSolver_FFT( const real Poi_Coeff, const int SaveSg );
void Cube_to_Slice( real *RhoK, real *SendBuf, real *RecvBuf );
void Slice_to_Cube( real *RhoK, real *SendBuf, real *RecvBuf, const int SaveSg );
void End_MemFree_PoissonGravity();
void Gra_AdvanceDt( const int lv, const double PrepTime, const double dt, const int SaveSg, const bool GraAcc,
                    const bool OverlapMPI, const bool Overlap_Sync );
void Gra_Close( const int lv, const int SaveSg, 
                const real h_Flu_Array_G[][GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE], 
                const int NPG, const int *PID0_List );
void Gra_Prepare_Flu( const int lv, real h_Flu_Array_G[][GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE],
                      const int NPG, const int *PID0_List );
void Gra_Prepare_Pot( const int lv, const double PrepTime, real h_Pot_Array_P_Out[][GRA_NXT][GRA_NXT][GRA_NXT], 
                      const int NPG, const int *PID0_List );
void End_FFTW();
void Init_FFTW();
void Init_MemAllocate_PoissonGravity( const int Pot_NPatchGroup );
void Init_Set_Default_MG_Parameter( int &Max_Iter, int &NPre_Smooth, int &NPost_Smooth, real &Tolerated_Error );
void Init_Set_Default_SOR_Parameter( real &SOR_Omega, int &SOR_Max_Iter, int &SOR_Min_Iter );
void Output_PreparedPatch_Poisson( const int TLv, const int TPID, const int TComp,
                                   const real h_Rho_Array_P   [][RHO_NXT][RHO_NXT][RHO_NXT],
                                   const real h_Pot_Array_P_In[][POT_NXT][POT_NXT][POT_NXT], 
                                   const int NPG, const int *PID0_List, const int CLv, const char *comment );
void Poi_Close( const int lv, const int SaveSg, const real h_Pot_Array_P_Out[][GRA_NXT][GRA_NXT][GRA_NXT], 
                const int NPG, const int *PID0_List );
void Poi_GetAverageDensity();
void Poi_Prepare_Pot( const int lv, const double PrepTime, real h_Pot_Array_P_In[][POT_NXT][POT_NXT][POT_NXT], 
                      const int NPG, const int *PID0_List );
void Poi_Prepare_Rho( const int lv, const double PrepTime, real h_Rho_Array_P[][RHO_NXT][RHO_NXT][RHO_NXT], 
                      const int NPG, const int *PID0_List );
#ifndef SERIAL
void MPI_CubeSlice( const bool Dir, real *SendBuf, real *RecvBuf );
#endif 
#endif // #ifdef GRAVITY


// Tables
int TABLE_01( const int SibIndex, const char dim, const int w0, const int w1, const int w2 );
int TABLE_02( const int LocalID, const char dim, const int w0, const int w1 );
int TABLE_03( const int SibID, const int Count );
int TABLE_04( const int SibID );
int TABLE_05( const int SibID );
int TABLE_06( const int SibID, const int FlagLayer );
int TABLE_07( const int SibID, const int Count );


// LoadBalance
#ifdef LOAD_BALANCE
void LB_AllocateBufferPatch_Father( const int SonLv, 
                                    const bool SearchAllSon, const int NInput, int* TargetSonPID0, 
                                    const bool RecordFaPID, int* NNewFaBuf0, int** NewFaBufPID0 );
void LB_AllocateBufferPatch_Sibling_Base();
void LB_AllocateBufferPatch_Sibling( const int lv );
void LB_AllocateFluxArray( const int FaLv );
void LB_ExchangeFlaggedBuffer( const int lv );
void LB_FindFather( const int SonLv, const bool SearchAllSon, const int NInput, int* TargetSonPID0 );
void LB_FindSonNotHome( const int FaLv, const bool SearchAllFa, const int NInput, int* TargetFaPID );
void LB_GetBufferData( const int lv, const int FluSg, const int PotSg, const GetBufMode_t GetBufMode, 
                       const int TVar, const int ParaBuf );
void LB_GrandsonCheck( const int lv );
void LB_Init_LoadBalance( const bool DuringRestart );
void LB_SetCutPoint( const int lv, long *CutPoint, const bool InputLBIdxList, long *LBIdx_AllRank );
void LB_Output_LBIdx( const int lv );
void LB_RecordExchangeDataPatchID( const int Lv, const bool AfterRefine );
void LB_RecordExchangeFixUpDataPatchID( const int Lv );
void LB_RecordExchangeRestrictDataPatchID( const int FaLv );
void LB_RecordOverlapMPIPatchID( const int Lv );
void LB_Refine( const int FaLv );
void LB_SFC_to_Slice( real *RhoK, real *SendBuf_Rho, real *RecvBuf_Rho, long *SendBuf_SIdx, long *RecvBuf_SIdx,
                      int **List_PID, int **List_k, int *List_NSend_Rho, int *List_NRecv_Rho, 
                      const int *List_z_start, const int local_nz );
void LB_SiblingSearch( const int lv, const bool SearchAllPID, const int NInput, int *TargetPID0 );
void LB_Slice_to_SFC( const real *RhoK, real *SendBuf, real *RecvBuf, const int SaveSg, const long *List_SIdx, 
                      int **List_PID, int **List_k, int *List_NSend, int *List_NRecv, const int local_nz );
long LB_Corner2Index( const int lv, const int Corner[], const Check_t Check );
void LB_Index2Corner( const int lv, const long Index, int Corner[], const Check_t Check );
int  LB_Index2Rank( const int lv, const long LB_Idx, const Check_t Check );
#endif // #ifdef LOAD_BALANCE


// Hydro model
#if    ( MODEL == HYDRO )
void Hydro_Aux_Check_Negative( const int lv, const int Mode, const char *comment );
void Hydro_GetTimeStep_Fluid( double &dt, double &dTime, int &MinDtLv, real MinDtVar[], const double dt_dTime );
void Hydro_GetTimeStep_Gravity( double &dt, double &dTime, int &MinDtLv, real &MinDtVar, const double dt_dTime );
void Hydro_GetMaxCFL( real MaxCFL[], real MinDtVar_AllLv[][NCOMP] );
void Hydro_GetMaxAcc( real MaxAcc[] );
void Hydro_Init_StartOver_AssignData( const int lv );
void Hydro_Init_UM_AssignData( const int lv, const real *UM_Data, const int NVar );


// MHD model
#elif ( MODEL == MHD )
#warning : WAIT MHD !!!


// ELBDM model
#elif ( MODEL == ELBDM )
void ELBDM_Init_StartOver_AssignData( const int lv );
void ELBDM_Init_UM_AssignData( const int lv, const real *UM_Data, const int NVar );
void ELBDM_GetTimeStep_Fluid( double &dt, double &dTime, int &MinDtLv, const double dt_dTime );
void ELBDM_GetTimeStep_Gravity( double &dt, double &dTime, int &MinDtLv, real &MinDtVar, const double dt_dTime );
void ELBDM_GetTimeStep_Phase( double &dt, double &dTime, int &MinDtLv, real *MinDtVar, const double dt_dTime );
void ELBDM_GetMaxPot( real MaxPot[] );
void ELBDM_GetMaxPhaseDerivative( real MaxdS_dt[], real MinDtVar_AllLv[][3] );
bool ELBDM_Flag_EngyDensity( const int i, const int j, const int k, const real Real_Array[], 
                             const real Imag_Array[], const double Angle_2pi, const double Eps );
real ELBDM_UnwrapPhase( const real Phase_Ref, const real Phase_Wrapped );


#else
#error : ERROR : unsupported MODEL !!
#endif 


// OOC
#ifdef OOC
void  Init_OOC();
void  OOC_SwapPointer( AMR_t * &patch1, AMR_t * &patch2 );
void  OOC_ExchangeBoundaryFlag( AMR_t *Tpatch, const int lv, const bool Dir, const int OOC_MyRank );
void  OOC_MPI_ExchangeBoundaryFlag( const int lv );
void  OOC_MPI_GetBufferData( const int lv, const int Package, const int ParaBuffer );
void  OOC_GetAverageDensity();
void  OOC_UpdateBuffer( const int lv, const int Mode, const int FluSg, const int PotSg, const bool Sorted );
void  OOC_SortPatch( const AMR_t *Tpatch, TargetP_t *TargetP, const int lv );
void  OOC_Init_SortPatch();
void  OOC_SetTargetP( TargetP_t *TargetP, const uint R_Start, const uint R_NPatch, uint *R_Table,
                     const uint B_Start = 0, const uint B_NPatch = 0, uint *B_Table = NULL );

int   OOC_Dump( const AMR_t *amr, const int rank, const int lv, const int mode, const int Sg, 
                const TargetP_t *list );
int   OOC_Load( AMR_t *amr, const int rank, const int lv, const int mode, const int Sg, 
                const TargetP_t *list );
int   OOC_Send( const void *ptr, const ulong size, const int S_rank, const int D_rank, const int tag );
int   OOC_Recv( void *ptr, const ulong size, const int S_rank, const int D_rank, const int tag );
void  OOC_Initialize( const int NRank );
void  Final_OOC();
void  OOC_partialDump( OOC_info * info, const void* ptr, const ulong start, const ulong size );
void  OOC_partialLoad( const OOC_info *info, void * ptr, const ulong start, const ulong size );
void  OOC_Release( OOC_info *&info );
void* OOC_Load( const OOC_info *info );
OOC_info* OOC_Dump( const void *ptr, const ulong size );

void OOC_Init_StartOver_ConstructAllLevels( const int lv );
void OOC_Init_StartOver_Restrict( const int lv );
void OOC_Init_Reload_LoadData( const int lv, const int PatchDataSize, const long int Offset, 
                               const bool DataOrder_xyzv, real (*InvData_Flu)[PATCH_SIZE][PATCH_SIZE][NCOMP],
                               const bool opt__output_pot );
void OOC_Init_Reload_FindFather();
void OOC_Init_Reload_ConstructAllLevels( const int lv );
void OOC_Init_Reload_Restrict( const int lv );
void OOC_Integration_IndiviTimeStep( const int lv, const double dt, const double dTime );
void OOC_Integration_IndiviTimeStep_FluAdvanceDt( const int lv, const double dt, const double dTime );
void OOC_Integration_IndiviTimeStep_FixUp( const int lv, const double dt );
void OOC_Integration_IndiviTimeStep_Flag( const int lv, const double dt );
void OOC_Integration_IndiviTimeStep_Refine( const int lv );
void OOC_Mis_GetMaxCFL( const int lv, real MaxCFL[], real MaxVal[] );
void OOC_Output_DumpData_Total( const int lv, FILE *File, real (*InvData_Flu)[PATCH_SIZE][PATCH_SIZE][NCOMP] );
void OOC_Output_DumpData_Part( const OptOutputPart_t Part, const bool BaseOnly, const real x, const real y, 
                               const real z, const int lv, FILE *File );
void OOC_GetBufferData( AMR_t *Tpatch, const int lv, const int Sg, const int Package, const int ParaBuffer, 
                        const int Mode, const int OOC_MyRank );
#ifdef GRAVITY
void OOC_Init_DAINO_GetPot( const int lv );
void OOC_Integration_IndiviTimeStep_GraAdvanceDt( const int lv, const double dt, const double dTime );
void OOC_Mis_GetMaxAcc( const int lv, const real Coeff, real (*Acc_Array)[GRA_NXT][GRA_NXT][GRA_NXT], 
                        real MaxAcc [] );
void OOC_Gra_AdvanceDt_InvokeSolver( const int lv, const double PrepTime, const double dt, const int SaveSg );
void OOC_CPU_PoissonSolver_FFT_Cube_to_Slice( real *SendBuf );
void OOC_CPU_PoissonSolver_FFT_Slice_to_Cube( real *RecvBuf, const int SaveSg );
#endif
#endif // ifdef OOC


// GPU API
#ifdef GPU
void CUAPI_Asyn_FluidSolver( real h_Flu_Array_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                             real h_Flu_Array_Out[][FLU_NOUT][ PS2*PS2*PS2 ],
                             real h_Flux_Array[][9][NCOMP   ][ PS2*PS2 ], 
                             real h_MinDtInfo_Array[],
                             const int NPatchGroup, const real dt, const real dh, const real Gamma, 
                             const bool StoreFlux, const bool XYZ, const LR_Limiter_t LR_Limiter, 
                             const real MinMod_Coeff, const real EP_Coeff, const WAF_Limiter_t WAF_Limiter,
                             const real Eta, const bool GetMinDtInfo, const int GPU_NStream );
void CUAPI_DiagnoseDevice();
void CUAPI_MemAllocate_Fluid( const int Flu_NPatchGroup, const int GPU_NStream );
void CUAPI_MemFree_Fluid( const int GPU_NStream );
void CUAPI_Set_Default_GPU_Parameter( int &GPU_NStream, int &FLU_GPU_NPGroup, int &POT_GPU_NPGroup );
void CUAPI_SetDevice( const int Mode );
void CUAPI_Synchronize();
#ifdef GRAVITY
void CUAPI_Asyn_PoissonGravitySolver( const real h_Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT], 
                                            real h_Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                            real h_Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                                            real h_Flu_Array    [][GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE], 
                                      const int NPatchGroup, const real dt, const real dh, const int SOR_Min_Iter,
                                      const int SOR_Max_Iter, const real SOR_Omega, const int MG_Max_Iter,
                                      const int MG_NPre_Smooth, const int MG_NPost_Smooth, 
                                      const real MG_Tolerated_Error, const real Poi_Coeff,
                                      const IntScheme_t IntScheme, const bool P5_Gradient, const real Eta,
                                      const bool Poisson, const bool GraAcc, const int GPU_NStream );
void CUAPI_MemAllocate_PoissonGravity( const int Pot_NPatchGroup );
void CUAPI_MemFree_PoissonGravity();
#endif // #ifdef GRAVITY
#endif // #ifdef GPU



#endif // __PROTOTYPE_H__
