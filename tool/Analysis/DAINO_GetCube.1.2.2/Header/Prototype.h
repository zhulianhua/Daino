#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__



// DAINO_GetSlice functions
void ReadOption( int argc, char **argv );
void Init_MPI( int *argc, char ***argv );
void CheckParameter();
void TakeNote();
void Refine2TargetLevel();
void PreparePatch( const int lv, const int PID, const int Buffer, real FData[], const real CData[] );
void Interpolation( const int CSize, const real CData[], const int FSize, real FData[] );
void AllocateOutputArray();
void StoreData( const int lv, const int PID, real FData[], const int Buffer, real *Out );
void Output();
void SumOverRanks();
void Init_TargetDomain();
void GetCandidateBox();
bool WithinCandidateBox( const int *Corner, const int Size, const int Buf );
void GetDivVel( real FData[], const int FSize, const int Buffer, const int NextIdx, const real dh );
void GetCurlVel( real FData[], const int FSize, const int Buffer, const int NextIdx, const real dh );


// DAINO functions
void LoadData( const char *FileName );
//void LoadData( const char *FileName, const bool OldDataFormat );
//void LoadData_Old( const char *FileName );
void Init_MemAllocate();
void Init_RecordBasePatch();
void End_MemFree();
void Buf_AllocateBufferPatch( const int lv );
void Buf_AllocateBufferPatch_Base();
void Buf_RecordBoundaryPatch( const int lv );
void Buf_RecordBoundaryPatch_Base();
void Buf_SortBoundaryPatch( const int NPatch, int *IDList, int *PosList );
void Buf_GetBufferData( const int lv, const int Package, const int ParaBuffer );
void Buf_RecordExchangeDataPatchID( const int lv );
void FindFather( const int lv );
void SiblingSearch( const int lv );
void SiblingSearch_Base();
void Flu_Restrict( const int lv, const bool GetAvePot );
void MPI_Exit();
void MPI_ExchangeBufferPosition( int NSend[26], int NRecv[26], int *Send_PosList[26], int *Recv_PosList[26] );
void MPI_ExchangeInfo( const int TargetRank[2], const int SendSize[2], const int RecvSize[2],
                       real *SendBuffer[2], real *RecvBuffer[2] );
int  TABLE_01( const int SibIndex, const char dim, const int w0, const int w1, const int w2 );
int  TABLE_02( const int LocalID, const char dim, const int w0, const int w1 );
int  TABLE_03( const int SibID, const int Count );
int  TABLE_04( const int SibID );
int  TABLE_05( const int SibID );
int  TABLE_07( const int SibID, const int Count );



#endif // __PROTOTYPE_H__
