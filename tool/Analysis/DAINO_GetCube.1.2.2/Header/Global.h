#ifndef __GLOBAL_H__
#define __GLOBAL_H__



extern int        CanMin[3][3], CanMax[3][3], CanBuf;                    
extern DAINO_t    patch;
extern ParaVar_t  ParaVar;
extern real       GAMMA;        
extern int        NX0_TOT[3], NX0[3], MyRank, MyRank_X[3], SibRank[26], NGPU, NGPU_X[3], TargetLevel, DumpID;
extern int        *BaseP, *BounP_IDMap[NLEVEL][26], SendP_NList[NLEVEL][26], *SendP_IDList[NLEVEL][26];  
extern int        RecvP_NList[NLEVEL][26], *RecvP_IDList[NLEVEL][26], NPatchComma[NLEVEL][28];
extern int        TargetX[3], X[3], NLoad, NOut;
extern double     Time[NLEVEL];
extern long int   Step;
extern bool       OutputPot;



#endif // #ifndef __GLOBAL_H__
