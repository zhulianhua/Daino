
#include "DAINO.h"
#include "CUFLU.h"

#ifdef GPU



#if   ( MODEL == HYDRO )
#if   ( FLU_SCHEME == RTVD )
__global__ void CUFLU_FluidSolver_RTVD( real g_Fluid_In [][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                        real g_Fluid_Out[][5][ PS2*PS2*PS2 ],
                                        real g_Flux[][9][5][ PS2*PS2 ], 
                                        const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                        const bool XYZ );
#elif ( FLU_SCHEME == WAF )
__global__ void CUFLU_FluidSolver_WAF( real g_Fluid_In [][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                       real g_Fluid_Out[][5][ PS2*PS2*PS2 ],
                                       real g_Flux[][9][5][ PS2*PS2 ], 
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const bool XYZ, const WAF_Limiter_t WAF_Limiter );
#elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
__global__ void CUFLU_FluidSolver_MHM( const real g_Fluid_In[][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                       real g_Fluid_Out  [][5][ PS2*PS2*PS2 ],
                                       real g_Flux    [][9][5][ PS2*PS2 ], 
                                       real g_PriVar     [][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Slope_PPM_x[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_y[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_z[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_FC_Var_xL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_xR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_zR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Flux_x  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_y  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_z  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, 
                                       const real EP_Coeff );
#elif ( FLU_SCHEME == CTU )
__global__ void CUFLU_FluidSolver_CTU( const real g_Fluid_In[][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                       real g_Fluid_Out  [][5][ PS2*PS2*PS2 ],
                                       real g_Flux    [][9][5][ PS2*PS2 ], 
                                       real g_PriVar     [][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Slope_PPM_x[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_y[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_z[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_FC_Var_xL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_xR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_zR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Flux_x  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_y  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_z  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, 
                                       const real EP_Coeff );
#endif
//__global__ void CUFLU_GetMaxCFL( real g_Fluid[][5][ PS2*PS2*PS2 ], real g_MaxCFL[], const real Gamma );
#elif ( MODEL == MHD )
#warning : WAIT MHD !!!

#elif ( MODEL == ELBDM )
__global__ void CUFLU_ELBDMSolver( real g_Fluid_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                   real g_Fluid_Out[][FLU_NOUT][ PS2*PS2*PS2 ],
                                   const real dt, const real _dh, const real Eta, const bool XYZ );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL


// device pointers
extern real (*d_Flu_Array_F_In )[FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ];
extern real (*d_Flu_Array_F_Out)[FLU_NOUT][ PS2*PS2*PS2 ];
extern real (*d_Flux_Array)[9][5][ PS2*PS2 ];
extern real  *d_MinDtInfo_Fluid_Array;
#if ( MODEL == HYDRO )
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
extern real (*d_PriVar)     [5][ FLU_NXT*FLU_NXT*FLU_NXT ];
extern real (*d_Slope_PPM_x)[5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ];
extern real (*d_Slope_PPM_y)[5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ];
extern real (*d_Slope_PPM_z)[5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ];
extern real (*d_FC_Var_xL)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_xR)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_yL)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_yR)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_zL)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_zR)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Flux_x)  [5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ];
extern real (*d_FC_Flux_y)  [5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ];
extern real (*d_FC_Flux_z)  [5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ];
#endif // FLU_SCHEME
#elif ( MODEL == MHD )
#warning : WAIT MHD !!!
#endif // MODEL

extern cudaStream_t *Stream;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Asyn_FluidSolver
// Description :  1. MODEL == HYDRO : use GPU to solve the Euler equations by different schemes
//                                    --> invoke the kernel "CUFLU_FluidSolver_XXX"
//                2. MODEL == ELBDM : use GPU to solve the kinematic operator in the Schrodinger's equations 
//                                    --> invoke the kernel "CUFLU_ELBDMSolver"
//           
//                ***********************************************************
//                **                Asynchronous Function                  **
//                **                                                       ** 
//                **  will return before the execution in GPU is complete  **
//                ***********************************************************
//
// Note        :  1. Use streams for the asychronous memory copy between device and host
//                2. Prefix "d" : for pointers pointing to the "Device" memory space
//                   Prefix "h" : for pointers pointing to the "Host"   memory space
//                3. Use the input pamameter "XYZ" to control the order of update for dimensional-splitting
//                   methods (RTVD/WAF)
//                4. Currently five hydro schemes are supported : 
//                   1. Relaxing TVD scheme                            (RTVD  ) -->   split
//                   2. Weighted-Average-Flux scheme                   (WAF   ) -->   split
//                   3. MUSCL-Hancock scheme                           (MHM   ) --> unsplit
//                   4. MUSCL-Hancock scheme with Riemann prediction   (MHM_RP) --> unsplit
//                   5. Corner-Transport-Upwind scheme                 (CTU   ) --> unsplit
//
// Parameter   :  h_Flu_Array_In    : Host array to store the input variables
//                h_Flu_Array_Out   : Host array to store the output variables
//                h_Flux_Array      : Host array to store the output fluxes
//                h_MinDtInfo_Array : Host array to store the minimum time-step information in each patch group
//                                    --> useful only if "GetMinDtInfo == true"
//                                    --> not supported yet
//                NPatchGroup       : Number of patch groups evaluated simultaneously by GPU 
//                dt                : Time interval to advance solution
//                dh                : Grid size
//                Gamma             : Ratio of specific heats
//                XYZ               : true  : x->y->z ( forward sweep)
//                                    false : z->y->x (backward sweep)
//                                    ~ useless in directionally unsplit schemes
//                LR_Limiter        : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                    (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                   vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff      : Coefficient of the generalized MinMod limiter
//                EP_Coeff          : Coefficient of the extrema-preserving limiter
//                WAF_Limiter       : Flux limiter for the WAF scheme
//                                    (0/1/2/3) = (SuperBee/vanLeer/vanAlbada/MinBee)
//                StoreFlux         : true --> store the coarse-fine fluxes
//                Eta               : Particle mass / Planck constant
//                GetMinDtInfo      : true --> Invoke the kernel "CUFLU_GetMinDtInfo_Fluid" to gather the
//                                             minimum time-step information (the CFL condition in HYDRO )
//                                             in each patch group
//                                         --> NOT supported yet
//                GPU_NStream       : Number of CUDA streams for the asynchronous memory copy
//
// Useless parameters in HYDRO : Eta
// Useless parameters in ELBDM : h_Flux_Array, Gamma, LR_Limiter, MinMod_Coeff, EP_Coeff, WAF_Limite
//-------------------------------------------------------------------------------------------------------
void CUAPI_Asyn_FluidSolver( real h_Flu_Array_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                             real h_Flu_Array_Out[][FLU_NOUT][ PS2*PS2*PS2 ],
                             real h_Flux_Array[][9][NCOMP   ][ PS2*PS2 ], 
                             real h_MinDtInfo_Array[],
                             const int NPatchGroup, const real dt, const real dh, const real Gamma, 
                             const bool StoreFlux, const bool XYZ, const LR_Limiter_t LR_Limiter, 
                             const real MinMod_Coeff, const real EP_Coeff, const WAF_Limiter_t WAF_Limiter,
                             const real Eta, const bool GetMinDtInfo, const int GPU_NStream )
{

// check
#  if   ( MODEL == HYDRO )
   if ( LR_Limiter != VANLEER  &&  LR_Limiter != GMINMOD  &&  LR_Limiter != ALBADA  &&  LR_Limiter != EXTPRE  &&
        LR_Limiter != VL_GMINMOD  &&  LR_Limiter != LR_LIMITER_NONE )
      Aux_Error( ERROR_INFO, "unsupported limiter (%d) !!\n", LR_Limiter );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   if ( StoreFlux )  
      Aux_Message( stderr, "WARNING : currently the option \"%s\" in \"%s\" for ELBDM is useless !!\n",
                   "StoreFlux", __FUNCTION__ ); 

#  else
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif


   const real _dh = (real)1.0/dh;
   const dim3 BlockDim_FluidSolver ( FLU_BLOCK_SIZE_X, FLU_BLOCK_SIZE_Y, 1 ); // for the fluidsolvers
#  if   ( MODEL == HYDRO )
// const dim3 BlockDim_GetMinDtInfo( PS2, PS2, 1 );                           // for the time-step estimation
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif // MODEL

   int *NPatch_per_Stream  = new int [GPU_NStream];
   int *UsedPatch          = new int [GPU_NStream];
   int *Flu_MemSize_In     = new int [GPU_NStream];
   int *Flu_MemSize_Out    = new int [GPU_NStream];
   int *Flux_MemSize       = new int [GPU_NStream];


// set the number of patches of each stream
   UsedPatch[0] = 0;

   if ( GPU_NStream == 1 )    NPatch_per_Stream[0] = NPatchGroup;
   else
   {
      for (int s=0; s<GPU_NStream-1; s++)    
      {
         NPatch_per_Stream[s] = NPatchGroup / GPU_NStream;
         UsedPatch[s+1] = UsedPatch[s] + NPatch_per_Stream[s];
      }

      NPatch_per_Stream[GPU_NStream-1] = NPatchGroup - UsedPatch[GPU_NStream-1];
   }


// set the size of data to be transferred into GPU in each stream
   for (int s=0; s<GPU_NStream; s++)
   {
      Flu_MemSize_In [s] = sizeof(real)*NPatch_per_Stream[s]*FLU_NIN *FLU_NXT*FLU_NXT*FLU_NXT;
      Flu_MemSize_Out[s] = sizeof(real)*NPatch_per_Stream[s]*FLU_NOUT*PS2*PS2*PS2;
      Flux_MemSize   [s] = sizeof(real)*NPatch_per_Stream[s]*NCOMP*9*PS2*PS2;
   }


// a. copy data from host to device
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Flu_Array_F_In + UsedPatch[s], h_Flu_Array_In + UsedPatch[s],  
                         Flu_MemSize_In[s], cudaMemcpyHostToDevice, Stream[s] )  );
   } // for (int s=0; s<GPU_NStream; s++)


// b. execute the kernel 
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

#     if   ( MODEL == HYDRO )

#        if   ( FLU_SCHEME == RTVD )

         CUFLU_FluidSolver_RTVD <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>> 
                                ( d_Flu_Array_F_In  + UsedPatch[s],
                                  d_Flu_Array_F_Out + UsedPatch[s],
                                  d_Flux_Array      + UsedPatch[s],
                                  dt, _dh, Gamma, StoreFlux, XYZ );

#        elif ( FLU_SCHEME == WAF )

         CUFLU_FluidSolver_WAF <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>> 
                               ( d_Flu_Array_F_In  + UsedPatch[s],
                                 d_Flu_Array_F_Out + UsedPatch[s],
                                 d_Flux_Array      + UsedPatch[s],
                                 dt, _dh, Gamma, StoreFlux, XYZ, WAF_Limiter );
         
#        elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )

         CUFLU_FluidSolver_MHM <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
                               ( d_Flu_Array_F_In  + UsedPatch[s], 
                                 d_Flu_Array_F_Out + UsedPatch[s], 
                                 d_Flux_Array      + UsedPatch[s], 
                                 d_PriVar          + UsedPatch[s],
                                 d_Slope_PPM_x     + UsedPatch[s],
                                 d_Slope_PPM_y     + UsedPatch[s],
                                 d_Slope_PPM_z     + UsedPatch[s],
                                 d_FC_Var_xL       + UsedPatch[s], 
                                 d_FC_Var_xR       + UsedPatch[s], 
                                 d_FC_Var_yL       + UsedPatch[s], 
                                 d_FC_Var_yR       + UsedPatch[s], 
                                 d_FC_Var_zL       + UsedPatch[s], 
                                 d_FC_Var_zR       + UsedPatch[s], 
                                 d_FC_Flux_x       + UsedPatch[s], 
                                 d_FC_Flux_y       + UsedPatch[s], 
                                 d_FC_Flux_z       + UsedPatch[s],
                                 dt, _dh, Gamma, StoreFlux, LR_Limiter, MinMod_Coeff, EP_Coeff );

#        elif ( FLU_SCHEME == CTU )

         CUFLU_FluidSolver_CTU <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
                               ( d_Flu_Array_F_In  + UsedPatch[s], 
                                 d_Flu_Array_F_Out + UsedPatch[s], 
                                 d_Flux_Array      + UsedPatch[s], 
                                 d_PriVar          + UsedPatch[s],
                                 d_Slope_PPM_x     + UsedPatch[s],
                                 d_Slope_PPM_y     + UsedPatch[s],
                                 d_Slope_PPM_z     + UsedPatch[s],
                                 d_FC_Var_xL       + UsedPatch[s], 
                                 d_FC_Var_xR       + UsedPatch[s], 
                                 d_FC_Var_yL       + UsedPatch[s], 
                                 d_FC_Var_yR       + UsedPatch[s], 
                                 d_FC_Var_zL       + UsedPatch[s], 
                                 d_FC_Var_zR       + UsedPatch[s], 
                                 d_FC_Flux_x       + UsedPatch[s], 
                                 d_FC_Flux_y       + UsedPatch[s], 
                                 d_FC_Flux_z       + UsedPatch[s],
                                 dt, _dh, Gamma, StoreFlux, LR_Limiter, MinMod_Coeff, EP_Coeff );

#        else

#        error : unsupported GPU hydro scheme

#        endif // FLU_SCHEME

         /*
         if ( GetMinDtInfo )
            CUFLU_GetMaxCFL <<< NPatch_per_Stream[s], BlockDim_GetMinDtInfo, 0, Stream[s] >>> 
                            ( d_Flu_Array_F_Out + UsedPatch[s], 
                              d_MinDtInfo_Array + UsedPatch[s], 
                              Gamma );
                              */

#     elif ( MODEL == MHD )
#     warning :: WAIT MHD !!!

#     elif ( MODEL == ELBDM )

         CUFLU_ELBDMSolver <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>> 
                               ( d_Flu_Array_F_In  + UsedPatch[s],
                                 d_Flu_Array_F_Out + UsedPatch[s],
                                 dt, _dh, Eta, XYZ );

#     else

#        error : unsupported MODEL !!

#     endif // MODEL


      CUDA_CHECK_ERROR( cudaGetLastError() );
   } // for (int s=0; s<GPU_NStream; s++)


// c. copy data from device to host
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Flu_Array_Out + UsedPatch[s], d_Flu_Array_F_Out + UsedPatch[s],  
                         Flu_MemSize_Out[s], cudaMemcpyDeviceToHost, Stream[s] )  );

      if ( StoreFlux )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Flux_Array + UsedPatch[s], d_Flux_Array + UsedPatch[s],  
                         Flux_MemSize[s], cudaMemcpyDeviceToHost, Stream[s] )  );

#     if   ( MODEL == HYDRO )
      /*
      if ( GetMinDtInfo )
      {
         CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_MinDtInfo_Array       + UsedPatch[s], 
                                             d_MinDtInfo_Fluid_Array + UsedPatch[s], 
                            NPatch_per_Stream[s]*sizeof(real), cudaMemcpyDeviceToHost, Stream[s] )  );
      }
      */

#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!
#     endif // MODEL
   } // for (int s=0; s<GPU_NStream; s++)


   delete [] NPatch_per_Stream;
   delete [] UsedPatch;
   delete [] Flu_MemSize_In;
   delete [] Flu_MemSize_Out;
   delete [] Flux_MemSize;

} // FUNCTION : CUAPI_Asyn_FluidSolver



#endif // #ifdef GPU
