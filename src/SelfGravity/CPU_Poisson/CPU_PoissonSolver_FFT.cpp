
#include "DAINO.h"

#ifdef GRAVITY

#ifndef DENS
#error : ERROR : variable "DENS" is NOT defined in the CPU FFT solver !!
#endif



static void FFT( real *RhoK, const real Poi_Coeff, const int j_start, const int dj, const int RhoK_Size );

#ifdef SERIAL
extern rfftwnd_plan     FFTW_Plan, FFTW_Plan_Inv;
#else
extern rfftwnd_mpi_plan FFTW_Plan, FFTW_Plan_Inv;
#endif

#ifdef OOC
extern Timer_t *Timer_Gra_Advance[NLEVEL];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Cube_to_Slice
// Description :  Rectangular domain decomposition --> slice domain decomposition (for density)
//
// Parameter   :  RhoK     : Array to store the FFT results
//                SendBuf  : Sending buffer
//                RecvBuf  : Receiving buffer
//-------------------------------------------------------------------------------------------------------
void Cube_to_Slice( real *RhoK, real *SendBuf, real *RecvBuf )
{

   int OOC_NRank_X[3];
#  ifdef OOC
   OOC_NRank_X[0] = ooc.NRank_X[0];
   OOC_NRank_X[1] = ooc.NRank_X[1];
   OOC_NRank_X[2] = ooc.NRank_X[2];
#  else
   OOC_NRank_X[0] = 1;
   OOC_NRank_X[1] = 1;
   OOC_NRank_X[2] = 1;
#  endif
   
   const int scale0 = patch->scale[0];       
   const int LS[3]  = { NX0[0]*OOC_NRank_X[0], 
                        NX0[1]*OOC_NRank_X[1], 
                        NX0_TOT[2]/MPI_NRank  };            // LS : Layer Size
   const int PS[2] = { 2*(NX0_TOT[0]/2+1), NX0_TOT[1] };    // PS : Padded Size for FFTW

   int start_i, start_j, start_k;                           // the starting (i,j,k) indices of each patch
   int ii, jj, kk, layer, temp, ID1, ID2, i2, j2;


// copy the density into the send buffer layer-by-layer
#ifndef OOC

   for (int PID=0; PID<patch->NPatchComma[0][1]; PID++)
   {
      start_i = patch->ptr[0][0][PID]->corner[0] / scale0 - MPI_Rank_X[0]*NX0[0]*OOC_NRank_X[0];
      start_j = patch->ptr[0][0][PID]->corner[1] / scale0 - MPI_Rank_X[1]*NX0[1]*OOC_NRank_X[1];
      start_k = patch->ptr[0][0][PID]->corner[2] / scale0 - MPI_Rank_X[2]*NX0[2]*OOC_NRank_X[2];

      for (int k=0; k<PATCH_SIZE; k++)    {  temp  = start_k + k;
                                             layer = temp / LS[2];
                                             kk    = temp % LS[2];
      for (int j=0; j<PATCH_SIZE; j++)    {  jj    = start_j + j;
      for (int i=0; i<PATCH_SIZE; i++)    {  ii    = start_i + i;
   
         ID1 = ((layer*LS[2] + kk)*LS[1] + jj)*LS[0] + ii;

         SendBuf[ID1] = patch->ptr[ patch->FluSg[0] ][0][PID]->fluid[DENS][k][j][i];

      }}}

   } // for (int P=0; P<patch->NPatchComma[0][1]; P++)

#else // OOC

   OOC_CPU_PoissonSolver_FFT_Cube_to_Slice( SendBuf );

#endif


#  if ( defined OOC  &&  defined TIMING )
   MPI_Barrier( MPI_COMM_WORLD );
   Timer_Gra_Advance[0]->Start();
#  endif


// send different layers to different ranks
   MPI_CubeSlice( true, SendBuf, RecvBuf );


// copy the density from the recv buffer to the RhoK array layer-by-layer
   for (int k=0; k<NX0_TOT[2]/MPI_NRank; k++)   {  kk = k;
   for (int j=0; j<NX0_TOT[1]; j++)             {  jj = j%(NX0[1]*OOC_NRank_X[1]); j2 = j/(NX0[1]*OOC_NRank_X[1]);
   for (int i=0; i<NX0_TOT[0]; i++)             {  ii = i%(NX0[0]*OOC_NRank_X[0]); i2 = i/(NX0[0]*OOC_NRank_X[0]);
                                                   
      layer = j2*MPI_NRank_X[0] + i2;
      ID1   = (k*PS[1] + j)*PS[0] + i;
      ID2   = ((layer*LS[2] + kk)*LS[1] + jj)*LS[0] + ii;

      RhoK[ID1] = RecvBuf[ID2]; 

   }}}

#  if ( defined OOC  &&  defined TIMING )
   MPI_Barrier( MPI_COMM_WORLD );
   Timer_Gra_Advance[0]->Stop( false );
#  endif

} // FUNCTION : Cube_to_Slice



//-------------------------------------------------------------------------------------------------------
// Function    :  Slice_to_Cube
// Description :  Slice domain decomposition --> rectangular domain decomposition (for density)
// 
// Parameter   :  RhoK     : Array to store the FFT results
//                SendBuf  : Sending buffer
//                RecvBuf  : Receiving buffer
//                SaveSg   : Sandglass to store the updated data 
//-------------------------------------------------------------------------------------------------------
void Slice_to_Cube( real *RhoK, real *SendBuf, real *RecvBuf, const int SaveSg )
{
   
   int OOC_NRank_X[3];
#  ifdef OOC
   OOC_NRank_X[0] = ooc.NRank_X[0];
   OOC_NRank_X[1] = ooc.NRank_X[1];
   OOC_NRank_X[2] = ooc.NRank_X[2];
#  else
   OOC_NRank_X[0] = 1;
   OOC_NRank_X[1] = 1;
   OOC_NRank_X[2] = 1;
#  endif

   const int scale0 = patch->scale[0];       
   const int LS[3]  = { NX0[0]*OOC_NRank_X[0], 
                        NX0[1]*OOC_NRank_X[1], 
                        NX0_TOT[2]/MPI_NRank  };                  // LS : Layer Size
   const int PS[2]  = { 2*(NX0_TOT[0]/2+1), NX0_TOT[1] };         // PS : Padded Size for FFTW

   int start_i, start_j, start_k;                                 // the starting (i,j,k) indices of each patch
   int i2, j2, ii, jj, kk, layer, temp, ID1, ID2;


// copy the potential in the RhoK array to the send buffer layer-by-layer
#  if ( defined OOC  &&  defined TIMING )
   MPI_Barrier( MPI_COMM_WORLD );
   Timer_Gra_Advance[0]->Start();
#  endif

   for (int k=0; k<NX0_TOT[2]/MPI_NRank; k++)   {  kk = k;
   for (int j=0; j<NX0_TOT[1]; j++)             {  jj = j%(NX0[1]*OOC_NRank_X[1]); j2 = j/(NX0[1]*OOC_NRank_X[1]);
   for (int i=0; i<NX0_TOT[0]; i++)             {  ii = i%(NX0[0]*OOC_NRank_X[0]); i2 = i/(NX0[0]*OOC_NRank_X[0]);
                                                   
      layer = j2*MPI_NRank_X[0] + i2;
      ID1   = (k*PS[1] + j)*PS[0] + i;
      ID2   = ((layer*LS[2] + kk)*LS[1] + jj)*LS[0] + ii;

      SendBuf[ID2] = RhoK[ID1];

   }}}


// send different layers to different ranks
   MPI_CubeSlice( false, SendBuf, RecvBuf );


#  if ( defined OOC  &&  defined TIMING )
   MPI_Barrier( MPI_COMM_WORLD );
   Timer_Gra_Advance[0]->Stop( false );
#  endif


// copy the potential from the recv buffer to the patch->ptr 
#ifndef OOC

   for (int PID=0; PID<patch->NPatchComma[0][1]; PID++)
   {
      start_i = patch->ptr[0][0][PID]->corner[0] / scale0 - MPI_Rank_X[0]*NX0[0]*OOC_NRank_X[0];
      start_j = patch->ptr[0][0][PID]->corner[1] / scale0 - MPI_Rank_X[1]*NX0[1]*OOC_NRank_X[1];
      start_k = patch->ptr[0][0][PID]->corner[2] / scale0 - MPI_Rank_X[2]*NX0[2]*OOC_NRank_X[2];

      for (int k=0; k<PATCH_SIZE; k++)    {  temp  = start_k + k;
                                             layer = temp / LS[2];
                                             kk    = temp % LS[2];   
      for (int j=0; j<PATCH_SIZE; j++)    {  jj    = start_j + j;    
      for (int i=0; i<PATCH_SIZE; i++)    {  ii    = start_i + i;    
      
         ID1 = ((layer*LS[2] + kk)*LS[1] + jj)*LS[0] + ii;

         patch->ptr[SaveSg][DENS][PID]->pot[k][j][i] = RecvBuf[ID1];
      
      }}}
   }

#else // OOC

   OOC_CPU_PoissonSolver_FFT_Slice_to_Cube( RecvBuf, SaveSg );

#endif

} // FUNCTION : Slice_to_Cube



//-------------------------------------------------------------------------------------------------------
// Function    :  FFT
// Description :  Evaluate the gravitational potential by FFT 
//
// Note        :  Work with the periodic B.C.
//
// Parameter   :  RhoK        : Array storing the input density and output potential
//                Poi_Coeff   : Poi_Coefficient in front of density in the Poisson equation (4*Pi*Newton_G*a)   
//                j_start     : Starting j index
//                dj          : Size of array in the j (y) direction after the forward FFT
//                RhoK_Size   : Size of the array "RhoK"
//-------------------------------------------------------------------------------------------------------
void FFT( real *RhoK, const real Poi_Coeff, const int j_start, const int dj, const int RhoK_Size )
{

   const int Nx        = NX0_TOT[0];
   const int Ny        = NX0_TOT[1];
   const int Nz        = NX0_TOT[2];
   const int Nx_Padded = Nx/2 + 1;
   const real dh       = patch->dh[0];
   real Deno;
   fftw_complex *cdata;


// forward FFT
#  ifdef SERIAL
   rfftwnd_one_real_to_complex( FFTW_Plan, RhoK, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif


// the data are now complex, so typecast a pointer
   cdata = (fftw_complex*) RhoK;


// set up the dimensionless wave number and the corresponding sin(k)^2 function
   real kx[Nx_Padded], ky[Ny], kz[Nz];
   real sinkx2[Nx_Padded], sinky2[Ny], sinkz2[Nz];

   for (int i=0; i<Nx_Padded; i++) {   kx    [i] = 2.0*M_PI/Nx*i;
                                       sinkx2[i] = SQR(  SIN( (real)0.5*kx[i] )  );    }
   for (int j=0; j<Ny;        j++) {   ky    [j] = ( j <= Ny/2 ) ? 2.0*M_PI/Ny*j : 2.0*M_PI/Ny*(j-Ny);
                                       sinky2[j] = SQR(  SIN( (real)0.5*ky[j] )  );    }
   for (int k=0; k<Nz;        k++) {   kz    [k] = ( k <= Nz/2 ) ? 2.0*M_PI/Nz*k : 2.0*M_PI/Nz*(k-Nz);
                                       sinkz2[k] = SQR(  SIN( (real)0.5*kz[k] )  );    }


// divide the Rho_K by -k^2
#  ifdef SERIAL // serial mode
   int ID;

   for (int k=0; k<Nz; k++)
   {
      for (int j=0; j<Ny; j++)   
      for (int i=0; i<Nx_Padded; i++)
      {
         ID = (k*Ny + j)*Nx_Padded + i;

#  else // parallel mode
   int j, ID;

   for (int jj=0; jj<dj; jj++)   
   {  
      j = j_start + jj;

      for (int k=0; k<Nz; k++)
      for (int i=0; i<Nx_Padded; i++)
      {
         ID = (jj*Nz + k)*Nx_Padded + i;

#  endif // #ifdef SERIAL ... else ...

      
//       this form is more consistent with the "second-order discrete" Laplacian operator
         Deno = -4.0 * ( sinkx2[i] + sinky2[j] + sinkz2[k] );
//       Deno = -( kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k] );


//       remove the DC mode
         if ( Deno == 0.0 )
         {
            cdata[ID].re = 0.0;
            cdata[ID].im = 0.0;
         }

         else
         {
            cdata[ID].re =  cdata[ID].re * Poi_Coeff / Deno;
            cdata[ID].im =  cdata[ID].im * Poi_Coeff / Deno;
         }
      } // i,j,k
   } // i,j,k


// backward FFT
#  ifdef SERIAL
   rfftwnd_one_complex_to_real( FFTW_Plan_Inv, cdata, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan_Inv, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif


// normalization
   const real norm = dh*dh / ( (real)Nx*Ny*Nz );

   for (int t=0; t<RhoK_Size; t++)  RhoK[t] *= norm;

} // FUNCTION : FFT



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_PoissonSolver_FFT 
// Description :  Evaluate the base-level potential by FFT 
//
// Parameter   :  Poi_Coeff   : Coefficient in front of the RHS in the Poisson eq.
//                SaveSg      : Sandglass to store the updated data 
//-------------------------------------------------------------------------------------------------------
void CPU_PoissonSolver_FFT( const real Poi_Coeff, const int SaveSg )
{

// get the array indices using by FFTW
   int local_nz, local_z_start, local_ny_after_transpose, local_y_start_after_transpose, total_local_size;

#  ifdef SERIAL
   local_ny_after_transpose      = NULL_INT;
   local_y_start_after_transpose = NULL_INT;
   total_local_size              = 2*(NX0_TOT[0]/2+1)*NX0_TOT[1]*NX0_TOT[2];
#  else
   rfftwnd_mpi_local_sizes( FFTW_Plan, &local_nz, &local_z_start, &local_ny_after_transpose,
                            &local_y_start_after_transpose, &total_local_size );
#  endif


// allocate memory
   real *RhoK    = new real [ total_local_size ];
   real *SendBuf = new real [ patch->NPatchComma[0][1]*PS1*PS1*PS1 ];
#  ifdef SERIAL
   real *RecvBuf = SendBuf;
#  else
   real *RecvBuf = new real [ NX0_TOT[0]*NX0_TOT[1]*local_nz ];
#  endif
#  ifdef LOAD_BALANCE
   long *SendBuf_SIdx = new long [ patch->NPatchComma[0][1]*PS1 ];   // Sending MPI buffer of 1D coordinate in slab
   long *RecvBuf_SIdx = new long [ NX0_TOT[0]*NX0_TOT[1]*local_nz/PS1/PS1 ];
   int  *List_PID    [MPI_NRank];   // PID of each patch slice sent to each rank
   int  *List_k      [MPI_NRank];   // local z coordinate of each patch slice sent to each rank
   int   List_NSend  [MPI_NRank];   // size of data (density/potential) sent to each rank
   int   List_NRecv  [MPI_NRank];   // size of data (density/potential) received from each rank
   int   List_nz     [MPI_NRank];   // slab thickness of each rank in the FFTW slab decomposition
   int   List_z_start[MPI_NRank+1]; // starting z coordinate of each rank in the FFTW slab decomposition

// collect "local_nz" from all ranks and set the corresponding list "List_z_start"
   MPI_Allgather( &local_nz, 1, MPI_INT, List_nz, 1, MPI_INT, MPI_COMM_WORLD ); 

   List_z_start[0] = 0;
   for (int r=0; r<MPI_NRank; r++)  List_z_start[r+1] = List_z_start[r] + List_nz[r];

   if ( List_z_start[MPI_NRank] != NX0_TOT[2] )    
      Aux_Error( ERROR_INFO, "List_z_start[%d] (%d) != expectation (%d) !!\n", 
                 MPI_NRank, List_z_start[MPI_NRank],  NX0_TOT[2] );
#  endif // #ifdef LOAD_BALANCE


// rearrange data from cube to slice (or space-filling curve decomposition if load-balance is enabled)
#  ifdef LOAD_BALANCE 
   LB_SFC_to_Slice( RhoK, SendBuf, RecvBuf, SendBuf_SIdx, RecvBuf_SIdx, List_PID, List_k, List_NSend, List_NRecv,
                    List_z_start, local_nz );
#  else
   Cube_to_Slice( RhoK, SendBuf, RecvBuf );
#  endif


// evaluate the potential by FFT
#  if ( defined OOC  &&  defined TIMING )
   MPI_Barrier( MPI_COMM_WORLD );
   Timer_Gra_Advance[0]->Start();
#  endif

   FFT( RhoK, Poi_Coeff, local_y_start_after_transpose, local_ny_after_transpose, total_local_size );

#  if ( defined OOC  &&  defined TIMING )
   MPI_Barrier( MPI_COMM_WORLD );
   Timer_Gra_Advance[0]->Stop( false );
#  endif


// rearrange data from slice back to cube (or space-filling curve decomposition if load-balance is enabled)
#  ifdef LOAD_BALANCE 
   LB_Slice_to_SFC( RhoK, RecvBuf, SendBuf, SaveSg, RecvBuf_SIdx, List_PID, List_k, List_NRecv, List_NSend, local_nz );
#  else
   Slice_to_Cube( RhoK, RecvBuf, SendBuf, SaveSg );
#  endif


   delete [] RhoK;
   delete [] SendBuf;
#  ifndef SERIAL
   delete [] RecvBuf;
#  endif
#  ifdef LOAD_BALANCE
   delete [] SendBuf_SIdx;
   delete [] RecvBuf_SIdx;
#  endif

} // FUNCTION : CPU_PoissonSolver_FFT



#endif // #ifdef GRAVITY
