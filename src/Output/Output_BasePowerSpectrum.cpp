
#include "DAINO.h"

#ifdef GRAVITY

#ifndef DENS
#error : ERROR : variable "DENS" is NOT defined in the function Output_BasePowerSpectrum !!
#endif


//output the dimensionless power spectrum
//#define DIMENSIONLESS_FORM


extern void GetBasePowerSpectrum( real *RhoK, const int j_start, const int dj, const int RhoK_Size, real *PS_total );

#ifdef SERIAL
extern rfftwnd_plan     FFTW_Plan, FFTW_Plan_Inv;
#else
extern rfftwnd_mpi_plan FFTW_Plan, FFTW_Plan_Inv;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_BasePowerSpectrum 
// Description :  Evaluate and output the base-level power spectrum by FFT 
//
// Parameter   :  FileName : Name of the output file 
//-------------------------------------------------------------------------------------------------------
void Output_BasePowerSpectrum( const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// 1. get the array indices using by FFTW
   const int Nx_Padded = NX0_TOT[0]/2+1;
   int local_nz, local_z_start, local_ny_after_transpose, local_y_start_after_transpose, total_local_size;

#  ifdef SERIAL
   local_ny_after_transpose      = NULL_INT;
   local_y_start_after_transpose = NULL_INT;
   total_local_size              = 2*Nx_Padded*NX0_TOT[1]*NX0_TOT[2];
#  else
   rfftwnd_mpi_local_sizes( FFTW_Plan, &local_nz, &local_z_start, &local_ny_after_transpose,
                            &local_y_start_after_transpose, &total_local_size);
#  endif


// 2. allocate memory
   real *PS_total = NULL;
   real *RhoK     = new real [ total_local_size ];
   real *SendBuf  = new real [ patch->NPatchComma[0][1]*PS1*PS1*PS1 ];
#  ifdef SERIAL
   real *RecvBuf  = SendBuf;
#  else
   real *RecvBuf  = new real [ NX0_TOT[0]*NX0_TOT[1]*local_nz ];
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

   if ( MPI_Rank == 0 )    PS_total = new real [Nx_Padded];


// 3. rearrange data from cube to slice (or space-filling curve decomposition if load-balance is enabled)
#  ifdef LOAD_BALANCE 
   LB_SFC_to_Slice( RhoK, SendBuf, RecvBuf, SendBuf_SIdx, RecvBuf_SIdx, List_PID, List_k, List_NSend, List_NRecv,
                    List_z_start, local_nz );
#  else
   Cube_to_Slice( RhoK, SendBuf, RecvBuf );
#  endif


// 4. evaluate the base-level power spectrum by FFT
   GetBasePowerSpectrum( RhoK, local_y_start_after_transpose, local_ny_after_transpose, total_local_size, PS_total );


// 5. output the power spectrum
   if ( MPI_Rank == 0 )
   {
//    check if the targeted file already exists
      FILE *File_Check = fopen( FileName, "r" );
      if ( File_Check != NULL )
      {
         Aux_Message( stderr, "WARNING : the file \"%s\" already exists and will be overwritten !!\n", FileName );
         fclose( File_Check );
      }


//    output the power spectrum
      const real WaveK0 = 2.0*M_PI/patch->BoxSize[0];
      FILE *File = fopen( FileName, "w" );

      fprintf( File, "%13s%4s%13s\n", "k", "", "Power" );

      for (int b=0; b<Nx_Padded; b++)     fprintf( File, "%13.6e%4s%13.6e\n", WaveK0*b, "", PS_total[b] );

      fclose( File );
   } // if ( MPI_Rank == 0 )


// 6. free memory
   delete [] RhoK;
   delete [] SendBuf;
#  ifndef SERIAL
   delete [] RecvBuf;
#  endif
#  ifdef LOAD_BALANCE
   delete [] SendBuf_SIdx;
   delete [] RecvBuf_SIdx;
#  endif
   if ( MPI_Rank == 0 )    delete [] PS_total; 


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_BasePowerSpectrum



//-------------------------------------------------------------------------------------------------------
// Function    :  GetBasePowerSpectrum
// Description :  Evaluate and base-level power spectrum by FFT 
//
// Note        :  Invoked by the function "Output_BasePowerSpectrum"
//
// Parameter   :  RhoK        : Array storing the input density and output potential
//                j_start     : Starting j index
//                dj          : Size of array in the j (y) direction after the forward FFT
//                RhoK_Size   : Size of the array "RhoK"
//                PS_total    : Power spectrum summed over all MPI ranks
//
// Return      :  PS_total
//-------------------------------------------------------------------------------------------------------
void GetBasePowerSpectrum( real *RhoK, const int j_start, const int dj, const int RhoK_Size, real *PS_total )
{

// check
   if ( MPI_Rank == 0  &&  PS_total == NULL )   Aux_Error( ERROR_INFO, "PS_total == NULL at the root rank !!\n" );


   const int Nx        = NX0_TOT[0];
   const int Ny        = NX0_TOT[1];
   const int Nz        = NX0_TOT[2];
   const int Nx_Padded = Nx/2 + 1;

   fftw_complex *cdata=NULL;
   real PS_local[Nx_Padded];
   long Count_local[Nx_Padded], Count_total[Nx_Padded];
   int  bin, bin_i[Nx_Padded], bin_j[Ny], bin_k[Nz]; 
   

// forward FFT
#  ifdef SERIAL
   rfftwnd_one_real_to_complex( FFTW_Plan, RhoK, NULL );
#  else
   rfftwnd_mpi( FFTW_Plan, 1, RhoK, NULL, FFTW_TRANSPOSED_ORDER );
#  endif


// the data are now complex, so typecast a pointer
   cdata = (fftw_complex*) RhoK;


// set up the dimensionless wave number coefficients according to the FFTW data format
   for (int i=0; i<Nx_Padded; i++)     bin_i[i] = i;
   for (int j=0; j<Ny;        j++)     bin_j[j] = ( j <= Ny/2 ) ? j : j-Ny;
   for (int k=0; k<Nz;        k++)     bin_k[k] = ( k <= Nz/2 ) ? k : k-Nz;


// estimate the power spectrum
   int Idx;

   for (int b=0; b<Nx_Padded; b++)  
   {
      PS_local   [b] = (real)0.0;
      Count_local[b] = 0;
   }

#  ifdef SERIAL // serial mode

   for (int k=0; k<Nz; k++)
   {
      for (int j=0; j<Ny; j++)   
      for (int i=0; i<Nx_Padded; i++)
      {
         Idx = (k*Ny + j)*Nx_Padded + i;

#  else // parallel mode
   int j;

   for (int jj=0; jj<dj; jj++)   
   {  
      j = j_start + jj;

      for (int k=0; k<Nz; k++)
      for (int i=0; i<Nx_Padded; i++)
      {
         Idx = (jj*Nz + k)*Nx_Padded + i;

#  endif // #ifdef SERIAL ... else ...

         bin = int(   SQRT(   real( SQR(bin_i[i]) + SQR(bin_j[j]) + SQR(bin_k[k]) )  )   );

         if ( bin < Nx_Padded )
         {
            PS_local   [bin] += SQR( cdata[Idx].re ) + SQR( cdata[Idx].im );
            Count_local[bin] ++;
         }
      
      } // i,j,k
   } // i,j,k


// sum over all ranks
#  ifdef FLOAT8
   MPI_Reduce( PS_local,    PS_total,    Nx_Padded, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
#  else
   MPI_Reduce( PS_local,    PS_total,    Nx_Padded, MPI_FLOAT,  MPI_SUM, 0, MPI_COMM_WORLD );
#  endif

   MPI_Reduce( Count_local, Count_total, Nx_Padded, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );


// normalization: SQR(AveRho) accounts for Delta=Rho/AveRho
   const real Coeff = patch->BoxSize[0]*patch->BoxSize[1]*patch->BoxSize[2] / SQR( (real)Nx*(real)Ny*(real)Nz*AveDensity );
   real Norm;

#  ifdef DIMENSIONLESS_FORM
   const real k0    = (real)2.0*M_PI/patch->BoxSize[0];  // assuming cubic box
   real WaveK;
#  endif

   if ( MPI_Rank == 0 )
   {
      for (int b=0; b<Nx_Padded; b++)
      {
//       average
         PS_total[b] /= (real)Count_total[b];

//       normalization
#        ifdef DIMENSIONLESS_FORM
         WaveK        = b*k0;
         Norm         = Coeff*CUBE(WaveK)/(2.0*M_PI*M_PI);     // dimensionless power spectrum
#        else
         Norm         = Coeff;                                 // dimensional power spectrum [Mpc^3/h^3]
#        endif
         PS_total[b] *= Norm;
      }
   }

} // FUNCTION : GetBasePowerSpectrum



#endif // #ifdef GRAVITY
