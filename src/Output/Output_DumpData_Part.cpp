
#include "DAINO.h"

static void WriteFile( FILE *File, const int lv, const int PID, const int i, const int j, const int k,
                       const int ii, const int jj, const int kk );




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpData_Part
// Description :  Output part of data in the ASCII form
//
// Parameter   :  Part     : OUTPUT_XY   : xy plane  
//                           OUTPUT_YZ   : yz plane
//                           OUTPUT_XZ   : xz plane
//                           OUTPUT_X    : x  line
//                           OUTPUT_Y    : y  line 
//                           OUTPUT_Z    : z  line
//                           OUTPUT_DIAG : diagonal along (+1,+1,+1)
//
//                BaseOnly : Only output the base-level data
//
//                x        : x coordinate
//                y        : y coordinate
//                z        : z coordinate
//
//                FileName : Name of the output file
//-------------------------------------------------------------------------------------------------------
void Output_DumpData_Part( const OptOutputPart_t Part, const bool BaseOnly, const real x, const real y, 
                           const real z, const char *FileName )
{  

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check the input parameters
   if ( Part != OUTPUT_XY  &&  Part != OUTPUT_YZ  &&  Part != OUTPUT_XZ  &&  
        Part != OUTPUT_X   &&  Part != OUTPUT_Y   &&  Part != OUTPUT_Z   &&  Part != OUTPUT_DIAG )
      Aux_Error( ERROR_INFO, "unsupported option \"Part = %d\" [0 ~ 6] !!\n", Part );

   if (  ( Part == OUTPUT_YZ  ||  Part == OUTPUT_Y  ||  Part == OUTPUT_Z )  &&
         ( x < 0.0  ||  x >= patch->BoxSize[0] )  )
      Aux_Error( ERROR_INFO, "incorrect x (out of range [0<=X<%f]) !!\n", patch->BoxSize[0] );

   if (  ( Part == OUTPUT_XZ  ||  Part == OUTPUT_X  ||  Part == OUTPUT_Z )  &&
         ( y < 0.0  ||  y >= patch->BoxSize[1] )  )
      Aux_Error( ERROR_INFO, "incorrect y (out of range [0<=Y<%f]) !!\n", patch->BoxSize[1] );
      
   if (  ( Part == OUTPUT_XY  ||  Part == OUTPUT_X  ||  Part == OUTPUT_Y )  &&
         ( z < 0.0  ||  z >= patch->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "incorrect z (out of range [0<=Z<%f]) !!\n", patch->BoxSize[2] );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_Check_Synchronization( Time[0], Time[lv], __FUNCTION__, true );


// check if the file already exists
   if ( MPI_Rank == 0 )  
   {
      FILE *File_Check = fopen( FileName, "r" );
      if ( File_Check != NULL )  
      {
         Aux_Message( stderr, "WARNING : the file \"%s\" already exists and will be overwritten !!\n", FileName );
         fclose( File_Check );

         FILE *Temp = fopen( FileName, "w" );
         fclose( Temp );
      }
   }


   const real dh_min = patch->dh[NLEVEL-1];;
   const int  NLv    = ( BaseOnly ) ? 1 : NLEVEL;

   int  ii, jj, kk, scale;
   real dh, PW;
   real xx, yy, zz;        // grid physical coordinates
   int *Corner  = NULL;    // corner grid ID 
   bool Check_x = false;
   bool Check_y = false;
   bool Check_z = false;

   switch ( Part )
   {
      case OUTPUT_XY :                                      Check_z = true;   break;
      case OUTPUT_YZ :  Check_x = true;                                       break;
      case OUTPUT_XZ :                    Check_y = true;                     break;
      case OUTPUT_X  :                    Check_y = true;   Check_z = true;   break;
      case OUTPUT_Y  :  Check_x = true;                     Check_z = true;   break;
      case OUTPUT_Z  :  Check_x = true;   Check_y = true;                     break;
   }


   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
         FILE *File = fopen( FileName, "a" );

//       output header
         if ( TargetMPIRank == 0 )  
         {
#           if   ( MODEL == HYDRO )
            fprintf( File, "%4s %4s %4s %10s %10s %10s %13s %13s %13s %13s %13s %13s %13s\n", 
                     "i", "j", "k", "x", "y", "z",
                     "Density", "Momentum x", "Momentum y", "Momentum z", "Energy", "Pressure", "Potential" );  

#           elif ( MODEL == MHD )
#           warning : WAIT MHD !!!

#           elif ( MODEL == ELBDM )
            fprintf( File, "%4s %4s %4s %10s %10s %10s %13s %13s %13s %13s\n", 
                     "i", "j", "k", "x", "y", "z",
                     "Density", "Real", "Imag", "Potential" );

#           else
#           error : ERROR : unsupported MODEL !!
#           endif // MODEL
         }


//       output data
         for (int lv=0; lv<NLv; lv++)                             
         {  

#ifndef OOC

            dh    = patch->dh   [lv];
            scale = patch->scale[lv];
            PW    = PATCH_SIZE*dh;

            for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)    
            {
//             output the patch data only if it has no son (if the option "BaseOnly" is turned off)
               if ( patch->ptr[0][lv][PID]->son == -1  ||  BaseOnly )
               {
                  Corner = patch->ptr[0][lv][PID]->corner;

                  if ( Part == OUTPUT_DIAG ) // (+1,+1,+1) diagonal
                  {
                     if ( Corner[0] == Corner[1]  &&  Corner[0] == Corner[2] )
                     {
                        for (int k=0; k<PS1; k++)  
                        {
                           kk = Corner[2] + k*scale;

                           WriteFile( File, lv, PID, k, k, k, kk, kk, kk );
                        }
                     }
                  } // if ( Part == OUTPUT_DIAG )


                  else // x/y/z lines || xy/yz/xz slices
                  {
//                   check whether the patch corner is within the targeted range
                     if (  !Check_x  ||  ( Corner[0]*dh_min <= x && Corner[0]*dh_min+PW > x )  )
                     if (  !Check_y  ||  ( Corner[1]*dh_min <= y && Corner[1]*dh_min+PW > y )  )
                     if (  !Check_z  ||  ( Corner[2]*dh_min <= z && Corner[2]*dh_min+PW > z )  )
                     {
//                      check whether the cell is within the targeted range
                        for (int k=0; k<PS1; k++)  {  kk = Corner[2] + k*scale;  zz = kk*dh_min;
                                                      if ( Check_z && ( zz>z || zz+dh<=z ) )    continue;

                        for (int j=0; j<PS1; j++)  {  jj = Corner[1] + j*scale;  yy = jj*dh_min;
                                                      if ( Check_y && ( yy>y || yy+dh<=y ) )    continue;

                        for (int i=0; i<PS1; i++)  {  ii = Corner[0] + i*scale;  xx = ii*dh_min;
                                                      if ( Check_x && ( xx>x || xx+dh<=x ) )    continue;

                           WriteFile( File, lv, PID, i, j, k, ii, jj, kk );

                        }}}
                     } // if patch corner is within the targeted range

                  } // if ( Part == OUTPUT_DIAG ... else ... )
               } // if ( patch->ptr[0][lv][PID]->son == -1 )
            } // for (int PID=0; PID<patch->NPatchComma[lv][1]; PID++)

#else // OOC

            OOC_Output_DumpData_Part( Part, BaseOnly, x, y, z, lv, File );

#endif
         } // for (int lv=0; lv<NLv; lv++)

         fclose( File );

      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_DumpData_Part



//-------------------------------------------------------------------------------------------------------
// Function    :  WriteFile 
// Description :  Output data to file
//
// Parameter   :  File     : File pointer 
//                lv       : Targeted refinement level
//                PID      : Patch ID
//                i/j/k    : Cell indices within the patch
//                ii/jj/kk : Cell scale indices in the simulation domain
//-------------------------------------------------------------------------------------------------------
void WriteFile( FILE *File, const int lv, const int PID, const int i, const int j, const int k,
                const int ii, const int jj, const int kk )
{

   const real dh_min  = patch->dh[NLEVEL-1];
   const real scale_2 = 0.5*patch->scale[lv];
   real u[NCOMP];

   for (int v=0; v<NCOMP; v++)   u[v] = patch->ptr[ patch->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

// output cell indices and coordinates
   fprintf( File, "%4d %4d %4d %10.4e %10.4e %10.4e", 
            ii, jj, kk, (ii+scale_2)*dh_min, (jj+scale_2)*dh_min, (kk+scale_2)*dh_min );

// output all variables in the fluid array
   for (int v=0; v<NCOMP; v++)   fprintf( File, " %13.6e", u[v] );

// output pressure in HYDRO
#  if   ( MODEL == HYDRO )
   fprintf( File, " %13.6e", ( u[ENGY]-0.5*(u[MOMX]*u[MOMX]+u[MOMY]*u[MOMY]+
                                            u[MOMZ]*u[MOMZ])/u[DENS] )*(GAMMA-1.0) );
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif // MODEL

// output potential
#  ifdef GRAVITY
   if ( OPT__OUTPUT_POT ) 
   fprintf( File, " %13.6e", patch->ptr[ patch->PotSg[lv] ][lv][PID]->pot[k][j][i] );
#  endif

   fprintf( File, "\n" );

} // FUNCTION : WriteFile
