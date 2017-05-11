
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Parallelization
// Description :  Initialize parameters for parallelization 
//
// Note        :  This function should be invoked even in the SERIAL mode
//-------------------------------------------------------------------------------------------------------
void Init_Parallelization()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Parallelization ... " );


// 1. number of coarse-grid cells in one rank for the rectangular domain decomposition
   for (int d=0; d<3; d++)    NX0[d] = NX0_TOT[d] / DAINO_NRANK_X(d);


// 2. mpi rank in different spatial directions
   MPI_Rank_X[0] =  MPI_Rank%MPI_NRank_X[0];
   MPI_Rank_X[1] = (MPI_Rank/MPI_NRank_X[0]) % MPI_NRank_X[1];
   MPI_Rank_X[2] = (MPI_Rank/MPI_NRank_X[0]) / MPI_NRank_X[1];


// 3. sibling mpi ranks
   MPI_SibRank[ 0] = ( MPI_Rank_X[2]                  )               *MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]                  )               *MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]+MPI_NRank_X[0]-1 )%MPI_NRank_X[0];

   MPI_SibRank[ 1] = ( MPI_Rank_X[2]                  )               *MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]                  )               *MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]               +1 )%MPI_NRank_X[0];

   MPI_SibRank[ 2] = ( MPI_Rank_X[2]                  )               *MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]+MPI_NRank_X[1]-1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]                  );

   MPI_SibRank[ 3] = ( MPI_Rank_X[2]                  )               *MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]               +1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                +
                     ( MPI_Rank_X[0]                  );

   MPI_SibRank[ 4] = ( MPI_Rank_X[2]+MPI_NRank_X[2]-1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]                  )               *MPI_NRank_X[0]                +
                     ( MPI_Rank_X[0]                  );

   MPI_SibRank[ 5] = ( MPI_Rank_X[2]               +1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]                  )               *MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]                  );

   MPI_SibRank[ 6] = ( MPI_Rank_X[2]                  )               *MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]+MPI_NRank_X[1]-1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]+MPI_NRank_X[0]-1 )%MPI_NRank_X[0];

   MPI_SibRank[ 7] = ( MPI_Rank_X[2]                  )               *MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]+MPI_NRank_X[1]-1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]               +1 )%MPI_NRank_X[0];

   MPI_SibRank[ 8] = ( MPI_Rank_X[2]                  )               *MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]               +1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]+MPI_NRank_X[0]-1 )%MPI_NRank_X[0];

   MPI_SibRank[ 9] = ( MPI_Rank_X[2]                  )               *MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]               +1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]               +1 )%MPI_NRank_X[0];
   
   MPI_SibRank[10] = ( MPI_Rank_X[2]+MPI_NRank_X[2]-1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]+MPI_NRank_X[1]-1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                +
                     ( MPI_Rank_X[0]                  );

   MPI_SibRank[11] = ( MPI_Rank_X[2]+MPI_NRank_X[2]-1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]               +1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]                  );

   MPI_SibRank[12] = ( MPI_Rank_X[2]               +1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]+MPI_NRank_X[1]-1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                +
                     ( MPI_Rank_X[0]                  );                  

   MPI_SibRank[13] = ( MPI_Rank_X[2]               +1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]               +1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                +
                     ( MPI_Rank_X[0]                  );

   MPI_SibRank[14] = ( MPI_Rank_X[2]+MPI_NRank_X[2]-1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] +
                     ( MPI_Rank_X[1]                  )               *MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]+MPI_NRank_X[0]-1 )%MPI_NRank_X[0];

   MPI_SibRank[15] = ( MPI_Rank_X[2]               +1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]                  )               *MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]+MPI_NRank_X[0]-1 )%MPI_NRank_X[0];
   
   MPI_SibRank[16] = ( MPI_Rank_X[2]+MPI_NRank_X[2]-1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]                  )               *MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]               +1 )%MPI_NRank_X[0];
   
   MPI_SibRank[17] = ( MPI_Rank_X[2]               +1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] +
                     ( MPI_Rank_X[1]                  )               *MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]               +1 )%MPI_NRank_X[0];
   
   MPI_SibRank[18] = ( MPI_Rank_X[2]+MPI_NRank_X[2]-1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]+MPI_NRank_X[1]-1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]+MPI_NRank_X[0]-1 )%MPI_NRank_X[0];
   
   MPI_SibRank[19] = ( MPI_Rank_X[2]+MPI_NRank_X[2]-1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]+MPI_NRank_X[1]-1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]               +1 )%MPI_NRank_X[0];
   
   MPI_SibRank[20] = ( MPI_Rank_X[2]+MPI_NRank_X[2]-1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]               +1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]+MPI_NRank_X[0]-1 )%MPI_NRank_X[0];
   
   MPI_SibRank[21] = ( MPI_Rank_X[2]+MPI_NRank_X[2]-1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]               +1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]               +1 )%MPI_NRank_X[0];
   
   MPI_SibRank[22] = ( MPI_Rank_X[2]               +1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]+MPI_NRank_X[1]-1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]+MPI_NRank_X[0]-1 )%MPI_NRank_X[0];
   
   MPI_SibRank[23] = ( MPI_Rank_X[2]               +1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] +
                     ( MPI_Rank_X[1]+MPI_NRank_X[1]-1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]               +1 )%MPI_NRank_X[0];
   
   MPI_SibRank[24] = ( MPI_Rank_X[2]               +1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]               +1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]+MPI_NRank_X[0]-1 )%MPI_NRank_X[0];
   
   MPI_SibRank[25] = ( MPI_Rank_X[2]               +1 )%MPI_NRank_X[2]*MPI_NRank_X[0]*MPI_NRank_X[1] + 
                     ( MPI_Rank_X[1]               +1 )%MPI_NRank_X[1]*MPI_NRank_X[0]                + 
                     ( MPI_Rank_X[0]               +1 )%MPI_NRank_X[0];


// 4. number of parallel buffer cells for setting boundary conditions by either data copy or interpolation
// 4.1 get the number of ghost zones required in different interpolations
   int NSide_Useless, NGhost_Flu, NGhost_RefFlu;
   Int_Table( OPT__FLU_INT_SCHEME,     NSide_Useless, NGhost_Flu );
   Int_Table( OPT__REF_FLU_INT_SCHEME, NSide_Useless, NGhost_RefFlu );

   int IntGhostSize_Flu;
   IntGhostSize_Flu = ( (FLU_GHOST_SIZE == 0) ? 0 : (FLU_GHOST_SIZE+1)/2 + NGhost_Flu );
   IntGhostSize_Flu = MAX( IntGhostSize_Flu, NGhost_RefFlu );

#  ifdef GRAVITY
   int NGhost_Pot, NGhost_Rho, NGhost_Gra, NGhost_RefPot;
   Int_Table( OPT__POT_INT_SCHEME,     NSide_Useless, NGhost_Pot );
   Int_Table( OPT__RHO_INT_SCHEME,     NSide_Useless, NGhost_Rho );
   Int_Table( OPT__GRA_INT_SCHEME,     NSide_Useless, NGhost_Gra );
   Int_Table( OPT__REF_POT_INT_SCHEME, NSide_Useless, NGhost_RefPot );

   int IntGhostSize_Pot, IntGhostSize_Rho;
   IntGhostSize_Pot = ( (POT_GHOST_SIZE == 0) ? 0 : (POT_GHOST_SIZE+1)/2 + NGhost_Pot );
   IntGhostSize_Pot = MAX( IntGhostSize_Pot, ( (GRA_GHOST_SIZE == 0) ? 0 : (GRA_GHOST_SIZE+1)/2 + NGhost_Gra ) );
   IntGhostSize_Pot = MAX( IntGhostSize_Pot, NGhost_RefPot );
   IntGhostSize_Rho = ( (RHO_GHOST_SIZE == 0) ? 0 : (RHO_GHOST_SIZE+1)/2 + NGhost_Rho );
#  endif


// 4.2 set the sizes of parallel buffers
   Flu_ParaBuf = MAX( FLU_GHOST_SIZE, IntGhostSize_Flu );
#  ifdef GRAVITY
   Pot_ParaBuf = MAX( GRA_GHOST_SIZE, IntGhostSize_Pot );
   Rho_ParaBuf = MAX( RHO_GHOST_SIZE, IntGhostSize_Rho );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : Init_Parallelization
