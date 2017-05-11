
#include "DAINO.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  SiblingSearch
// Description :  Construct the sibling relation for level "lv" 
//
// Parameter   :  Targeted refinement level
//-------------------------------------------------------------------------------------------------------
void SiblingSearch( const int lv )
{

// check
   if ( lv < 0  ||  lv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );


// lv == 0 : invoke the function "SiblingSearch_Base" for the root level
   if ( lv == 0 )
   {
      SiblingSearch_Base();

      return;
   }

   
// lv > 0 :
#  pragma omp parallel for
   for (int PID=0; PID<patch->num[lv]; PID++)
   {
      int sibson;
      patch_t *fa = patch->ptr[0][lv-1][ patch->ptr[0][lv][PID]->father ];

      // initialize all siblings as -1
      for (int s=0; s<26; s++)      patch->ptr[0][lv][PID]->sibling[s] = -1;  

      switch ( PID%8 )
      {
         case 0:
            patch->ptr[0][lv][PID]->sibling[ 1] = PID+1;
            patch->ptr[0][lv][PID]->sibling[ 3] = PID+2;
            patch->ptr[0][lv][PID]->sibling[ 5] = PID+3;
            patch->ptr[0][lv][PID]->sibling[ 9] = PID+4;
            patch->ptr[0][lv][PID]->sibling[13] = PID+5;
            patch->ptr[0][lv][PID]->sibling[17] = PID+6;
            patch->ptr[0][lv][PID]->sibling[25] = PID+7;

            if ( fa->sibling[0] != -1  &&  patch->ptr[0][lv-1][fa->sibling[0]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[0]]->son;

               patch->ptr[0][lv][PID]->sibling[ 0] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[ 8] = sibson + 4;
               patch->ptr[0][lv][PID]->sibling[15] = sibson + 6;
               patch->ptr[0][lv][PID]->sibling[24] = sibson + 7;
            }

            if ( fa->sibling[2] != -1  &&  patch->ptr[0][lv-1][fa->sibling[2]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[2]]->son;
               
               patch->ptr[0][lv][PID]->sibling[ 2] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[ 7] = sibson + 4;
               patch->ptr[0][lv][PID]->sibling[12] = sibson + 5;
               patch->ptr[0][lv][PID]->sibling[23] = sibson + 7;
            }

            if ( fa->sibling[4] != -1  &&  patch->ptr[0][lv-1][fa->sibling[4]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[4]]->son;

               patch->ptr[0][lv][PID]->sibling[ 4] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[11] = sibson + 5;
               patch->ptr[0][lv][PID]->sibling[16] = sibson + 6;
               patch->ptr[0][lv][PID]->sibling[21] = sibson + 7;
            }

            if ( fa->sibling[6] != -1  &&  patch->ptr[0][lv-1][fa->sibling[6]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[6]]->son;
               
               patch->ptr[0][lv][PID]->sibling[ 6] = sibson + 4;
               patch->ptr[0][lv][PID]->sibling[22] = sibson + 7;
            }

            if ( fa->sibling[10] != -1  &&  patch->ptr[0][lv-1][fa->sibling[10]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[10]]->son;

               patch->ptr[0][lv][PID]->sibling[10] = sibson + 5;
               patch->ptr[0][lv][PID]->sibling[19] = sibson + 7;
            }

            if ( fa->sibling[14] != -1  &&  patch->ptr[0][lv-1][fa->sibling[14]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[14]]->son;
               
               patch->ptr[0][lv][PID]->sibling[14] = sibson + 6;
               patch->ptr[0][lv][PID]->sibling[20] = sibson + 7;
            }

            if ( fa->sibling[18] != -1  &&  patch->ptr[0][lv-1][fa->sibling[18]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[18]]->son;

               patch->ptr[0][lv][PID]->sibling[18] = sibson + 7;
            }

            break;


         case 1:
            patch->ptr[0][lv][PID]->sibling[ 8] = PID+1;
            patch->ptr[0][lv][PID]->sibling[15] = PID+2;
            patch->ptr[0][lv][PID]->sibling[ 3] = PID+3;
            patch->ptr[0][lv][PID]->sibling[24] = PID+4;
            patch->ptr[0][lv][PID]->sibling[ 5] = PID+5;
            patch->ptr[0][lv][PID]->sibling[13] = PID+6;
            patch->ptr[0][lv][PID]->sibling[ 0] = PID-1;
 
            if ( fa->sibling[1] != -1  &&  patch->ptr[0][lv-1][fa->sibling[1]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[1]]->son;
               
               patch->ptr[0][lv][PID]->sibling[ 1] = sibson ;
               patch->ptr[0][lv][PID]->sibling[ 9] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[17] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[25] = sibson + 5;
            }

            if ( fa->sibling[2] != -1  &&  patch->ptr[0][lv-1][fa->sibling[2]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[2]]->son;
               
               patch->ptr[0][lv][PID]->sibling[ 6] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[ 2] = sibson + 4;
               patch->ptr[0][lv][PID]->sibling[22] = sibson + 5;
               patch->ptr[0][lv][PID]->sibling[12] = sibson + 7;
            }

            if ( fa->sibling[4] != -1  &&  patch->ptr[0][lv-1][fa->sibling[4]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[4]]->son;
               
               patch->ptr[0][lv][PID]->sibling[14] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[20] = sibson + 5;
               patch->ptr[0][lv][PID]->sibling[ 4] = sibson + 6;
               patch->ptr[0][lv][PID]->sibling[11] = sibson + 7;
            }

            if ( fa->sibling[7] != -1  &&  patch->ptr[0][lv-1][fa->sibling[7]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[7]]->son;
               
               patch->ptr[0][lv][PID]->sibling[ 7] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[23] = sibson + 5;
            }

            if ( fa->sibling[10] != -1  &&  patch->ptr[0][lv-1][fa->sibling[10]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[10]]->son;
               
               patch->ptr[0][lv][PID]->sibling[18] = sibson + 5;
               patch->ptr[0][lv][PID]->sibling[10] = sibson + 7;
            }

            if ( fa->sibling[16] != -1  &&  patch->ptr[0][lv-1][fa->sibling[16]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[16]]->son;
               
               patch->ptr[0][lv][PID]->sibling[16] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[21] = sibson + 5;
            }

            if ( fa->sibling[19] != -1  &&  patch->ptr[0][lv-1][fa->sibling[19]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[19]]->son;
               
               patch->ptr[0][lv][PID]->sibling[19] = sibson + 5;
            }

           break;


         case 2:
            patch->ptr[0][lv][PID]->sibling[12] = PID+1;
            patch->ptr[0][lv][PID]->sibling[ 1] = PID+2;
            patch->ptr[0][lv][PID]->sibling[ 5] = PID+3;
            patch->ptr[0][lv][PID]->sibling[23] = PID+4;
            patch->ptr[0][lv][PID]->sibling[17] = PID+5;
            patch->ptr[0][lv][PID]->sibling[ 7] = PID-1;
            patch->ptr[0][lv][PID]->sibling[ 2] = PID-2;
 
            if ( fa->sibling[0] != -1  &&  patch->ptr[0][lv-1][fa->sibling[0]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[0]]->son;
               
               patch->ptr[0][lv][PID]->sibling[ 6] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[ 0] = sibson + 4;
               patch->ptr[0][lv][PID]->sibling[22] = sibson + 6;
               patch->ptr[0][lv][PID]->sibling[15] = sibson + 7;
            }

            if ( fa->sibling[3] != -1  &&  patch->ptr[0][lv-1][fa->sibling[3]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[3]]->son;
               
               patch->ptr[0][lv][PID]->sibling[ 3] = sibson;
               patch->ptr[0][lv][PID]->sibling[ 9] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[13] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[25] = sibson + 6;
            }

            if ( fa->sibling[4] != -1  &&  patch->ptr[0][lv-1][fa->sibling[4]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[4]]->son;
               
               patch->ptr[0][lv][PID]->sibling[10] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[ 4] = sibson + 5;
               patch->ptr[0][lv][PID]->sibling[19] = sibson + 6;
               patch->ptr[0][lv][PID]->sibling[16] = sibson + 7;
            }

            if( fa->sibling[8] != -1  &&  patch->ptr[0][lv-1][fa->sibling[8]]->son != - 1)
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[8]]->son;
               
               patch->ptr[0][lv][PID]->sibling[ 8] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[24] = sibson + 6;
            }

            if ( fa->sibling[11] != -1  &&  patch->ptr[0][lv-1][fa->sibling[11]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[11]]->son;
               
               patch->ptr[0][lv][PID]->sibling[11] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[21] = sibson + 6;
            }

            if ( fa->sibling[14] != -1  &&  patch->ptr[0][lv-1][fa->sibling[14]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[14]]->son;
               
               patch->ptr[0][lv][PID]->sibling[18] = sibson + 6;
               patch->ptr[0][lv][PID]->sibling[14] = sibson + 7;
            }

            if ( fa->sibling[20] != -1  &&  patch->ptr[0][lv-1][fa->sibling[20]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[20]]->son;

               patch->ptr[0][lv][PID]->sibling[20] = sibson + 6;
            }
            
           break;
      

         case 3:
            patch->ptr[0][lv][PID]->sibling[21] = PID+1;
            patch->ptr[0][lv][PID]->sibling[ 3] = PID+2;
            patch->ptr[0][lv][PID]->sibling[ 1] = PID+3;
            patch->ptr[0][lv][PID]->sibling[ 9] = PID+4;
            patch->ptr[0][lv][PID]->sibling[11] = PID-1;
            patch->ptr[0][lv][PID]->sibling[16] = PID-2;
            patch->ptr[0][lv][PID]->sibling[ 4] = PID-3;
 
            if ( fa->sibling[0] != -1  &&  patch->ptr[0][lv-1][fa->sibling[0]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[0]]->son;

               patch->ptr[0][lv][PID]->sibling[14] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[20] = sibson + 4;
               patch->ptr[0][lv][PID]->sibling[ 0] = sibson + 6;
               patch->ptr[0][lv][PID]->sibling[ 8] = sibson + 7;
            }

            if ( fa->sibling[2] != -1  &&  patch->ptr[0][lv-1][fa->sibling[2]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[2]]->son;
               
               patch->ptr[0][lv][PID]->sibling[10] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[19] = sibson + 4;
               patch->ptr[0][lv][PID]->sibling[ 2] = sibson + 5;
               patch->ptr[0][lv][PID]->sibling[ 7] = sibson + 7;
            }

            if ( fa->sibling[5] != -1  &&  patch->ptr[0][lv-1][fa->sibling[5]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[5]]->son;
               
               patch->ptr[0][lv][PID]->sibling[ 5] = sibson;
               patch->ptr[0][lv][PID]->sibling[17] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[13] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[25] = sibson + 4;
            }

            if ( fa->sibling[6] != -1  &&  patch->ptr[0][lv-1][fa->sibling[6]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[6]]->son;
               
               patch->ptr[0][lv][PID]->sibling[18] = sibson + 4;
               patch->ptr[0][lv][PID]->sibling[ 6] = sibson + 7;
            }

            if ( fa->sibling[12] != -1  &&  patch->ptr[0][lv-1][fa->sibling[12]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[12]]->son;
               
               patch->ptr[0][lv][PID]->sibling[12] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[23] = sibson + 4;
            }

            if ( fa->sibling[15] != -1  &&  patch->ptr[0][lv-1][fa->sibling[15]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[15]]->son;

               patch->ptr[0][lv][PID]->sibling[15] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[24] = sibson + 4;
            }

            if ( fa->sibling[22] != -1  &&  patch->ptr[0][lv-1][fa->sibling[22]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[22]]->son;
               
               patch->ptr[0][lv][PID]->sibling[22] = sibson + 4;
            }

           break;
      

         case 4:
            patch->ptr[0][lv][PID]->sibling[15] = PID+1;
            patch->ptr[0][lv][PID]->sibling[12] = PID+2;
            patch->ptr[0][lv][PID]->sibling[ 5] = PID+3;
            patch->ptr[0][lv][PID]->sibling[22] = PID-1;
            patch->ptr[0][lv][PID]->sibling[ 0] = PID-2;
            patch->ptr[0][lv][PID]->sibling[ 2] = PID-3;
            patch->ptr[0][lv][PID]->sibling[ 6] = PID-4;
  
            if ( fa->sibling[1] != -1  &&  patch->ptr[0][lv-1][fa->sibling[1]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[1]]->son;
               
               patch->ptr[0][lv][PID]->sibling[ 7] = sibson ;
               patch->ptr[0][lv][PID]->sibling[ 1] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[23] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[17] = sibson + 5;
            }

            if ( fa->sibling[3] != -1  &&  patch->ptr[0][lv-1][fa->sibling[3]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[3]]->son;

               patch->ptr[0][lv][PID]->sibling[ 8] = sibson;
               patch->ptr[0][lv][PID]->sibling[ 3] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[24] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[13] = sibson + 6;
            }

            if ( fa->sibling[4] != -1  &&  patch->ptr[0][lv-1][fa->sibling[4]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[4]]->son;

               patch->ptr[0][lv][PID]->sibling[18] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[14] = sibson + 5;
               patch->ptr[0][lv][PID]->sibling[10] = sibson + 6;
               patch->ptr[0][lv][PID]->sibling[ 4] = sibson + 7;
            }

            if ( fa->sibling[9] != -1  &&  patch->ptr[0][lv-1][fa->sibling[9]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[9]]->son;

               patch->ptr[0][lv][PID]->sibling[ 9] = sibson;
               patch->ptr[0][lv][PID]->sibling[25] = sibson + 3;
            }

            if ( fa->sibling[11] != -1  &&  patch->ptr[0][lv-1][fa->sibling[11]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[11]]->son;
               
               patch->ptr[0][lv][PID]->sibling[20] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[11] = sibson + 6;
            }

            if ( fa->sibling[16] != -1  &&  patch->ptr[0][lv-1][fa->sibling[16]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[16]]->son;
               
               patch->ptr[0][lv][PID]->sibling[19] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[16] = sibson + 5;
            }

            if ( fa->sibling[21] != -1  &&  patch->ptr[0][lv-1][fa->sibling[21]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[21]]->son;
               
               patch->ptr[0][lv][PID]->sibling[21] = sibson + 3;
            }
            
           break;
      
           
         case 5:
            patch->ptr[0][lv][PID]->sibling[ 7] = PID+1;
            patch->ptr[0][lv][PID]->sibling[ 1] = PID+2;
            patch->ptr[0][lv][PID]->sibling[16] = PID-1;
            patch->ptr[0][lv][PID]->sibling[ 2] = PID-2;
            patch->ptr[0][lv][PID]->sibling[ 4] = PID-3;
            patch->ptr[0][lv][PID]->sibling[19] = PID-4;
            patch->ptr[0][lv][PID]->sibling[10] = PID-5;
  
            if ( fa->sibling[0] != -1  &&  patch->ptr[0][lv-1][fa->sibling[0]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[0]]->son;
               
               patch->ptr[0][lv][PID]->sibling[18] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[14] = sibson + 4;
               patch->ptr[0][lv][PID]->sibling[ 6] = sibson + 6;
               patch->ptr[0][lv][PID]->sibling[ 0] = sibson + 7;
            }

            if ( fa->sibling[3] != -1  &&  patch->ptr[0][lv-1][fa->sibling[3]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[3]]->son;
               
               patch->ptr[0][lv][PID]->sibling[11] = sibson;
               patch->ptr[0][lv][PID]->sibling[21] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[ 3] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[ 9] = sibson + 6;
            }

            if ( fa->sibling[5] != -1  &&  patch->ptr[0][lv-1][fa->sibling[5]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[5]]->son;
               
               patch->ptr[0][lv][PID]->sibling[12] = sibson;
               patch->ptr[0][lv][PID]->sibling[23] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[ 5] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[17] = sibson + 4;
            }

            if ( fa->sibling[8] != -1  &&  patch->ptr[0][lv-1][fa->sibling[8]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[8]]->son;
               
               patch->ptr[0][lv][PID]->sibling[20] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[ 8] = sibson + 6;
            }

            if ( fa->sibling[13] != -1  &&  patch->ptr[0][lv-1][fa->sibling[13]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[13]]->son;
               
               patch->ptr[0][lv][PID]->sibling[13] = sibson;
               patch->ptr[0][lv][PID]->sibling[25] = sibson + 1;
            }

            if ( fa->sibling[15] != -1  &&  patch->ptr[0][lv-1][fa->sibling[15]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[15]]->son;
               
               patch->ptr[0][lv][PID]->sibling[22] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[15] = sibson + 4;
            }

            if ( fa->sibling[24] != -1  &&  patch->ptr[0][lv-1][fa->sibling[24]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[24]]->son;
               
               patch->ptr[0][lv][PID]->sibling[24] = sibson + 1;
            }
            
           break;
            
           
         case 6:
            patch->ptr[0][lv][PID]->sibling[ 3] = PID+1;
            patch->ptr[0][lv][PID]->sibling[ 8] = PID-1;
            patch->ptr[0][lv][PID]->sibling[11] = PID-2;
            patch->ptr[0][lv][PID]->sibling[ 0] = PID-3;
            patch->ptr[0][lv][PID]->sibling[20] = PID-4;
            patch->ptr[0][lv][PID]->sibling[ 4] = PID-5;
            patch->ptr[0][lv][PID]->sibling[14] = PID-6;
  
            if ( fa->sibling[1] != -1  &&  patch->ptr[0][lv-1][fa->sibling[1]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[1]]->son;
               
               patch->ptr[0][lv][PID]->sibling[16] = sibson;
               patch->ptr[0][lv][PID]->sibling[21] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[ 1] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[ 9] = sibson + 5;
            }

            if ( fa->sibling[2] != -1  &&  patch->ptr[0][lv-1][fa->sibling[2]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[2]]->son;
               
               patch->ptr[0][lv][PID]->sibling[18] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[10] = sibson + 4;
               patch->ptr[0][lv][PID]->sibling[ 6] = sibson + 5;
               patch->ptr[0][lv][PID]->sibling[ 2] = sibson + 7;
            }

            if ( fa->sibling[5] != -1  &&  patch->ptr[0][lv-1][fa->sibling[5]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[5]]->son;
               
               patch->ptr[0][lv][PID]->sibling[15] = sibson;
               patch->ptr[0][lv][PID]->sibling[ 5] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[24] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[13] = sibson + 4;
            }

            if ( fa->sibling[7] != -1  &&  patch->ptr[0][lv-1][fa->sibling[7]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[7]]->son;
               
               patch->ptr[0][lv][PID]->sibling[19] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[ 7] = sibson + 5;
            }

            if ( fa->sibling[12] != -1  &&  patch->ptr[0][lv-1][fa->sibling[12]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[12]]->son;
               
               patch->ptr[0][lv][PID]->sibling[22] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[12] = sibson + 4;
            }

            if ( fa->sibling[17] != -1  &&  patch->ptr[0][lv-1][fa->sibling[17]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[17]]->son;
               
               patch->ptr[0][lv][PID]->sibling[17] = sibson;
               patch->ptr[0][lv][PID]->sibling[25] = sibson + 2;
            }

            if ( fa->sibling[23] != -1  &&  patch->ptr[0][lv-1][fa->sibling[23]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[23]]->son;
               
               patch->ptr[0][lv][PID]->sibling[23] = sibson + 2;
            }
            
           break;
            
           
         case 7:
            patch->ptr[0][lv][PID]->sibling[ 2] = PID-1;
            patch->ptr[0][lv][PID]->sibling[ 0] = PID-2;
            patch->ptr[0][lv][PID]->sibling[ 4] = PID-3;
            patch->ptr[0][lv][PID]->sibling[ 6] = PID-4;
            patch->ptr[0][lv][PID]->sibling[14] = PID-5;
            patch->ptr[0][lv][PID]->sibling[10] = PID-6;
            patch->ptr[0][lv][PID]->sibling[18] = PID-7;
   
            if ( fa->sibling[1] != -1  &&  patch->ptr[0][lv-1][fa->sibling[1]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[1]]->son;
               
               patch->ptr[0][lv][PID]->sibling[19] = sibson;
               patch->ptr[0][lv][PID]->sibling[16] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[ 7] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[ 1] = sibson + 5;
            }

            if ( fa->sibling[3] != -1  &&  patch->ptr[0][lv-1][fa->sibling[3]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[3]]->son;
               
               patch->ptr[0][lv][PID]->sibling[20] = sibson;
               patch->ptr[0][lv][PID]->sibling[11] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[ 8] = sibson + 3;
               patch->ptr[0][lv][PID]->sibling[ 3] = sibson + 6;
            }

            if ( fa->sibling[5] != -1  &&  patch->ptr[0][lv-1][fa->sibling[5]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[5]]->son;
               
               patch->ptr[0][lv][PID]->sibling[22] = sibson;
               patch->ptr[0][lv][PID]->sibling[12] = sibson + 1;
               patch->ptr[0][lv][PID]->sibling[15] = sibson + 2;
               patch->ptr[0][lv][PID]->sibling[ 5] = sibson + 4;
            }

            if ( fa->sibling[9] != -1  &&  patch->ptr[0][lv-1][fa->sibling[9]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[9]]->son;
               
               patch->ptr[0][lv][PID]->sibling[21] = sibson;
               patch->ptr[0][lv][PID]->sibling[ 9] = sibson + 3;
            }

            if ( fa->sibling[13] != -1  &&  patch->ptr[0][lv-1][fa->sibling[13]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[13]]->son;
               
               patch->ptr[0][lv][PID]->sibling[24] = sibson;
               patch->ptr[0][lv][PID]->sibling[13] = sibson + 1;
            }

            if ( fa->sibling[17] != -1  &&  patch->ptr[0][lv-1][fa->sibling[17]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[17]]->son;
               
               patch->ptr[0][lv][PID]->sibling[23] = sibson;
               patch->ptr[0][lv][PID]->sibling[17] = sibson + 2;
            }

            if ( fa->sibling[25] != -1  &&  patch->ptr[0][lv-1][fa->sibling[25]]->son != -1 )
            {
               sibson = patch->ptr[0][lv-1][fa->sibling[25]]->son;
               
               patch->ptr[0][lv][PID]->sibling[25] = sibson;
            }

           break;

      }  // switch ( PID%8 )
   }  // for (int PID=0; PID<patch->num[lv]; PID++)

} // FUNCTION : SiblingSearch
