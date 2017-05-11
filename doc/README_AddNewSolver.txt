------------------------------------------------------
Basic steps for adding a new hydro solver into DAINO :
------------------------------------------------------

1. In the header file "include/Macro.h", add the new symbolic constant
   and the number of ghost cells required for the new hydro solver.
   For example:

      #define NEW_SOLVER_NAME    6
      #define FLU_GHOST_SIZE     4

   Enable the new fluid solver in the Makefile:

      FLU_SCHEME=NEW_SOLVER_NAME


2. Put the source code of the new GPU solver in the directory

      "src/Model_Hydro/GPU_Hydro/",

   and add the new filename in the Makefile. For example,

      CUDA_FILE += CUFLU_FluidSolver_NewSolverName.cu


3. Invoke the new GPU kernel in

      "src/GPU_API/CUAPI_Asyn_FluidSolver.cu"


4. Set the default CUDA parameters for the new GPU solver in

      "src/GPU_API/CUAPI_Set_Default_GPU_Parameter.cu",

   and set the number of CUDA threads for the new GPU solver in 
   
      "include/CUFLU.h"

   For example,

      #define FLU_BLOCK_SIZE_X    128
      #define FLU_BLOCK_SIZE_Y    1


********************* NOTICE *********************************************** 
For all steps given above, one could use the "CTU" scheme as a good example.
**************************************************************************** 

