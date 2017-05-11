
   ********************************
   ** HYDRO Riemann problem test ** 
   ********************************

Procedure to run the test problem:
-------------------------------------------------------------------------------
1. Execute the script "Copy.sh" to copy files to their corresponding
   directories 

      Init_TestProb.cpp      --> DAINO/src/Init
      Makefile               --> DAINO/src
      Input__*               --> DAINO/bin/Run

2. Compile DAINO by typing "make" in the directory "DAINO/src"
   ("make clean" first if the program has been compiled previously)

3. Run DAINO by executing the file "Dizzy" in the directory "DAINO/bin/Run"


Note:
-------------------------------------------------------------------------------
1. Test problem parameters can be set in the file "Input__TestProb"
2. By default, the maximum refinement level = 2
3. By default, the refinement threshold of the Lohner error estimator = 0.50
4. Periodic boundary condition is adopted
5. Reference solutions have been put in the directory "ReferenceSolution". For
   example, one could use gnuplot and try "plot 'Xline_y0.000_z0.000_000011' 
   u 4:7 w p, '../../test_problem/Model_Hydro/Riemann/ReferenceSolution/Gamma_1.67/Sod_Shock_Tube'
   u 1:2 w l" 
