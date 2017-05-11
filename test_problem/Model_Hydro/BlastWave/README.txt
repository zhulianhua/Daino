
   ***************************
   ** HYDRO blast wave test **
   ***************************

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
1. Lohner's error estimator is applied to the energy density
2. By default, the maximum refinement level = 2
3. By default, the refinement threshold of the Lohner error estimator = 0.80
4. To plot the density profile, one could use gnuplot and try
   "plot 'Xline_y0.500_z0.500_000010' u 4:7 w p"
5. The tool "DAINO_SphereAnalysis" can be used to compute the density profile
