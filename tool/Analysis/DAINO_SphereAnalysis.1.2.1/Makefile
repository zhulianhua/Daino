


# file names
#######################################################################################################
PROGRAM    := DAINO_SphereAnalysis
EXECUTABLE := DAINO_SphereAnalysis



# siimulation options
#######################################################################################################
# model: HYDRO/MHD/ELBDM
SIMU_OPTION += -DMODEL=HYDRO

# double precision
#SIMU_OPTION += -DFLOAT8



# siimulation parameters
#######################################################################################################
NLEVEL      := 6        # level : 0 ~ NLEVEL-1 
MAX_PATCH   := 400000   # the maximum number of patches in each level

NLEVEL      := $(strip $(NLEVEL))
MAX_PATCH   := $(strip $(MAX_PATCH))

SIMU_PARA = -DNLEVEL=$(NLEVEL) -DMAX_PATCH=$(MAX_PATCH)



# rules and targets
#######################################################################################################
CC    := icpc
CFLAG := -O3 -mp1 -w1
#CFLAG += -debug -g
#CC    := g++ 
#CFLAG := -O3 -Wall


$(EXECUTABLE): $(PROGRAM).o
	 $(CC) $(CFLAG) -o $@ $< 
	 cp $(EXECUTABLE) ./Run/

$(PROGRAM).o: $(PROGRAM).cpp
	 $(CC) $(CFLAG) $(SIMU_PARA) $(SIMU_OPTION) -o $@ -c $<

clean:
	 rm -f *.o
	 rm -f $(EXECUTABLE)

