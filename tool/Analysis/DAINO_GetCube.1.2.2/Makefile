


# name of the executable file 
#######################################################################################################
EXECUTABLE  = DAINO_GetCube



# siimulation options
#######################################################################################################
# model: HYDRO/MHD/ELBDM
SIMU_OPTION += -DMODEL=HYDRO

# double precision
#SIMU_OPTION += -DFLOAT8

# enable OpenMP parallelization
SIMU_OPTION += -DOPENMP



# siimulation parameters
#######################################################################################################
NLEVEL         = 6         # level : 0 ~ NLEVEL-1 
MAX_PATCH      = 400000    # the maximum number of patches in each level
BUF_SIZE       = 2         # buffer size for interpolation 

NLEVEL         := $(strip $(NLEVEL))
MAX_PATCH      := $(strip $(MAX_PATCH))
BUF_SIZE       := $(strip $(BUF_SIZE))

SIMU_PARA = -DNLEVEL=$(NLEVEL) -DMAX_PATCH=$(MAX_PATCH) -DBUF_SIZE=$(BUF_SIZE)



# sources
#######################################################################################################
SOURCE      = Main.cpp  LoadData.cpp  Buf_AllocateBufferPatch_Base.cpp \
              Buf_AllocateBufferPatch.cpp  Buf_GetBufferData.cpp  Buf_RecordBoundaryPatch_Base.cpp \
              Buf_RecordBoundaryPatch.cpp  Buf_RecordExchangeDataPatchID.cpp  End_MemFree.cpp  FindFather.cpp \
              Flu_Restrict.cpp  Init_MemAllocate.cpp  Init_RecordBasePatch.cpp  MPI_Exit.cpp \
              SiblingSearch_Base.cpp  SiblingSearch.cpp  MPI_ExchangeBufferPosition.cpp  MPI_ExchangeInfo.cpp \
              Table_01.cpp  Table_02.cpp  Table_03.cpp  Table_04.cpp  Table_05.cpp  Table_07.cpp \
              Buf_SortBoundaryPatch.cpp 
#LoadData_Old.cpp



# rules and targets
#######################################################################################################
CC       := mpiCC
CFLAG    := -O3 -mp1 -w1
#CFLAG   := -O3 -Wall
#CFLAG   += -g -debug
LIB      := -limf

ifeq "$(findstring OPENMP, $(SIMU_OPTION))" "OPENMP"
CFLAG    += -openmp
#CFLAG    += -fopenmp
else
CFLAG    += -Wno-unknown-pragmas
endif

OBJ_DIR  := ./Object
INCLUDE  := -I./Header/

vpath %.cpp ./ DAINO_Functions

LIB := $(strip $(LIB))
OBJ := $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SOURCE))


$(OBJ_DIR)/%.o : %.cpp
	$(CC) $(CFLAG) $(SIMU_PARA) $(SIMU_OPTION) $(INCLUDE) -o $@ -c $<

$(EXECUTABLE): $(OBJ)
	$(CC) $(CFLAG) -o $@ $^ $(LIB)
	cp $(EXECUTABLE) ./Run/

clean:
	rm -f $(OBJ)
	rm -f $(EXECUTABLE)
