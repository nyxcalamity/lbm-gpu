# --------------------------------------------------------------------------------------------------
SRC_DIR=src
BUILD_DIR=build

SOURCES=$(SRC_DIR)/initialization.c $(SRC_DIR)/visualization.c $(SRC_DIR)/boundary.c \
	$(SRC_DIR)/collision.c $(SRC_DIR)/streaming.c $(SRC_DIR)/cell_computation.c \
	$(SRC_DIR)/utils.c $(SRC_DIR)/main.c $(SRC_DIR)/lbm_model.c
SOURCES_CU=$(SRC_DIR)/lbm_solver_gpu.cu $(SRC_DIR)/utils_gpu.cu $(SRC_DIR)/initialization_gpu.cu #$(SRC_DIR)/cell_computation_gpu.cu
EXECUTABLE=$(BUILD_DIR)/lbm-sim
#CUDA_LINKER=$(BUILD_DIR)/cuda-linker.o

#cuda settings
COMPUTE_CAPABILITY=20
# --------------------------------------------------------------------------------------------------
#Compiler command
CC=g++
CC_CU=nvcc
#Compiler flags
CFLAGS=-g -Wall -ofast -funroll-loops #-Werror -pedantic
CFLAGS_CU=-g --ptxas-options=-v -arch=sm_$(COMPUTE_CAPABILITY)
# Linker flags
LDFLAGS=-lm
LDFLAGS_CU=-lcudart #-lcuda
# --------------------------------------------------------------------------------------------------
OBJECTS=$(SOURCES:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
OBJECTS_CU=$(SOURCES_CU:$(SRC_DIR)/%.cu=$(BUILD_DIR)/%.o)

all: clean $(EXECUTABLE)

clean:
	rm -f $(OBJECTS) $(OBJECTS_CU) $(EXECUTABLE)
	rm -f img/*.vtk

$(EXECUTABLE): $(OBJECTS) $(OBJECTS_CU) #$(CUDA_LINKER)
	$(CC) $(OBJECTS_CU) $(OBJECTS) -o $@ $(LDFLAGS) $(LDFLAGS_CU)
	#$(CC) $(CUDA_LINKER) $(OBJECTS) -o $@ $(LDFLAGS) $(LDFLAGS_CU)
# --------------------------------------------------------------------------------------------------
$(OBJECTS): $(BUILD_DIR)/%.o : $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)
# --------------------------------------------------------------------------------------------------
# CUDA compilation
# --------------------------------------------------------------------------------------------------
#$(CUDA_LINKER): $(OBJECTS_CU)
#	$(CC_CU) $(CFLAGS_CU) -dlink $(OBJECTS_CU) -o $@ $(LDFLAGS_CU)

$(OBJECTS_CU): $(BUILD_DIR)/%.o : $(SRC_DIR)/%.cu
	$(CC_CU) $(CFLAGS_CU) -c $< -o $@ $(LDFLAGS_CU)