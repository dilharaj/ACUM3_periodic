#####################
# Makefile          #
#                   #
# Author      :     #
# Version     :     #
# Date        :     #
# Discription :     #
#####################

TARGET            := acum_cuda


CUDA_INSTALL_PATH ?= /usr/local/cuda

CXX = gcc
NVCC ?= $(CUDA_INSTALL_PATH)/bin/nvcc
INCD = -I"$(CUDA_INSTALL_PATH)/include" -I"./"
LIBS = -lgfortran -L"$(CUDA_INSTALL_PATH)/lib64" -lcudart -lcublas -lcufft
NVCCFLAGS := -arch=sm_35 # --ptxas-options=-v -G -g 
CXXFLAGS := -O3 


# files
CU_SOURCES        := variables.cu io.cu mathfunc.cu preprocess.cu core.cu bb_core.cu bb_core_gpu.cu bb_utils.cu main.cu
C_SOURCES         := #fft.c 
#HEADERS           := fft.h
OBJS              := $(patsubst %.c, %.o, $(C_SOURCES)) $(patsubst %.cu, %.o, $(CU_SOURCES))

%.o : %.c
	$(CXX) $(CXXFLAGS) -c $(INCD) -o $@ $<

%.o : %.cu 
	$(NVCC) $(NVCCFLAGS) -c $(INCD) -o $@ $<

$(TARGET): $(OBJS)
	$(NVCC) -o $(TARGET) $(OBJS) $(LDFLAGS) $(INCD) $(LIBS) 

clean:
	rm -f $(TARGET) *.o


