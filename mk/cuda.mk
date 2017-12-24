# Redefine any make variables you want to override
# from mk/MxApp.mk here.


CUDA_PATH ?= /usr/local/cuda
#TARGET_ARCH ?= x86_64
#TARGET_SIZE ?= 64
HOST_COMPILER ?= g++

NVCC := $(CUDA_PATH)/bin/nvcc -ccbin $(HOST_COMPILER)

#NVCCFLAGS := -m${TARGET_SIZE}
NVCCFLAGS += $(addprefix -Xcompiler ,$(CXXFLAGS))
NVCCFLAGS += --expt-relaxed-constexpr
NVCCFLAGS += --std=c++14 -O3 

CUDA_LDFLAGS :=
CUDA_LDFLAGS += $(addprefix -Xlinker ,$(LDFLAGS))
CUDA_LDFLAGS += $(addprefix -Xlinker ,$(EXTRA_LDFLAGS))

# Common includes and paths for CUDA
INCLUDES  := -I/home/jrmales/cuda/samples/NVIDIA_CUDA-9.0_Samples/common/inc/
LIBRARIES :=

# Gencode arguments
SMS ?= 30 35 37 50 52 60 61 70

ifeq ($(GENCODE_FLAGS),)
# Generate SASS code for each SM architecture listed in $(SMS)
$(foreach sm,$(SMS),$(eval GENCODE_FLAGS += -gencode arch=compute_$(sm),code=sm_$(sm)))


# Generate PTX code from the highest SM architecture in $(SMS) to guarantee forward-compatibility
HIGHEST_SM := $(lastword $(sort $(SMS)))
ifneq ($(HIGHEST_SM),)
GENCODE_FLAGS += -gencode arch=compute_$(HIGHEST_SM),code=compute_$(HIGHEST_SM)
endif
endif

CUDA_LIBRARIES += -L/usr/local/cuda/lib64 -lcudart -lcufft

NVCC_INCLUDES := -I/home/jrmales/include -DMXLIB_MKL

%.o : %.cu
	$(EXEC) $(NVCC) $(NVCC_INCLUDES) $(NVCCFLAGS) $(GENCODE_FLAGS) -c $<

LDFLAGS += $(CUDA_LDFLAGS) 
LDLIBS += $(CUDA_LIBRARIES)

#fhi: fhi.o
#	g++ -o fhi fhi.o $(CUDA_LDFLAGS) $(CUDA_LIBRARIES) -L/home/jrmales/lib -lfftw3_threads -lfftw3f_threads -lfftw3l_threads -lfftw3 -lfftw3f  -lfftw3l -lmxlib -lsofa_c -llevmar -lcfitsio -lboost_system -lboost_filesystem -lgsl

