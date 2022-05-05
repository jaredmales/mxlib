-include local/MxLib.mk

# Must include path to the include directory, and to sofa
INCLUDES += -I../include

# SOFA
CURRENT_SOFA ?= 20210125
SOFA_PATH ?= $(abspath $(SELF_DIR)/../source/vendor/sofa/$(CURRENT_SOFA)/c/src/)
EXTRA_LDFLAGS += -L$(SOFA_PATH)
INCLUDES += -I$(SOFA_PATH)
SOFA_LIB = -lsofa_c