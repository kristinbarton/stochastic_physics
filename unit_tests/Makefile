DEBUG ?= NO
BUILD_DIR = build

# Set flags based on debug mode
ifeq ($(DEBUG), YES)
	FLAGS = $(DEBUG_FFLAGS)
else
	FLAGS = $(FFLAGS)
endif

# Look in parent directory and src for files
VPATH = src:..

# ---------------------------------------------------
# Cellular Automata (standalone_ca.x)
# ---------------------------------------------------
CA_SRCS = kinddef.F90 \
		  mpi_wrapper.F90 \
		  fv_arrays_stub.F90 \
		  fv_mp_stub_mod.F90 \
          fv_control_stub.F90 \
		  atmosphere_stub.F90 \
		  halo_exchange.fv3.F90 \
		  random_numbers.F90 \
		  plumes.F90 \
		  update_ca.F90 \
          cellular_automata_sgs.F90 \
		  cellular_automata_global.F90

CA_OBJS = $(addprefix $(BUILD_DIR)/, $(CA_SRCS:.F90=.o))

$(BUILD_DIR)/libca.a: $(CA_OBJS)
	ar rv $@ $^

ca: $(BUILD_DIR)/libca.a $(BUILD_DIR)/standalone_ca.o
	$(FC) $(FLAGS) -o $(BUILD_DIR)/standalone_ca.x $(BUILD_DIR)/standalone_ca.o -L$(BUILD_DIR) -lca $(INCS) $(LIBS)

# ---------------------------------------------------
# Gaussian Grid Cellular Automata (standalone_ggca.x)
# ---------------------------------------------------
GGCA_SRCS = kinddef.F90 \
		    mpi_wrapper.F90 \
			stochy_namelist_def.F90 \
			stochy_internal_state_mod.F90 \
			spectral_transforms.F90 \
		    fv_arrays_stub.F90 \
		    fv_mp_stub_mod.F90 \
            fv_control_stub.F90 \
		    atmosphere_stub.F90 \
		    halo_exchange.fv3.F90 \
		    random_numbers.F90 \
		    plumes.F90 \
		    update_ca.F90 \
		    cellular_automata.F90

GGCA_OBJS = $(addprefix $(BUILD_DIR)/, $(GGCA_SRCS:.F90=.o))

$(BUILD_DIR)/libggca.a: $(GGCA_OBJS)
	ar rv $@ $^

ggca: $(BUILD_DIR)/libggca.a $(BUILD_DIR)/standalone_ggca.o
	$(FC) $(FLAGS) -o $(BUILD_DIR)/standalone_ggca.x $(BUILD_DIR)/standalone_ggca.o -L$(BUILD_DIR) -lggca $(INCS) $(LIBS)

# ---------------------------------------------------
# General rules
# ---------------------------------------------------
# Make the build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# This makes sure BUILD_DIR is created before compiling
$(BUILD_DIR)/%.o: %.F90 | $(BUILD_DIR)
	$(FC) $(FLAGS) $(MOD_FLAG) -I$(BUILD_DIR) $(INCS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) err out
