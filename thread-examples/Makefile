FLINT_DIR=..

openmp_run_all:

include $(FLINT_DIR)/Makefile

CFLAGS += -fopenmp
INCS += -I$(FLINT_DIR)
LIBS += -L$(FLINT_DIR)

OMP_BUILD_DIR = build

OMP_PROGS = t-tshift-single.c t-tshift-pthreads.c t-tshift-openmp.c

OMP_FUNCS = taylor_shift_multi_mod_openmp.c

OMP_EXES = $(OMP_PROGS:%.c=$(OMP_BUILD_DIR)/%)

OMP_OBJS = $(OMP_FUNCS:%.c=$(OMP_BUILD_DIR)/%.o)

OMP_RUNS = $(OMP_PROGS:%.c=%)

$(OMP_OBJS): $(OMP_BUILD_DIR)/%.o: %.c
	$(AT)mkdir -p $(OMP_BUILD_DIR)
	$(CC) -c $(CFLAGS) $(INCS) $^ -o $@

$(OMP_EXES): $(OMP_BUILD_DIR)/%: %.c $(OMP_OBJS)
	$(AT)mkdir -p $(OMP_BUILD_DIR)
	$(CC) $(CFLAGS) $(INCS) $^ -o $@ $(LIBS)

$(OMP_RUNS): %: $(OMP_BUILD_DIR)/%
	$(AT)echo "Running $<..."
	$(AT)LD_LIBRARY_PATH=$(FLINT_DIR) $<
	$(AT)echo

openmp_run_all: $(OMP_RUNS)

.PHONY: openmp_run_all $(OMP_RUNS)
