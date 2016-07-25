FLINT_DIR=..

paramul_run_all:

include $(FLINT_DIR)/Makefile

#CFLAGS += -O0 -ggdb -DWANT_ASSERT=1
INCS += -I$(FLINT_DIR)
LIBS += -L$(FLINT_DIR)

PARAMUL_BUILD_DIR = build

PARAMUL_PROGS = t-pmul.c

PARAMUL_FUNCS = fmpz_poly_pmul.c

PARAMUL_EXES = $(PARAMUL_PROGS:%.c=$(PARAMUL_BUILD_DIR)/%)

PARAMUL_OBJS = $(PARAMUL_FUNCS:%.c=$(PARAMUL_BUILD_DIR)/%.o)

PARAMUL_RUNS = $(PARAMUL_PROGS:%.c=%)

$(PARAMUL_OBJS): $(PARAMUL_BUILD_DIR)/%.o: %.c
	$(AT)mkdir -p $(PARAMUL_BUILD_DIR)
	$(CC) -c $(CFLAGS) $(INCS) $^ -o $@

$(PARAMUL_EXES): $(PARAMUL_BUILD_DIR)/%: %.c $(PARAMUL_OBJS)
	$(AT)mkdir -p $(PARAMUL_BUILD_DIR)
	$(CC) $(CFLAGS) $(INCS) $^ -o $@ $(LIBS)

$(PARAMUL_RUNS): %: $(PARAMUL_BUILD_DIR)/%
	$(AT)echo "Running $<..."
	$(AT)LD_LIBRARY_PATH=$(FLINT_DIR) $<
	$(AT)echo

paramul_run_all: $(PARAMUL_RUNS)

.PHONY: paramul_run_all $(PARAMUL_RUNS)
