SHELL=/bin/sh
LD_PRE=LD_LIBRARY_PATH=..

CC=gcc

WARNS=-ansi -pedantic -Wall -Wextra -Werror -Wno-unused-parameter -fmax-errors=4 #-fopenmp
CFLAGS=$(WARNS) -DWANT_ASSERT=1 -mpopcnt -ggdb -O0 -gdwarf
CFLAGS_PROF=$(WARNS) -Wno-unused-parameter -mpopcnt -funroll-loops -O3 -ggdb -gdwarf -fno-omit-frame-pointer
LIBS=-lflint -lmpfr -lgmp -lm

HEADERS=../*.h
SOURCES=$(wildcard ./*.c)
TEST_SOURCES=$(wildcard test/t-*.c)
OURTEST_SOURCES=$(wildcard ourtest/ot-*.c)
PROF_SOURCES=$(wildcard profile/p-*.c)
TESTS=$(patsubst test/t-%.c,%,$(TEST_SOURCES))
OURTESTS=$(patsubst ourtest/ot-%.c,%,$(OURTEST_SOURCES))
PROFS=$(patsubst profile/p-%.c,%,$(PROF_SOURCES))
TEST_RUNS=$(patsubst %,t-%,$(TESTS))
OURTEST_RUNS=$(patsubst %,ot-%,$(OURTESTS))
PROF_RUNS=$(patsubst %,p-%,$(PROFS))

BD=build

TESTS_EXTRA_SRC=../test_helpers.c
PROFS_EXTRA_SRC=../profiler.c

TEST_EXES=$(patsubst %,$(BD)/t-%,$(TESTS))
OURTEST_EXES=$(patsubst %,$(BD)/ot-%,$(OURTESTS))
PROF_EXES=$(patsubst %,$(BD)/p-%,$(PROFS))

all: tests profile ourtests

$(TEST_RUNS): %: build/%
	@echo
	@echo "RUNNING TEST $@"
	@$(LD_PRE) $<

tests: $(TEST_RUNS)

$(TEST_EXES): $(BD)/%: test/%.c $(SOURCES) $(HEADERS) $(TESTS_EXTRA_SRC)
	@echo
	@echo -n "COMPILING $@..."
	$(CC) $(CFLAGS) -I.. $< $(SOURCES) $(TESTS_EXTRA_SRC) -L.. $(LIBS) -o $@
	@echo done

$(OURTEST_RUNS): %: build/%
	@echo
	@echo "RUNNING TEST $@"
	@$(LD_PRE) $<

ourtests: $(OURTEST_RUNS)

$(OURTEST_EXES): $(BD)/%: ourtest/%.c $(SOURCES) $(HEADERS) $(TESTS_EXTRA_SRC)
	@echo
	@echo -n "COMPILING $@..."
	$(CC) $(CFLAGS) -I.. $< $(SOURCES) $(TESTS_EXTRA_SRC) -L.. $(LIBS) -o $@
	@echo done

$(PROF_RUNS): %: build/%
	@echo
	@echo "RUNNING PROFILE $@"
	@$(LD_PRE) $<

profile: $(PROF_RUNS)

$(PROF_EXES): $(BD)/%: profile/%.c $(SOURCES) $(HEADERS) $(PROFS_EXTRA_SRC)
	@echo
	@echo -n "COMPILING $@..."
	@$(CC) $(CFLAGS_PROF) -I.. $< $(SOURCES) $(PROFS_EXTRA_SRC) -L.. $(LIBS) -o $@
	@echo done

clean:
	rm -f $(BD)/*

.PHONY: clean all tests profile $(TEST_RUNS) $(PROF_RUNS)
