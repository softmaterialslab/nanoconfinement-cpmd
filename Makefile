#This make file builds the sub folder make files
PROG = confined_ions
JOBSCR = iu_cluster_job_script.pbs
TESTSCR = test.pbs
BIN = bin
BASE = src
SCRIPT = scripts

all:
	@echo "Starting build of the $(BASE) directory";
ifeq ($(CCF),BigRed2)	
	+$(MAKE) -C $(BASE) cluster-install
else ifeq ($(CCF),nanoHUB)
	+$(MAKE) -C $(BASE) nanoHUB-install
else
	+$(MAKE) -C $(BASE) install
endif
	@echo "Ending the build of the $(BASE) directory";
	@echo "installing the $(PROG) into $(BIN) directory"; cp -f $(BASE)/$(PROG) $(BIN)

local-install:  create-dirs
	make all

cluster-install: create-dirs
	make CCF=BigRed2 all


nanoHUB-install: create-dirs
	make CCF=nanoHUB all

create-dirs:
	@echo "Checking and creating needed sub-directories in the $(BIN) directory"
	if ! [ -d $(BIN) ]; then mkdir $(BIN); fi
	if ! [ -d $(BIN)/outfiles ]; then mkdir $(BIN)/outfiles; fi
	@echo "Directory creation is over."

cluster-submit:
	@echo "Installing jobscript into $(BIN) directory"
	cp -f $(SCRIPT)/$(JOBSCR) $(BIN)
	+$(MAKE) -C $(BIN) submit

cluster-test-submit:
	@echo "Installing test jobscript into $(BIN) directory"
	cp -f $(SCRIPT)/$(TESTSCR) $(BIN)
	+$(MAKE) -C $(BIN) test

clean:
	rm -f $(BASE)/*.o
	rm -f $(BASE)/$(PROG)
	rm -f $(BIN)/$(PROG)

dataclean:
	rm -f $(BIN)/outfiles/*.dat $(BIN)/outfiles/*.xyz  $(BIN)/outfiles/*.lammpstrj
	rm -f $(BIN)/*.log
	rm -f $(BIN)/*.pbs

.PHONY: all clean
