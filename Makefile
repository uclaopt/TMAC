# Makefile for TMAC_SOCP

CC := g++ # This is the main compiler
# CC := clang --analyze # and comment out the linker last line for sanity
# directory of objective files
BUILDDIR := build
# directory of binary files
BINDIR := bin
# directory of source code
SRCDIR := src
# extension of source file
SRCEXT := cc
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

DEPENDENCY := $(shell find $(BUILDDIR) -type f -name *.d 2>/dev/null)

# app for PRS solving SOCP
TMAC_PRS_SOCP := $(BINDIR)/tmac_prs_socp

CFLAGS := -g -std=c++0x -MMD -w


LIB := -lblas -lgfortran -lpthread -llapack
INC := -I include

all: $(TMAC_PRS_SOCP)

$(TMAC_PRS_SOCP): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/tmac_prs_socp.o
	@echo " $(CC) $^ -o $(TMAC_PRS_SOCP) $(LIB)"; $(CC) $^ -o $(TMAC_PRS_SOCP) $(LIB)
	@echo " $(TMAC_PRS_SOCP) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

# Compile code to objective files
###################################
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) $(BINDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

run:
	./$(TMAC_PRS_SOCP) -problem_size 400 -nthread 4 -epoch 100 -A data/A.mtx -b data/b.mtx -c data/c.mtx -alpha 0.5 -gamma 2 -lambda 1.8 -cone_num 2 -cone_dim 200 200 
	
# clean up the executables and objective fils
##############################################
clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(BINDIR)"; $(RM) -r $(BUILDDIR) $(BINDIR)

-include $(DEPENDENCY)

