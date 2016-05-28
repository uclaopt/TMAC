# Makefile for MOTAC

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
# directory of test files
TESTDIR := test
# directory of application files
APPSDIR := apps

DEPENDENCY := $(shell find $(BUILDDIR) -type f -name *.d 2>/dev/null)

# app for applying gradient descent for solving least squares
MOTAC_GD_LS := $(BINDIR)/motac_gd_ls
# app for applying gradient descent for minimizing huber loss
MOTAC_GD_HB := $(BINDIR)/motac_gd_huber
# app for jacobi method for solving linear equations
MOTAC_JACOBI := $(BINDIR)/motac_jacobi
# app for applying FBS for solving lasso
MOTAC_FBS_LASSO := $(BINDIR)/motac_fbs_lasso
# app for applying FBS for solving l2_svm
MOTAC_FBS_L2_SVM := $(BINDIR)/motac_fbs_l2_svm
# app for applying FBS for solving l2_svm
MOTAC_FBS_DUAL_SVM := $(BINDIR)/motac_fbs_dual_svm
# app for applying FBS for solving l1 regularized logistic regression
MOTAC_FBS_L1LOG := $(BINDIR)/motac_fbs_l1_log
# app for applying BFS for solving l2 ball constrained quadratic programming
MOTAC_BFS_L2BALL_QP := $(BINDIR)/motac_bfs_l2_ball_qp
# app for demo PRS for finding the intersection of two sets
MOTAC_PRS_DEMO := $(BINDIR)/motac_prs_demo

MOTAC_3S_PORTFOLIO := $(BINDIR)/motac_3s_portfolio

TEST_LIBSVM := $(BINDIR)/test_libsvm
TEST_BLAS := $(BINDIR)/test_blas
TEST_BENCHMARK := $(BINDIR)/test_benchmark

CFLAGS := -g -std=c++0x -MMD -w


LIB := -lblas -lgfortran -lpthread
INC := -I include

all:  $(MOTAC_GD_LS) $(MOTAC_GD_HB) $(MOTAC_FBS_LASSO) $(MOTAC_FBS_L2_SVM) $(MOTAC_FBS_L1LOG) $(MOTAC_BFS_L2BALL_QP) $(MOTAC_PRS_DEMO) $(MOTAC_JACOBI) $(MOTAC_3S_PORTFOLIO) $(MOTAC_FBS_DUAL_SVM)


# APPS
##########################################################
# App1: gradient descent method for least square problem
$(MOTAC_GD_LS): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/motac_gd_least_square.o
	@echo " $(CC) $^ -o $(MOTAC_GD_LS) $(LIB)"; $(CC) $^ -o $(MOTAC_GD_LS) $(LIB)
	@echo " $(MOTAC_GD_LS) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'


$(MOTAC_GD_HB): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/motac_gd_huber_loss.o
	@echo " $(CC) $^ -o $(MOTAC_GD_HB) $(LIB)"; $(CC) $^ -o $(MOTAC_GD_HB) $(LIB)
	@echo " $(MOTAC_GD_HB) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'


$(MOTAC_JACOBI): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/motac_jacobi_linear_eqn.o
	@echo " $(CC) $^ -o $(MOTAC_JACOBI) $(LIB)"; $(CC) $^ -o $(MOTAC_JACOBI) $(LIB)
	@echo " $(MOTAC_JACOBI) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'


$(MOTAC_FBS_LASSO): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/motac_fbs_lasso.o
	@echo " $(CC) $^ -o $(MOTAC_FBS_LASSO) $(LIB)"; $(CC) $^ -o $(MOTAC_FBS_LASSO) $(LIB)
	@echo " $(MOTAC_FBS_LASSO) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

$(MOTAC_FBS_L2_SVM): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/motac_fbs_l2_svm.o
	@echo " $(CC) $^ -o $(MOTAC_FBS_L2_SVM) $(LIB)"; $(CC) $^ -o $(MOTAC_FBS_L2_SVM) $(LIB)
	@echo " $(MOTAC_FBS_L2_SVM) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

$(MOTAC_FBS_DUAL_SVM): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/motac_fbs_dual_svm.o
	@echo " $(CC) $^ -o $(MOTAC_FBS_DUAL_SVM) $(LIB)"; $(CC) $^ -o $(MOTAC_FBS_DUAL_SVM) $(LIB)
	@echo " $(MOTAC_FBS_DUAL_SVM) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'


$(MOTAC_FBS_L1LOG): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/motac_fbs_l1_log.o
	@echo " $(CC) $^ -o $(MOTAC_FBS_L1LOG) $(LIB)"; $(CC) $^ -o $(MOTAC_FBS_L1LOG) $(LIB)
	@echo " $(MOTAC_FBS_L1LOG) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

$(MOTAC_BFS_L2BALL_QP): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/motac_bfs_l2_ball_qp.o
	@echo " $(CC) $^ -o $(MOTAC_BFS_L2BALL_QP) $(LIB)"; $(CC) $^ -o $(MOTAC_BFS_L2BALL_QP) $(LIB)
	@echo " $(MOTAC_BFS_L2BALL_QP) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

$(MOTAC_PRS_DEMO): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/motac_prs_demo.o
	@echo " $(CC) $^ -o $(MOTAC_PRS_DEMO) $(LIB)"; $(CC) $^ -o $(MOTAC_PRS_DEMO) $(LIB)
	@echo " $(MOTAC_PRS_DEMO) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

$(MOTAC_3S_PORTFOLIO): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/motac_3s_portfolio.o
	@echo " $(CC) $^ -o $(MOTAC_3S_PORTFOLIO) $(LIB)"; $(CC) $^ -o $(MOTAC_3S_PORTFOLIO) $(LIB)
	@echo " $(MOTAC_3S_PORTFOLIO) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'


# Compile code to objective files
###################################
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) $(BINDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(BUILDDIR)/%.o: $(APPSDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) $(BINDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<


run:
	./$(MOTAC_PRS_DEMO) -problem_size 1000 -nthread 2 -epoch 100
	./$(MOTAC_FBS_L1LOG) -data ./data/rcv1_train.svm -epoch 100 -lambda 1. -nthread 2
	./$(MOTAC_GD_LS) -data ./data/ds_large_A.mtx -label ./data/ds_large_b.mtx -epoch 100 -nthread 2
	./$(MOTAC_FBS_LASSO) -data ./data/ds_large_A.mtx -label ./data/ds_large_b.mtx -epoch 100 -nthread 2 -lambda 1.
	./$(MOTAC_FBS_L2_SVM) -data ./data/ds_large_A.mtx -label ./data/ds_large_b.mtx -epoch 100 -nthread 2 -lambda 1.
	./$(MOTAC_FBS_DUAL_SVM) -data ./data/ds_large_A.mtx -label ./data/ds_large_b.mtx -epoch 100 -nthread 2 -lambda 1.
	./$(MOTAC_BFS_L2BALL_QP) -problem_size 1000 -nthread 2 -epoch 100


# clean up the executables and objective fils
##############################################
clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(BINDIR)"; $(RM) -r $(BUILDDIR) $(BINDIR)

-include $(DEPENDENCY)

