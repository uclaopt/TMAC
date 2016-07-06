# Makefile for TMAC

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
TMAC_GD_LS := $(BINDIR)/tmac_gd_ls
# app for applying gradient descent for minimizing huber loss
TMAC_GD_HB := $(BINDIR)/tmac_gd_huber
# app for jacobi method for solving linear equations
TMAC_JACOBI := $(BINDIR)/tmac_jacobi
# app for applying FBS for solving lasso
TMAC_FBS_LASSO := $(BINDIR)/tmac_fbs_lasso
# app for applying FBS for solving l2_svm
TMAC_FBS_L2_SVM := $(BINDIR)/tmac_fbs_l2_svm
# app for applying FBS for solving l2_svm
TMAC_FBS_DUAL_SVM := $(BINDIR)/tmac_fbs_dual_svm
# app for applying FBS for solving l1 regularized logistic regression
TMAC_FBS_L1LOG := $(BINDIR)/tmac_fbs_l1_log
# app for applying BFS for solving l2 ball constrained quadratic programming
TMAC_BFS_L2BALL_QP := $(BINDIR)/tmac_bfs_l2_ball_qp
# app for demo PRS for finding the intersection of two sets
TMAC_PRS_DEMO := $(BINDIR)/tmac_prs_demo

TMAC_3S_PORTFOLIO := $(BINDIR)/tmac_3s_portfolio

TEST_LIBSVM := $(BINDIR)/test_libsvm
TEST_BLAS := $(BINDIR)/test_blas
TEST_BENCHMARK := $(BINDIR)/test_benchmark

CFLAGS := -g -std=c++0x -MMD -w


LIB := -lblas -lpthread
INC := -I include

all:  $(TMAC_GD_LS) $(TMAC_GD_HB) $(TMAC_FBS_LASSO) $(TMAC_FBS_L2_SVM) $(TMAC_FBS_L1LOG) $(TMAC_BFS_L2BALL_QP) $(TMAC_PRS_DEMO) $(TMAC_JACOBI) $(TMAC_3S_PORTFOLIO) $(TMAC_FBS_DUAL_SVM)


# APPS
##########################################################
# App1: gradient descent method for least square problem
$(TMAC_GD_LS): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/tmac_gd_least_square.o
	@echo " $(CC) $^ -o $(TMAC_GD_LS) $(LIB)"; $(CC) $^ -o $(TMAC_GD_LS) $(LIB)
	@echo " $(TMAC_GD_LS) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'


$(TMAC_GD_HB): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/tmac_gd_huber_loss.o
	@echo " $(CC) $^ -o $(TMAC_GD_HB) $(LIB)"; $(CC) $^ -o $(TMAC_GD_HB) $(LIB)
	@echo " $(TMAC_GD_HB) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'


$(TMAC_JACOBI): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/tmac_jacobi_linear_eqn.o
	@echo " $(CC) $^ -o $(TMAC_JACOBI) $(LIB)"; $(CC) $^ -o $(TMAC_JACOBI) $(LIB)
	@echo " $(TMAC_JACOBI) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'


$(TMAC_FBS_LASSO): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/tmac_fbs_lasso.o
	@echo " $(CC) $^ -o $(TMAC_FBS_LASSO) $(LIB)"; $(CC) $^ -o $(TMAC_FBS_LASSO) $(LIB)
	@echo " $(TMAC_FBS_LASSO) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

$(TMAC_FBS_L2_SVM): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/tmac_fbs_l2_svm.o
	@echo " $(CC) $^ -o $(TMAC_FBS_L2_SVM) $(LIB)"; $(CC) $^ -o $(TMAC_FBS_L2_SVM) $(LIB)
	@echo " $(TMAC_FBS_L2_SVM) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

$(TMAC_FBS_DUAL_SVM): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/tmac_fbs_dual_svm.o
	@echo " $(CC) $^ -o $(TMAC_FBS_DUAL_SVM) $(LIB)"; $(CC) $^ -o $(TMAC_FBS_DUAL_SVM) $(LIB)
	@echo " $(TMAC_FBS_DUAL_SVM) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'


$(TMAC_FBS_L1LOG): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/tmac_fbs_l1_log.o
	@echo " $(CC) $^ -o $(TMAC_FBS_L1LOG) $(LIB)"; $(CC) $^ -o $(TMAC_FBS_L1LOG) $(LIB)
	@echo " $(TMAC_FBS_L1LOG) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

$(TMAC_BFS_L2BALL_QP): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/tmac_bfs_l2_ball_qp.o
	@echo " $(CC) $^ -o $(TMAC_BFS_L2BALL_QP) $(LIB)"; $(CC) $^ -o $(TMAC_BFS_L2BALL_QP) $(LIB)
	@echo " $(TMAC_BFS_L2BALL_QP) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

$(TMAC_PRS_DEMO): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/tmac_prs_demo.o
	@echo " $(CC) $^ -o $(TMAC_PRS_DEMO) $(LIB)"; $(CC) $^ -o $(TMAC_PRS_DEMO) $(LIB)
	@echo " $(TMAC_PRS_DEMO) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

$(TMAC_3S_PORTFOLIO): build/algebra.o build/util.o build/matrices.o build/nist_spblas.o build/tmac_3s_portfolio.o
	@echo " $(CC) $^ -o $(TMAC_3S_PORTFOLIO) $(LIB)"; $(CC) $^ -o $(TMAC_3S_PORTFOLIO) $(LIB)
	@echo " $(TMAC_3S_PORTFOLIO) is successfully built."
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
	./$(TMAC_PRS_DEMO) -problem_size 1000 -nthread 2 -epoch 100
	./$(TMAC_FBS_L1LOG) -data ./data/rcv1_train.svm -epoch 100 -lambda 1. -nthread 2
	./$(TMAC_GD_LS) -data ./data/ds_large_A.mtx -label ./data/ds_large_b.mtx -epoch 100 -nthread 2
	./$(TMAC_FBS_LASSO) -data ./data/ds_large_A.mtx -label ./data/ds_large_b.mtx -epoch 100 -nthread 2 -lambda 1.
	./$(TMAC_FBS_L2_SVM) -data ./data/ds_large_A.mtx -label ./data/ds_large_b.mtx -epoch 100 -nthread 2 -lambda 1.
	./$(TMAC_FBS_DUAL_SVM) -data ./data/ds_large_A.mtx -label ./data/ds_large_b.mtx -epoch 100 -nthread 2 -lambda 1.
	./$(TMAC_BFS_L2BALL_QP) -problem_size 1000 -nthread 2 -epoch 100


# clean up the executables and objective fils
##############################################
clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(BINDIR)"; $(RM) -r $(BUILDDIR) $(BINDIR)

-include $(DEPENDENCY)

