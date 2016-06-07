#
# Unittests for MOTAC related functions
# 

# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR = ../lib/googletest

# Where to find test code
USER_DIR = .

SRC_DIR = ../src

INC := -I ../include

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include

# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wextra -pthread -lm -lblas -std=c++0x -w

BLASFLAGS += -lm -lblas

CC := g++
# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.


TESTS = operators_unittest.exe algebra_unittest.exe splitting_schemes_unittest.exe test_libsvm_loader.exe benchmark_test.exe objective_function_unittest.exe cmp_sync_lasso.exe


# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# House-keeping build targets.

all : $(TESTS)

lasso : cmp_sync_lasso.exe

clean :
	rm -f $(TESTS) gtest.a gtest_main.a *.o

# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CC) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CC) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Builds a sample test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.

algebra.o : $(SRC_DIR)/algebra.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(INC) -c $(SRC_DIR)/algebra.cc

nist_spblas.o : $(SRC_DIR)/nist_spblas.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(INC) -c $(SRC_DIR)/nist_spblas.cc

util.o: $(SRC_DIR)/util.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(INC) -c $(SRC_DIR)/util.cc


operators_unittest.o : $(USER_DIR)/operators_unittest.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(INC) -c $(USER_DIR)/operators_unittest.cc
operators_unittest.exe : operators_unittest.o algebra.o nist_spblas.o gtest_main.a util.o
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@ $(BLASFLAGS)


algebra_unittest.o : $(USER_DIR)/algebra_unittest.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(INC) -c $(USER_DIR)/algebra_unittest.cc
algebra_unittest.exe : algebra_unittest.o algebra.o nist_spblas.o gtest_main.a
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@ $(BLASFLAGS)


splitting_schemes_unittest.o : $(USER_DIR)/splitting_schemes_unittest.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(INC) -c $(USER_DIR)/splitting_schemes_unittest.cc
splitting_schemes_unittest.exe : splitting_schemes_unittest.o algebra.o nist_spblas.o gtest_main.a
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@ $(BLASFLAGS)


test_libsvm_loader.o : $(USER_DIR)/test_libsvm_loader.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(INC) -c $(USER_DIR)/test_libsvm_loader.cc
test_libsvm_loader.exe: test_libsvm_loader.o util.o algebra.o nist_spblas.o
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@ $(BLASFLAGS)


benchmark_test.o : $(USER_DIR)/benchmark_test.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(INC) -c $(USER_DIR)/benchmark_test.cc
benchmark_test.exe: benchmark_test.o util.o algebra.o nist_spblas.o
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@ $(BLASFLAGS)

objective_function_unittest.o : $(USER_DIR)/objective_function_unittest.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(INC) -c $(USER_DIR)/objective_function_unittest.cc
objective_function_unittest.exe : objective_function_unittest.o algebra.o util.o nist_spblas.o gtest_main.a
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@ $(BLASFLAGS)

cmp_sync_lasso.o : $(USER_DIR)/cmp_sync_lasso.cc
	$(CC) $(CPPFLAGS) $(CXXFLAGS) $(INC) -c $(USER_DIR)/cmp_sync_lasso.cc

cmp_sync_lasso.exe: cmp_sync_lasso.o util.o algebra.o nist_spblas.o
	$(CC) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@ $(BLASFLAGS)
