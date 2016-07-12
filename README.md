# TMAC

Welcome to **TMAC**: A toolbox that implements a set of modern methods:

   * First-order methods: gradient descent, proximal-point, and prox-gradient algorithms
   * Operator splitting methods based on forward-backward, backward-forward, Douglas-Rachford, Peaceman-Rachford splittings
   * (Block) coordinate update: cyclic, random, parallel Gauss-Seidel index rules
   * Parallel and asynchronous parallel coordinate updates
   * Stochastic gradient methods (coming soon)

You can easily apply these methods to your application by just plugging in the functions or operators specific to your problem.

## Components

   * A rich set of operators: ...
   * Operator splitting schemes: forward-backward, backward-forward, Douglas-Rachford, Peaceman-Rachford
   * A shared-memory (asynchronous) parallel driver
   * Examples: (sparse) logistic regression, LASSO, portfolio optimization, nonnegative matrix factorization, insection of two sets
   * Supported dataset formats: matrix market, LIBSVM

## Support platforms

Linux, Mac OS X, and Windows (32 and 64 bits)

## Installation

   * Linux [(guide)](https://github.com/uclaopt/TMAC/blob/master/doc/source/build.rst#requirements-for-linux)
   * Mac OS X [(guide)](https://github.com/uclaopt/TMAC/blob/master/doc/source/build.rst#requirements-for-mac)
   * Windows [(MINGW32 (32-bit) guide, ](http://www.math.ucla.edu/~wotaoyin/software/tmac_windows_installation_mingw32.html) [MINGW64 (64-bit) guide,](http://www.math.ucla.edu/~wotaoyin/software/tmac_windows_installation_mingw64.html)  [Cygwin (32/64-bit) guide)](http://www.math.ucla.edu/~wotaoyin/software/tmac_windows_installation_cygwin.html)

## Examples

We build a user-friendly interface to run TMAC for different applications through shell. You can run it through the following command:

```
./run_me.sh
```

## Full documentation

We use [sphinx](http://www.sphinx-doc.org/en/stable/) for documentation. To use it, you need to install [sphinx](http://www.sphinx-doc.org/en/stable/install.html)
 and [sphinx_rtd_theme](https://github.com/snide/sphinx_rtd_theme). To compile the source code for the documentation, you can run the following commands::

```bash
cd doc
make html
```

Then the html files are in the ```build/html``` folder.


## Contributing Code
We welcome patches. If you plan to contribute a patch, you need to write tests for any new code. Changes should be verified to not
break existing tests before they are submitted for review. We use Google Test for unit testing. Our unit tests are in the test folder.
If you are adding tests to a new file, you will need to modify the Makefile to compile the source code, otherwise, nothing needs to be
changed in the Makefile. The following commands should build and run the unit tests.

The details for contributing to TMAC is [here](https://github.com/uclaopt/TMAC/blob/master/doc/source/contribute.rst).

The contributors' list is [here](https://github.com/uclaopt/TMAC/blob/master/CONTRIBUTORS.md).



## Getting the Data

So far, TMAC support the following data format:
   * [Matrix Market Format](http://math.nist.gov/MatrixMarket/formats.html#MMformat)
   * [LIBSVM format](https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/)

You can download some example datasets from [here](https://github.com/uclaopt/datasets)
You can also obtain the regression and classification datasets from the LIBSVM [website](https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/).



## Troubleshooting

* You can report bugs through the [issue tab](https://github.com/uclaopt/TMAC/issues/new)

* fatal error: ```<omp.h>``` file not found
  * First make sure you are using the correct compiler, i.e., GCC >= 4.7 from the GNU Project. For OSX,
  please follow the instruction [here](http://stackoverflow.com/questions/20340117/omp-h-library-isnt-found-in-the-gcc-version-4-2-1-in-mavericks) to download and install gcc-4.9.
  For Ubuntu, you can follow this [link](http://askubuntu.com/questions/428198/getting-installing-gcc-g-4-9-on-ubuntu).

* install g++ without root privileges. This [link](http://luiarthur.github.io/gccinstall) is very helpful.

## Acknowledgement

We would like to acknowledge the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library, the [sparse BLAS](http://math.nist.gov/spblas/) library, and
[Google Test](https://github.com/google/googletest).


## Related Papers

```
@article{edmunds2016tmac,
  title={TMAC: A Toolbox of Modern Async-Parallel, Coordinate, Splitting, and Stochastic Methods},
  author={Edmunds, Brent and Peng, Zhimin and Yin, Wotao},
  journal={arXiv preprint arXiv:1606.04551},
  year={2016}
}

@article{peng2015arock,
  title = {ARock: an Algorithmic Framework for Asynchronous Parallel Coordinate Updates},
  author = {Peng, Zhimin and Xu, Yangyang and Yan, Ming and Yin, Wotao},
  journal = {arXiv:1506.02396},
  year = {2015},
  publisher = {http://arxiv.org/abs/1506.02396}
}

@article{PengWuXuYanYin2016_coordinate,
  title = {Coordinate friendly structures, algorithms and applications},
  volume = {1},
  number = {1},
  journal = {Annals of Mathematical Sciences and Applications},
  author = {Peng, Zhimin and Wu, Tianyu and Xu, Yangyang and Yan, Ming and Yin, Wotao},
  month = Jan,
  year = {2016},
  pages = {59--119}
}

```
