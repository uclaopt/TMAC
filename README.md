# TMAC

Welcome to **TMAC**: A toolbox that implements a set of modern methods:

   * First-order methods: gradient descent, proximal-point, and prox-gradient algorithms
   * Operator splitting methods based on forward-backward, backward-forward, Douglas-Rachford, Peaceman-Rachford splittings
   * (Block) coordinate update: cyclic, random, parallel Gauss-Seidel index rules
   * Parallel and asynchronous parallel coordinate updates
   * Stochastic gradient methods (coming soon)

You can easily apply these methods to your application by just plugging in the functions or operators specific to your problem.

## Components

   * A rich set of operators: proximal operators, projection operators, and gradient operators
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

```bash
./run_me.sh
```

## Full documentation
You can view the detailed documentations [here](http://uclaopt.github.io/TMAC/).


## Develop a new algorithm with the TMAC features
Once you understand TMAC's architecture and interface, it is easy to develop a new algorithm for your problem that inherits many TMAC methods and features. A tutorial is here.


## Open collaborations
TMAC includes contributions from a number of people around the world. They are listed [here](https://github.com/uclaopt/TMAC/blob/master/CONTRIBUTORS.md).

Guidelines for contributions: [here](http://uclaopt.github.io/TMAC/contribute.html).

## Getting the data

So far, TMAC support the following data format:
   * [Matrix Market Format](http://math.nist.gov/MatrixMarket/formats.html#MMformat)
   * [LIBSVM format](https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/)

You can download some example datasets from [here](https://github.com/uclaopt/datasets).
You can also obtain the regression and classification datasets from the LIBSVM [website](https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/).


## Troubleshooting
Please refer to the [FAQ](uclaopt.github.io/TMAC/faq.html) page for troubleshooting. You can also report bugs through the [issue tab](https://github.com/uclaopt/TMAC/issues/new)


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
