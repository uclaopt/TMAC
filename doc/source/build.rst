Download and install
==========================
Latest version
-----------------------
The current version of the entire TMAC package can be downloaded at `https://github.com/uclaopt/TMAC <https://github.com/uclaopt/TMAC>`_. The package includes all source codes, help documents, and demos. Some demos require external datasets at `our dataset center <https://github.com/uclaopt/datasets>`_ and `the LIBSVM website <https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/>`_.

In the future, precompiled executables for Windows users will be available for download.

Required and optional software for installation
-------------------------------------------------
- A C++11 compiler (e.g., gcc version 4.7 or higher)
- BLAS
- (recommended) `GNU make <https://www.gnu.org/software/make/>`_

Linux installation
-----------------------
1. Install the packages of coding essentials, BLAS, and LAPACK

.. code-block:: bash

  sudo apt-get install build-essential
  sudo apt-get install libblas-dev liblapack-dev

2. Download TMAC and build it:

.. code-block:: bash

  git clone https://github.com/uclaopt/tmac.git
  cd tmac
  make

3. You can test run a Peaceman-Rachford Splitting example, first with 1 thread and then with 2 threads:

.. code-block:: bash

  ./bin/tmac_prs_demo -problem_size 1500 -nthread 1
  ./bin/tmac_prs_demo -problem_size 1500 -nthread 2

Mac inistallation
--------------------
1. Install XCode: in App Store, update Xcode; then, launch Xcode and accept its license terms
2. Install XCode's command line tools:

.. code-block:: bash

  xcode-select --install
  gcc -v    # verify you have version >= 4.7

3. Download TMAC and build it::

.. code-block:: bash

  git clone https://github.com/uclaopt/tmac.git
  cd tmac
  make

4. Once TMAC is successfully built, you can test run a Peaceman-Rachford Splitting example, first with 1 thread and then with 2 threads::

.. code-block:: bash

  ./bin/tmac_prs_demo -problem_size 1500 -nthread 1
  ./bin/tmac_prs_demo -problem_size 1500 -nthread 2

Windows installation
-------------------------
Please choose one of the following approaches. Any one of them will set up a coding environment to install the required software and build TMAC.

- Cygwin: `32/64-bit step-by-step installation <http://www.math.ucla.edu/~wotaoyin/software/tmac_windows_installation_cygwin.html>`_;
- MinGW:  `32-bit installation <http://www.math.ucla.edu/~wotaoyin/software/tmac_windows_installation_mingw32.html>`_ and `64-bit installation <http://www.math.ucla.edu/~wotaoyin/software/tmac_windows_installation_mingw64.html>`_.
- Visual Studio: (under construction)

What to do next?
----------------
- Run the examples: `least squares <uclaopt.github.io/opt/least_squares.html>`_, `regularized regression <uclaopt.github.io/opt/regression.html>`_, and `more <uclaopt.github.io/optimization.html>`.
- Build your own algorithms.

