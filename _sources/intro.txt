Introduction
**************


Overview
===========
TMAC is a toolbox written in C++11 that implements algorithms based on a set of modern methods for large-scale optimization. It covers a variety of optimization problems, which can be both smooth and nonsmooth, convex and nonconvex, as well as constrained and unconstrained.

The algorithms implemented in TMAC, such as the coordinate update method and operator splitting method, are scalable as they decompose a problem into simple subproblems. These algorithms can run in a multi-threaded fashion, either synchronously or asynchronously, to take advantages of all the cores available.

TMAC is separated into several layers, and this architecture mimics how a scientist writes down an optimization algorithm. Therefore, it is easy for one to obtain a new algorithm by making simple modifications such as adding a new operator and adding a new splitting, while maintaining the multicore parallelism and other features.

  * A rich set of predefined operators.
  * A few predefined operator splitting schemes.
  * An asynchronous coordinate update framework.
  * Options to load datasets in various format.
  * A rich set of applications build on top of TMAC.
  * Support for Matlab, Python and Julia (coming soon)

  
TMAC is written in C++. Python, Julia and Matlab interfaces are under development. The following graph is an overview of the design.

.. image:: arch.png
    :width: 400px
    :align: center




Dependencies
==============
* Functioning C++ compilers (`gcc <https://www.gnu.org/software/gcc/releases.html>`_)
* BLAS library  
* `GNU make <https://www.gnu.org/software/make/>`_
* or Microsoft Visual C++


Related Papers
==============
.. [TMAC2016] Brent Edmunds, Zhimin Peng, Wotao Yin, *TMAC: A Toolbox of Modern Async-Parallel, Coordinate, Splitting, and Stochastic Methods*, UCLA-CAM report 16-38 (2016)
.. [AROCK2015] Zhimin Peng, Yangyang Xu, Ming Yan, Wotao Yin, *ARock: an Algorithmic Framework for Asynchronous Parallel Coordinate Updates*,  arXiv preprint arXiv:1506.02396 (2015)
.. [CF2016] Zhimin Peng, Tianyu Wu, Yangyang Xu, Ming Yan, Wotao Yin, *Coordinate Friendly Structures, Algorithms and applications*, Annals of Mathematical Sciences and Applications, 1(1), pp. 57–119, 2016. 
  

License and copyright
=====================
Except for the Eigen library header files, all files distributed with TMAC are made available under the `New BSD license <http://www.opensource.org/licenses/bsd-license.php>`_,
which states::

    TMAC is a toolbox of modern async-parallel, coordinate, splitting, and stochastic methods.
    Copyright (C) 2016 Brent Edmunds, Zhimin Peng, Wotao Yin 

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.


Note that Eigen is `MPL2-licensed <https://www.mozilla.org/MPL/2.0/>`_.
