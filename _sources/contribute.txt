Contribute to TMAC
********************
We welcome patches. If you plan to contribute a new module, you can consider to add testcases in `tests <https://github.com/uclaopt/TMAC/tree/master/test>`_. Changes should be verified to not break existing tests before they are submitted for review. We use Google test for unit testing.

Setup the repository
=====================
1. Fork the TMAC `repository <https://github.com/uclaopt/TMAC>`_ reporsitory;
2. Open a Terminal to clone the TMAC project to your machine::

     git clone https://github.com/YOUR-USERNAME/TMAC

3. Config the remote repository::

     git remote add upstream https://github.com/uclaopt/TMAC

4. Sync with the up-to-date code::

     git pull upstream master

Create pull request (PR)
========================
1. You first want to do code development on a branch::

     git branch YOUR-BRANCH-NAME
     git checkout YOUR-BRANCH-NAME

2. After you have done some coding, you can stage and commit your code locally::

     git add FILE_NAME
     git commit -m "MESSAGE TO DESCRIBE YOUR CHANGES"

3. Now you can push the update to GitHub::
     
     git push upstream YOUR-BRANCH-NAME

4. Go to the remote `repo <https://github.com/uclaopt/TMAC>`_ page, you will see a message that ask you to create a pull request (PR).

5. Create the PR, and add sufficient information to describe your modifications.


Combine multiple commits into one
=================================
Sometimes we want to combine multiple small commits, especially when later commits are only fixes to previous ones.
Assume we want to merge the last 3 commits, you can do it with the following command::

  git rebase -i HEAD~3

You will see a pop up editor, then you can replace "pick" on the second and subsequent commits with "squash".

Unit tests
==========
All the unit tests are in test folder. We use gtest for C++ unittests. We will use Python nose for Python test cases. If you are adding tests to a new file, you will need to modify the Makefile to compile the source code, otherwise, nothing needs to be changed in the Makefile. The following commands should build and run the unit tests::

  make
  ./your_unittest

  
Update docs
============
If you add new features to the codebase, you are also required to update the documentation to reflect the new feature. We use `sphinx <http://www.sphinx-doc.org/en/stable/>`_
and `sphinx_rtd_theme <https://github.com/snide/sphinx_rtd_theme>`_ for documentation. You can find the installation instructions on their websites. You can build document locally with the following commands::

  cd doc
  make html


Applications
=============
Applications and demonstration examples are in the `apps <https://github.com/uclaopt/TMAC/tree/master/apps>`_ folder. Source code for new applications can be added here.
We would love to hear your story. If you have interesting applications that can use TMAC, please feel free to add it to the toolbox.


Code style
==========
We follow `Google's C++ style guide <https://google.github.io/styleguide/cppguide.html>`_ on C++ code.


Language independent features
=============================
If you want to implement TMAC in a different programming language or add interfaces for other programming languages, the following features should be taken into consideration:

    1. Fast-indexing. The major data structures for TMAC are dense vector, dense matrix, and sparse matrix. Fast-indexing means that either the rows or the cols in a matrix should be accessed quickly.
    2. Fast level 1 BLAS. Level 1 BLAS operations are the main computations in each iteration. Fast implementation should be used.
    3. Multi-threading. The programming language should provide a threading library that has the following features:
        * spawn a given number of threads;
	* create atomic variables;
	* achieve inter-threads communication through mutex and conditional variables.
    4. Shared variables. Threads should have fast read and write access to shared variables (such as data, unknown variable, maintained variables).
    5. Profile tools. The programming language should have multi-threading profiling tools so that we can track the execution time, identify which functions are consuming the most time so that we can evaluate them for possible performance improvements. 
