Download, build and run
==========================
TMAC's build system relies on `GNU make <https://www.gnu.org/software/make/>`_. It can be easily build on Linux and Unix environments, and various versions of Microsoft Windows. A relative up-to-date C++ compiler (e.g., gcc >= 4.7) is required in all cases.

Requirements for Linux
-----------------------
You can use the following commands to install g++, and BLA::

  sudo apt-get install build-essential
  sudo apt-get install libblas-dev liblapack-dev


Requirements for Mac
--------------------
Step 1. Install XCode and the Command Line Tools for XCode.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can install them with the following commands:

  1. To install XCode, You can go to the App Store application and check "Updates". After updating Xcode, be sure to launch the Xcode application and accept the Apple license terms.

  2. To install the command line tools, you can use the following command::
       
       xcode-select --install

Step 2. Install GNU gcc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To install GNU gcc, you can follow the following steps.

   1. download gcc-4.9-bin.tar.gz download or newer from `here <http://prdownloads.sourceforge.net/hpc/gcc-4.9-bin.tar.gz>`_;

   2. cd to your downloads folder and un-gzip the archive gunzip gcc-4.9-bin.tar.gz;
	
   3. in the same folder run::

	sudo tar -xvf gcc-4.9-bin.tar -C /

      this will place new executable to /usr/local/bin

   4. add the following to ~/.bash_profile::

	export PATH=/usr/local/bin:$PATH;

   5. open new terminal and run::

	which gcc.

      This should point to /usr/local/bin/gcc.


Requirements for Windows
-------------------------
You can either use MinGW or Cygwin to setup a coding environment to run TMAC. Here we include the installation instructions.

Cygwin
^^^^^^^
You can find the instructions for installing Cygwin `here <http://www.math.ucla.edu/~wotaoyin/windows_coding_cygwin.html>`_;

MinGW
^^^^^^
You can find the instructions for installing MinGW `here <http://www.math.ucla.edu/~wotaoyin/windows_coding.html>`_;


Download TMAC
----------------
The TMAC package can be downloaded from the following link::

  https://github.com/uclaopt/TMAC
  
  
Build TMAC
----------------
On Linux or Unix machines with g++ and GNU make installed in standard locations, building TMAC can be as simple as::

  cd TMAC
  make

The executable files are in the bin folder. The following is a sample of the executable files::

  tmac_fbs_l1_log [solver for regularized logistic regression with forward backward splitting]



Test TMAC
-------------------------
Once TMAC has been successfully compiled, it is a good idea to verify that the executable files are functioning properly. To test TMAC for l1 regularized logistic regression problem, run the following commands::

  ./bin/tmac_fbs_l1_log -data ./data/rcv1_train.svm -epoch 10 -nthread 2 -lambda 1.

