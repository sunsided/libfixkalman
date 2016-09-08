.. image:: https://raw.githubusercontent.com/sunsided/libfixkalman/static/kalman.png
   :align: right

Fixed point Kalman filter library
=================================

libfixkalman is a Kalman filter computation library for microcontrollers.
It is based on the libfixmatrix_ and libfixmath_ libraries, which use 16.16 bit fixed point values.
The main focus is processors without an FPU, such as ARM Cortex-M0 or M3.

Matrix inversion in the correction step is implemented using Cholesky decomposition and an optimized
inversion algorithm ported from EJML_.

See `function reference`_ for further details and `example_gravity.c`_ for example code.

.. _libfixmath: http://code.google.com/p/libfixmath/
.. _libfixmatrix: https://github.com/PetteriAimonen/libfixmatrix
.. _EJML: https://code.google.com/p/efficient-java-matrix-library/
.. _function reference: https://github.com/sunsided/libfixkalman/blob/master/FUNCTIONS.rst
.. _`example_gravity.c`: https://github.com/sunsided/libfixkalman/blob/master/example_gravity.c

conan.io
--------

This library now has experimental support for the `conan.io`_ package manager and is aimed at CMake. Both ``libfixmath`` and ``libfixmatrix`` dependencies are available on conan.io and you should be able to verify the package building process by calling::

    conan test_package --build

In general, to reference the library you'd provide a ``conanfile.txt`` with the following content::

    [requires]
    libfixkalman/20161008@sunside/stable

where ``20161008`` would be replaced with the latest version as found via ``conan search -v libfixkalman*` -r=conan.io`. You can then just::

    conan install

or

    conan install --build

to obtain all required references.

.. _`conan.io`: https://conan.io/