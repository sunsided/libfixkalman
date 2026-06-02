.. image:: https://raw.githubusercontent.com/sunsided/libfixkalman/static/kalman.png
   :align: right

Fixed point Kalman filter library
=================================

.. image:: https://github.com/sunsided/libfixkalman/actions/workflows/ci.yml/badge.svg?branch=main
   :target: https://github.com/sunsided/libfixkalman/actions/workflows/ci.yml
   :alt: CI

libfixkalman is a Kalman filter computation library for microcontrollers.
It is based on the libfixmatrix_ and libfixmath_ libraries, which use 16.16 bit fixed point values.
The main focus is processors without an FPU, such as ARM Cortex-M0 or M3.

A 🦀 Rust variant, albeit not a direct port, is available at `sunsided/minikalman-rs`_.

----

Matrix inversion in the correction step is implemented using Cholesky decomposition and an optimized
inversion algorithm ported from EJML_.

See `function reference`_ for further details and `example_gravity.c`_ for example code.


.. _libfixmath: http://code.google.com/p/libfixmath/
.. _libfixmatrix: https://github.com/PetteriAimonen/libfixmatrix
.. _EJML: https://code.google.com/p/efficient-java-matrix-library/
.. _function reference: https://github.com/sunsided/libfixkalman/blob/master/FUNCTIONS.rst
.. _`example_gravity.c`: https://github.com/sunsided/libfixkalman/blob/master/example_gravity.c
.. `sunsided/minikalman-rs`: https://github.com/sunsided/minikalman-rs

Building
--------

The library depends on libfixmath_ and libfixmatrix_. These are not published on
any current package remote, so the CMake build fetches them straight from GitHub
(pinned commits) via ``FetchContent`` - no package manager required.

With CMake::

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --parallel
    ctest --test-dir build --output-on-failure

Or with go-task_ (see ``Taskfile.dist.yaml``)::

    task build      # configure + compile
    task example    # run the gravity example
    task test       # run the ctest suite

Compile-time options (pass as ``-D<option>=ON``, or via ``task build -- -D...``):

============================  =========================================================
Option                        Effect
============================  =========================================================
``FIXMATRIX_MAX_SIZE`` (=8)   Max matrix dimension; >= #states / #inputs / #measurements
``KALMAN_JOSEPH_FORM``        Joseph-form covariance update
``KALMAN_TIME_VARYING``       Covariance prediction uses ``B*Q*B'``
``KALMAN_DISABLE_C``          Drop functions for systems with control inputs
``KALMAN_DISABLE_UC``         Drop functions for systems without control inputs
``KALMAN_DISABLE_LAMBDA``     Drop certainty (lambda) tuning
============================  =========================================================

Continuous integration builds and runs the example across these configurations;
see ``.github/workflows/ci.yml``.

.. _go-task: https://taskfile.dev/
