<img align="right" src="https://raw.githubusercontent.com/sunsided/libfixkalman/static/kalman.png" alt="" />

# Fixed point Kalman filter library

[![CI](https://github.com/sunsided/libfixkalman/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/sunsided/libfixkalman/actions/workflows/ci.yml)

libfixkalman is a Kalman filter computation library for microcontrollers.
It is based on the [libfixmatrix][libfixmatrix] and [libfixmath][libfixmath] libraries, which use 16.16 bit fixed point values.
The main focus is processors without an FPU, such as ARM Cortex-M0 or M3.

A 🦀 Rust variant, albeit not a direct port, is available at [sunsided/minikalman-rs](https://github.com/sunsided/minikalman-rs).

---

Matrix inversion in the correction step is implemented using Cholesky decomposition and an optimized
inversion algorithm ported from [EJML][ejml].

See the [function reference][function-reference] for further details and [`example_gravity.c`][example] for example code.

## Building

The library depends on [libfixmath][libfixmath] and [libfixmatrix][libfixmatrix]. These are not published on
any current package remote, so the CMake build fetches them straight from GitHub
(pinned commits) via `FetchContent` - no package manager required.

With CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Or with [go-task](https://taskfile.dev/) (see `Taskfile.dist.yaml`):

```sh
task build      # configure + compile
task example    # run the gravity example
task test       # run the ctest suite
```

Compile-time options (pass as `-D<option>=ON`, or via `task build -- -D...`):

| Option | Effect |
| ------ | ------ |
| `FIXMATRIX_MAX_SIZE` (=8) | Max matrix dimension; >= #states / #inputs / #measurements |
| `KALMAN_JOSEPH_FORM` | Joseph-form covariance update |
| `KALMAN_TIME_VARYING` | Covariance prediction uses `B*Q*B'` |
| `KALMAN_DISABLE_C` | Drop functions for systems with control inputs |
| `KALMAN_DISABLE_UC` | Drop functions for systems without control inputs |
| `KALMAN_DISABLE_LAMBDA` | Drop certainty (lambda) tuning |

Continuous integration builds and runs the example across these configurations;
see [`.github/workflows/ci.yml`](.github/workflows/ci.yml).

[libfixmath]: http://code.google.com/p/libfixmath/
[libfixmatrix]: https://github.com/PetteriAimonen/libfixmatrix
[ejml]: https://code.google.com/p/efficient-java-matrix-library/
[function-reference]: FUNCTIONS.md
[example]: https://github.com/sunsided/libfixkalman/blob/main/example_gravity.c
