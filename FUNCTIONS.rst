================================
libfixkalman: Function reference
================================

.. contents ::

Requirements
============

libfixkalman relies on `libfixmath <https://code.google.com/p/libfixmath/>`_ and `libfixmatrix <https://github.com/PetteriAimonen/libfixmatrix>`_ for all calculations.

libfixmatrix requires ``FIXMATRIX_MAX_SIZE`` to be defined on the project to the largest matrix size. libfixkalman adds the constraint that ``FIXMATRIX_MAX_SIZE`` must be greater than or equal to
the number of system states, inputs or observations (i.e. the largest of these values). If, for example, the system has 4 states, 1 input and 2 measured outputs (observations), ``FIXMATRIX_MAX_SIZE`` must be defined
to at least 4.

Preprocessor control
====================

In order to reduce code size, some features can be disabled by globally defining one ore multiple of the following:

:KALMAN_DISABLE_UC:         Disable filter functions for systems without control input. Set this if all your systems use control input or if you use certainty tuning.
:KALMAN_DISABLE_C:          Disable filter functions for systems with control inputs. Set this if none of your systems use control input.
:KALMAN_DISABLE_LAMBDA:     Disable filter functions with certainty tuning (lambda parameter). Set this if you do not need to control filter convergence speed.

Definitions
===========

There are two kinds of filters that can be processed with this library: Regular filters described by state transition and control input and such systems without control inputs (dubbed 'uncontrolled').
For each of these types a separate set of functions exists: Regular functions and functions postfixed with ``_uc`` (such as `kalman_predict_uc`_) for uncontrolled systems.

kalman16_t
----------
Data structure for the filter state. ::

    typedef struct {
        mf16 x; // S x 1
        mf16 A; // S x S
        mf16 P; // S x S
        mf16 u; // C x 1
        mf16 B; // S x C
        mf16 Q; // C x C
    } kalman16_t;

:x:         System state vector. Number of rows in the vector, 1 <= rows <= FIXMATRIX_MAX_SIZE.
:A:         Square system state transition model matrix. Number of rows and columns in the matrix is equal to the number of rows in the state vector ``x``.
:P:         Square system state covariance matrix. Number of rows and columns is identical to ``A``.
:u:         Input vector. Number of rows in the vector, 1 <= rows <= FIXMATRIX_MAX_SIZE.
:B:         Control input model matrix. Number of rows in the matrix is equal to the number of rows in the state vector ``x``. Number of columns in the matrix is equal to the number of rows in the input vector ``u``.
:Q:         Square contol input covariance matrix. Number of rows and columns in the matrix is equal to the number of rows in the input vector ``u``.

The filter structure can be initialized by calling

    ``kalman_filter_initialize(&filter, NUM_STATES, NUM_INPUTS);``

kalman16_uc_t
-------------
Data structure for the uncontrolled (inputless) filter state. ::

    typedef struct {
        mf16 x; // S x 1
        mf16 A; // S x S
        mf16 P; // S x S
        mf16 Q; // S x S
    } kalman16_t;

:x:         System state vector. Number of rows in the vector, 1 <= rows <= FIXMATRIX_MAX_SIZE.
:A:         Square system state transition model matrix. Number of rows and columns in the matrix is equal to the number of rows in the state vector ``x``.
:P:         Square system state covariance matrix. Number of rows and columns is identical to ``A``.
:Q:         Square system process noise matrix. Number of rows and columns is identical to ``A``.

The filter structure can be initialized by calling

    ``kalman_filter_initialize_uc(&filter, NUM_STATES);``
    
kalman16_observation_t
----------------------
Data structure for the measurement updates. ::

    typedef struct {
        mf16 z; // Z x 1
        mf16 H; // Z x S
        mf16 R; // Z x Z
    } kalman16_t;

:z:         Observation vector. Number of rows in the vector, 1 <= rows <= FIXMATRIX_MAX_SIZE.
:H:         Observation model matrix. Number of rows in the matrix is equal to the number of rows in the measurement vector ``z``. Number of columns in the matrix is equal to the number of rows in the state vector ``x``.
:R:         Square observation covariance matrix. Number of rows and columns is identical to the number of rows in the measurement vector ``z``.

The filter structure can be initialized by calling

    ``kalman_observation_initialize(&filter, NUM_STATES, NUM_OBSERVATIONS);`` 

Initialization Functions
========================

kalman_filter_initialize
------------------------
Initializes a *kalman16_t* structure::

    void kalman_filter_initialize(kalman16_t *const kf, uint_fast8_t num_states, uint_fast8_t num_inputs);

:kf:          The filter structure to initialize.
:num_states:  The number of system states.
:num_inputs:  The number of system inputs.

kalman_filter_initialize_uc
----------------------------
Initializes a *kalman16_uc_t* structure::

    void kalman_filter_initialize_uc(kalman16_uc_t *const kf, uint_fast8_t num_states);

:kf:          The filter structure to initialize.
:num_states:  The number of system states.

kalman_observation_initialize
-----------------------------
Initializes a *kalman16_observation_t* structure::

    void kalman_observation_initialize(kalman16_observation_t *const kfm, uint_fast8_t num_states, uint_fast8_t num_observations);

:kf:                The observation structure to initialize.
:num_states:        The number of system states.
:num_observations:  The number of observations.

Prediction and Update Functions (regular systems)
=================================================

kalman_predict
--------------
Kalman filter prediction (time update) step::
    
    void kalman_predict(kalman16_t *kf);

:kf:        The filter to update.

This performs a state and covariance update according to the state transition model *A* and the input model *B*. If *B* has zero dimensions, only the state transition model will be used.

This function is a thin wrapper around `kalman_predict_x`_ and `kalman_predict_P`_.
It is often more efficient to perform the state update manually instead of relying on the matrix multiplication algorithm. In this case, `kalman_predict_P`_ can be used to update the system covariance
matrix afterwards.

If input values are used, the user is required to set the values in ``kfm.u`` prior to calling this function.

kalman_predict_tuned
--------------------
Kalman filter prediction (time update) step with applied certainty tuning::
    
    void kalman_predict_tuned(kalman16_t *kf, fix16_t lambda);

:kf:        The filter to update.
:lambda:    The estimation certainty tuning factor. 0.0 < lambda <= 1.0;

This performs a state and covariance update according to the state transition model *A* and the input model *B*. If *B* has zero dimensions, only the state transition model will be used.
In addition, the system covariance matrix will be scaled by the factor 1/lambda^2. This can be used to artificially increase prediction uncertainty to prevent convergence.

If input values are used, the user is required to set the values in ``kfm.u`` prior to calling this function.

Similar to *kalman_predict()*, this function is a thin wrapper around `kalman_predict_x`_ and `kalman_predict_P_tuned`_.
It is often more efficient to perform the state update manually instead of relying on the matrix multiplication algorithm. In this case, `kalman_predict_P_tuned`_ can be used to update the system covariance
matrix afterwards.

kalman_predict_x
----------------
Kalman filter state-only prediction (time update) step::
    
    void kalman_predict_x(kalman16_t *kf);

:kf:        The filter to update.

This performs a state-only (i.e. no covariance) update according to the state transition model *A* and the input model *B*. If *B* has zero dimensions, only the state transition model will be used.

If input values are used, the user is required to set the values in ``kfm.u`` prior to calling this function.

kalman_predict_P
----------------
Kalman filter covariance-only prediction (time update) step::
    
    void kalman_predict_P(kalman16_t *kf);

:kf:        The filter to update.

This performs a covariance-only (i.e. no state) update according to the state transition model *A* and the input model *B*. If *B* has zero dimensions, only the state transition model will be used.

In cases where it is more efficient to calculate the state update manually (i.e. by not calling `kalman_predict`_), *kalman_predict_P* can be used to update the covariance matrix.

kalman_predict_P_tuned
----------------------
Kalman filter covariance-only prediction (time update) step with certainty tuning::
    
    void kalman_predict_P_tuned(kalman16_t *kf, fix16_t lambda);

:kf:        The filter to update.
:lambda:    The estimation certainty tuning factor. 0.0 < lambda <= 1.0

Similar to ``kalman_predict_P()``, this function performs a covariance-only (i.e. no state) update according to the state transition model *A* and the input model *B*. If *B* has zero dimensions, only the state transition model will be used.
In addition, the system covariance matrix will be scaled by the factor 1/lambda^2. This can be used to artificially increase prediction uncertainty to prevent convergence.

In cases where it is more efficient to calculate the state update manually (i.e. by not calling `kalman_predict_tuned`_), *kalman_predict_P_tuned* can be used to update the covariance matrix.

kalman_correct
--------------
Kalman filter correction (measurement update) step::

    void kalman_correct(kalman16_t *kf, kalman16_observation_t *kfm);

:kf:        The filter to update.
:kfm:       The observation used to update the filter.

This updates the state estimation as retrieved from the prediction functions and corrects the estimate using the observation in *kfm*.

The user is required to set the values in ``kfm.z`` (and ``kfm.R`` if required) prior to calling this function.

Prediction and Update Functions (systems without control input)
===============================================================

kalman_predict_uc
------------------
Kalman filter prediction (time update) step::
    
    void kalman_predict_uc(kalman16_uc_t *kf);

:kf:        The filter to update.

This performs a state and covariance update according to the state transition model *A*.

This function is a thin wrapper around `kalman_predict_x_uc`_ and `kalman_predict_P_uc`_.
It is often more efficient to perform the state update manually instead of relying on the matrix multiplication algorithm. In this case, `kalman_predict_P_uc`_ can be used to update the system covariance
matrix afterwards.

If input values are used, the user is required to set the values in ``kfm.u`` prior to calling this function.

kalman_cpredict_uc
------------------
Kalman filter continuous-time prediction (time update) and integration step::
    
    void kalman_predict_uc(kalman16_uc_t *kf, register fix16_t deltaT);

:kf:        The filter to update.
:deltaT:    The time differential in seconds.

This performs a state and covariance update according to the state transition model *A*.

This function is a thin wrapper around `kalman_cpredict_x_uc`_ and `kalman_cpredict_P_uc`_.
It is often more efficient to perform the state update manually instead of relying on the matrix multiplication algorithm. In this case, `kalman_cpredict_P_uc`_ can be used to update the system covariance
matrix afterwards.

If input values are used, the user is required to set the values in ``kfm.u`` prior to calling this function.

kalman_predict_tuned_uc
-----------------------
Kalman filter prediction (time update) step with applied certainty tuning::
    
    void kalman_predict_tuned_uc(kalman16_uc_t *kf, fix16_t lambda);

:kf:        The filter to update.
:lambda:    The estimation certainty tuning factor. 0.0 < lambda <= 1.0;

This performs a state and covariance update according to the state transition model *A*.
In addition, the system covariance matrix will be scaled by the factor 1/lambda^2. This can be used to artificially increase prediction uncertainty to prevent convergence.

If input values are used, the user is required to set the values in ``kfm.u`` prior to calling this function.

Similar to *kalman_predict_uc()*, this function is a thin wrapper around `kalman_predict_x_uc`_ and `kalman_predict_P_tuned_uc`_.
It is often more efficient to perform the state update manually instead of relying on the matrix multiplication algorithm. In this case, `kalman_predict_P_tuned_uc`_ can be used to update the system covariance
matrix afterwards.

kalman_predict_x_uc
-------------------
Kalman filter state-only prediction (time update) step::
    
    void kalman_predict_x_uc(kalman16_uc_t *kf);

:kf:        The filter to update.

This performs a state-only (i.e. no covariance) update according to the state transition model *A*.

If input values are used, the user is required to set the values in ``kfm.u`` prior to calling this function.

kalman_cpredict_x_uc
-------------------
Kalman filter continuous-time state-only prediction (time update) and integration step::
    
    void kalman_predict_x_uc(kalman16_uc_t *kf, register fix16_t deltaT);

:kf:        The filter to update.
:deltaT:	The time differential in seconds,

This performs a state-only (i.e. no covariance) update according to the state transition model *A*.

If input values are used, the user is required to set the values in ``kfm.u`` prior to calling this function.

kalman_predict_P_uc
-------------------
Kalman filter covariance-only prediction (time update) step::
    
    void kalman_predict_P_uc(kalman16_uc_t *kf);

:kf:        The filter to update.

This performs a covariance-only (i.e. no state) update according to the state transition model *A*.

In cases where it is more efficient to calculate the state update manually (i.e. by not calling `kalman_predict`_), *kalman_predict_P* can be used to update the covariance matrix.

kalman_cpredict_P_uc
-------------------
Kalman filter cintinuous-time covariance-only prediction (time update) and integration step::
    
    void kalman_predict_P_uc(kalman16_uc_t *kf, register fix16_t deltaT);

:kf:        The filter to update.
:deltaT:	The time differential.

This performs a covariance-only (i.e. no state) update according to the state transition model *A*.

In cases where it is more efficient to calculate the state update manually (i.e. by not calling `kalman_cpredict`_), *kalman_cpredict_P* can be used to update the covariance matrix.

kalman_predict_P_tuned_uc
--------------------------
Kalman filter covariance-only prediction (time update) step with certainty tuning::
    
    void kalman_predict_P_tuned_uc(kalman16_uc_t *kf, fix16_t lambda);

:kf:        The filter to update.
:lambda:    The estimation certainty tuning factor. 0.0 < lambda <= 1.0

Similar to ``kalman_predict_P()``, this function performs a covariance-only (i.e. no state) update according to the state transition model *A*.
In addition, the system covariance matrix will be scaled by the factor 1/lambda^2. This can be used to artificially increase prediction uncertainty to prevent convergence.

In cases where it is more efficient to calculate the state update manually (i.e. by not calling `kalman_predict_tuned_uc`_), *kalman_predict_P_tuned_uc* can be used to update the covariance matrix.

kalman_correct_uc
-----------------
Kalman filter correction (measurement update) step::

    void kalman_correct_uc(kalman16_uc_t *kf, kalman16_observation_t *kfm);

:kf:        The filter to update.
:kfm:       The observation used to update the filter.

This updates the state estimation as retrieved from the prediction functions and corrects the estimate using the observation in *kfm*.

The user is required to set the values in ``kfm.z`` (and ``kfm.R`` if required) prior to calling this function.


Helper Functions (regular systems)
==================================

kalman_get_state_vector
-----------------------
Retrieves a pointer to the state vector *x*::

    mf16* kalman_get_state_vector(kalman16_t *kf);

:kf:        The filter.

kalman_get_state_transition
---------------------------
Retrieves a pointer to the state transition model *A*::

    mf16* kalman_get_state_transition(kalman16_t *kf);

:kf:        The filter.

kalman_get_system_covariance
----------------------------
Retrieves a pointer to the system covariance matrix *P*::

    mf16* kalman_get_system_covariance(kalman16_t *kf);

:kf:        The filter.

kalman_get_input_vector
-----------------------
Retrieves a pointer to the control input vector *u*::

    mf16* kalman_get_input_vector(kalman16_t *kf);

:kf:        The filter.

kalman_get_input_transition
---------------------------
Retrieves a pointer to the control input transition model *B*::

    mf16* kalman_get_input_transition(kalman16_t *kf)

:kf:        The filter.

kalman_get_input_covariance
---------------------------
Retrieves a pointer to the control input covariance matrix *Q*::

    mf16* kalman_get_input_covariance(kalman16_t *kf)

:kf:        The filter.

Helper Functions (systems without control input)
================================================

kalman_get_state_vector_uc
--------------------------
Retrieves a pointer to the state vector *x*::

    mf16* kalman_get_state_vector_uc(kalman16_uc_t *kf);

:kf:        The filter.

kalman_get_state_transition_uc
------------------------------
Retrieves a pointer to the state transition model *A*::

    mf16* kalman_get_state_transition_uc(kalman16_uc_t *kf);

:kf:        The filter.

kalman_get_system_covariance_uc
--------------------------------
Retrieves a pointer to the system covariance matrix *P*::

    mf16* kalman_get_system_covariance_uc(kalman16_uc_t *kf);

:kf:        The filter.

kalman_get_system_process_noise_uc
----------------------------------
Retrieves a pointer to the system process noise matrix *Q*::

    mf16* kalman_get_system_process_noise_uc(kalman16_t *kf)

:kf:        The filter.

Helper Functions (observations)
===============================

kalman_get_observation_vector
-----------------------------
Retrieves a pointer to the observation vector *z*::

    mf16* kalman_get_observation_transformation(kalman16_observation_t *kfm)

:kfm:        The measurement.

kalman_get_observation_process_noise
------------------------------------
Retrieves a pointer to the process noise matrix *R*::

    mf16* kalman_get_observation_process_noise(kalman16_observation_t *kfm)

:kfm:        The measurement.