/*!
* \brief Example of the Kalman filter
*
* In this example, the gravity constant (~9.81 m/s^2) will be estimated using only
* measurements of the position. These measurements have a variance of var(s) = 0.5m.
*
* The formulas used are:
* s = s + v*T + g*0.5*T^2
* v = v + g*T
* g = g
*
* The time constant is set to T = 1s.
*
* The initial estimation of the gravity constant is set to 6 m/s^2.
*/

#include <assert.h>

#include "kalman.h"

/**
* \brief Enables usage of the functions specialized for systems without control inputs
*
* Define to 0 in order to use the regular functions
*/
#define USE_UNCONTROLLED 1

// create the filter structure
#define KALMAN_NAME gravity
#define KALMAN_NUM_STATES 3
#define KALMAN_NUM_INPUTS 0

#if USE_UNCONTROLLED
kalman16_uc_t kf;
#else
kalman16_t kf;
#endif


// create the measurement structure
#define KALMAN_MEASUREMENT_NAME position
#define KALMAN_NUM_MEASUREMENTS 1
kalman16_observation_t kfm;

#define matrix_set(matrix, row, column, value) \
    matrix->data[row][column] = value

#define matrix_set_symmetric(matrix, row, column, value) \
    matrix->data[row][column] = value; \
    matrix->data[column][row] = value

#ifndef FIXMATRIX_MAX_SIZE
#error FIXMATRIX_MAX_SIZE must be defined and greater or equal to the number of states, inputs and measurements.
#endif

#if (FIXMATRIX_MAX_SIZE < KALMAN_NUM_STATES) || (FIXMATRIX_MAX_SIZE < KALMAN_NUM_INPUTS) || (FIXMATRIX_MAX_SIZE < KALMAN_NUM_MEASUREMENTS)
#error FIXMATRIX_MAX_SIZE must be greater or equal to the number of states, inputs and measurements.
#endif

/*!
* \brief Initializes the gravity Kalman filter
*/
static void kalman_gravity_init()
{
    /************************************************************************/
    /* initialize the filter structures                                     */
    /************************************************************************/
#if USE_UNCONTROLLED
    kalman_filter_initialize_uc(&kf, KALMAN_NUM_STATES);
#else
    kalman_filter_initialize(&kf, KALMAN_NUM_STATES, KALMAN_NUM_INPUTS);
#endif
    kalman_observation_initialize(&kfm, KALMAN_NUM_STATES, KALMAN_NUM_MEASUREMENTS);

    /************************************************************************/
    /* set initial state                                                    */
    /************************************************************************/
#if USE_UNCONTROLLED
    mf16 *x = kalman_get_state_vector_uc(&kf);
#else
    mf16 *x = kalman_get_state_vector(&kf);
#endif
    x->data[0][0] = 0; // s_i
    x->data[1][0] = 0; // v_i
    x->data[2][0] = fix16_from_float(6); // g_i

    /************************************************************************/
    /* set state transition                                                 */
    /************************************************************************/
#if USE_UNCONTROLLED
    mf16 *A = kalman_get_state_transition_uc(&kf);
#else
    mf16 *A = kalman_get_state_transition(&kf);
#endif
    
    // set time constant
    const fix16_t T = fix16_one;
    const fix16_t Tsquare = fix16_sq(T);

    // helper
    const fix16_t fix16_half = fix16_from_float(0.5);

    // transition of x to s
    matrix_set(A, 0, 0, fix16_one);   // 1
    matrix_set(A, 0, 1, T);   // T
    matrix_set(A, 0, 2, fix16_mul(fix16_half, Tsquare)); // 0.5 * T^2
    
    // transition of x to v
    matrix_set(A, 1, 0, 0);   // 0
    matrix_set(A, 1, 1, fix16_one);   // 1
    matrix_set(A, 1, 2, T);   // T

    // transition of x to g
    matrix_set(A, 2, 0, 0);   // 0
    matrix_set(A, 2, 1, 0);   // 0
    matrix_set(A, 2, 2, fix16_one);   // 1

    /************************************************************************/
    /* set covariance                                                       */
    /************************************************************************/
#if USE_UNCONTROLLED
    mf16 *P = kalman_get_system_covariance_uc(&kf);
#else
    mf16 *P = kalman_get_system_covariance(&kf);
#endif

    matrix_set_symmetric(P, 0, 0, fix16_half);   // var(s)
    matrix_set_symmetric(P, 0, 1, 0);   // cov(s,v)
    matrix_set_symmetric(P, 0, 2, 0);   // cov(s,g)

    matrix_set_symmetric(P, 1, 1, fix16_one);   // var(v)
    matrix_set_symmetric(P, 1, 2, 0);   // cov(v,g)

    matrix_set_symmetric(P, 2, 2, fix16_one);   // var(g)

    /************************************************************************/
    /* set system process noise                                             */
    /************************************************************************/
#if USE_UNCONTROLLED
    mf16 *Q = kalman_get_system_process_noise_uc(&kf);
    mf16_fill(Q, F16(0.0001));
#endif

    /************************************************************************/
    /* set measurement transformation                                       */
    /************************************************************************/
    mf16 *H = kalman_get_observation_transformation(&kfm);

    matrix_set(H, 0, 0, fix16_one);     // z = 1*s 
    matrix_set(H, 0, 1, 0);     //   + 0*v
    matrix_set(H, 0, 2, 0);     //   + 0*g

    /************************************************************************/
    /* set process noise                                                    */
    /************************************************************************/
    mf16 *R = kalman_get_observation_process_noise(&kfm);

    matrix_set(R, 0, 0, fix16_half);     // var(s)
}

// define measurements.
//
// MATLAB source
// -------------
// s = s + v*T + g*0.5*T^2; 
// v = v + g*T;
#define MEAS_COUNT (15)
static fix16_t real_distance[MEAS_COUNT] = {
    F16(4.905),
    F16(19.62),
    F16(44.145),
    F16(78.48),
    F16(122.63),
    F16(176.58),
    F16(240.35),
    F16(313.92),
    F16(397.31),
    F16(490.5),
    F16(593.51),
    F16(706.32),
    F16(828.94),
    F16(961.38),
    F16(1103.6) };

// define measurement noise with variance 0.5
//
// MATLAB source
// -------------
// noise = 0.5^2*randn(15,1);
static fix16_t measurement_error[MEAS_COUNT] = {
    F16(0.13442),
    F16(0.45847),
    F16(-0.56471),
    F16(0.21554),
    F16(0.079691),
    F16(-0.32692),
    F16(-0.1084),
    F16(0.085656),
    F16(0.8946),
    F16(0.69236),
    F16(-0.33747),
    F16(0.75873),
    F16(0.18135),
    F16(-0.015764),
    F16(0.17869) };

/*!
* \brief Runs the gravity Kalman filter.
*/
#if USE_UNCONTROLLED
void kalman_gravity_demo()
{
    // initialize the filter
    kalman_gravity_init();

    mf16 *x = kalman_get_state_vector_uc(&kf);
    mf16 *z = kalman_get_observation_vector(&kfm);
    
    // filter!
    uint_fast16_t i;
    for (i = 0; i < MEAS_COUNT; ++i)
    {
        // prediction.
        kalman_predict_uc(&kf);

        // measure ...
        fix16_t measurement = fix16_add(real_distance[i], measurement_error[i]);
        matrix_set(z, 0, 0, measurement);

        // update
        kalman_correct_uc(&kf, &kfm);
    }

    // fetch estimated g
    const fix16_t g_estimated = x->data[2][0];
    const float value = fix16_to_float(g_estimated);
    assert(value > 9.7 && value < 10);
}
#else
void kalman_gravity_demo()
{
    // initialize the filter
    kalman_gravity_init();

    mf16 *x = kalman_get_state_vector(&kf);
    mf16 *z = kalman_get_observation_vector(&kfm);

    // filter!
    uint_fast16_t i;
    for (i = 0; i < MEAS_COUNT; ++i)
    {
        // prediction.
        kalman_predict(&kf);

        // measure ...
        fix16_t measurement = fix16_add(real_distance[i], measurement_error[i]);
        matrix_set(z, 0, 0, measurement);

        // update
        kalman_correct(&kf, &kfm);
    }

    // fetch estimated g
    const fix16_t g_estimated = x->data[2][0];
    const float value = fix16_to_float(g_estimated);
    assert(value > 9.7 && value < 10);
}
#endif


/*!
* \brief Runs the gravity Kalman filter with lambda tuning.
*/
#if USE_UNCONTROLLED
void kalman_gravity_demo_lambda()
{
    // initialize the filter
    kalman_gravity_init();

    mf16 *x = kalman_get_state_vector_uc(&kf);
    mf16 *z = kalman_get_observation_vector(&kfm);

    // forcibly increase uncertainty in every prediction step by ~20% (1/lambda^2)
    const fix16_t lambda = F16(0.9);

    // filter!
    uint_fast16_t i;
    for (i = 0; i < MEAS_COUNT; ++i)
    {
        // prediction.
        kalman_predict_tuned_uc(&kf, lambda);

        // measure ...
        fix16_t measurement = fix16_add(real_distance[i], measurement_error[i]);
        matrix_set(z, 0, 0, measurement);

        // update
        kalman_correct_uc(&kf, &kfm);
    }

    // fetch estimated g
    const fix16_t g_estimated = x->data[2][0];
    const float value = fix16_to_float(g_estimated);
    assert(value > 9.7 && value < 10);
}
#else
void kalman_gravity_demo_lambda()
{
    // initialize the filter
    kalman_gravity_init();

    mf16 *x = kalman_get_state_vector(&kf);
    mf16 *z = kalman_get_observation_vector(&kfm);

    // forcibly increase uncertainty in every prediction step by ~20% (1/lambda^2)
    const fix16_t lambda = F16(0.9);

    // filter!
    uint_fast16_t i;
    for (i = 0; i < MEAS_COUNT; ++i)
    {
        // prediction.
        kalman_predict_tuned(&kf, lambda);

        // measure ...
        fix16_t measurement = fix16_add(real_distance[i], measurement_error[i]);
        matrix_set(z, 0, 0, measurement);

        // update
        kalman_correct(&kf, &kfm);
    }

    // fetch estimated g
    const fix16_t g_estimated = x->data[2][0];
    const float value = fix16_to_float(g_estimated);
    assert(value > 9.7 && value < 10);
}
#endif


void main()
{
    kalman_gravity_demo();
    kalman_gravity_demo_lambda();
}