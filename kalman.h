#ifndef _KALMAN_H_
#define _KALMAN_H_

#include <stdint.h>
#include "compiler.h"
#include "fixmatrix.h"

/*!
* \def EXTERN_INLINE_KALMAN Helper inline to switch from local inline to extern inline
*/
#ifndef EXTERN_INLINE_KALMAN
#define EXTERN_INLINE_KALMAN EXTERN_INLINE
#endif

/*!
* \brief Kalman Filter structure
* \see kalman_measurement_t
*/
typedef struct
{
    /*!
    * \brief State vector The state transition matrix (#states x #states)
    */
    mf16 x;

    /*!
    * \brief System matrix (#states x 1)
    * \see P
    */
    mf16 A;

    /*!
    * \brief System covariance matrix (#states x #states)
    * \see A
    */
    mf16 P;
    
    /*!
    * \brief Input vector (#inputs x 1)
    */
    mf16 u;

    /*!
    * \brief Input matrix (#inputs x #inputs)
    * \see Q
    */
    mf16 B;

    /*!
    * \brief Input covariance/uncertainty matrix (#inputs x #inputs)
    * \see B
    */
    mf16 Q;

} kalman_t;

/*!
* \brief Kalman Filter measurement structure
* \see kalman_t
*/
typedef struct
{
    /*!
    * \brief Measurement vector (#measurements x 1)
    */
    mf16 z;

    /*!
    * \brief Measurement transformation matrix (#measurements x #measurements)
    * \see R
    */
    mf16 H;

    /*!
    * \brief Process noise covariance matrix (#measurements x #measurements)
    * \see H
    */
    mf16 R;

    /*!
    * \brief Innovation vector (#measurements x 1)
    */
    mf16 y;

    /*!
    * \brief Residual covariance matrix (#measurements x #measurements)
    * \see y
    */
    mf16 S;

    /*!
    * \brief Kalman gain matrix (#states x #measurements)
    */
    mf16 K;

} kalman_measurement_t;

/*!
* \brief Initializes the Kalman Filter
* \param[in] kf The Kalman Filter structure to initialize
* \param[in] num_states The number of state variables
* \param[in] num_inputs The number of input variables
*/
COLD LEAF NONNULL
void kalman_filter_initialize(kalman_t *kf, uint_fast8_t num_states, uint_fast8_t num_inputs);

/*!
* \brief Sets the measurement vector
* \param[in] kfm The Kalman Filter measurement structure to initialize
* \param[in] num_states The number of states
* \param[in] num_measurements The number of measurements
*/
COLD LEAF NONNULL
void kalman_measurement_initialize(kalman_measurement_t *kfm, uint_fast8_t num_states, uint_fast8_t num_measurements);

/*!
* \brief Performs the time update / prediction step of only the state vector
* \param[in] kf The Kalman Filter structure to predict with.
*
* \see kalman_predict
* \see kalman_predict_tuned
*/
HOT LEAF NONNULL
void kalman_predict_x(register kalman_t *const kf);

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*
* \see kalman_predict
* \see kalman_predict_Q_tuned
*/
HOT LEAF NONNULL
void kalman_predict_Q(register kalman_t *const kf);

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*
* \see kalman_predict_tuned
* \see kalman_predict_Q
*/
HOT LEAF NONNULL
void kalman_predict_Q_tuned(register kalman_t *const kf, fix16_t lambda);

/*!
* \brief Performs the time update / prediction step.
* \param[in] kf The Kalman Filter structure to predict with.
* \param[in] lambda Lambda factor (\c 0 < {\ref lambda} <= \c 1) to forcibly reduce prediction certainty. Smaller values mean larger uncertainty.
*
* This call assumes that the input covariance and variables are already set in the filter structure.
*
* \see kalman_predict_x
* \see kalman_predict_Q
*/
LEAF NONNULL
EXTERN_INLINE_KALMAN void kalman_predict(kalman_t *kf)
{
    /************************************************************************/
    /* Predict next state using system dynamics                             */
    /* x = A*x                                                              */
    /************************************************************************/

    kalman_predict_x(kf);

    /************************************************************************/
    /* Predict next covariance using system dynamics and input              */
    /* P = A*P*A' + B*Q*B'                                                  */
    /************************************************************************/

    kalman_predict_Q(kf);
}

/*!
* \brief Performs the time update / prediction step.
* \param[in] kf The Kalman Filter structure to predict with.
* \param[in] lambda Lambda factor (\c 0 < {\ref lambda} <= \c 1) to forcibly reduce prediction certainty. Smaller values mean larger uncertainty.
*
* This call assumes that the input covariance and variables are already set in the filter structure.
*
* \see kalman_predict_x
* \see kalman_predict_Q_tuned
*/
LEAF HOT NONNULL
EXTERN_INLINE_KALMAN void kalman_predict_tuned(kalman_t *kf, fix16_t lambda)
{
    /************************************************************************/
    /* Predict next state using system dynamics                             */
    /* x = A*x                                                              */
    /************************************************************************/

    kalman_predict_x(kf);

    /************************************************************************/
    /* Predict next covariance using system dynamics and input              */
    /* P = A*P*A' * 1/lambda^2 + B*Q*B'                                     */
    /************************************************************************/

    kalman_predict_Q_tuned(kf, lambda);
}

/*!
* \brief Performs the measurement update step.
* \param[in] kf The Kalman Filter structure to correct.
*/
HOT LEAF NONNULL
void kalman_correct(kalman_t *kf, kalman_measurement_t *kfm);

/*!
* \brief Gets a pointer to the state vector x.
* \param[in] kf The Kalman Filter structure
* \return The state vector x.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_state_vector(kalman_t *kf)
{
    return &(kf->x);
}

/*!
* \brief Gets a pointer to the state transition matrix A.
* \param[in] kf The Kalman Filter structure
* \return The state transition matrix A.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_state_transition(kalman_t *kf)
{
    return &(kf->A);
}

/*!
* \brief Gets a pointer to the system covariance matrix P.
* \param[in] kf The Kalman Filter structure
* \return The system covariance matrix.
*/
LEAF RETURNS_NONNULL NONNULL PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_system_covariance(kalman_t *kf)
{
    return &(kf->P);
}

/*!
* \brief Gets a pointer to the input vector u.
* \param[in] kf The Kalman Filter structure
* \return The input vector u.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_input_vector(kalman_t *kf)
{
    return &(kf->u);
}

/*!
* \brief Gets a pointer to the input transition matrix B.
* \param[in] kf The Kalman Filter structure
* \return The input transition matrix B.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_input_transition(kalman_t *kf)
{
    return &(kf->B);
}

/*!
* \brief Gets a pointer to the input covariance matrix P.
* \param[in] kf The Kalman Filter structure
* \return The input covariance matrix.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_input_covariance(kalman_t *kf)
{
    return &(kf->Q);
}

/*!
* \brief Gets a pointer to the measurement vector z.
* \param[in] kfm The Kalman Filter measurement structure.
* \return The measurement vector z.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_measurement_vector(kalman_measurement_t *kfm)
{
    return &(kfm->z);
}

/*!
* \brief Gets a pointer to the measurement transformation matrix H.
* \param[in] kfm The Kalman Filter measurement structure.
* \return The measurement transformation matrix H.
*/
RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_measurement_transformation(kalman_measurement_t *kfm)
{
    return &(kfm->H);
}

/*!
* \brief Gets a pointer to the process noise matrix R.
* \param[in] kfm The Kalman Filter measurement structure.
* \return The process noise matrix R.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_process_noise(kalman_measurement_t *kfm)
{
    return &(kfm->R);
}

#undef EXTERN_INLINE_KALMAN
#endif