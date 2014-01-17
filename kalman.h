#ifndef _KALMAN_H_
#define _KALMAN_H_

#include <stdint.h>
#include "compiler.h"
#include "fixmatrix.h"

/*!
* \def KALMAN_DISABLE_UC Global define to disable functions for systems without control inputs
*/
#ifdef KALMAN_DISABLE_UC
#pragma message("KALMAN_DISABLE_UC defined. Disabling Kalman filter functions for systems without control inputs.")
#endif

/*!
* \def KALMAN_DISABLE_C Global define to disable functions for systems with control inputs
*/
#ifdef KALMAN_DISABLE_C
#pragma message("KALMAN_DISABLE_C defined. Disabling Kalman filter functions for systems with control inputs.")
#endif

/*!
* \def KALMAN_DISABLE_LAMBDA Global define to disable certainty tuning. 
*/
#ifdef KALMAN_DISABLE_LAMBDA
#pragma message("KALMAN_DISABLE_LAMBDA defined. Disabling Kalman filter functions with certainty tuning.")
#endif

#if defined(KALMAN_DISABLE_C) && defined(KALMAN_DISABLE_UC)
#error KALMAN_DISABLE_C and KALMAN_DISABLE_UC defined; This removes all Kalman filter functions. Remove one of the defines and rebuild.
#endif

#ifndef KALMAN_DISABLE_UC

/*!
* \brief Kalman Filter structure for a system without control inputs
* \see kalman16_observation_t
* \see kalman16_t
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
    * \brief System process noise matrix (#states x #states)
    * \see P
    */
    mf16 Q;

} kalman16_uc_t;

#endif // KALMAN_DISABLE_UC

#if !defined(KALMAN_DISABLE_C) || !defined(KALMAN_DISABLE_UC) // some UC methods require this struct

/*!
* \brief Kalman Filter structure
* \see kalman16_observation_t
* \see kalman16_uc_t
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

} kalman16_t;

#endif // KALMAN_DISABLE_C

/*!
* \brief Kalman Filter measurement structure
* \see kalman16_t
*/
typedef struct
{
    /*!
    * \brief Measurement vector (#measurements x 1)
    */
    mf16 z;

    /*!
    * \brief Measurement transformation matrix (#measurements x #states)
    * \see R
    */
    mf16 H;

    /*!
    * \brief Observation process noise covariance matrix (#measurements x #measurements)
    * \see H
    */
    mf16 R;

} kalman16_observation_t;

/************************************************************************/
/* Helper defines                                                       */
/************************************************************************/

/*!
* \def EXTERN_INLINE_KALMAN Helper inline to switch from local inline to extern inline
*/
#ifndef EXTERN_INLINE_KALMAN
#define EXTERN_INLINE_KALMAN EXTERN_INLINE
#endif

/************************************************************************/
/* Functions for systems with control inputs                            */
/************************************************************************/

#ifndef KALMAN_DISABLE_C

/*!
* \brief Initializes the Kalman Filter
* \param[in] kf The Kalman Filter structure to initialize
* \param[in] num_states The number of state variables
* \param[in] num_inputs The number of input variables
*/
COLD LEAF NONNULL
void kalman_filter_initialize(kalman16_t *const kf, uint_fast8_t num_states, uint_fast8_t num_inputs);

/*!
* \brief Performs the time update / prediction step of only the state vector
* \param[in] kf The Kalman Filter structure to predict with.
*
* \see kalman_predict
* \see kalman_predict_tuned
*/
HOT LEAF NONNULL
void kalman_predict_x(register kalman16_t *const kf);

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*
* \see kalman_predict
* \see kalman_predict_P_tuned
*/
HOT LEAF NONNULL
void kalman_predict_P(register kalman16_t *const kf);

#ifndef KALMAN_DISABLE_LAMBDA

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*
* \see kalman_predict_tuned
* \see kalman_predict_P
*/
HOT LEAF NONNULL
void kalman_predict_P_tuned(register kalman16_t *const kf, fix16_t lambda);

#endif // KALMAN_DISABLE_LAMBDA

/*!
* \brief Performs the time update / prediction step.
* \param[in] kf The Kalman Filter structure to predict with.
* \param[in] lambda Lambda factor (\c 0 < {\ref lambda} <= \c 1) to forcibly reduce prediction certainty. Smaller values mean larger uncertainty.
*
* This call assumes that the input covariance and variables are already set in the filter structure.
*
* \see kalman_predict_x
* \see kalman_predict_P
*/
LEAF NONNULL
EXTERN_INLINE_KALMAN void kalman_predict(kalman16_t *kf)
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

    kalman_predict_P(kf);
}

#ifndef KALMAN_DISABLE_LAMBDA

/*!
* \brief Performs the time update / prediction step.
* \param[in] kf The Kalman Filter structure to predict with.
* \param[in] lambda Lambda factor (\c 0 < {\ref lambda} <= \c 1) to forcibly reduce prediction certainty. Smaller values mean larger uncertainty.
*
* This call assumes that the input covariance and variables are already set in the filter structure.
*
* \see kalman_predict_x
* \see kalman_predict_P_tuned
*/
LEAF HOT NONNULL
EXTERN_INLINE_KALMAN void kalman_predict_tuned(kalman16_t *kf, fix16_t lambda)
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

    kalman_predict_P_tuned(kf, lambda);
}

#endif // KALMAN_DISABLE_LAMBDA

/*!
* \brief Performs the measurement update step.
* \param[in] kf The Kalman Filter structure to correct.
*/
HOT LEAF NONNULL
void kalman_correct(kalman16_t *kf, kalman16_observation_t *kfm);

/*!
* \brief Gets a pointer to the state vector x.
* \param[in] kf The Kalman Filter structure
* \return The state vector x.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_state_vector(kalman16_t *kf)
{
    return &(kf->x);
}

/*!
* \brief Gets a pointer to the state transition matrix A.
* \param[in] kf The Kalman Filter structure
* \return The state transition matrix A.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_state_transition(kalman16_t *kf)
{
    return &(kf->A);
}

/*!
* \brief Gets a pointer to the system covariance matrix P.
* \param[in] kf The Kalman Filter structure
* \return The system covariance matrix.
*/
LEAF RETURNS_NONNULL NONNULL PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_system_covariance(kalman16_t *kf)
{
    return &(kf->P);
}

/*!
* \brief Gets a pointer to the input vector u.
* \param[in] kf The Kalman Filter structure
* \return The input vector u.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_input_vector(kalman16_t *kf)
{
    return &(kf->u);
}

/*!
* \brief Gets a pointer to the input transition matrix B.
* \param[in] kf The Kalman Filter structure
* \return The input transition matrix B.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_input_transition(kalman16_t *kf)
{
    return &(kf->B);
}

/*!
* \brief Gets a pointer to the input covariance matrix P.
* \param[in] kf The Kalman Filter structure
* \return The input covariance matrix.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_input_covariance(kalman16_t *kf)
{
    return &(kf->Q);
}

#endif // KALMAN_DISABLE_C

/************************************************************************/
/* Functions for observation handling                                   */
/************************************************************************/

/*!
* \brief Sets the measurement vector
* \param[in] kfm The Kalman Filter measurement structure to initialize
* \param[in] num_states The number of states
* \param[in] num_observations The number of observations
*/
COLD LEAF NONNULL
void kalman_observation_initialize(kalman16_observation_t *const kfm, uint_fast8_t num_states, uint_fast8_t num_observations);


/*!
* \brief Gets a pointer to the measurement vector z.
* \param[in] kfm The Kalman Filter measurement structure.
* \return The measurement vector z.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_observation_vector(kalman16_observation_t *kfm)
{
    return &(kfm->z);
}

/*!
* \brief Gets a pointer to the measurement transformation matrix H.
* \param[in] kfm The Kalman Filter measurement structure.
* \return The measurement transformation matrix H.
*/
RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_observation_transformation(kalman16_observation_t *kfm)
{
    return &(kfm->H);
}

/*!
* \brief Gets a pointer to the process noise matrix R.
* \param[in] kfm The Kalman Filter measurement structure.
* \return The process noise matrix R.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE 
EXTERN_INLINE_KALMAN mf16* kalman_get_observation_process_noise(kalman16_observation_t *kfm)
{
    return &(kfm->R);
}

/************************************************************************/
/* Functions for systems without control inputs                         */
/************************************************************************/

#ifndef KALMAN_DISABLE_UC

/*!
* \brief Initializes the Kalman Filter
* \param[in] kf The Kalman Filter structure to initialize
* \param[in] num_states The number of state variables
*/
COLD LEAF NONNULL
void kalman_filter_initialize_uc(kalman16_uc_t *const kf, uint_fast8_t num_states);

/*!
* \brief Performs the time update / prediction step of only the state vector
* \param[in] kf The Kalman Filter structure to predict with.
*
* \see kalman_predict_uc
* \see kalman_predict_tuned_uc
*/
HOT LEAF NONNULL
void kalman_predict_x_uc(register kalman16_uc_t *const kf);

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*
* \see kalman_predict_uc
* \see kalman_predict_P_tuned_uc
*/
HOT LEAF NONNULL
void kalman_predict_P_uc(register kalman16_uc_t *const kf);

#ifndef KALMAN_DISABLE_LAMBDA

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*
* \see kalman_predict_tuned_uc
* \see kalman_predict_P_uc
*/
HOT LEAF NONNULL
void kalman_predict_P_tuned_uc(register kalman16_uc_t *const kf, fix16_t lambda);

#endif // KALMAN_DISABLE_LAMBDA

/*!
* \brief Performs the time update / prediction step.
* \param[in] kf The Kalman Filter structure to predict with.
* \param[in] lambda Lambda factor (\c 0 < {\ref lambda} <= \c 1) to forcibly reduce prediction certainty. Smaller values mean larger uncertainty.
*
* This call assumes that the input covariance and variables are already set in the filter structure.
*
* \see kalman_predict_x
* \see kalman_predict_P
*/
LEAF NONNULL
EXTERN_INLINE_KALMAN void kalman_predict_uc(kalman16_uc_t *kf)
{
    /************************************************************************/
    /* Predict next state using system dynamics                             */
    /* x = A*x                                                              */
    /************************************************************************/

    kalman_predict_x_uc(kf);

    /************************************************************************/
    /* Predict next covariance using system dynamics and input              */
    /* P = A*P*A' + B*Q*B'                                                  */
    /************************************************************************/

    kalman_predict_P_uc(kf);
}

#ifndef KALMAN_DISABLE_LAMBDA

/*!
* \brief Performs the time update / prediction step.
* \param[in] kf The Kalman Filter structure to predict with.
* \param[in] lambda Lambda factor (\c 0 < {\ref lambda} <= \c 1) to forcibly reduce prediction certainty. Smaller values mean larger uncertainty.
*
* This call assumes that the input covariance and variables are already set in the filter structure.
*
* \see kalman_predict_x
* \see kalman_predict_P_tuned
*/
LEAF HOT NONNULL
EXTERN_INLINE_KALMAN void kalman_predict_tuned_uc(kalman16_uc_t *kf, fix16_t lambda)
{
    /************************************************************************/
    /* Predict next state using system dynamics                             */
    /* x = A*x                                                              */
    /************************************************************************/

    kalman_predict_x_uc(kf);

    /************************************************************************/
    /* Predict next covariance using system dynamics and input              */
    /* P = A*P*A' * 1/lambda^2 + B*Q*B'                                     */
    /************************************************************************/

    kalman_predict_P_tuned_uc(kf, lambda);
}

#endif // KALMAN_DISABLE_LAMBDA

/*!
* \brief Performs the measurement update step.
* \param[in] kf The Kalman Filter structure to correct.
*/
HOT LEAF NONNULL
void kalman_correct_uc(kalman16_uc_t *kf, kalman16_observation_t *kfm);

/*!
* \brief Gets a pointer to the state vector x.
* \param[in] kf The Kalman Filter structure
* \return The state vector x.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE
EXTERN_INLINE_KALMAN mf16* kalman_get_state_vector_uc(kalman16_uc_t *kf)
{
    return &(kf->x);
}

/*!
* \brief Gets a pointer to the state transition matrix A.
* \param[in] kf The Kalman Filter structure
* \return The state transition matrix A.
*/
LEAF RETURNS_NONNULL NONNULL HOT PURE
EXTERN_INLINE_KALMAN mf16* kalman_get_state_transition_uc(kalman16_uc_t *kf)
{
    return &(kf->A);
}

/*!
* \brief Gets a pointer to the system covariance matrix P.
* \param[in] kf The Kalman Filter structure
* \return The system covariance matrix.
*/
LEAF RETURNS_NONNULL NONNULL PURE
EXTERN_INLINE_KALMAN mf16* kalman_get_system_covariance_uc(kalman16_uc_t *kf)
{
    return &(kf->P);
}

/*!
* \brief Gets a pointer to the system process noise vector.
* \param[in] kf The Kalman Filter structure
* \return The process noise vector.
*/
LEAF RETURNS_NONNULL NONNULL PURE
EXTERN_INLINE_KALMAN mf16* kalman_get_system_process_noise_uc(kalman16_uc_t *kf)
{
    return &(kf->Q);
}

#endif // KALMAN_DISABLE_UC

#undef EXTERN_INLINE_KALMAN
#endif
