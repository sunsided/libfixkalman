#include <stdint.h>
#include <assert.h>

#include "fixmatrixadditions.h"

#define EXTERN_INLINE_KALMAN INLINE
#include "kalman.h"

/*!
* \brief Initializes the Kalman Filter
* \param[in] kf The Kalman Filter structure to initialize
* \param[in] num_states The number of state variables
* \param[in] num_inputs The number of input variables
*/
void kalman_filter_initialize(kalman_t *kf, uint_fast8_t num_states, uint_fast8_t num_inputs)
{
    kf->A.rows = num_states;
    kf->A.columns = num_states;
    
    kf->P.rows = num_states;
    kf->P.columns = num_states;

    kf->B.rows = num_inputs;
    kf->B.columns = num_inputs;

    kf->Q.rows = num_inputs;
    kf->Q.columns = num_inputs;

    kf->x.rows = num_states;
    kf->x.columns = 1;

    kf->u.rows = num_inputs;
    kf->u.columns = 1;
}


/*!
* \brief Sets the measurement vector
* \param[in] kfm The Kalman Filter measurement structure to initialize
* \param[in] num_states The number of states
* \param[in] num_measurements The number of measurements
*/
void kalman_measurement_initialize(kalman_measurement_t *kfm, uint_fast8_t num_states, uint_fast8_t num_measurements)
{
    kfm->H.rows = num_measurements;
    kfm->H.columns = num_states;

    kfm->R.rows = num_measurements;
    kfm->R.columns = num_measurements;

    kfm->z.rows = num_measurements;
    kfm->z.columns = 1;

    kfm->K.rows = num_states;
    kfm->K.columns = num_measurements;

    kfm->S.rows = num_measurements;
    kfm->S.columns = num_measurements;

    kfm->y.rows = num_measurements;
    kfm->y.columns = 1;
}

/*!
* \brief Performs the time update / prediction step of only the state vector
* \param[in] kf The Kalman Filter structure to predict with.
*/
void kalman_predict_x(register kalman_t *const kf)
{
    // matrices and vectors
    const mf16 *RESTRICT const A = &kf->A;
    mf16 *RESTRICT const x = &kf->x;
    
    /************************************************************************/
    /* Predict next state using system dynamics                             */
    /* x = A*x                                                              */
    /************************************************************************/

    // x = A*x
    mf16_mul(x, A, x);
    
    // TODO: Implement more efficient matrix/row vector multiply
}

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*/
void kalman_predict_Q(register kalman_t *const kf)
{
    // matrices and vectors
    const mf16 *RESTRICT const A = &kf->A;
    const mf16 *RESTRICT const B = &kf->B;
    mf16 *RESTRICT const P = &kf->P;
    mf16 *RESTRICT const Q = &kf->Q;

    mf16 P_temp = { P->rows, P->columns };
    mf16 BQ_temp = { A->rows, B->columns };

    /************************************************************************/
    /* Predict next covariance using system dynamics and input              */
    /* P = A*P*A' + B*Q*B'                                                  */
    /************************************************************************/

    // P = A*P*A'
    mf16_mul_abat(P, A, P);                 // P = A*P*A'

    // TODO: this should work without P_temp
    // TODO: extract the code block below

    // P = P + B*Q*B'
    if (B->rows > 0)
    {
        mf16_mul_abat(&BQ_temp, B, Q);       // temp = B*Q*B'
        mf16_add(P, P, &BQ_temp);             // P += P_temp

        // TODO: Implement matrix-matrix multiply-and-add
        // TODO: Implement matrix-matrix add-in-place
    }
}

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*/
void kalman_predict_Q_tuned(register kalman_t *const kf, fix16_t lambda)
{
    // matrices and vectors
    const mf16 *RESTRICT const A = &kf->A;
    const mf16 *RESTRICT const B = &kf->B;
    mf16 *RESTRICT const P = &kf->P;
    
    mf16 P_temp = { P->rows, P->columns };
    mf16 BQ_temp = { A->rows, B->columns };

    /************************************************************************/
    /* Predict next covariance using system dynamics and input              */
    /* P = A*P*A' * 1/lambda^2 + B*Q*B'                                     */
    /************************************************************************/

    // lambda = 1/lambda^2
    lambda = fix16_div(fix16_one, fix16_sq(lambda)); // TODO: This should be precalculated, e.g. using kalman_set_lambda(...);

    // P = A*P*A'
    mf16_mul_abat(P, A, P);                         // temp = A*P*A'
    mf16_mul_s(P, P, lambda);                  // P *= 1/(lambda^2)

    // TODO: an a-b-ct multiplication might come in handy
    // TODO: extract the code block below

    // P = P + B*Q*B'
    if (B->rows > 0)
    {
        mf16_mul_abat(&BQ_temp, B, &kf->Q);         // temp = B*Q*B'
        mf16_add(P, &BQ_temp, P);               // P += temp

        // TODO: Implement matrix-matrix multiply-and-add
        // TODO: Implement matrix-matrix add-in-place
    }
}

/*!
* \brief Performs the measurement update step.
* \param[in] kf The Kalman Filter structure to correct.
*/
void kalman_correct(kalman_t *kf, kalman_measurement_t *kfm)
{
    mf16 *RESTRICT const P = &kf->P;
    const mf16 *RESTRICT const H = &kfm->H;
    mf16 *RESTRICT const K = &kfm->K;
    mf16 *RESTRICT const S = &kfm->S;
    mf16 *RESTRICT const y = &kfm->y;
    mf16 *RESTRICT const x = &kf->x;
    
    mf16 temp_HP = { H->rows, H->columns };
    mf16 temp_PHt = { P->rows, H->columns };

    /************************************************************************/
    /* Calculate innovation and residual covariance                         */
    /* y = z - H*x                                                          */
    /* S = H*P*H' + R                                                       */
    /************************************************************************/

    // y = z - H*x
    mf16_mul(y, H, x);
    mf16_sub(y, &kfm->z, y);

    // S = H*P*H' + R
    mf16_mul_abat(S, H, P);               // temp = H*P*H'
    mf16_add(S, &kfm->R, S);                // S += R 

    /************************************************************************/
    /* Calculate Kalman gain                                                */
    /* K = P*H' * S^-1                                                      */
    /************************************************************************/

    // K = P*H' * S^-1
    mf16_cholesky(S, S);
    mf16_invert_lt(S, S);               // Sinv = S^-1
    mf16_mul_bt(&temp_PHt, P, H);           // temp = P*H'
    mf16_mul(K, &temp_PHt, S);          // K = temp*Sinv

    /************************************************************************/
    /* Correct state prediction                                             */
    /* x = x + K*y                                                          */
    /************************************************************************/

    // x = x + K*y 
    mf16_mul(y, K, y);
    mf16_add(x, x, y);

    /************************************************************************/
    /* Correct state covariances                                            */
    /* P = (I-K*H) * P                                                      */
    /*   = P - K*(H*P)                                                      */
    /************************************************************************/

    // P = P - K*(H*P)
    mf16_mul(&temp_HP, H, P);            // temp_HP = H*P
    mf16_mul(K, K, &temp_HP);     // temp_KHP = K*temp_HP
    mf16_sub(P, P, K);           // P -= temp_KHP 
}