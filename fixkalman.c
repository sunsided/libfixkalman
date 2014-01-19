#include <stdint.h>
#include <assert.h>

#include "fixarray.h"

#define EXTERN_INLINE_KALMAN INLINE
#include "fixkalman.h"

typedef enum {
    ADD = 0,
    SUB = 1,
} addsub_t;

// Calculates dest = dest +/- A * B
HOT PURE NONNULL void mf16_mul_addsub(mf16 *dest, const mf16 *a, const mf16 *b, addsub_t mode)
{
    int row, column;
    const int 
        acolumns = a->columns,
        drows = dest->rows, 
        dcolumns = dest->columns;

    // If dest and input matrices alias, we have to use a temp matrix.
    mf16 tmp;
    fa16_unalias(dest, (void**)&a, (void**)&b, &tmp, sizeof(tmp));

    dest->errors = a->errors | b->errors;

    if (a->columns != b->rows)
        dest->errors |= FIXMATRIX_DIMERR;

    dest->rows = a->rows;
    dest->columns = b->columns;

    for (row = drows-1; row >= 0; --row)
    {
        const fix16_t *aptr = &a->data[row][0];
        const fix16_t *bptr = &b->data[0][0];
        fix16_t *rowptr = &dest->data[row][0];

        for (column = dcolumns-1; column >= 0; --column)
        {
            fix16_t value = fa16_dot(
                aptr, 1,
                bptr + column, FIXMATRIX_MAX_SIZE,
                acolumns);
            
            if (mode == ADD)
                rowptr[column] += value;
            else
                rowptr[column] -= value;

            if (rowptr[column] == fix16_overflow)
                dest->errors |= FIXMATRIX_OVERFLOW;
        }
    }
}

/*!
* \brief Calculates A*B*A'
*/
NONNULL
STATIC_INLINE void mf16_mul_abat(mf16 *dest, const mf16 *a, const mf16 *b)
{
    mf16_mul(dest, a, b);
    mf16_mul_bt(dest, dest, a);
}

/*!
* \brief Calculates dest += A*B*A'
*/
NONNULL
STATIC_INLINE void mf16_mul_abat_add(mf16 *dest, const mf16 *a, const mf16 *b)
{
    mf16_mul_bt(dest, b, a);                // dest = B * A'
    mf16_mul_addsub(dest, a, dest, ADD);    // dest += A * dest

}


#ifndef KALMAN_DISABLE_C

/*!
* \brief Initializes the Kalman Filter
* \param[in] kf The Kalman Filter structure to initialize
* \param[in] num_states The number of state variables
* \param[in] num_inputs The number of input variables
*/
void kalman_filter_initialize(kalman16_t *const kf, uint_fast8_t num_states, uint_fast8_t num_inputs)
{
    kf->A.rows = num_states;
    kf->A.columns = num_states;
    kf->A.errors = 0;
    mf16_fill(&kf->A, 0);
    
    kf->P.rows = num_states;
    kf->P.columns = num_states;
    kf->P.errors = 0;
    mf16_fill(&kf->P, 0);

    kf->B.rows = num_inputs;
    kf->B.columns = num_inputs;
    kf->B.errors = 0;
    mf16_fill(&kf->B, 0);

    kf->Q.rows = num_inputs;
    kf->Q.columns = num_inputs;
    kf->Q.errors = 0;
    mf16_fill(&kf->Q, 0);

    kf->x.rows = num_states;
    kf->x.columns = 1;
    kf->x.errors = 0;
    mf16_fill(&kf->x, 0);

    kf->u.rows = num_inputs;
    kf->u.columns = 1;
    kf->u.errors = 0;
    mf16_fill(&kf->u, 0);
}

#endif

#ifndef KALMAN_DISABLE_UC

/*!
* \brief Initializes the Kalman Filter
* \param[in] kf The Kalman Filter structure to initialize
* \param[in] num_states The number of state variables
*/
void kalman_filter_initialize_uc(kalman16_uc_t *const kf, uint_fast8_t num_states)
{
    kf->A.rows = num_states;
    kf->A.columns = num_states;
    kf->A.errors = 0;
    mf16_fill(&kf->A, 0);

    kf->P.rows = num_states;
    kf->P.columns = num_states;
    kf->P.errors = 0;
    mf16_fill(&kf->P, 0);

    kf->x.rows = num_states;
    kf->x.columns = 1;
    kf->x.errors = 0;
    mf16_fill(&kf->x, 0);

    kf->Q.rows = num_states;
    kf->Q.columns = num_states;
    kf->Q.errors = 0;
    mf16_fill(&kf->Q, 0);
}

#endif // KALMAN_DISABLE_UC

/*!
* \brief Sets the measurement vector
* \param[in] kfm The Kalman Filter measurement structure to initialize
* \param[in] num_states The number of states
* \param[in] num_measurements The number of observations
*/
void kalman_observation_initialize(kalman16_observation_t *const kfm, uint_fast8_t num_states, uint_fast8_t num_observations)
{
    kfm->H.rows = num_observations;
    kfm->H.columns = num_states;
    kfm->H.errors = 0;
    mf16_fill(&kfm->H, 0);

    kfm->R.rows = num_observations;
    kfm->R.columns = num_observations;
    kfm->R.errors = 0;
    mf16_fill(&kfm->R, 0);

    kfm->z.rows = num_observations;
    kfm->z.columns = 1;
    kfm->z.errors = 0;
    mf16_fill(&kfm->z, 0);
}

#ifndef KALMAN_DISABLE_C

/*!
* \brief Performs the time update / prediction step of only the state vector
* \param[in] kf The Kalman Filter structure to predict with.
*/
void kalman_predict_x(register kalman16_t *const kf)
{
    // matrices and vectors
    const mf16 *RESTRICT const A = &kf->A;
    const mf16 *RESTRICT const B = &kf->B;
    const mf16 *RESTRICT const u = &kf->u;
    mf16 *RESTRICT const x = &kf->x;

    /************************************************************************/
    /* Predict next state using system dynamics                             */
    /* x = A*x                                                              */
    /************************************************************************/

    // x = A*x
    mf16_mul(x, A, x);
    
    // TODO: Implement more efficient matrix/row vector multiply

    if (B->rows > 0)
    {
        mf16_mul_addsub(x, B, u, ADD);       // x += B*u
    }
}

#endif

#ifndef KALMAN_DISABLE_UC

/*!
* \brief Performs the time update / prediction step of only the state vector
* \param[in] kf The Kalman Filter structure to predict with.
*/
void kalman_predict_x_uc(register kalman16_uc_t *const kf)
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
}

#endif // KALMAN_DISABLE_UC

#ifndef KALMAN_DISABLE_C

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*/
void kalman_predict_P(register kalman16_t *const kf)
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
    mf16_mul_abat(P, A, P);                     // P = A*P*A'

    // TODO: this should work without P_temp
    // TODO: extract the code block below

    // P = P + B*Q*B'
    if (B->rows > 0)
    {
        mf16_mul_abat_add(P, B, Q);             // P += B*Q*B'
    }
}

#endif

#ifndef KALMAN_DISABLE_UC

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*/
void kalman_predict_P_uc(register kalman16_uc_t *const kf)
{
    // matrices and vectors
    const mf16 *RESTRICT const A = &kf->A;
    const mf16 *RESTRICT const Q = &kf->Q;
    mf16 *RESTRICT const P = &kf->P;

    /************************************************************************/
    /* Predict next covariance using system dynamics and input              */
    /* P = A*P*A' + B*Q*B'                                                  */
    /************************************************************************/

    // P = A*P*A'
    mf16_mul_abat(P, A, P);                 // P = A*P*A'
    mf16_add(P, P, Q);                      // P += Q
}

#endif // KALMAN_DISABLE_UC

#ifndef KALMAN_DISABLE_C
#ifndef KALMAN_DISABLE_LAMBDA

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*/
void kalman_predict_P_tuned(register kalman16_t *const kf, fix16_t lambda)
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
    mf16_mul_abat(P, A, P);                     // temp = A*P*A'
    mf16_mul_s(P, P, lambda);                   // P *= 1/(lambda^2)

    // TODO: extract the code block below

    // P = P + B*Q*B'
    if (B->rows > 0)
    {
        mf16_mul_abat_add(P, B, &kf->Q);     // temp = B*Q*B'
    }
}

#endif // KALMAN_DISABLE_LAMBDA
#endif

#ifndef KALMAN_DISABLE_UC
#ifndef KALMAN_DISABLE_LAMBDA

/*!
* \brief Performs the time update / prediction step of only the state covariance matrix
* \param[in] kf The Kalman Filter structure to predict with.
*/
void kalman_predict_P_tuned_uc(register kalman16_uc_t *const kf, fix16_t lambda)
{
    // matrices and vectors
    const mf16 *RESTRICT const A = &kf->A;
    const mf16 *RESTRICT const Q = &kf->Q;
    mf16 *RESTRICT const P = &kf->P;

    /************************************************************************/
    /* Predict next covariance using system dynamics and input              */
    /* P = A*P*A' * 1/lambda^2 + B*Q*B'                                     */
    /************************************************************************/

    // lambda = 1/lambda^2
    lambda = fix16_div(fix16_one, fix16_sq(lambda)); // TODO: This should be precalculated, e.g. using kalman_set_lambda(...);

    // P = A*P*A'
    mf16_mul_abat(P, A, P);                 // temp = A*P*A'
    mf16_mul_s(P, P, lambda);               // P *= 1/(lambda^2)
    mf16_add(P, P, Q);                      // P += Q
}

#endif // KALMAN_DISABLE_LAMBDA
#endif // KALMAN_DISABLE_UC

#if !defined(KALMAN_DISABLE_C) || !defined(KALMAN_DISABLE_UC)

/*!
* \brief Performs the measurement update step.
* \param[in] kf The Kalman Filter structure to correct.
*/
void kalman_correct(kalman16_t *kf, kalman16_observation_t *kfm)
{
    mf16 *RESTRICT const        x = &kf->x;
    mf16 *RESTRICT const        P = &kf->P;
    const mf16 *RESTRICT const  H = &kfm->H;

    static mf16 K, S, y,
        temp_HP,    // { H->rows, H->columns };
        temp_PHt;   // { P->rows, H->columns };

    /************************************************************************/
    /* Calculate innovation and residual covariance                         */
    /* y = z - H*x                                                          */
    /* S = H*P*H' + R                                                       */
    /************************************************************************/

    // y = z - H*x
    mf16_mul(&y, H, x);
    mf16_sub(&y, &kfm->z, &y);

    // S = H*P*H' + R
    mf16_mul_abat(&S, H, P);               // temp = H*P*H'
    mf16_add(&S, &kfm->R, &S);                // S += R 

    /************************************************************************/
    /* Calculate Kalman gain                                                */
    /* K = P*H' * S^-1                                                      */
    /************************************************************************/

    // K = P*H' * S^-1
    mf16_cholesky(&S, &S);
    mf16_invert_lt(&S, &S);               // Sinv = S^-1
    mf16_mul_bt(&temp_PHt, P, H);         // temp = P*H'
    mf16_mul(&K, &temp_PHt, &S);          // K = temp*Sinv

    /************************************************************************/
    /* Correct state prediction                                             */
    /* x = x + K*y                                                          */
    /************************************************************************/

    // x = x + K*y 
    mf16_mul_addsub(x, &K, &y, ADD);

    /************************************************************************/
    /* Correct state covariances                                            */
    /* P = (I-K*H) * P                                                      */
    /*   = P - K*(H*P)                                                      */
    /************************************************************************/

    // P = P - K*(H*P)
    mf16_mul(&temp_HP, H, P);                  // temp_HP = H*P
    mf16_mul_addsub(P, &K, &temp_HP, SUB);       // P -= K*temp_HP
}

#endif

#ifndef KALMAN_DISABLE_UC

/*!
* \brief Performs the measurement update step.
* \param[in] kf The Kalman Filter structure to correct.
*/
void kalman_correct_uc(kalman16_uc_t *kf, kalman16_observation_t *kfm)
{
    // just be careful, kid
    kalman_correct((kalman16_t*)kf, kfm);
}

#endif // KALMAN_DISABLE_UC
