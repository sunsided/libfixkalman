// Settings header for libfixkalman

#ifndef _SETTINGS_H_
#define _SETTINGS_H_

/* ----------------------------------------------------------------------- */
/* Replace of the covariance update equation with the Joseph form.         */
/* Set this if Kalman gain is non-optimal or your filter has problems with */ 
/* numerical stability. Please be aware that this form of the covariance   */
/* update equation is computationally more expensive.                      */
/* ----------------------------------------------------------------------- */
//#define KALMAN_JOSEPH_FORM

/* ----------------------------------------------------------------------- */
/* Transforms the square system process noise matrix into square contol    */
/* input covariance matrix for Kalman filter with control input.           */
/*                                                                         */
/* The square system process noise matrix has as many of rows and columns  */
/* as the matrix A.                                                        */
/* In this case the covariance prediction step looks like: P = A*P*A' + Q  */
/*                                                                         */
/* By the square contol input covariance matrix the number of rows and     */
/* columns are equal to the number of rows in the input vector u.          */
/* In this case the covariance prediction step looks like:                 */
/* P = A*P*A' + B*Q*B'                                                     */  
/* ----------------------------------------------------------------------- */
//#define KALMAN_TIME_VARYING

#endif
