// Settings header for libfixkalman

#ifndef _SETTINGS_H_
#define _SETTINGS_H_

/* ----------------------------------------------------------------------- */
/* Replace of the covariance update equation with the Joseph form.         */
/* Set this if Kalman gain is non-optimal or your filter has problems with */ 
/* numerical stability. Please be aware that this form of the covariance   */
/* update equation is computationally more expensive.                      */
/* ----------------------------------------------------------------------- */
#define KALMAN_JOSEPH_FORM

#endif
