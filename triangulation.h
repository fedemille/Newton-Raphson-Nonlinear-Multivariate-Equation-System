/*
 * triangulation.h
 *
 *  Created on: Jul 28, 2022
 *      Author: fedemille
 */

#ifndef INC_TRIANGULATION_H_
#define INC_TRIANGULATION_H_


#include "main.h"


//#define PRINT_ERRORS
//#define PRINT_FLOAT_VALUES





typedef double (*funs)(double*);
typedef funs* jacobian;


double* mul(double* result, double* m1, double* m2, uint16_t r1, uint16_t c1, uint16_t r2, uint16_t c2);
double* sub(double* result, double* m1, double* m2, int r1, int c1, int r2, int c2);
double* transpose(double *trans, double *m, uint16_t r, uint16_t c);
void toReducedRowEchelonForm(double* m, uint16_t row, uint16_t col);
double* inverse(double* inv, double* m, uint16_t r, uint16_t c);
void solve(double *gu2, funs *fs, int nfuns, jacobian *jacob, double* guesses, int glen);



///////// EQUATIONS //////////

static double functo1(double x[]) {
	return ((x[0] - 10) * (x[0] - 10)) + ((x[1] - 0) * (x[1] - 0)) + ((x[2] - 10) * (x[2] - 10)) - 144;		// equation 1
}
static double functo2(double x[]) {
	return ((x[0] - 0) * (x[0] - 0)) + ((x[1] - 10) * (x[1] - 10)) + ((x[2] - 10) * (x[2] - 10)) - 144;		// equation 2
}
static double functo3(double x[]) {
	return ((x[0] - 10) * (x[0] - 10)) + ((x[1] - 10) * (x[1] - 10)) + ((x[2] - 10) * (x[2] - 10)) - 144;	// equation 3
}







///////// DERIVATES //////////

static double func1(double x[]) { return 2 * (x[0] - 10); }	// derivata di x[0] nella prima eq
static double func2(double x[]) { return 2 * (x[1] - 0); }		// derivata di x[1] nella prima eq
static double func3(double x[]) { return 2 * (x[2] - 10); }

static double func4(double x[]) { return 2 * (x[0] - 0); }		// derivata di x[0] nella seconda eq
static double func5(double x[]) { return 2 * (x[1] - 10); }	// derivata di x[1] nella seconda eq
static double func6(double x[]) { return 2 * (x[2] - 10); }

static double func7(double x[]) { return 2 * (x[0] - 10); }
static double func8(double x[]) { return 2 * (x[1] - 10); }
static double func9(double x[]) { return 2 * (x[2] - 10); }





#endif /* INC_TRIANGULATION_H_ */
