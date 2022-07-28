/*
 * triangulation.c
 *
 *  Created on: 28 lug 2022
 *      Author: fedemille
 */
#include "triangulation.h"





double* mul(double* result, double* m1, double* m2, uint16_t r1, uint16_t c1, uint16_t r2, uint16_t c2){
	if (c1 != r2) {
		#ifdef PRINT_ERRORS
		printf("Errore molt");
		#endif
		return NULL;
	}


	for (int i = 0; i < r1; i++) {
		for (int j = 0; j < c2; j++) {
			for (int k = 0; k < r2; k++) {
				*((result+i*c2) + j) += *((m1+i*c1) + k) * *((m2+k*c2) + j);
			}
		}
	}
	return result;
}



double* sub(double* result, double* m1, double* m2, int r1, int c1, int r2, int c2){
	if (r1 != r2 || c1 != c2) {
		#ifdef PRINT_ERRORS
		printf("Errore sottr %d %d %d %d\n", r1, c1, r2, c2);
		#endif
		return NULL;
	}


	for (int i = 0; i < r1; i++) {
		for (int j = 0; j < c1; j++){
			*((result+i*c1) + j) = *((m1+i*c1) + j) - *((m2+i*c2) + j);
		}
	}
	return result;
}





double* transpose(double *trans, double *m, uint16_t r, uint16_t c){
	for (int i = 0; i < c; i++) {
		for (int j = 0; j < r; j++) {
  			*((trans+j*c) + i) = *((m+j*c) + i);
  		}
  	}
  	return trans;
}



void toReducedRowEchelonForm(double* m, uint16_t row, uint16_t col){
	int lead = 0;
	int i;
	for (int r = 0; r < row; r++) {
		if (col <= lead) {
			return;
		}
		i = r;

		for (;*((m+i*col) + lead) == 0;) {
			i++;
			if (row == i) {
				i = r;
				lead++;
				if (col == lead) {
					return;
				}
			}
		}


		double temp;
		for(int d=0; d<col; d++){
			temp = *((m+i*col) + d);
			*((m+i*col) + d) = *((m+r*col) + d);
			*((m+r*col) + d) = temp;
		}

		double div;
		if ((div = *((m+r*col) + lead)) != 0) {
			for (int j = 0; j < col; j++) {
				*((m+r*col) + j) /= div;
			}
		}

		double mult;
		for (int k = 0; k < row; k++) {
			if (k != r) {
				mult = *((m+k*col) + lead);
				for (int j = 0; j < col; j++) {
					*((m+k*col) + j) -= *((m+r*col) + j) * mult;
				}
			}
		}
		lead++;
	}
}



double* inverse(double* inv, double* m, uint16_t r, uint16_t c){

	if(r != c){
		#ifdef PRINT_ERRORS
		printf("Errore Deve essere una matrice quadrata per invertirla");
		#endif
		return NULL;
	}

	double aug[r][2*r];

	for(int i=0;i<r;i++){
		for(int j=0; j<r*2; j++){
			aug[i][j] = 0;
		}
	}

	for (int i = 0; i < r; i++) {
		for(int j=0; j<r; j++){
			aug[i][j] = *((m + i*c) + j);
		}

		aug[i][i+r] = 1;
	}


	toReducedRowEchelonForm(aug, r, 2*r);


	for(int i = 0; i < r; i++) {
		for(int j=0; j<r; j++){
		   *((inv + i*c) + j)  = aug[i][j+r];
		}
	}

	return inv;
}



void printfMatrix(double *matr, int a, int b){
#ifdef PRINT_FLOAT_VALUES
    for(int i=0; i<a; i++){
        for(int j=0; j<b; j++){

            printf("%f - ", *((matr+i*b) + j));
        }
        printf("\n");

    }
#endif
}


void solve(double *gu2, funs *fs, int nfuns, jacobian *jacob, double* guesses, int glen){
	int size = nfuns;
	double gu1[glen];

	memcpy(gu2, guesses, glen*sizeof(double));

	double jac[size][size];
	double tol = 1e-8;
	int maxIter = 12;
	int iter = 0;


	while (1){
		for(int i=0; i<glen; i++)
			gu1[i] = gu2[i];

		double gtemp[1][glen];
		for(int i=0; i<glen; i++)
			gtemp[0][i] = gu1[i];



		double g[glen][1];
		transpose(g, gtemp, 1, glen);


		double t[size];
		for (int i = 0; i < size; i++) {
			t[i] = fs[i](gu1);
		}


		double f[size][1];

		for(int i=0; i<size; i++)
			f[i][0] = t[i];


		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				jac[i][j] = jacob[i][j](gu1);
			}
		}

		double g1[glen][1];
		double tempInv[size][1];
		double tempM[size][size];

		for(int i=0; i<glen; i++)
			g1[i][0] = 0;
		for(int i=0; i<size; i++)
			tempInv[i][0] = 0;
		for(int i=0; i<size; i++)
			for(int j=0; j<size; j++)
			tempM[i][j] = 0;


		inverse(tempM, jac, size, size);
		mul( tempInv, tempM , f , size, size, size , 1);
		sub( g1 , g , tempInv , glen, 1, size, 1 );


		for (int i = 0; i < size; i++) {
			gu2[i] = g1[i][0];
		}
		iter++;
		uint8_t any = 0;


		for(int i=0; i<size; i++){
			double v = gu2[i];
			if (((v<0?-v:v)-gu1[i]) > tol) {
				any = 1;
				break;
			}
		}
		if (!any || iter >= maxIter) {
			break;
		}
	}

	return gu2;
}




