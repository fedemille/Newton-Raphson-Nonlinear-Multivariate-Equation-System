#include <stdio.h>
#include "triangulation.h"

int main()
{
    
    funs fs[] = {&functo1, &functo2, &functo3};

  	funs f1[] = {&func1, &func2, &func3};
  	funs f2[] = {&func4, &func5, &func6};
  	funs f3[] = {&func7, &func8, &func9};
  	jacobian jacob[] = {f1, f2, f3};



  	double guesses[] = {1, 1, 0};		// initial values

  	double sol[3];

    solve(sol, fs, 3, jacob, guesses, 3);

  	printf("The approximate solutions are x = %.8f - y = %.8f - z = %.8f\n", sol[0], sol[1], sol[2]);

    return 0;
}
