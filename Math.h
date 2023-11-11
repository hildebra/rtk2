#pragma once
#include "options.h"


//const bool verbose=1;


int getRand(int until);

void swap(int &x,int &y);



/*
int fac(int n) { // calculates factorial
	int ret;
	for (ret = 1; n > 0; --n) { ret *= n; }
	return ret;
}
float hypergeometricProb(int a, int b, int c, int d) {
	int num = fac(a + b) * fac(c + d) * fac(a + c) * fac(b + d);
	int den = fac(a) * fac(b) * fac(c) * fac(d) * fac(a + b + c + d);
	return (float)num / (float)den;
}

float fisher_test(int a, int b, int c, int d) {
	int n = a + b + c + d;
	// find cutoff probability
	float pCutoff = hypergeometricProb(a, b, c, d);
	float pValue = 0;
	// sum over probability smaller than the cutoff
	for (int x = 0; x <= n; ++x) { // among all possible x
		if (a + b - x >= 0 && a + c - x >= 0 && d - a + x >= 0) { // consider valid x
			float p = hypergeometricProb(x, a + b - x, a + c - x, d - a + x);
			if (p <= pCutoff) pValue += p;
		}
	}
	std::cout << "Two-sided p-value is " << pValue << std::endl;
	return(pValue);
}
*/
