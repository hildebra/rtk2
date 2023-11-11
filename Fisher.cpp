#include "Fisher.h"

double FET::ulogHypergeometricProb(int a, int b, int c, int d) {
	return logFacs[a + b] + logFacs[c + d] + logFacs[a + c] + logFacs[b + d]
		- logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a + b + c + d];
}


void FET::uinitLogFacs(uint n) {
	unsigned int curLFs = logFacs.size();
	if ((n + 1) > curLFs) {
		lock_guard<mutex> guard(protect);
		logFacs.resize(n + 1);
	}
	else {
		return;
	}
	for (uint i = curLFs; i < (n + 1); ++i) {
		if (i == 0) {
			logFacs[0] = 0;
		}
		else {
			logFacs[i] = logFacs[i - 1] + log((double)i);
		}
	}
}


double FET::ufet(int a, int b, int c, int d, int n) {
	double f = 10000000;
	double logpCutoff = round(ulogHypergeometricProb( a, b, c, d) * f) / f;
	double pFraction = 0;
	for (int x = 0; x <= n; ++x) {
		if (a + b - x >= 0 && a + c - x >= 0 && d - a + x >= 0) {
			double l = round(ulogHypergeometricProb( x, a + b - x, a + c - x, d - a + x) * f) / f;
			if (l <= logpCutoff) {
				pFraction += exp(l - logpCutoff);
			}
		}
	}
	double logpValue = logpCutoff + log(pFraction); // normalization: p_i = p_i / sum_i(p_i)
	double pval = exp(logpValue);
	pval = fmax(fmin(pval, 1.0), 0.0);

	return pval;
}


float FET::ultrafastfet(int a, int b, int c, int d) {

	float pval(0);

	int m = a + b + c + d;

	uinitLogFacs((uint) m);

	pval = (float)ufet(a, b, c, d, m);

	return (pval);
}