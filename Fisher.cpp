#include "Fisher.h"

double FET::ulogHypergeometricProb(int a, int b, int c, int d) {
	return logFacs[a + b] + logFacs[c + d] + logFacs[a + c] + logFacs[b + d]
		- logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a + b + c + d];
}


void FET::uinitLogFacs(uint n) {
	uint curLFs = static_cast<uint>(logFacs.size());
	if ((n + 1) <= curLFs) { return; }
	logFacs.resize(n + 1);
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

	if (a < 0 || b < 0 || c < 0 || d < 0) {
		throw invalid_argument("Fisher exact test received a negative table cell");
	}
	long long total = static_cast<long long>(a) + b + c + d;
	if (total > numeric_limits<int>::max()) {
		throw overflow_error("Fisher exact test table is too large");
	}
	int m = static_cast<int>(total);

	// Keep readers protected as well as vector growth: resize invalidates storage used by ufet.
	lock_guard<mutex> guard(protect);
	uinitLogFacs((uint) m);

	pval = (float)ufet(a, b, c, d, m);

	return (pval);
}
