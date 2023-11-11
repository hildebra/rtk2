#pragma once
#include "options.h"
class FET {
public:
	FET():logFacs(0){}
	float ultrafastfet(int, int, int, int);
private:
	double ufet(int a, int b, int c, int d, int n);
	void uinitLogFacs(uint n);
	double ulogHypergeometricProb(int a, int b, int c, int d);

	//storage
	std::vector<double> logFacs;
	mutex protect;
};