#pragma once
//#include "IO.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
//#include <iterator>
#include <cstring>
#include <map>
#include <list>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <time.h>
#include <random>
#include <assert.h>
#include <unordered_map>
#include <numeric>
#include <future>
#include <mutex>
#include <chrono>
#include <random>
#include "robin_hood.h"



#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#define _gziprea//d
#pragma warning(disable:4996)
#else
#define _gzipread
#endif
#define notRpackage


#ifdef _gzipread
#include "gzstream.h"
#endif
typedef double mat_fl;
typedef float smat_fl;

bool isGZfile(const std::string fi);

using namespace std;
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unordered_map <uint, uint> rare_map;


typedef std::map<std::string, int> GeneIDidx;
//contains gene ID, taxa
//typedef robin_hood::map<std::string, vector<string>> LvlUp;
typedef robin_hood::unordered_map<std::string, vector<string>> LvlUp;
typedef std::unordered_map<string, smat_fl>::iterator SmplAbunIT;
typedef std::unordered_map<string, smat_fl> SmplAbun;
typedef std::unordered_map<string, vector<int> >::iterator SmplOccurITmult;
typedef std::unordered_map<string, vector<int> > SmplOccurMult;
typedef std::unordered_map<string, string > string2string;
typedef std::unordered_map<string, int >::iterator SmplOccurIT;
typedef std::unordered_map<string, int> SmplOccur;
typedef std::unordered_map<string, int> ModOccur;



struct options
{
public:
	options(int argc, char** argv);
	//	options(std::string, std::string , int repeats, std::vector<double> depth, 
	//		int NoOfMatrices, 
	//		bool verbose, unsigned int threads);
	void print_rare_details();
	//~options();

	//vars
	std::string input = "";
	std::string output = "";
	std::string mode  = "";
	std::string referenceDir = "";
	std::string referenceFile = "";
	std::string map = "";
	std::vector<double> depth;
	long depthMin;
	char sepChar;
	unsigned int repeats = 10;
	unsigned int write = 0;
	unsigned int threads = 1;
	unsigned int occPerSmpl = 0;
	unsigned int occMin = 0;
	bool writeSwap = true;
	bool verbose = false;
	bool oldMapStyle = true;
	bool sparse = true;
	float pval = (float)1e-5;//for fisher test in decluter

	std::string modDB;
	int modRedund;
	float modEnzCompl;
	float modModCompl;
	bool modWrXtraInfo;
	bool modCollapse;
	bool calcCoverage;
	bool calcCovMedian;
	bool extendHierachy = false;
	bool mean;
	bool median;
	bool check4idxMatch;//assummes tab separated row name, that is number (idx); used in "lineExtr", "extractRows"
	bool gzOut;
	bool header;

	std::string modDescr;
	std::string modHiera;
	std::string xtra;
	string funcHieraSep; // = ";" #in sumMat, denotes hierachical levels of annotation string
	string funcHAnnoAND; // = ","  #in sumMat, denotes annotations that should be merged (summed across genes)
	string funcAnnoOR; // = "|"  #in sumMat, denotes annotations that are equal but undecided (averaged across genes)
};

void cerr2(const std::string x, int ex = 0);