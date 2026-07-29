#pragma once
//#include "IO.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
//#include <iterator>
#include <cstring>
#include <map>
#include <list>
#include <deque>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <time.h>
#include <random>
#include <assert.h>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <future>
#include <mutex>
#include <chrono>
#include <random>
#include <sys/stat.h>
#include <cstdint>
#include <cerrno>
#include <limits>
#include <stdexcept>
#include <cctype>
#include "robin_hood.h"




#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#define _gziprea//d
#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif
#else
#define _gzipread
#endif
#ifndef notRpackage
#define notRpackage
#endif


#ifdef _gzipread
#include "gzstream.h"
#endif
typedef double mat_fl;
typedef float smat_fl;

inline bool fileExists(const std::string fi) {
	//struct stat buffer;
	//return (stat(fi.c_str(), &buffer) == 0);
	//return exists(fi);
	std::ifstream file(fi);
	if (!file.is_open()) {
		return false;
	}
	file.close();
	return true;
}

bool isGZfile(const std::string fi);

using namespace std;
typedef unsigned int uint;
typedef unsigned long ulong;
using rtk_count = std::uint64_t;
using global_count_cube = vector<vector<vector<rtk_count>>>;
typedef unordered_map <uint, uint> rare_map;


typedef std::map<std::string, int> GeneIDidx;
//contains gene ID, taxa
//typedef robin_hood::map<std::string, vector<string>> LvlUp;
typedef robin_hood::unordered_map<std::string, vector<string>> LvlUp;
typedef std::unordered_map<string, mat_fl>::iterator SmplAbunIT;
typedef std::unordered_map<string, mat_fl> SmplAbun;
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
	std::string pafFile = "";
	std::string groupFile = "";
	std::vector<double> depth;
	rtk_count depthMin;
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
	double mosaicMinIdentity = 0.90;
	double mosaicMinQueryCoverage = 0.80;
	double mosaicMinTargetCoverage = 0.80;
	double mosaicMinScoreMargin = 0.02;
	double mosaicMinLengthRatio = 0.80;
	double mosaicMinBaseline = 1e-12;
	double mosaicDetectionFraction = 0.10;
	double mosaicMinMergedCopy = 0.70;
	double mosaicMaxMergedCopy = 1.30;
	double mosaicMaxMergedError = 0.30;
	double mosaicMinFitGain = 0.15;
	double mosaicMaxSampleOverlapFraction = 0.15;
	double mosaicMaxAbundanceOverlapFraction = 0.15;
	unsigned int mosaicMinSamples = 5;
	unsigned int mosaicMinVariantSamples = 2;
	unsigned int mosaicMinBaselineGenes = 5;
	bool mosaicPreferMarkerBaseline = true;
	bool mosaicRequireCigar = false;

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
	bool highFromLow;
	bool hieraSrtDown;
	bool median;
	bool check4idxMatch;//assummes tab separated row name, that is number (idx); used in "lineExtr", "extractRows"
	bool gzOut;
	bool header;
	bool calcSupplCov;//for "geneMat" -> calc also abundance for suppl coverage?
	bool normalizeReadPass;


	std::string modDescr;
	std::string modHiera;
	std::string xtra;
	string funcHieraSep; // = ";" #in sumMat, denotes hierachical levels of annotation string
	string funcAnnoAND; // = ","  #in sumMat, denotes annotations that should be merged (summed across genes)
	string funcAnnoOR; // = "|"  #in sumMat, denotes annotations that are equal but undecided (averaged across genes)
};

void cerr2(const std::string x, int ex = 0);
