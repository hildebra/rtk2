#pragma once
#include <iomanip>
#include "IO.h"
#include "options.h"
#include "Math.h"
#include "Fisher.h"


void vecPurge(vector<vector<mat_fl>>& vec, mat_fl val);//removes val from each entry in vec<vec<>>
string join(const vector<string>& in, const string &X);
string join(const vector<int>& in, const string &X);
bool writeChunk(ostream* o, string  s);

typedef robin_hood::unordered_map<uint, mat_fl> sparseMatV;
//typedef std::unordered_map<uint, mat_fl> sparseMatV;

class column{
	public:
		column() :colsum(0.f),id("") {}
		~column() {}
		double colsum;
		string id;

};

class HMat
{
public:
	HMat(string L, vector<string> Samples, vector<string> Features);
	~HMat(){}
	void setDelims(options* opts) { funcAnnoAND = opts->funcAnnoAND; funcAnnoOR = opts->funcAnnoOR; }
	//get
	//unsigned long operator [](int i) const    { return registers[i]; }
	//add single matrix entry (v) to at position j to hierachy level kk
	void set(string& kk, int j, mat_fl v);
	void set(string& kk, vector<mat_fl>& v);
	void prep();

	void print(ofstream&,bool,uint, uint);

private:
	GeneIDidx Feat2mat;
	string LvlName;

	vector<string> FeatureNs, SampleNs;
	vector< mat_fl > empty;
	vector< vector< mat_fl > > mat;
	vector<vector<uint>> catCnt;
	vector<uint> maxCatCnt;
	uint hiTaNAcnt;
	string funcAnnoAND; // = ","  #in sumMat, denotes annotations that should be merged (summed across genes)
	string funcAnnoOR; // = "|"  #in sumMat, denotes annotations that are equal but undecided (averaged across genes)
};

class SparseMatrix
{
public:
	SparseMatrix();
	~SparseMatrix();
	void addCount(string smpl, int row, smat_fl abund);
	//void newRow(string x) { mat.push_back(vector<SmplAbun>(0)); rowIDs.push_back(x); }
	void newRow(void) { mat.push_back(SmplAbun(0)); }
private:
	vector<SmplAbun> mat;
	SmplOccur colNames;
	vector<string> rowIDs;
};

class VecFiles
{
public:
	VecFiles(const string inF, const string outF, const string xtra);
	~VecFiles(){}
private:
	void readVecFile(const string inF);
	int getIdx(const string&);

	//VARIABLES
	vector<string> infiles;
};



struct ModStep
{
public:
	ModStep() :alternates(0), redundancy(0){}
	ModStep(const string&, bool&, vector<string>&);
	void getAllKOs(list<string>&);
	void setRedund(ModOccur& m);
	void abundParts(const vector<mat_fl>& v, const unordered_map<string, int>& IDX,
		vector<mat_fl>&, vector<bool>&, vector<string>&,
		float hitComplRatio =0.8, int redund=0);



	//e.g. alt[0][0] = KO001 and requires alt[0][1] =KO002 or alternatively only alt[1][0]=K0003
	vector< vector< string > > alternates;
	//how redundant is each KO and the different steps (basically occurence of KOs across DB)
	vector< vector <int> > redundancy;
};

class Module
{
public:
	Module(vector<string>& n);
	void getAllKOs(list<string>& r) { for (size_t i = 0; i < steps.size(); i++) { steps[i].getAllKOs(r); } }
	void setReddundancy(ModOccur& m) { for (size_t i = 0; i < steps.size(); i++) { steps[i].setRedund(m); } }
	mat_fl pathAbundance(const vector<mat_fl>& v, const unordered_map<string, int>& IDX,
		const int redun, const float PathwCompl, const float enzymCompl, string & savedNmsKO,float & modScoreOut);
	string name;
	string description;
	vector<ModStep> steps;
//handles module recurrency
	vector<string> submods;
	bool containsMods,usedInOtherMods;
};


struct declutRes {
	declutRes(int nsize) :smpls(0),
		idxs(0), curClus(nsize, ""), totalMergedID(0), totalMergedCE(0) {}
	vector<vector<string>> smpls;
	vector<vector<int>> idxs;
	vector<string> curClus;
	int totalMergedID;
	int totalMergedCE;

};
//internal structures
struct job2 {
	std::future <declutRes*> fut;
	bool inUse = false;
};

struct rowRes {
	rowRes(vector<mat_fl>& x) :rowV(x), NAcnt(0), curRowIdx(0), rowID(""), taxa(0),Occ(0) {}
	rowRes(int s) :rowV(s, (mat_fl)0), NAcnt(0), curRowIdx(0), rowID(""), taxa(0),Occ(0){}
	vector<mat_fl> rowV;
	vector<string> taxa;
	int NAcnt;
	int curRowIdx;
	int Occ;
	string rowID;
};

struct job3 {
	std::future <rowRes*> fut;
	bool inUse = false;
};


struct job4 {
	std::future <string> fut;
	bool inUse = false;
};

class Matrix
{
// convention: mat[smpl][feature]
public:
	//Matrix(const string inF);
	//read and write
	Matrix(options*, const string xtra,vector<string>& outFName, bool highLvl = false, bool NumericRowId = false, bool writeTmpFiles = true);
	//read to mem
	Matrix(options*, const string xtra, bool highLvl = false);
	//module abundance matrix
	Matrix(const vector<string>& rnms, const vector<string>& cnms);
	//empty opbject
	Matrix(void);
	//normalize on the fly on vector colSums
	Matrix(const string inF, const string outF, vector< double> colsums, 
		vector<string>colNmds, const char sep, options* opts);
	~Matrix(void);
	void addTtlSmpl(vector<mat_fl> x, int idx) { mat[idx] = x; }
	void splitOnHDD(string out_seed);
	void writeSums(string);
	void normalize();
	void transpose();
	//requires clStr object from
	declutRes* declutCalc(vector<string>* txt, options* opts);
	declutRes* declutCalc2(vector<string>* txt, options* opts);
	void decluter(options*);//CD-HIT version
	void decluter2(options*);//mmseqs2 version
	void writeMatrix(const string ofile,bool onlyFilled=false,int threads=1);
	size_t smplNum(){ return colIDs.size(); }
	int rowNum(){ return rowIDs.size(); }

	smplVec* getSampleVec(uint which){ return new smplVec(mat[which],1); }
	string getSampleName(uint which){ return colIDs[which]; }

	int SmplNum() { return (int)mat.size(); }
	int FtNum() {
		if (mat.size() >= 1) { return (int)mat[0].size(); }
		else { return 0; }
	}
	void estimateModuleAbund(char ** args, int argc);
	void estimateModuleAbund(options*);
	void resizeMatRows(uint x,mat_fl def=(mat_fl)0);
	//for the R module, all used for rarefactions only
	void addRow(string rname, vector<mat_fl>);//idea is that a single row is added on to the matrix
	//for multicore file reader
	void addRow(rowRes* ret);
	void setSampleNames(vector<string> in) { colIDs = in; }
	void setRowNames(vector<string> in) { rowIDs = in; }
	vector < string > getSampleNames(){ return(colIDs); }
	vector < string > getRowNames(){ return(rowIDs); }
	//void addCount(string, int, mat_fl);

	double getMinColSum();
	column getMinColumn(uint offset = 0);
	vector< pair <double, string>> getColSums(bool sorted = false);
	vector<double> getCSum() { return colSum; }
	void writeColSums(string outF);
	void writeRowOcc(options*);
	void prepCoex(options* opt);
protected:

	rowRes* processLine2Mat1(string  line, int, int, int, const char);
	string processLine2line(string  line, int, int, int, const char);
	string prepRow(rowRes* ret, string SEP, bool normalize = false);
	string row2line(int i);

	//test a list of genes for complimentarity, e.g. actually same gene just clutered weird
	void complimentarity(vector<string>&, vector<float>&, options*,
		vector<vector<int>> &, 
		vector<vector<string>>&, int&, int&);
	bool check_compli(int x, int y);
	//reads the number of columns and checks in first few lines
	void readColNms(istream* in, const char sep);
	int iniCols(istream* in, const char sep);
	void read_subset_genes(const string);
	void read_hierachy(const string, options*);// bool);
	void addColumn(string);
	vector<mat_fl> getRow(uint idx);
	//merge two rows in matrix, into first entry, setting subsequent rows to 0
	void mergeRows(vector<int>&);
	void readModuleFile(const string&);
	vector<mat_fl> getRowSums();
	void ini_mat();

	//storage
	//convention is: mat[columnIdx][rowIDx]
	vector< vector< mat_fl > > mat;
	vector< sparseMatV > matSp;
	vector< string > rowIDs,colIDs;
	unordered_map<string, int> colID_hash, rowID_hash;
	vector<HMat*> HI;
	LvlUp LUp; //contains hierachy info (gene -> tax levels/function levels..
	int maxLvl;
	vector<string> LvlNms;
	string sampleNameSep;
	GeneIDidx subset;
	bool doSubsets, doHigh;
	
	vector<double> colSum;
	vector<int> rowOcc;
	int maxCols;//number samples
	int coexcl_thr;
	float CoEx_thresh; 
	float CoEx_minp;
	
	FET* fet;
	mutex protect;

	vector< pair <double, string>> colsums;

	bool sparse; //sparse matrix?
	bool loadPrimaryMat; //save matrix in ram? or just high etc?

};



/*
class SigMatrix :public Matrix {
public:
	SigMatrix(const string& inf) { Matrix(inf, ""); }
private:
	void estimateBinModel() { 
#ifdef notRpackage
cerr << "todo estimateBinModel"; exit(33); 
#endif
}
	//functions to determine parameters
	//stored model parameters

};
*/
class Modules : public Matrix
{
public:
	Modules(const string&, vector<string>);
	~Modules() {}

	void addDescription(const string&);
	void addHierachy(const string&);
	void writeModDescr(const string&,bool onlyUsed);


	void setRedund(int x) { redund = x; }
	void setPathwCompl(float x) { PathwCompl = x; }
	void setEnzymCompl(float x) { enzymCompl = x; }

	void calcModAbund(vector<mat_fl>&, const int pos, const unordered_map<string, int>&,
		vector<string>&, vector<float>&);


	vector<string> & modNms() { return rowIDs; }
	vector<string> modNms_numbered();//puts a number behind double modules
	vector<string> & modDescr() { return moduleDescriptions; }
	vector<string> & getRedundantUsedMods() { return redundantUsedMods; }
//special implementation to collapse rows..
	void writeMatrix(const string ofile, bool onlyFilled = false, bool collapseDblFeats = false);


private:
	void calc_redund();
	//contains the modules in the DB, each entry being one module
	vector<Module> mods;
	vector<string>  moduleDescriptions, redundantUsedMods; 
	vector<vector<string>> hierachy;
	//moduleNames = rowIDs
	//list of KOs used in DB, and how often they occur
	ModOccur MO;
	//in case of double entries, track these
	unordered_map<string, vector<int>> ModPos;
	vector<int> recurrentMods;
	vector<bool> ModUsed;

	//list of options
	int redund; // max redundancy of KOs used
	float PathwCompl; //corresponds to -c
	float enzymCompl; //single enzymes complexes - how much needs to be present to trigger present

};



