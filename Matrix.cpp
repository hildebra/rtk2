#include "Matrix.h"

/*
Matrix::Matrix(const string inF):rowIDs(0),colIDs(0),sampleNameSep("")
{
	//reads matrix from HDD
	string line;
	ifstream in(inF.c_str());
	int ini_ColPerRow(0),cnt(0);

	//check MAP format
	while(getline(in,line,'\n')) {
		if(line.substr(0,1) == "#" || line.length()<2){continue;}
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss,segments,'\t')) {
			ColsPerRow++;
		}

		if (cnt==0){
			ini_ColPerRow = ColsPerRow;
		} else {
			if (ColsPerRow != ini_ColPerRow){

#ifdef notRpackage
cerr<<"Number of columns on line "<<cnt<<" is "<<ColsPerRow<<". Expected "<<ini_ColPerRow<<" columns.\n";
				std::exit(6);
			}
		}
		cnt++;
	}
	//vector<mat_fl> ini_vec(ini_ColPerRow-1,0.f);
	mat.resize(cnt-1,vector<mat_fl>(ini_ColPerRow-1,0.f));
	colIDs.resize(ini_ColPerRow,"");
	rowIDs.resize(cnt-1,"");
	int lineCnt= cnt;
	//reset input stream
	in.clear();
	in.seekg(0, ios::beg);
	cnt=-2;
	string segments;

	while(getline(in,line,'\n')) {
		cnt++;
		if(line.substr(0,1) == "#"){continue;}
		if (line.length()<10){continue;}
		stringstream ss;
		ss << line;
		int cnt2(-2);
		if (cnt==-1){//read header
			cnt2++;
			while (getline(ss,segments,'\t')) {
				cnt2++;
				colIDs[cnt2] = segments;
			}
			continue;
		}
		while (getline(ss,segments,'\t')) {
			cnt2++;
			if (cnt2==-1){
				rowIDs[cnt] = segments;
				continue;
			}
			mat[cnt][cnt2] = (mat_fl)atof(segments.c_str());
		}

	}
	in.close();

}
*/

inline mat_fl median(std::vector<mat_fl> vec, bool ignoreZeros)
{
	if (vec.size() == 0) { return (mat_fl)0; }
	//std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, vec.end());
	sort(vec.begin(), vec.end());
	size_t i(0);
	if (ignoreZeros) {
		for (; i < vec.size(); i++) {
			if (vec[i] > 0) {
				break;
			}
		}
		if (vec.size() == i) { return (mat_fl)0; }
	}
	size_t size = vec.size() - i;

	if (size % 2 == 0) {
		return (vec[(size / 2) - 1 + i] + vec[size / 2 + i]) / 2;
	}
	return  vec[size / 2 + i];
}
void vecPurge(vector<vector<mat_fl>>& vec, mat_fl val) {
	for (size_t i = 0; i < vec.size(); i++) {
		for (size_t j = 0; j < vec[i].size(); j++) {
			vec[i][j] -= val;
		}
	}
}
string join(const vector<string>& in, const string &X) {
	string ret(in[0]); for (size_t i = 1; i < in.size(); i++) { ret += X + in[i]; } return ret;
}
string join(const vector<int>& in, const string &X) {
	string ret(itos(in[0])); for (size_t i = 1; i < in.size(); i++) { ret += X + itos(in[i]); } return ret;
}
bool writeChunk(ostream* o, string  s) { (*o) << s; return true; }




//*********************************************************
ModStep::ModStep(const string & s, bool & recMod, vector<string>& subMod) :alternates(0), redundancy(0) {
	istringstream ss(s);
	string token(""), tok2("");
	//std::istream_iterator<std::string> beg(ss), end;
	//std::vector<std::string> tokens(beg, end); // done!
	//for (auto& s : tokens) { std::cout << s; }
	while (std::getline(ss, token, '\t')) {
		stringstream buff(token);
		vector<string> tmp(0);
		while (std::getline(buff, tok2, ',')) {
			tmp.push_back(tok2);
			if (tok2[0] == 'M') {
				recMod = true;
				subMod.push_back(tok2);
			}
		}
		alternates.push_back(tmp);
	}
}

//how often does the KO used occur in total dataset?
void ModStep::setRedund(ModOccur& m) {
	redundancy.resize(alternates.size());
	for (size_t i = 0; i < alternates.size(); i++) {
		vector<int> tmp(alternates[i].size(), 0);
		for (size_t j = 0; j < alternates[i].size(); j++) {
			tmp[j] = m[alternates[i][j]];
		}
		redundancy[i] = tmp;
	}
}
void ModStep::getAllKOs(list<string>& ret) {

	for (size_t i = 0; i < alternates.size(); i++) {
		for (size_t j = 0; j < alternates[i].size(); j++) {
			ret.push_back(alternates[i][j]);
		}
	}
}
void ModStep::abundParts(const vector<mat_fl>& v, const unordered_map<string, int>& IDX,
	vector<mat_fl>& abund, vector<bool>& active, vector<string>& KOdescr,
	float hitComplRatio, int redund) {
	//some params, should be fine tuned if possible
	//float hitComplRatio(0.8f);

	active.resize(alternates.size(), false);
	abund.resize(alternates.size(), (mat_fl)0);
	KOdescr.resize(alternates.size(), "");
	//actual deep routine to determine if KOs in this step satisfy presence conditions
	for (size_t i = 0; i < alternates.size(); i++) {
		size_t altS = alternates[i].size(); float hits(0);
		vector<mat_fl> tmpAB(altS, (mat_fl)0);
		for (size_t j = 0; j < altS; j++) {
			//redundant???
			if (redundancy[i][j] > redund) {
				altS--; //reduce size of this set

				continue;
			}
			//check for actual abundance in matrix-vector subpart
			auto fn = IDX.find(alternates[i][j]);
			if (fn == IDX.end()) {
				tmpAB[j] = 0;
			}
			else {
				tmpAB[j] = v[fn->second];
				if (tmpAB[j] > 0) {
					hits++;
					KOdescr[i] += alternates[i][j] + ",";
				}
			}
		}
		if (altS == 0) {
			//r[i] = (mat_fl)-1;//signal that removed due to redundancy
			continue;
		}
		/*		if (hits > 0) {
		int x = 0;
		}*/
		if (hits / (float)altS >= hitComplRatio) {
			abund[i] = median(tmpAB);
			active[i] = true;
		}
	}

	//return (active);
}



//*********************************************************
Module::Module(vector<string>& n) :name(""), description(""), steps(0), 
		submods(0),containsMods(false), usedInOtherMods(false){
	string token("");
	for (size_t i = 0; i < n.size(); i++) {
		if (i == 0) {//module name & description
			istringstream ss(n[i]);
			std::getline(ss, token, '\t'); name = token;
			std::getline(ss, token, '\t'); description = token;
		} else { //actual modules components (e.g. KOs)
			steps.push_back( ModStep(n[i], containsMods, submods));
		}
	}
}
mat_fl Module::pathAbundance(const vector<mat_fl>& v, const unordered_map<string, int>& IDX,
	const int redund, const float PathwCompl, const float enzymCompl, string & savedNmsKO, float& modScoreOut) {
	//initial parameters
	//float PathwCompl(0.6f); //corresponds to -c 
	//float enzymCompl(0.8f); 
	//int redund(0);

	vector< vector< mat_fl >> abunds(steps.size(), vector<mat_fl>(0));//contains abundance
	vector< vector< bool >> active(steps.size());//contains info if path was even active
	vector<mat_fl> preMed(steps.size(), (mat_fl)0), postMed(steps.size(), (mat_fl)0);
	vector<vector<string>> altKOs(steps.size(), vector<string>(0));//just for saving which KO's were exactly active

	//auto t = IDX.find("xx");
	for (size_t i = 0; i < steps.size(); i++) {
		steps[i].abundParts(v, IDX, abunds[i], active[i], altKOs[i], enzymCompl, redund);
		//determine median overall value
		preMed[i] = median(abunds[i]);
	}
	mat_fl pm = median(preMed);
	mat_fl retval(0);

	if (0) {
		//VAR 1
		//select one median value per step only, for the final pathway median abundance
		for (size_t i = 0; i < steps.size(); i++) {
			vector<mat_fl> tmp(abunds[i].size(), (mat_fl)0);
			for (size_t j = 0; j < abunds[i].size(); j++) {
				tmp[j] = abunds[i][j] - pm;
				exit(99);
			}
		}
	}
	else {
		//VAR 2
		//calc abundance for each possible (active) path and substract from others, but add remainders up
		//vector<mat_fl> tmp(steps.size(), (mat_fl)0);
		vector<mat_fl> curP(steps.size(), (mat_fl)0);
		float act(0.f), shldAct((float)steps.size());
		vector<int> decIdx(steps.size(), 0);
		//ini decIdx
		for (size_t i = 0; i < steps.size(); i++) {
			size_t dI = 0; double maxAB = 0;
			while (dI < uint(active[i].size())) {
				if (abunds[i][dI] > maxAB && active[i][dI]) {
					decIdx[i] = dI;
					maxAB = abunds[i][decIdx[i]];
				}
				dI++;
			}
		}
		bool saveKOnames(true);// save names of KOs used in extra file?

		while (1) { // this loop goes over every possible path combination
			for (size_t i = 0; i < steps.size(); i++) {
				if (abunds[i][decIdx[i]] > 0 && active[i][decIdx[i]]) {
					curP[i] = abunds[i][decIdx[i]]; act++;
					if (saveKOnames) {
						savedNmsKO += altKOs[i][decIdx[i]];// +",";
					}
				} /*else if (!active[i][decIdx[i]]) {
				  shldAct--;
				  } */
			}
			if ((act / shldAct) >= PathwCompl) { //shldAct > 0 &&
												 //this part parses out the KOs that are actually active
				mat_fl curM = median(curP, true);
				retval += curM;
				vecPurge(abunds, curM);
				modScoreOut = act / shldAct;
			}
			else {
				savedNmsKO = "";
				modScoreOut = 0.f;
			}
			break;
		}

	}
	return retval;
}



//*********************************************************
//read in module file
Modules::Modules(const string& inF, vector<string> cns) :
	Matrix(),
	moduleDescriptions(0), redundantUsedMods(0),
	recurrentMods(0), ModUsed(0),redund(1), PathwCompl(0.6f), enzymCompl(0.8f)
{
	//ini matrix base class members
	colIDs = cns; 
	maxCols = ((int)colIDs.size());


	ifstream is(inF.c_str());
	string line(""); vector<string> buffer(0);
	string ModToken = "M";
	while (safeGetline(is, line)) {
		//comment
		if (line[0] == '#' || line.substr(0,7)=="Mod Des") { continue; }
		//new module opens, create old module
		if (buffer.size()>0 && line.find(ModToken) == 0 && line.find("\t") != string::npos) {
			if (buffer.size() > 0) {
				mods.push_back(Module(buffer));
			}
			buffer.resize(0);
		}
		if (line[0] == ' ') {
			line.erase(0, 1);
		}
		if (line.size() > 3) {
			buffer.push_back(line);
		}
	}
	//create last module
	mods.push_back(Module(buffer));

	buffer.resize(0);
	is.close();
	//redundancy of KOs
	calc_redund();

	//set up names
	rowIDs.resize(mods.size(), "");
	for (size_t i = 0; i < mods.size(); i++) {
		rowIDs[i] = mods[i].name;
		//and track position
		if (ModPos.find(rowIDs[i]) == ModPos.end()) {
			ModPos[rowIDs[i]] = vector<int>(1, i);
		}
		else {
			ModPos[rowIDs[i]].push_back(i);
		}
	}
	//finished off all ModPos entries..
	for (size_t i = 0; i < mods.size(); i++) {
		if (mods[i].containsMods) {
			recurrentMods.push_back(i);
			for (size_t j = 0; j < mods[i].submods.size(); j++) {
				bool doAddThMod(true);
				for (size_t k = 0; k < redundantUsedMods.size(); k++) {
					if (redundantUsedMods[k] == mods[i].submods[j]) {
						doAddThMod = false; break;
					}
				}
				if (doAddThMod) {
					redundantUsedMods.push_back(mods[i].submods[j]);
					for (size_t kk = 0; kk < ModPos[mods[i].submods[j]].size(); kk++){
						mods[ModPos[mods[i].submods[j]][kk]].usedInOtherMods = true;
					}
				}
			}
		}
	}
	//set up descriptions
	moduleDescriptions.resize(mods.size(), "");
	for (size_t i = 0; i < mods.size(); i++) {
		moduleDescriptions[i] = mods[i].description;
	}
    #ifdef notRpackage
	std::cout << "Read " << mods.size() << " modules\n";
    #endif
	//ini base class
	//Matrix(this->modNms(), colIDs);
	this->ini_mat();


}

void Modules::addDescription(const string& inF) {
	if (inF == "") { return; }
	ifstream in(inF.c_str());
	if (!in) { 
		cerr2 ("Couldn't open module description infile " +inF +"\n",329); 
	}
	string line("");
	moduleDescriptions.resize(mods.size());
	while (safeGetline(in, line)) {
		stringstream ss;
		ss << line;
		std::map<std::string, vector<string>>::iterator fnd;
		string modN(""), modD("");
		getline(ss, modN, '\t');
		if (modN == "Mod" || modN == "") { continue; }
		getline(ss, modD, '\t');
		auto x = ModPos.find(modN);
		if (x == ModPos.end()) {
			cerr2 ("Couldn't find module " +modN +"\n",743);
		}	else {
			moduleDescriptions[ x->second[0] ] = modD;
		}

	}
	in.close();
}
void Modules::addHierachy(const string& inF){
	if (inF == "") { return; }
	ifstream in(inF.c_str());
	if (!in) { 
		cerr2 ("Couldn't open module hierachy infile " +inF +"\n",329); 
	}
	string line("");
	hierachy.resize(mods.size(),vector<string>(0));
	int cnt (-1); int modDef(-1);
	while (safeGetline(in, line)) {
		stringstream ss;
		ss << line;
		cnt++;
		string modN(""), modD("");
		if (cnt == 0) {
			int yy = 0;
			while (getline(ss, modD, '\t')) {
				if (modD == "Mod") {
					modDef = yy; break;
				}
				yy++;
			}
			if (modDef == -1) {
				cerr2 ("Error: couldn not find \"Mod\" tag in module hierachy file\n");
			}
			continue;
		}
		vector<string> tmp;
		int yy = 0;
		while (getline(ss, modD, '\t')) {			
			if (yy == modDef) { modN = modD;
			}	else {
				tmp.push_back(modD);
			}
			yy++;
		}
		if (modN == "") { continue; }

		auto x = ModPos.find(modN);
		if (x == ModPos.end()) {
			cerr2 ("Couldn't find module " +modN +"\n",743);
		}

		hierachy[x->second[0]] = tmp;
	}
	in.close();

}

void Modules::calc_redund() {
	list<string> fL; //full list of KOs
	for (size_t i = 0; i < mods.size(); i++) {
		mods[i].getAllKOs(fL);
	}
	MO.clear();
	//calc abundance
	for (auto& s : fL) {
		auto fnd = MO.find(s);
		if (fnd == MO.end()) {
			MO[s] = 1;
		}
		else {
			fnd->second++;
		}
	}
	//stats on KO redundancy (print out only)
	vector<int>statKOr = vector<int>(0, 0);
	int maxRed = 0;
	for (auto kor : MO) {
		if (kor.second > maxRed) {
			statKOr.resize((kor.second + 1), 0);
			maxRed = kor.second;
		}
		statKOr[kor.second] ++;
	}
	//std::cout << "stats on DB KO redundancy (redundancy : occurence):\n";
	for (size_t i = 0; i < statKOr.size(); i++) {
		if (statKOr[i] < 1) { continue; }
		//std::cout << i << " : " << statKOr[i] << endl;
	}
	for (size_t i = 0; i < mods.size(); i++) {
		mods[i].setReddundancy(MO);
	}

}
void Modules::writeMatrix(const string of, bool onlyFilled, bool collapseDblFeats) {
	ofstream out;
	out.open(of.c_str(), ios_base::out);
	out.precision(8); out << "Gene";
	for (size_t smpl = 0; smpl < (colIDs.size()); smpl++) {
		out << "\t" << colIDs[smpl];
	}
	out << endl;
	unordered_map<string, int> ModCnt;
	vector<mat_fl> rowSums;
	uint writeCnt(0);
	size_t cidS(colIDs.size());
	ModUsed.resize(rowIDs.size(), false);
	//if (onlyFilled) { 
	rowSums = getRowSums();
	for (size_t i = 0; i < rowIDs.size(); i++) {
		if (onlyFilled && rowSums[i] == 0) {
			continue;
		}
		if (rowSums[i] > 0) {
			ModUsed[i] = true;
		}
		writeCnt++;
		vector<mat_fl> wrVec(cidS, 0.f);
		if (collapseDblFeats && ModPos[rowIDs[i]].size() > 1) {//this is collapseable
			if (ModCnt.find(rowIDs[i]) != ModCnt.end()) {
				continue;
			}
			ModCnt[rowIDs[i]] = 1;
			//go over all alt representations of this mod and sum them up
			for (size_t x = 0; x < ModPos[rowIDs[i]].size(); x++) {
				int jj = ModPos[rowIDs[i]][x];
				for (size_t smpl = 0; smpl < cidS; smpl++) {
					wrVec[smpl] += mat[smpl][jj];
				}
			}

		}
		else {
			for (size_t smpl = 0; smpl < cidS; smpl++) {
				wrVec[smpl] = mat[smpl][i];
			}

		}
		out << rowIDs[i];
		for (size_t smpl = 0; smpl < cidS; smpl++) {
			out << "\t" << wrVec[smpl];
		}

		out << endl;
	}
	out.close();
	#ifdef notRpackage
	std::cout << "Wrote " << writeCnt << " modules in final matrix\n";
	#endif
}

vector<string> Modules::modNms_numbered() {
	vector<string> out = rowIDs;
	for (auto it = ModPos.begin(); it != ModPos.end(); ++it) {
		vector<int> pps = it->second;
		for (size_t i = 0; i < pps.size(); i++) {
			out[pps[i]] += "." + itos(i);
		}
	}
	return out;
}

void Modules::writeModDescr(const string& nos, bool onlyUsed) {
	if (onlyUsed && ModUsed.size() == 0) { return; }
	ofstream of(nos.c_str());
	bool useHie = false;
	if (hierachy.size() == moduleDescriptions.size()) { useHie = true; }
	unordered_map<string, int> usedMN;
	for (size_t i = 0; i < rowIDs.size(); i++) {
		if (usedMN.find(rowIDs[i]) != usedMN.end()) { continue; }
		if (onlyUsed && !ModUsed[i]) { continue; }
		of << rowIDs[i] ;
		if (useHie) {
			for (uint u = 0; u < hierachy[i].size(); u++){
				of << "\t" << hierachy[i][u];
			}
		}
		of << "\t" << moduleDescriptions[i];
		of << endl;
		usedMN[rowIDs[i]] = 1;
	}
	of.close();
}

void Modules::calcModAbund( vector<mat_fl>& v, const int pos, const unordered_map<string,
	int>& IDX, vector<string> &retStr, vector<float> &retScore) {
	vector<mat_fl> ret(mods.size(), (mat_fl)0);
	retStr.resize(mods.size(), ""); retScore.resize(mods.size(), 0.f);
	//mat_fl unass_cnt(0.f);//TOGO
						  //vector<bool> usedKOs(v.size()); //initial idea to save KOs used to estimated unassigned fraction - better to scale by seq depth external
	for (size_t i = 0; i < mods.size(); i++) {
		//if (mods[i].name == "M00022") {			int x = 0;		} /**/
		if (mods[i].containsMods) {
			if (mods[i].usedInOtherMods) {
    			// everyone knows, this is not the best fix, but what can we do then?
    			// we need to not have exit commands in the code for R
    			// termination should happen more pretty somehow
				cerr2("usedInOtherMods && containsMods - fatal error\n",823 );
			}
			continue; 
		}//these need to be estimated at the end
		ret[i] = mods[i].pathAbundance(v, IDX, redund, PathwCompl, enzymCompl, retStr[i], retScore[i]);
		if (mods[i].usedInOtherMods) { // add abundance to vec
			auto fnd = IDX.find(mods[i].name);
			if (fnd == IDX.end()) {
				cerr2 ("Could not find module " +mods[i].name +" but should be in AB matrix\n", 487); 
			}
			v[fnd->second] = ret[i];
		}
	}

	//estimate modules that contain other modules
	for (size_t j = 0; j < recurrentMods.size(); j++) {
		uint i = recurrentMods[j];
		ret[i] = mods[i].pathAbundance(v, IDX, redund, PathwCompl, enzymCompl, retStr[i], retScore[i]);
	}


	//add unkowns
	//ret.push_back(unass_cnt);
	//return ret;
	this->addTtlSmpl(ret, pos);
}


//*********************************************************
Matrix::Matrix(void)
	:rowIDs(0), colIDs(0), maxCols(0), HI(0), maxLvl(0), sampleNameSep(""), 
	doSubsets(false), doHigh(false),CoEx_minp(0),CoEx_thresh(0)
{
}

void Matrix::readColNms(istream* in,const char sep) {
	string segments; string line;
	safeGetline((*in), line);
	//getline(in, line, '\n');
	while (line.substr(0, 1) == "#") {
		safeGetline((*in), line);
	}
	stringstream sso;
	int cnt2(-2);
	sso << line;
	while (getline(sso, segments, sep)) {
		cnt2++;
		if (segments.length() > 150) {
			cerr2 (segments +" error!\n",5);
		}
		if (cnt2 == -1) { continue; }
		colIDs[cnt2] = segments;
	}
}

int Matrix::iniCols(istream* in, const char sep) {
	int ini_ColPerRow = 0;
	int cnt = 0;
	string line;
	while (safeGetline((*in), line)) {
//	while (getline((*in), line, '\n')) {
		if (line.substr(0, 1) == "#" || line.length() < 2) { continue; }
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss, segments, sep)) {
			ColsPerRow++;
		}

		if (cnt == 0) {
			ini_ColPerRow = ColsPerRow;
		}
		else {
			if (ColsPerRow != ini_ColPerRow) {
				cerr2 ("C1: Number of columns on line " + itos(cnt) +" is " + itos(ColsPerRow) +". Expected " + itos(ini_ColPerRow) + " columns.\n" + line +"\n",6);
			}
		}
		cnt++;
		if (cnt > 10) { break; }
	}
	if (ini_ColPerRow == 0) {
		cerr2 ("Could not find valid columns in matrix.. exiting\n",432);
	}
	//reset input stream
	(*in).clear();
	(*in).seekg(0, ios::beg);
	colIDs.resize(ini_ColPerRow - 1, "");
	colSum.resize(ini_ColPerRow - 1, 0.0);

	return ini_ColPerRow;
}
Matrix::Matrix(const string inF, const string outF, vector<double> colsums, vector<string> colNms,
	const char sep, options* opts)
	: rowIDs(0), colIDs(0), maxCols(0), HI(0), maxLvl(0), sampleNameSep(""), 
	doSubsets(false), doHigh(false), sparse(false),loadPrimaryMat(false)
{
	//reads matrix from HDD
	//and writes it simultaneously to single files
	colSum = colsums;

	
	istream* in;
	if (isGZfile(inF)) {
#ifdef _gzipread
		in = new igzstream(inF.c_str(), ios::in); cout << "Straming gzip input on the fly\n";
#else
		cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
	}
	else { in = new ifstream(inF.c_str()); }

	if (!(*in)) {cerr2 ("Cant open file " +inF +"\n", 11);}


	ostream * out;
	if (isGZfile(outF)) {
#ifdef _gzipread
		out = new ogzstream(outF.c_str(), ios::out);
		cout << "Writing gzip'd matrix " << outF << endl;
#else
		cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
	}
	else { out = new ofstream(outF); }
	//(*out).precision(9);


	if (!(*out)) {cerr2 ("Can't open out file " +outF +"\n", 11);}
	int ini_ColPerRow = iniCols(in,sep);
	//reset istream..
	delete in;
	if (isGZfile(inF)) {
#ifdef _gzipread
		in = new igzstream(inF.c_str(), ios::in); 
#else
		cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
	}
	else { in = new ifstream(inF.c_str()); }
	readColNms(in,sep);

	//check that colNames and colSums are in same order..
	for (uint i = 0; i < colIDs.size(); i++) {
		if (colNms[i] != colIDs[i]) {
			cout << "Unequal order!\n";
			exit(339);
		}
	}


	int cnt(-1), geneCnt(0);
	string SEP = "\t";
	stringstream ss;


	//create header of normalized matrix..
	*out << "Norm_rtk";
	for (uint i = 0; i < colIDs.size(); i++) {
		*out << SEP << colIDs[i];
	}
	*out << endl;

	//heavy reading


	int num_threads = opts->threads;
	vector < job4 > slots(num_threads);
	int j(0), cntNA(0);
	string line("");// = safeGetline2(in);
	//thread for ostream
	future<bool> wrStr; string wrTmp("");
	wrStr = async(std::launch::async, writeChunk, out, wrTmp);
	int rowCnt(0);

	while (true) {
		if (j >= num_threads) { j = 0; }
		if (slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(1)) == std::future_status::ready) {
			slots[j].inUse = false;
			wrTmp += slots[j].fut.get();
			rowCnt++;
		}
		if (slots[j].inUse == false) {
			//string line = readJ.get();readJ = async(std::launch::async, safeGetline2, in);
			if (line.substr(0, 1) != "#" && line.length() >= 5) {
				slots[j].fut = async(std::launch::async, &Matrix::processLine2line, this, line, cntNA, geneCnt, ini_ColPerRow, sep);
				//wrTmp += processLine2line( line,cntNA, geneCnt, ini_ColPerRow, sep);
				geneCnt++;
				slots[j].inUse = true;
				j++;
			}
			if (!(*in)) { break; }
			line = safeGetline2(in);
			if (wrTmp.length() > 1e5) {//outstream handling..
				bool x = wrStr.get();		
				wrStr = async(std::launch::async, writeChunk, out, wrTmp);
				//writeChunk(out, wrTmp);
				wrTmp = "";
			}

		}
	}

	//last jobs to collect..
	for (int j = 0; j < num_threads; j++) {
		if (slots[j].inUse == true) {
			slots[j].inUse = false;
			wrTmp += slots[j].fut.get();
			rowCnt++;
		}
	}
	bool x = wrStr.get();
	(*out) << wrTmp;

	/*
	string rowID,segments;
	while (safeGetline((*in), line)) {
		//while (getline(in, line, '\n')) {
		cnt++;
		if (line.substr(0, 1) == "#") { *out << line; continue; }
		if (line.length()<3) { continue; }
		int cnt2(-2);
		stringstream ss;
		ss << line;
		bool breaker(false);
		while (getline(ss, segments, sep)) {
			cnt2++;
			if (cnt2 == -1) {
				rowID = segments;

				//maybe add functionality later: only normalize subset
				if (doSubsets && subset.find(rowID) == subset.end()) {
					breaker = true;
					break;
				}
				
				geneCnt++;
				*out << rowID;
				continue;
			}
			mat_fl tmp = atof(segments.c_str());
			colSum[cnt2] += (double)tmp;
			if (tmp == 0) {
				*out << SEP << "0";
			}
			else {
				*out << SEP << (tmp / colsums[cnt2]);
			}
		}
		if (breaker) {
			continue;
		}
		if (cnt2 + 2 != ini_ColPerRow) {
			cerr2 ("C2: Number of columns on line " + itos( cnt ) + " is " +(itos(cnt2 + 2)) +". Expected " + itos(ini_ColPerRow) +" columns.\n", 62);
		}
		*out << endl;

	}
	*/
	cout << "Normalized matrix, writing " << rowCnt << " rows and " << colNms.size() << " columns\n";
	for (uint i = 0; i < colSum.size(); i++) {
		if (colsums[i] != colSum[i]) {
			cout << "Unequal colSum!\n";
			exit(339);
		}
	}



	delete in; //in.close();
	delete out; // out.close();
}



Matrix::Matrix(options* opts, const string xtra, vector<string>& outFName,
			bool highLvl, bool NumericRowId, bool writeTmpFiles)
	: rowIDs(0), colIDs(0), maxCols(0), HI(0), maxLvl(0), sampleNameSep(""), 
	doSubsets(false), doHigh(highLvl), colSum(0), rowOcc(0), sparse(true), loadPrimaryMat(false)
{
	//this implementation doesn't load the matrix into memory, just creates searches and sum stats
	const string inF = opts->input;
	const string outF = opts->output;
	const char sep = opts->sepChar;

	//bool medianSum = opts->median; 
	bool meanSum = opts->mean;
	//number of hits for a category (by different genes) per samples or overall
	uint occPerSmpl = opts->occPerSmpl;
	uint occMin = opts->occMin;
	//reads matrix from HDD
	//and writes it simultaneously to single files
	if (doHigh){
		read_hierachy(xtra,opts->extendHierachy);
	} else if (xtra.length() > 2){
		read_subset_genes(xtra);
	}

	istream* in;
	if (isGZfile(inF)) {
#ifdef _gzipread
		in = new igzstream(inF.c_str(), ios::in); cout << "Reading gzip input "<< inF<<endl;
#else
		cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
	}
	else { in = new ifstream(inF.c_str()); }
//ifstream in(inF.c_str());
	if (!(*in)){
		cerr2 ("Cant open file " + inF + "\n",11);
	}
	int ini_ColPerRow = iniCols(in,sep);
	delete in;
	if (isGZfile(inF)) {
#ifdef _gzipread
		in = new igzstream(inF.c_str(), ios::in); 
#else
		cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
	}
	else { in = new ifstream(inF.c_str()); }

	int cnt(0);
	cnt=-1;
	readColNms(in,sep);
	//set up levels to calc sum stat for
	vector<ofstream> outFs(ini_ColPerRow-1);
	vector<string> outStr(ini_ColPerRow-1);
	if (doHigh){
		for (int i = 0; i < maxLvl; i++){
			//set up empty higher level hierachies
			HI.push_back(new HMat(LvlNms[i], colIDs, vector<string>(0)));
		}
	}
	//set up tmp empty files
	if (!doHigh && writeTmpFiles) {
		for (uint i = 0; i < colIDs.size(); i++) {
			string oF2 = outF + sampleNameSep + colIDs[i];
			outFName.push_back(oF2);
			outFs[i].open(oF2.c_str(), ios_base::out);
			outFs[i].precision(12);
			outFs[i].close();
		}
	}
	int geneCnt(0);
	int cntNA(0);

	//future<string> readJ;
	//readJ = async(std::launch::async, safeGetline2, in);
	int num_threads = opts->threads;
	vector < job3 > slots(num_threads);
	int j(0);
	string line = safeGetline2(in);
	//while (getline_n((*in),line,100000)){
	//while (safeGetline((*in), line)) {
	while (true) {
		cnt++;
		if (j >= num_threads) { j = 0; }
		if (slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
			slots[j].inUse = false;
			rowRes* ret = slots[j].fut.get();
			if (ret->NAcnt) {cntNA++;}
			addRow(ret);
			rowOcc.push_back(ret->Occ);

			if (writeTmpFiles) {vector<mat_fl>& rowV = ret->rowV;//write to tmp File
				for (size_t i = 0; i < rowV.size();i++) {
					if (rowV[i] > 0) {
						// if id is numeric (number of row) or the actual id as string
						if (NumericRowId == true) {outStr[i] += std::to_string(ret->curRowIdx) + "\t" + to_string(rowV[i]) + "\n";
						} else { outStr[i] += ret->rowID + "\t" + to_string(rowV[i]) + "\n";
						}
					}
				}
			}
			delete ret;
		}
		if (slots[j].inUse == false) {
			if (line.substr(0, 1) != "#" &&  line.length() >= 5) {
				slots[j].fut = async(std::launch::async, &Matrix::processLine2Mat1, this, line,
					cntNA, geneCnt, ini_ColPerRow, sep);
				geneCnt++;
				slots[j].inUse = true;
			}
			//string line = readJ.get();
			//readJ = async(std::launch::async, safeGetline2, in);
			if (!(*in)) { break; }
			line = safeGetline2(in);
		}

		if (cnt % 1000 == 0 && writeTmpFiles){
			// every 1000 lines, write to file. The rest will be written later
			for (size_t cnt2=0;cnt2<outStr.size();cnt2++){
				// we open the file, and close it again, as we dont want to be limited in the number of
				// files to have open
				outFs[cnt2].open(outFName[cnt2], ios_base::out | ios::app);
				outFs[cnt2].precision(12);
				outFs[cnt2] << outStr[cnt2];outStr[cnt2] = "";
				outFs[cnt2].close();
			}
		}
		j++;

	}
	delete in;//in.close();


	for (int j = 0; j < num_threads; j++) {
		if (slots[j].inUse == true) {
			rowRes* ret = slots[j].fut.get();
			slots[j].inUse = false;

			addRow(ret);
			rowOcc.push_back(ret->Occ);

			if (writeTmpFiles) {
				vector<mat_fl>& rowV = ret->rowV;//write to tmp File
				for (size_t i = 0; i < rowV.size(); i++) {
					if (rowV[i] > 0) {
						// if id is numeric (number of row) or the actual id as string
						if (NumericRowId == true) {
							outStr[i] += std::to_string(ret->curRowIdx) + "\t" + to_string(rowV[i]) + "\n";
						}
						else {
							outStr[i] += ret->rowID + "\t" + to_string(rowV[i]) + "\n";
						}
					}
				}
			}
			delete ret;

		}
	}

	//output of matrices..
	ofstream out;
	if (doHigh ){//write out high lvl mats
		for (int i = 0; i < maxLvl; i++){
			HI[i]->prep();
			string oF2 = outF + LvlNms[i] + ".txt";
			out.open(oF2.c_str(), ios_base::out);
			HI[i]->print(out,meanSum, occMin, occPerSmpl);
			out.close();
		}
	}
	else if(writeTmpFiles){//close filestreams to single sample files
		// write the overlapp of the read in (since the last 1000)
		for (size_t i = 0; i < outFs.size(); i++){
			// we open the file, and close it again, as we dont want to be limited in the number of
			// files to have open
			outFs[i].open(outFName[i], ios_base::out | ios::app);
			outFs[i].precision(12);
			outFs[i] << outStr[i];
			outFs[i].close();
		}
	}
	//write colSums
	/*
	string oF2 = outF + sampleNameSep + "sums.txt";

	out.open(oF2.c_str(),ios_base::out);
	out.precision(12);
	for (size_t smpl=0;smpl<(colIDs.size()); smpl++){
		out<<colIDs[smpl]<<"\t"<<colSum[smpl]<<endl;
	}
	out.close();*/
	#ifdef notRpackage
	std::cout << "Read " << geneCnt << " rows and " << ini_ColPerRow << " columns\n";
	#endif
}

Matrix::Matrix(options* opts, const string xtra, bool highLvl)
	: rowIDs(0), colIDs(0), maxCols(0), HI(0), maxLvl(0), sampleNameSep(""), 
	doSubsets(false), doHigh(highLvl), sparse(opts->sparse), loadPrimaryMat(true)
{
	const string inF = opts->input;
	const char sep = opts->sepChar;
	//reads matrix from HDD
	//into mem.. careful with big files!
	if (doHigh){
		read_hierachy(xtra,opts->extendHierachy);
	}
	else if (xtra.length() > 2){
		read_subset_genes(xtra);
	}
//	ifstream in(inF.c_str());
	istream* in;
	if (isGZfile(inF)) {
#ifdef _gzipread
		in = new igzstream(inF.c_str(), ios::in);
		cout << "Reading gzip input\n";
#else
		cout << "gzip not supported in your rtk build\n"; exit(50);
#endif

	}
	else {
		in = new ifstream(inF.c_str());
	}
	if (!(*in)){
		cerr2 ("Cant open file " +inF +"\n",11);
	}
	int ini_ColPerRow = iniCols(in,sep);
	int cnt(0);
	int cntNA(0);
	//cout << "cols: " << ini_ColPerRow << endl;

	if (ini_ColPerRow == 0) {

		cerr2 ("Empty matrix provided\n");
		return;
	}
	cnt = -1;

	delete in;
	if (isGZfile(inF)) {
#ifdef _gzipread
		in = new igzstream(inF.c_str(), ios::in);
#else
		cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
	}
	else { in = new ifstream(inF.c_str()); }
	
	readColNms(in,sep);
	if (doHigh){
		for (int i = 0; i < maxLvl; i++){
			HI.push_back(new HMat(LvlNms[i], colIDs, vector<string>(0)));
		}
	}
	int geneCnt(0);
	//vector<mat_fl> emptyVec(ini_ColPerRow, (mat_fl)0);
	if (!sparse) {
		mat.resize(ini_ColPerRow - 1, vector<mat_fl>(0));
	}
	else {
		matSp.resize(ini_ColPerRow - 1);
	}
	//future<bool> job = async(std::launch::async, safeGetline, ref((*in)), ref(line));
	//bool usedThred(false);
	//future<string> readJ;
	//readJ = async(std::launch::async, safeGetline2, in);

	int num_threads = opts->threads;
	vector < job3 > slots(num_threads);
	int j(0); 
	string line = safeGetline2(in);

	while (true) {
		if (j >= num_threads) { j = 0; }
		if (slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(1)) == std::future_status::ready) {
			slots[j].inUse = false;
			rowRes* ret = slots[j].fut.get();
			if (ret->NAcnt) { cntNA++; }
			addRow(ret);
			delete ret;
		}
		if (slots[j].inUse == false) {
			//string line = readJ.get();readJ = async(std::launch::async, safeGetline2, in);
			if (line.substr(0, 1) != "#" && line.length() >= 5) {
				slots[j].fut = async(std::launch::async, &Matrix::processLine2Mat1, this, line,
					cntNA, geneCnt, ini_ColPerRow, sep);
				geneCnt++;
				slots[j].inUse = true;
			}
			if (!(*in)) { break; }
			line = safeGetline2(in);
		}
		j++;
	}

	//last jobs to collect..
	for (int j = 0; j < num_threads; j++) {
		if (slots[j].inUse == true) {
			rowRes* ret = slots[j].fut.get();
			addRow(ret);
			slots[j].inUse = false;
			delete ret;
		}
	}

	//in.close();
	delete in;
	if (sparse) {
		maxCols = (int)matSp.size();
		std::cout << "Sparse Matrix: Read " << maxCols << " colmuns and " << cnt << " rows\n";
	} else {
		maxCols = (int)mat.size();
#ifdef notRpackage
		std::cout << "Read " << maxCols << " colmuns and " << mat[0].size() << " rows\n";
#endif
	}


	if (geneCnt == 0) {
		cerr2("No genes read.. aborting\n", 73);
	}
}
Matrix::Matrix(const vector<string>& rnms, const vector<string>& cnms) :
	rowIDs(rnms), colIDs(cnms), maxCols((int)cnms.size()), loadPrimaryMat(true)
{
	ini_mat();
}
void Matrix::ini_mat() {
	if (maxCols != (int)colIDs.size()) {
		maxCols = (int)colIDs.size();
	}
	vector<mat_fl> iniV = vector<mat_fl>(rowIDs.size(), (mat_fl)0);
	mat.resize(maxCols, iniV);
}


string Matrix::processLine2line(string  line, int cntNA,
	int geneCnt, int ini_ColPerRow, const char sep) 
{

	rowRes* dres = this->processLine2Mat1( line,  cntNA, geneCnt,  ini_ColPerRow,  sep);
	 string str = prepRow(dres, string(1,sep), true);
	 delete dres;
	 return str;

}

rowRes* Matrix::processLine2Mat1(string  line, int cntNA,
	int geneCnt, int ini_ColPerRow, const char sep) {
	string segments;
	string rowID = "";
	stringstream ss;
	ss << line;
	vector<string> taxa(0);
	bool breaker(false);
	int cnt2(-2);
	//tmp objects to reduce mutex use
	//unordered_map<int, mat_fl> tmpSparse;
	rowRes* ret = new rowRes(ini_ColPerRow - 1);
	//vector<mat_fl> tmpV(ini_ColPerRow-1,(mat_fl)0);
	bool NAinrow(false);

	LvlUp::iterator fnd;
	while (getline(ss, segments, sep)) {
		cnt2++;
		if (cnt2 == -1) {//this is the row ID
			rowID = segments;
			if (rowID == "mapped") {//to deal with mocat files
				breaker = true;
				break;
			}
			if (doSubsets && subset.find(rowID) == subset.end()) {
				breaker = true;
				break;
			}
			//set index for rows..
			//lock_guard<mutex> guard(protect);
			ret->rowID = rowID;
			ret->curRowIdx = geneCnt;

			if (doHigh) {//check if present in hierachy already..
				fnd = LUp.find(rowID);
				if (fnd == LUp.end()) {//needs to be added to HMat
					ret->taxa = vector<string>(maxLvl, "-1");
					NAinrow = true;
					if (cntNA < 50) {
#ifdef notRpackage
						std::cout << "Row ID " << rowID << " is not in hierachy.\n";// \nAborting..\n"; std::exit(24);
						if (cntNA == 49) { std::cout << " ..\n"; }
#endif
					}
				}
				else {
					ret->taxa = (*fnd).second;
				}
				/*
				string lngTax = "";
				for (int tt = 0; tt < maxLvl; tt++){
				lngTax += taxa[tt] ;
				taxa[tt] = lngTax;
				lngTax += +";";
				}
				*/
			}
			continue;
		}
		mat_fl tmp = atof(segments.c_str());
		/*if (doHigh) {//1:finds relevant rowID, extracts taxa; 2:add on all HighLvl mats
			for (int tt = 0; tt < maxLvl; tt++) {
				HI[tt]->set(taxa[tt], cnt2, tmp);
			}
			colSum[cnt2] += (double)tmp;
			*/
		ret->rowV[cnt2] = tmp;
		if (tmp > 0) {
			ret->Occ++;
		}
	}
	if (cnt2 + 2 != ini_ColPerRow) {
		cerr2("C2: Number of columns on line " + itos(geneCnt) + " is " + itos(cnt2 + 2) + ". Expected " + itos(ini_ColPerRow) + " columns.\n", 64);
	}
	if (NAinrow) {
		ret->NAcnt = 1;
	}
	return ret;

}


string Matrix::prepRow(rowRes* ret, string SEP, bool normalize ) {
	//string retStr(ret->rowID);
	//precision(9)

	std::stringstream retStr;
	retStr << setprecision(9) << ret->rowID;

	for (uint i = 0; i < ret->rowV.size(); i++) {
		if (colSum[i] == 0) {
			retStr << SEP << "0";
		} else if (normalize) {
//			retStr += SEP + to_string(ret->rowV[i] / colSum[i]);
			retStr <<SEP << ret->rowV[i] / colSum[i];
		}else {
			retStr << SEP <<ret->rowV[i];
		}
		//not for multithreads
//		colSum[i] += (double)ret->rowV[i];
	}
	retStr << "\n";
	return retStr.str();
}

void Matrix::addRow(rowRes* ret){
	//set up indexing of current row
	size_t curRowIdx = ret->curRowIdx;//rowIDs.size();
	
	//single core'd
	lock_guard<mutex> lock(protect);
	rowID_hash[ret->rowID] = curRowIdx;
	if (rowIDs.size() <= curRowIdx) {
		rowIDs.resize((curRowIdx+1), "");
	}
	rowIDs[curRowIdx] = ret->rowID;
	//rowIDs.push_back(ret->rowID);
	if (doHigh) {//1:finds relevant rowID, extracts taxa; 2:add on all HighLvl mats
		for (int tt = 0; tt < maxLvl; tt++) {
			HI[tt]->set(ret->taxa[tt], ret->rowV);
		}
	}

	for (uint i = 0; i < ret->rowV.size(); i++) {
		if (loadPrimaryMat) {
			if (sparse) {
				if (ret->rowV[i] > 1e-15) {
					matSp[i][curRowIdx] = ret->rowV[i];
				}
			}
			else {
				if (mat[i].size() <= curRowIdx) {
					mat[i].resize(curRowIdx + 1, 0);
				}

				mat[i][curRowIdx] = ret->rowV[i];//.push_back(ret->rowV[i]);
			}
		}
		colSum[i] += (double)ret->rowV[i];
	}


}


//ClStr related functions
bool Matrix::check_compli(int x, int y) {
	int dblOcc(0), singlOcc(0), noOcc(0);
	int v2s(0), v1s(0);

	if (sparse) {
		for (int i = 0; i < maxCols; i++) {
			auto End = matSp[i].end();	
			auto f1 = matSp[i].find(x);
			auto f2 = matSp[i].find(y);
			if (f1 != End) {
				if (f2 != End) {
					dblOcc++;
				}
				else {
					v1s++;
				}
			}
			else if (f2 != End) {
				v2s++;
			}
			else {
				noOcc++;
			}
		}
	}
	else {
		vector<mat_fl> v1 = getRow(x);
		vector<mat_fl> v2 = getRow(y);
		for (uint i = 0; i < v1.size(); i++) {
			if (v1[i] != 0) {
				if (v2[i] != 0) {
					dblOcc++;
				}
				else {
					v1s++;
				}
			}
			else if (v2[i] != 0) {
				v2s++;
			}
			else {
				noOcc++;
			}
		}
	}
	singlOcc = v2s + v1s;
	//fisher_test();
	float ret = (float)singlOcc / ((float)singlOcc + (float)dblOcc);
	if (ret < CoEx_thresh || singlOcc < 3) {
		return(false);
	}
	float p = fet->ultrafastfet(dblOcc, v1s, v2s,noOcc);
	if (p > CoEx_minp) {	return(false);}
/*	if (v1occ < coexcl_thr || v2occ < coexcl_thr) {
		ret = 0.f;
	}
	*/
	//DEBUG
	//cerr << "Occ: " << ret << " S1:" << v1s << " S2:" << v2s << " D:" << dblOcc << " NO:" << noOcc << " p:"<< p<< endl;
	return (true);
}
void Matrix::prepCoex(options* opt) {
	coexcl_thr = int((float)maxCols / 10 + 0.5);
	CoEx_thresh = 0.9f;  CoEx_minp = opt->pval;
}
void Matrix::complimentarity(vector<string>& cand, vector<float>& id, options* opt,
	vector<vector<int>>& retI, vector<vector<string>>& retS, int& tcID, int& tcCE) {
	vector<int> idx1(cand.size(),0);
	if (coexcl_thr > 50) {
		coexcl_thr = 50;
	}	else if (coexcl_thr < 4) { coexcl_thr = 4; }
	for (uint i = 0; i < cand.size(); i++) {
		//stoi(cand[i]);
		auto idx = rowID_hash.find(cand[i]);
		if (idx == rowID_hash.end()) {
			cerr << cand.size() << " X " << i << " Y " << cand[0]  << "L "<<cand.back() <<endl;
			cerr2("Can't find matrix row \"" + cand[i] + "\"\n", 942);
		}
		else {
			idx1[i] = idx->second;
		}
	}
	//vector<float > tmp(cand.size());
	//vector<vector<float>> score(cand.size(),tmp);
	vector<bool> clust(cand.size(), false);
	for (uint i = 0; i < cand.size(); i++) {
		if (clust[i]) {continue;}
		//string curClus(cand[i]);
		int added = 1; clust[i] = true;
		vector<string> maybeS(1,cand[i]);
		vector<int> maybeI(1,idx1[i]);
		for (uint j = i+1; j < cand.size(); j++) {
			if (clust[j]) { continue; }//already clustered..
			//check co-occurrence pattern
			bool isCompl = check_compli(idx1[i], idx1[j]);
			if ( (i==0 && id[j]>0.99)
				|| isCompl) {
				//curClus += "," + cand[j]; 
				if (!isCompl) {tcID++;} else {tcCE++;}
				added++;
				maybeS.push_back(cand[j]);
				maybeI.push_back(idx1[j]);
				clust[j] = true;
			}
		}
		if (added > 1) {
			retI.push_back(maybeI);
			retS.push_back(maybeS);
		}
	}
}
void Matrix::mergeRows(vector<int>& idx) {
	if (idx.size() < 2) {
		return;
	}
	int tar = idx[0];
	if (sparse) {
		for (uint j = 0; j < matSp.size(); j++) {//sample wise addition
			auto fTar = matSp[j].find(tar);
			auto Etake = matSp[j].end();
			lock_guard<mutex> lock(protect);
			for (uint i = 1; i < idx.size(); i++) {
				int take = idx[i];
				auto fTake = matSp[j].find(take);
				if (fTake != Etake) {
					if (fTar != Etake) {
						fTar->second += fTake->second;
					}
					else {
						matSp[j][tar] = fTake->second;
					}
					fTake->second = 0;
				}
			}
		}
	}
	else {
		for (uint i = 1; i < idx.size(); i++) {
			int take = idx[i];
			lock_guard<mutex> lock(protect);
			for (uint j = 0; j < mat.size(); j++) {//sample wise addition
				mat[j][tar] += mat[j][take];
				mat[j][take] = 0;
			}
		}
	}
}

declutRes* Matrix::declutCalc(vector<string>* txt, options* opts) {
	int nsize = txt->size() - 1;
	//vector<string> curClus(nsize, "");
	vector<float> clusIDs(nsize, 0.f);
	declutRes* dRes = new declutRes(nsize);
	getGenesInClstr((*txt), dRes->curClus, clusIDs);
	if (dRes->curClus.size() < 2) {
		delete txt;
		return dRes;
	}
	//test if these genes have co-exclusion profiles
	//vector<vector<string>> smpls;
	//vector<vector<int>> idxs;
	complimentarity(dRes->curClus, clusIDs, opts,
		dRes->idxs, dRes->smpls, dRes->totalMergedID, dRes->totalMergedCE);//return csv list

	delete txt;
	return dRes;
}

declutRes* Matrix::declutCalc2(vector<string>* txt, options* opts) {
	int nsize = txt->size() - 1;
	//vector<string> curClus(nsize, "");
	declutRes* dRes = new declutRes(nsize);
	vector<float> clusIDs(nsize, 0.f);
	if (txt->size() > clusIDs.size()) {
		clusIDs.resize(txt->size() + 10, 0);
	}

	if (dRes->curClus.size() < 2) {
		delete txt;
		return dRes;
	}
	//test if these genes have co-exclusion profiles
	//vector<vector<string>> smpls;
	//vector<vector<int>> idxs;
	complimentarity((*txt), clusIDs, opts,
		dRes->idxs, dRes->smpls, dRes->totalMergedID, dRes->totalMergedCE);//return csv list

	delete txt;
	return dRes;
}

void Matrix::decluter(options* opts) {
	this->prepCoex(opts);
	string inF = opts->referenceFile; //should be .clstr from cd-hit
	const string outP = opts->output;
	const bool doGZ = opts->gzOut;
	fet = new FET();
	//set up path to all outfile
	string matCmpr = "";
	if (doGZ) { matCmpr = ".gz"; }
	string outF2 = outP + "mat.decl.mat" + matCmpr;
	string outFlist = outP + "concat.list";

	string lastline = "";
	textBlock* clbl;
	FILE* incl;
	int totalMergedID(0), totalMergedCE(0), cases(0);
	incl = fopen(inF.c_str(), "r");
	if (incl == NULL) {
		cerr2("Couldn't open clustering file " + inF + "\n", 55);
	}
	//first get number of total entries..
	int clusBlocks(0);
	bool cont(true);
	while (cont) {
		clbl = getClusBlock(incl, lastline);
		clusBlocks++;
		cont = clbl->cont;
		delete clbl;
	}
	fclose(incl);
	lastline = "";
	cout << "Found " << clusBlocks << "total cluster blocks\n";
	incl = fopen(inF.c_str(), "r");



	ofstream list(outFlist.c_str(), ios_base::out);
	int num_threads = opts->threads;
	vector<job2> slots(num_threads);
	int j(0); int cnt(0);
	while (true) {
		//work through clstr file
		if (j >= num_threads) { j = 0; }
		if (slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
			slots[j].inUse = false;
			declutRes* dR = slots[j].fut.get();
			totalMergedID += dR->totalMergedID;	totalMergedCE += dR->totalMergedCE;
			for (uint j = 0; j < dR->idxs.size(); j++) {//print merged row names
				list << join(dR->smpls[j], ",") << endl;
				this->mergeRows(dR->idxs[j]);
			}
			cases += dR->curClus.size() - dR->idxs.size();
			delete dR;
		}
		if (slots[j].inUse == false) {
			clbl = getClusBlock(incl, lastline);
			bool conti = clbl->cont;
			cnt++;
			vector<string>* trans = new vector<string>(clbl->txt);
			slots[j].fut = async(std::launch::async, &Matrix::declutCalc, this, trans, opts);
			slots[j].inUse = true;
			if (!conti) { break; }
		}
		if (cnt % 50000 == 0) {
			cerr << cnt << "...  ";
		}
		j++;

	}
	//finish jobs..
	for (int j = 0; j < num_threads; j++) {
		if (slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
			slots[j].inUse = false;
			declutRes* dR = slots[j].fut.get();
			totalMergedID += dR->totalMergedID;		totalMergedCE += dR->totalMergedCE;
			for (uint j = 0; j < dR->idxs.size(); j++) {
				list << join(dR->smpls[j], ",") << endl;
				this->mergeRows(dR->idxs[j]);
			}
			cases += dR->curClus.size() - dR->idxs.size();
			delete dR;
		}
	}


	list.close();
	fclose(incl);
	delete fet;
	std::cout << "Merged " << totalMergedID << "/" << totalMergedCE << "/" << cases << " cases\nWriting new matrix to " << outF2 << endl;

	//everything calculated, time to write results out..
	this->writeMatrix(outF2, false);

}


void Matrix::decluter2(options* opts) {
	
	this->prepCoex(opts);
	string inF = opts->referenceFile; //should be .clstr from cd-hit
	const string outP = opts->output;
	const bool doGZ = opts->gzOut;
	fet = new FET();
	cerr << "Starting declutering\n";
	//set up path to all outfile
	string matCmpr = "";
	if (doGZ) { matCmpr = ".gz"; }
	string outF2 = outP + "mat.decl.mat" + matCmpr;
	string outFlist = outP + "concat.list";

	string lastline = "";
	FILE* incl;
	bool cont(true);
	int totalMergedID(0), totalMergedCE(0), cases(0);
	incl = fopen(inF.c_str(), "r");
	if (incl == NULL) {
		cerr2("Couldn't open clustering file " + inF + "\n", 55);
	}

	int clusBlocks(0);
	while (cont) {
		vector<string> curClus = getMMSeqsClus(incl, lastline, cont);
		if (curClus.size() < 2) { continue; }
		clusBlocks++;
	}
	fclose(incl);
	cont = true; lastline = "";
	cout << "Found " << clusBlocks << " total cluster blocks\n";
	incl = fopen(inF.c_str(), "r");



	ofstream list(outFlist.c_str(), ios_base::out);
	vector<float> clusIDs(1000, 0);
	int num_threads = opts->threads;
	vector<job2> slots(num_threads);
	int cnt (0),j(0);
	while (cont) {
		//work through mmseqs file
		if (j >= num_threads) { j = 0; }
		if (slots[j].inUse == true && 
				slots[j].fut.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
			slots[j].inUse = false;
			declutRes* dR;
			dR = slots[j].fut.get();
			totalMergedID += dR->totalMergedID;	totalMergedCE += dR->totalMergedCE;
			for (uint j = 0; j < dR->idxs.size(); j++) {//print merged row names
				list << join(dR->smpls[j], ",") << endl;
				this->mergeRows(dR->idxs[j]);
			}
			cases += dR->curClus.size() - dR->idxs.size();
			delete dR;
		}
		if (slots[j].inUse == false) {
			vector<string> curClus = getMMSeqsClus(incl, lastline, cont);
			if (curClus.size() > 1) {
				cnt++;
				vector<string>* trans = new vector<string>(curClus);
//				declutRes*dR = declutCalc2(trans, opts);
//				delete dR;
				cerr << join(curClus, ",") << endl;
				slots[j].fut = async(std::launch::async, &Matrix::declutCalc2, this, trans, opts);
				slots[j].inUse = true;
				if (cnt % 50000 == 0) {
					cerr << cnt << "/" << clusBlocks << "...  " << endl;
				}
			}
			if (!incl) { cont = false; }
		}
		j++;

		if (!cont) { break; }

/*
		//test if these genes have co-exclusion profiles
		vector<vector<string>> smpls;
		vector<vector<int>> idxs;
		if (curClus.size() > clusIDs.size()) {
			clusIDs.resize(curClus.size() + 10, 0);
		}
		complimentarity(curClus, clusIDs, opts, idxs, smpls, totalMergedID, totalMergedCE);//return csv list
		for (uint j = 0; j < idxs.size(); j++) {
			//print merged row names
			list << join(smpls[j], ",") << endl;
			this->mergeRows(idxs[j]);
			//DEBUG
			//std::cout << "Merge: " << join(smpls[j], ",") << endl;
			//std::cout << "M2: " << join(idxs[j], " ") << endl;
			//totalMerged += idxs[j].size();
			//exit(12);
		}
		cases += curClus.size() - idxs.size();
		*/
	}
	std::cout << "Decluter2:\n";
	//finish jobs..
	for (int j = 0; j < num_threads; j++) {
		if (slots[j].inUse == true ) {
			slots[j].inUse = false;
			declutRes* dR = slots[j].fut.get();
			totalMergedID += dR->totalMergedID;		totalMergedCE += dR->totalMergedCE;
			for (uint j = 0; j < dR->idxs.size(); j++) {
				list << join(dR->smpls[j], ",") << endl;
				this->mergeRows(dR->idxs[j]);
			}
			cases += dR->curClus.size() - dR->idxs.size();
			delete dR;
		}
	}

	cout << "Merged " << totalMergedID << " / " << totalMergedCE << " / " << cases << " cases\n";
	cout << "Writing new matrix to " << outF2 << endl;

	list.close();
	fclose(incl);
	delete fet;

	//everything calculated, time to write results out..
	this->writeMatrix(outF2, false,num_threads);

}
void Matrix::estimateModuleAbund(char ** argv, int argc) {
	char* argX[1];
	options* psOpt = new options(0, argX);
	psOpt->modDB = argv[4];
	psOpt->output = argv[3];
	//string doModKOest = argv[4];
	//read module DB

	//ini options
	psOpt->modRedund = atoi(argv[5]);
	psOpt->modModCompl = (float)atof(argv[6]);
	psOpt->modEnzCompl = (float)atof(argv[7]);
	psOpt->modWrXtraInfo = false;
	//std::cout << argv[8] << endl;
	if (argc > 8 && strcmp(argv[8], "1") == 0) { psOpt->modWrXtraInfo = true; }
}

//minimal sum per column, used for rarefaction min raredep
double Matrix::getMinColSum() {
	if (colSum.size() > 0) {
		double minE = colSum[0];
		for (uint i = 0; i < colSum.size(); i++) {
			if (minE > colSum[i]) {
				minE = colSum[i];
			}
		}
		return minE;
	}
	else {
		return 0;
	}
}

column Matrix::getMinColumn(uint offset ) {
	column* minimalColumn = new column();
	if (colSum.size() > 0) {
		double minE = colSum[0];
		string ID;
		for (uint i = 0; i < colSum.size(); i++) {
			if (minE > colSum[i]) {
				minE = colSum[i];
				ID = colIDs[i];
			}
		}
		minimalColumn->id = ID;
		minimalColumn->colsum = minE;
		return *minimalColumn;
	}
	else {
		return *minimalColumn;
	}
}

void Matrix::resizeMatRows(uint x, mat_fl def ) {
	for (int i = 0; i < maxCols; i++) {
		mat[i].resize(x,def);
	}

}


void Matrix::estimateModuleAbund(options* opts) {
	string outFile = opts->output;
	//ini options
	int redundancy = opts->modRedund;
	float pathCompl = opts->modModCompl;
	float enzyCompl = opts->modEnzCompl;
	bool writeKOused = opts->modWrXtraInfo;

	//read module DB
	Modules* modDB = new Modules(opts->modDB, colIDs);
	modDB->addDescription(opts->modDescr);
	modDB->addHierachy(opts->modHiera);

	
	//get modules that are used in red mods and add them to inKOmatrix )(this matrix)
	vector<string> recUsedMods= modDB->getRedundantUsedMods();
	for (size_t i = 0; i < recUsedMods.size(); i++) {
		rowID_hash[recUsedMods[i]] = rowIDs.size() ;
		rowIDs.push_back(recUsedMods[i]);
	}
	//and resize mat
	resizeMatRows(rowIDs.size());

	//modWrXtraInfo
	
	modDB->setEnzymCompl(enzyCompl);
	modDB->setPathwCompl(pathCompl);
	modDB->setRedund(redundancy);
	//matrix vector of module abudnance
	//vector<vector<mat_fl>> modMat(maxCols);
//	Matrix modMat = Matrix(modDB->modNms(), colIDs);
	vector<vector<string>> modStr(maxCols);
	vector<vector<float>> modScore(maxCols);

	for (int i = 0; i < maxCols; i++) {
		//mat[i] contains KO abundances for 1 sample (i)
		modDB->calcModAbund(mat[i], i, rowID_hash, modStr[i], modScore[i]);			
	}
	//write description
	ofstream of; ofstream of2; string nos;
	of.precision(9);
	of2.precision(9);

	vector<string> moD = modDB->modDescr(); 
	vector<string> moN = modDB->modNms_numbered();
/*	string nos = outFile+".descr";
	of.open(nos.c_str());
	for (size_t i = 0; i < moD.size(); i++) {
		of << moN[i] << "\t" << moD[i] << endl;
	}
	of.close();
	*/
	//write KOs used
	if (writeKOused) {
		nos = outFile + ".KOused";
		of.open(nos.c_str());
		nos = outFile + ".MODscore";
		of2.open(nos.c_str());

		//write SmplIDs
		for (size_t i = 0; i < colIDs.size(); i++) {
			of << "\t" << colIDs[i]; of2 << "\t" << colIDs[i];
		}
		of << endl; of2 << endl;

		for (size_t i = 0; i < moD.size(); i++) {
			bool hasKOUse = false;
			for (size_t k = 0; k < modStr.size(); k++) {
				if (modStr[k][i] != "") { hasKOUse=true; break; }
			}
			if (!hasKOUse) { continue; }
			of << moN[i]; of2 << moN[i];
			for (size_t k = 0; k < modStr.size(); k++) {
				of << "\t" << modStr[k][i] ;
				of2 << "\t" << modScore[k][i] ;
			}
			of << endl; of2 << endl;
		}
		of.close(); of2.close();
	}

	//write module matrix
	modDB->writeMatrix(outFile+".mat",true,opts->modCollapse);//ModPos
	//after writeMatrix to get the ModUsed vec filled
	modDB->writeModDescr(outFile + ".descr", true);

	delete modDB;

}
void Matrix::addColumn(string cname) {
	//increase maxLvl cnt
	maxCols++;
	colIDs.push_back(cname);
	colID_hash[cname] = maxCols - 1;
	mat.push_back(vector<mat_fl>(mat[0].size()));
}
void Matrix::addRow(string rname,vector<mat_fl> in) {
	rowIDs.push_back(rname);
	int curRows = mat[0].size();
	rowID_hash[rname] = curRows;
	for (uint i = 0; i < mat.size(); i++) {
		mat[i].push_back(in[i]);
	}
}
vector<mat_fl> Matrix::getRow(uint idx) {
	if (idx > mat[0].size()) {
		cerr2("Requested index outside of matrix size: " + itos(idx) + "\n", 723);
	}
	vector<mat_fl> ret(mat.size());
	for (uint i = 0; i < mat.size(); i++) {
		ret[i] = mat[i][idx];
	}
	return ret;
}

Matrix::~Matrix(void)
{
	for (unsigned int i = 0; i < HI.size(); i++){
		delete HI[i];
	}

}

void Matrix::normalize() {
	vector<mat_fl> allSums(colIDs.size(), (mat_fl)0);
	for (size_t smpl = 0; smpl<(colIDs.size()); smpl++) {
		mat_fl sums(0);
		for (size_t i = 0; i<rowIDs.size(); i++) {
			sums += mat[smpl][i];
		}
		allSums[smpl] = sums;
	}
	for (size_t smpl = 0; smpl < (colIDs.size()); smpl++) {
		for (size_t i = 0; i < rowIDs.size(); i++) {
			mat[smpl][i] /= allSums[smpl];
		}
	}
	/* //DEBUG
	for (size_t smpl = 0; smpl < (colIDs.size() - 1); smpl++) {
		mat_fl sums(0);
		for (size_t i = 0; i < rowIDs.size(); i++) {
			sums += mat[smpl][i];
		}
		sums++;
	}
	*/
}

void Matrix::transpose(){
	// takes the matrix and transposes it
	// column ID and row ID have to be swapped as well
	vector< vector< mat_fl > >  transpMat(mat[0].size(), vector< mat_fl >(mat.size()));
	vector<double> TmpcolSum(transpMat.size(),0);
	for(uint i = 0; i < mat.size(); i++){
		for(uint j = 0; j < mat[i].size(); j++){
			transpMat[j][i] = mat[i][j];
			// build colsum on the go
			TmpcolSum[j] += (double)mat[i][j];
		}
	}

	// switch column and row names
	vector< string >  rowIDst 	= rowIDs;
	rowIDs 						= colIDs;
	colIDs						= rowIDst;

	// redo colSum so 0.95 * 0 works still for rarefaction
	colSum = TmpcolSum;


	// swap the matrices
	mat = transpMat;
}

vector<mat_fl> Matrix::getRowSums() {
	vector<mat_fl> rowSums(rowIDs.size(), 0);
	for (size_t i = 0; i < rowIDs.size(); i++) {
		for (size_t smpl = 0; smpl < (colIDs.size() ); smpl++) {
			rowSums[i] += mat[smpl][i];
		}
	}
	return rowSums;
}

void Matrix::writeMatrix(const string of, bool onlyFilled, int threads ) {
	ostream * out;
	if (isGZfile(of)) {
#ifdef _gzipread
		out = new ogzstream(of.c_str(), ios::out);
		cout << "Writing gzip'd matrix\n";
#else
		cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
	}
	else {
		out = new ofstream(of);
	}


	if (!(*out)) {
		cerr2("Couldn't open matrix output " + of + "\n", 57);
	}
	//out.open(of.c_str(), ios_base::out);
	out->precision(9); *out << "Gene";
	for (size_t smpl = 0; smpl < (colIDs.size() ); smpl++) {
		*out << "\t" << colIDs[smpl ];
	}
	*out << endl;
	vector<mat_fl> rowSums;
	//size_t cidS(colIDs.size());
	vector<future<string>> job (threads);
	future<bool> wrStr; string tmpSt("");
	//ini empty job
	//wrStr = async(std::launch::async, writeChunk, out, tmpSt);

	string tmp(""), tmp2("");
	if (onlyFilled) { rowSums = getRowSums(); }
	int i = -1;
	while (i<(int)rowIDs.size()) {
		uint j = 0;
		for (; j < threads; j++) {
			i++;
			if (onlyFilled && rowSums[i] == 0) {
				j--; continue;
			}
			if (i >= (int)rowIDs.size()) {
				break;
			}
			job[j] = async(std::launch::async, &Matrix::row2line, this, i);
		}
		//bool x= wrStr.get();		wrStr = async(std::launch::async, writeChunk, out,tmp);
		*out << tmp;
		tmp = "";
		for (uint jk = 0; jk < j; jk++) {
			tmp += job[jk].get();
		}
	}
	*out << tmp;
	
	delete out;// .close();
}
string Matrix::row2line(int i) {
	//*out << rowIDs[i] ;
	string tmp("");
	tmp += rowIDs[i];
	if (sparse) {
		for (size_t smpl = 0; smpl < colIDs.size(); smpl++) {
			auto fnd = matSp[smpl].find(i);
			if (fnd != matSp[smpl].end()) {
				tmp += "\t" + to_string(fnd->second);
			}
			else {
				tmp += "\t0";
			}
		}
	}
	else {
		for (size_t smpl = 0; smpl < colIDs.size(); smpl++) {
			tmp += "\t" + to_string(mat[smpl][i]);
		}
	}
	tmp += "\n";
	return tmp;
}
void Matrix::read_subset_genes(const string xtra){
	string line;
	ifstream in(xtra.c_str());
	string ID = ""; int cnt = 0;
	if (!in){
		cerr2("Can't open geneID file " +xtra +"\n",13);
	 }
	while (getline(in, ID, '\n')) {
		subset[ID] = 1;
		cnt++;
	}
	in.close();
	//just check that something was read
	if (cnt >= 1){
		doSubsets = true;
		#ifdef notRpackage
		std::cout << "Read " << cnt << " gene subsets\n";
		#endif
	}
}
void Matrix::read_hierachy(const string xtra, bool xtdHiera){
	int maxHir = 9;
	vector<string> features;
	string line;
	ifstream in(xtra.c_str());
	int cnt = 0;
	if (!in){
		cerr2 ("Can't open hierachy file " +xtra +"\n",13);
	}
	for (int k = 0; k < maxHir; k++){
		//string tmp = ;
		LvlNms.push_back("L" + stringify((double)k));
	}

	while (getline(in, line, '\n')) {
		if (line.substr(0, 1) == "#"){ continue; }
		vector<string> pseudo(maxHir, "?");
		cnt++;
		string segs;
		string segs2;
		stringstream ss;
		ss << line;
		getline(ss, segs, '\t');
		getline(ss, segs2, '\t');
		string spl;
		stringstream hir; hir << segs2;
		int i = 0;
		//string lngTax="";
		string prevH = "";
		while (getline(hir, spl, ';')){
			//add previous levels
			pseudo[i] = prevH+spl;
			if (xtdHiera) {
				prevH = pseudo[i] + ";";
			}
			i++;
			//lngTax += spl + ";";
			//index possible features
			if (i >= maxHir){ break; }
		}
		if (i > maxLvl){ maxLvl = i; }
		LUp[segs] = pseudo;
	}
	in.close();
	#ifdef notRpackage
	std::cout << "Read hierachy. Found " << maxLvl << " hierachical levels.\n";
	#endif
}
void Matrix::splitOnHDD(string out_seed){
	for (size_t smpl=0;smpl<(colIDs.size()-1); smpl++){
		string oF2 = out_seed + sampleNameSep + colIDs[smpl+1];
		ofstream out;
		out.open(oF2.c_str(),ios_base::out);
		out.precision(12);
		for (size_t i=0;i<rowIDs.size();i++){
			if (mat[i][smpl]==0){continue;}
			out<< rowIDs[i]<<"\t"<<mat[i][smpl]<<endl;
		}
		out.close();
		//if (smpl>20){std::exit(9);}
	}
}


void Matrix::writeSums(string out_seed){
	string oF2 = out_seed + sampleNameSep + "sums.txt";
	ofstream out;
	out.open(oF2.c_str(),ios_base::out);
	out.precision(12);
	//for (OTUid = colID_hash.begin(); OTUid != colID_hash.end(); OTUid++) {
	for (size_t smpl = 0; smpl<(colIDs.size() - 1); smpl++) {
		mat_fl sums(0);
		for (size_t i=0;i<rowIDs.size();i++){
			sums+=mat[i][smpl];
		}
		out<<colIDs[smpl+1]<<"\t"<<sums<<endl;
	}
	out.close();
}


bool sortPair(const pair<double, string>& left, const pair<double, string>& right){
	return left.second < right.second;
}

vector< pair <double, string>> Matrix::getColSums(bool sorted){
	// fill colsums stepwise
	if(colsums.size() == 0){
		for(uint i = 0; i < colIDs.size(); i++){
			pair<double, string> p(colSum[i], colIDs[i]);
			colsums.push_back(p);
		}
	}
	if(sorted == false){
		return colsums;
	}else{
		// now we sort this vector
		std::sort(colsums.begin(), colsums.end(), sortPair);
		return colsums;
	}
}
void Matrix::writeColSums(string outF){
	vector< pair <double, string>> colsumsUnsort 	= getColSums(false);
	vector< pair<double, string>> colsumsSort 		= getColSums(true);
	string outF2 = outF + "colSums.txt";
	ofstream out;
	out.open(outF2.c_str(),ios_base::out);
	out.precision(12);
	for( auto it = colsumsUnsort.begin(); it != colsumsUnsort.end(); it++ ){
		out << it->second<<"\t"<< it->first<<std::endl;
	}
	out.close();

	// sorted
	outF2 = outF + "colSums_sorted.txt";
	out.open(outF2.c_str(),ios_base::out);
	out.precision(12);
	for( auto it = colsumsSort.begin(); it != colsumsSort.end(); it++ ){
		out << it->second<<"\t"<< it->first<<std::endl;
	}
	out.close();
}

void Matrix::writeRowOcc(options* opt) {

	string outF = opt->output;
	//prevent overwritting...
	if (outF == opt->input) { outF += "rowOcc"; }
	ofstream out(outF);
	if (rowOcc.size() != rowIDs.size()) {
		cerr << "vector mismatch rowOcc,rowIDs:" << rowOcc.size() << ", " << rowIDs.size() << endl;
		exit(45);
	}
	for (uint i = 0; i < rowOcc.size(); i++) {
		out << rowIDs[i] << "\t" << rowOcc[i] << "\n";
	}
	out.close();
}
//Saprse Matrix class
SparseMatrix::SparseMatrix() :mat(0), colNames(0), rowIDs(0){}

void SparseMatrix::addCount(string smpl, int row, smat_fl abund) {
	//find sample
	SmplAbunIT tar = mat[row].find(smpl);
	if (tar == mat[row].end()) {//create entry & expand matrix
		mat[row][smpl] = abund;
	} else {
		(*tar).second += abund;
	}

//keep track of listed samples
	SmplOccurIT tar2 = colNames.find(smpl);
	if (tar2 == colNames.end()) {
		colNames[smpl] = 1;
	} else {
		(*tar2).second++;
	}
}


HMat::HMat(string L, vector<string> Samples, vector<string> Features)
:LvlName(L), FeatureNs(Features), SampleNs(Samples),mat(0), catCnt(0), maxCatCnt(0), hiTaNAcnt(0){
	empty = vector<mat_fl>(SampleNs.size(), 0);
	mat.resize(FeatureNs.size(), empty);
	catCnt.resize(FeatureNs.size(), vector<uint>(SampleNs.size(), 0));
	for (unsigned int i = 0; i < FeatureNs.size(); i++){
		Feat2mat[FeatureNs[i]] = i;
	}
}

void HMat::set(string& kk, int j, mat_fl v) {

	//no reason to continue..
	if (v <= 0.f) { return; }

	mat_fl div(1);
	//check if key is single key, or multiple keys (split by | )
	size_t pos(kk.find("|", 0)), npos(0);
	vector<string> subkk(0);
	while (pos != string::npos) {
		subkk.push_back(kk.substr(npos, pos - npos));
		npos = pos + 1;
		pos = kk.find("|", npos);
		div += 1.f;
	}
	subkk.push_back(kk.substr(npos));

	for (uint t = 0; t < subkk.size(); t++) {
		string yy = subkk[t];
		//find key
		std::map<string, int>::iterator i = Feat2mat.find(yy);
		if (i == Feat2mat.end()) {
			Feat2mat[yy] = (int)mat.size();
			mat.push_back(vector<mat_fl>(SampleNs.size(), 0));
			catCnt.push_back(vector<uint>(SampleNs.size(), 0));
			FeatureNs.push_back(yy);
			i = Feat2mat.find(yy);
			hiTaNAcnt++;
			//
			if (hiTaNAcnt < 100) {
				cerr2("Could not find entry " + yy + " in registered subset\n");
			}
			//std::exit(23);

		}
		if ((*i).second > (int)mat.size()) {
			cerr2("implied row index larger than high level mat!\nAborting..", 25);
		}
		mat[(*i).second][j] += v / div;
		catCnt[(*i).second][j]++;
		//if (j == 0) {			catCnt[(*i).second]++;		}
	}
}
void HMat::set(string& kk, vector<mat_fl>& v) {

	mat_fl div(1);
	//check if key is single key, or multiple keys (split by | )
	size_t pos(kk.find("|", 0)), npos(0);
	vector<string> subkk(0);
	while (pos != string::npos) {
		subkk.push_back(kk.substr(npos, pos - npos));
		npos = pos + 1;
		pos = kk.find("|", npos);
		div += 1.f;
	}
	div = 1/div;
	subkk.push_back(kk.substr(npos));

	for (uint t = 0; t < subkk.size(); t++) {
		string yy = subkk[t];
		//find key
		std::map<string, int>::iterator i = Feat2mat.find(yy);
		if (i == Feat2mat.end()) {
			Feat2mat[yy] = (int)mat.size();
			mat.push_back(vector<mat_fl>(SampleNs.size(), 0));
			catCnt.push_back(vector<uint>(SampleNs.size(), 0));
			FeatureNs.push_back(yy);
			i = Feat2mat.find(yy);
			hiTaNAcnt++;
			//
			if (hiTaNAcnt < 100) {
				cerr2("Could not find entry " + yy + " in registered subset\n");
			}
			//std::exit(23);

		}
		int idx = (*i).second;
		if (idx > (int)mat.size()) {
			cerr2("implied row index larger than high level mat!\nAborting..", 25);
		}
		for (size_t j = 0; j < v.size(); j++) {
			if (v[j] > 0) {
				mat[idx][j] += v[j] * div;
				catCnt[(*i).second][j]++;
			}
		}
		//if (j == 0) {			catCnt[(*i).second]++;		}
	}
}

void HMat::prep() {
	maxCatCnt.resize(catCnt.size());
	for (size_t i = 0; i < catCnt.size(); i++) {
		uint mxaOcc = 0;
		for (size_t j = 0; j < catCnt[i].size(); j++) {
			if (mxaOcc < catCnt[i][j]) { mxaOcc = catCnt[i][j]; }

		}
		maxCatCnt[i] = mxaOcc;
	}
}
void HMat::print(ofstream& O,bool mean,uint occMin ,uint occPerSmpl){
	O << LvlName ;
	for (size_t i = 0; i < SampleNs.size(); i++){
		O << "\t" << SampleNs[i];
	}
	//header done

	for (size_t i = 0; i < FeatureNs.size(); i++){
		//bool doMean = mean;
		if (FeatureNs[i] == "-1") {
			O << "\n" << FeatureNs[i];
			//doMean = false; //-1 does not get mean'd
			for (size_t j = 0; j < SampleNs.size(); j++) {
				O << "\t" << mat[i][j];
			}
			continue;
		}
		//not enough overall support for category?
		if (maxCatCnt[i] <= occMin) { continue; }
		O << "\n" << FeatureNs[i];
		for (size_t j = 0; j < SampleNs.size(); j++){
			if (catCnt[i][j] <= occPerSmpl) {
				O << "\t0";
			} else if (!mean) {
				O << "\t" << mat[i][j];
			} else {
				O << "\t" << (mat[i][j] / (mat_fl)catCnt[i][j]);
			}
		}
	}
}



VecFiles::VecFiles(const string inF, const string outF, const string xtra) :
infiles(0){
}

int VecFiles::getIdx(const string&){
	int idx(-1);
	return idx;
}


void VecFiles::readVecFile(const string inF){
	string line;
	ifstream in(inF.c_str());

	//int ini_ColPerRow(0),cnt(0);

	//string rowID="";
	while(getline(in,line,'\n')) {
		if(line.substr(0,1) == "#" || line.length()<2){continue;}
		string segments;
		//int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		int cnt2(-1);
		//int CurIdx(-1);
		while (getline(ss,segments,'\t')) {
			cnt2++;
			if (cnt2 == -1){
				//rowID = segments;
				//CurIdx = this->getIdx(segments);
				continue;
			}
			mat_fl tmp =  atof(segments.c_str());
			if (tmp==0){continue;}

		}
	}

}
