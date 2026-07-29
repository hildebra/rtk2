#include "ClStr2Mat.h"

static string trimWhitespace(const string& value) {
	size_t begin = 0;
	size_t end = value.size();
	while (begin < end && isspace(static_cast<unsigned char>(value[begin]))) { begin++; }
	while (end > begin && isspace(static_cast<unsigned char>(value[end - 1]))) { end--; }
	return value.substr(begin, end - begin);
}


void printVec(clusWrk * curClus, CentralOutputWriter* mO, ofstream* gN) {//, const vector<bool>& useSmpl) {
	if (curClus == nullptr) {
		cerr2("printVec received null cluster pointer\n", 68);
	}
	int64_t CLidx = curClus->Clnum;
	(*gN) << CLidx << curClus->geneNamesStr;
	mO->write(to_string(CLidx) + curClus->matStr);
	if (!(*gN)) {
		delete curClus;
		cerr2("I/O error while writing geneMat output rows\n", 69);
	}
	delete curClus;
}
clusWrk* ClStr2Mat::workClusBlock(textBlock* inVs,
	 int64_t CLidx, bool doSups) {
	clusWrk* ret(nullptr);
	if (!doSups) {
		ret = new clusWrk(static_cast<int>(smplN), CLidx);
	} else {
		ret = new clusWrk(static_cast<int>(smplNsup), CLidx);
	}
	auto reportMalformed = [&](const string& reason, const string& entry) {
		constexpr rtk_count maxMalformedWarnings = 5;
		lock_guard<mutex> guard(reportMtx);
		malformedCDHit++;
		if (malformedCDHitWarned < maxMalformedWarnings) {
			malformedCDHitWarned++;
			cerr << "Warning: malformed CD-HIT entry in cluster " << CLidx << " (" << reason << "): " << entry << "\n";
		}
	};
	bool repFound = false;
	for (uint i = 0; i < inVs->txt.size(); i++) {
		string &line = inVs->txt[i];
		if (i == 0) {//new cluster, add headerInfo
			if (line.empty() || line[0] != '>') { reportMalformed("invalid cluster header", line); }
			ret->geneNamesStr = "\t" + line;
			continue;
		}


		//1 get gene, sample (deparse)
		size_t pos = line.find("nt, >");
		if (pos == string::npos) {
			reportMalformed("missing 'nt, >' marker", line);
			continue;
		}
		size_t pos2 = line.find("...", pos + 4);
		if (pos2 == string::npos || pos2 <= pos + 5) {
			reportMalformed("missing or empty gene identifier", line);
			continue;
		}
		string gene = line.substr(pos + 5, pos2 - pos - 5);

		const size_t finalMarker = line.find_last_not_of(" \t\r");
		if (!repFound && finalMarker != string::npos && line[finalMarker] == '*') {//report representative gene
			ret->geneNamesStr += "\t" + gene + "\n";
			repFound = true;
		}

		//bool geneInAssembl(true);
		pos = gene.find(sampleSep);
		if (pos == string::npos || pos == 0) {
			reportMalformed("gene identifier is missing a sample separator", line);
			continue;
		}
		string sample;
		//SmplOccurITmult smNum;
		pos2 = gene.find("_L", pos + sampleSep.size());
		if (pos2 == string::npos) {
			reportMalformed("gene identifier is missing a contig locus marker", line);
			continue;
		}
		{ //has the characteristic "__" sample separator
			sample = gene.substr(0, pos);
			auto smNum = smpls.find(sample);
			if (smNum == smpls.end()) {
#ifdef notRpackage
				lock_guard<mutex> guard(reportMtx);
				cerr << "Warning: unknown geneMat sample name: " << sample << "::" << gene << endl;
#endif
				continue;
			}
			//2 get abundance
			//Contig + Sample Info will allow to create contig linkage between samples (CCH)
			//int contig = atoi(gene.substr(pos + 3, pos2-pos-3).c_str());
			//can be several samples (combined assembly)
			const vector<int>& smplLocs = (*smNum).second;
			for (size_t jj = 0; jj < smplLocs.size(); jj++) { //the loop takes account for multiple samples being grouped together in map
															  //CCH->addHit(smplLocs[jj], contig);
				int idxM = smplLocs[jj];
				//don't use sample, if not last in map group!
				//if (!useSmpl[idxM]) { continue; }
				//smat_fl abundance = this->GAabundance(idxM, gene);//GAs[idxM]->getAbundance(gene);
				//3 add to matrix / output vector
				ret->Vec[idxM] += this->GAabundance(idxM, gene);//abundance;
				//SmplSum[idxM] += abundance;
			}
			if (doSups) {
				for (size_t jj = 0; jj < smplLocs.size(); jj++) { //the loop takes account for multiple samples being grouped together in map
					int idxM = smplLocs[jj];
					if (smplSupIdx[idxM] != -1) {
						ret->Vec[smplSupIdx[idxM]] += this->GASupAbundance(idxM, gene);
					}
				}
			}
		}

	}


	if (!repFound) {
		ret->geneNamesStr += "\t?\n";
		lock_guard<mutex> guard(reportMtx);
		cerr << "Warning: no representative sequence found for cluster " << CLidx << "\n";
	}
	delete inVs;

	this->addSums(ret);

	ret->generateMatStr();


	return (ret);

}







ClStr2Mat::ClStr2Mat(options* opts):lastClIdWr(1),GAs(0), GAsSup(0),
	smplLoc(0), //baseP(0), 
	mapGr(0),SmplSum(0),
	smplN(0), smplNsup(0),
	sampleSep("__"), malformedCDHit(0), malformedCDHitWarned(0)
{

	std::cout << "Gene clustering matrix creation\n";

	ifstream incl2;
	//FILE* incl;
	const string inF = opts->input;
	const string outF = opts->output;
	const string mapF = opts->map;
	if (mapF.empty()) { cerr2("geneMat requires at least one mapping file via -map.\n", 56); }
	//const string basePX = opts->referenceDir;
	const bool doGZ = opts->gzOut;
	//flag whether supplementary abundances will be included in output matrix
	const bool calcSupplCov = opts->calcSupplCov;
	//set up baseP
	//stringstream ss(basePX); string segments;
	//while (getline(ss, segments,',') ) { baseP.push_back(segments); }

	/*incl = fopen(inF.c_str(), "r");
	if (incl == NULL) {
		cerr2("Couldn't open clustering file " + inF + "\n", 55);
	}*/

	incl2.open(inF.c_str());
	if (!incl2) {
		cerr2("Couldn't open clustering file " + inF + "\n", 55);
	}

	string matCmpr = "";
	if (doGZ) { matCmpr = ".gz"; }
	string outF2 = outF + ".mat" + matCmpr;
	matO = new CentralOutputWriter(
		outF2, static_cast<size_t>(std::max(1u, opts->threads)));
	geneNames = new ofstream(outF + ".genes2rows.txt");
	if (!(*geneNames)) {
		cerr2( "Couldn't open report file " + outF + ".genes2rows.txt" +"\n", 56);
	}
	


	//read map(s) and check that
	stringstream ss2(mapF);  string segments;
	while (getline(ss2, segments, ',')) {
		segments = trimWhitespace(segments);
		if (segments.empty()) { cerr2("geneMat received an empty mapping-file entry.\n", 56); }
		string baseP = ""; //SmplOccurMult CntMapGrps ;
		//also reads abundances
		read_map(segments,  baseP, opts);
		//read_abundances(CntMapGrps, baseP, opts->calcCoverage, opts->calcCovMedian, opts->oldMapStyle);
		
	}

	//smplnames in out matrix
	vector<string> SmplNmsVec(smplN, "");
	for (auto it = smpls.begin(); it != smpls.end(); it++) {
		vector<int> XX = (*it).second;
		if (XX.size() == 0) { continue; }
		SmplNmsVec[ XX[XX.size()-1] ] = smplRid[(*it).first];
	}
	if (GAs.size() != SmplNmsVec.size()) { cerr2("Primary abundance/sample count mismatch while building gene matrix.\n", 68); }
	for (size_t i = 0; i < SmplNmsVec.size(); i++) {
		if (SmplNmsVec[i].empty()) { cerr2("Could not resolve output name for geneMat sample index " + to_string(i) + ".\n", 68); }
	}
	std::cout << "Calculating abundance matrix on "<< smplN <<" samples (this may take awhile)\n";
	std::ostringstream matrixHeader;
	matrixHeader << "Genes";
	for (size_t i = 0; i < SmplNmsVec.size(); i++) {
		//if (!useSmpl[i]) { continue; }
		matrixHeader << "\t" << SmplNmsVec[i];
	}
	SmplSum.resize(smplN, 0.f);
	if (calcSupplCov) {
		//some basic assertions..
		if (SmplNmsVec.size() != GAsSup.size() || GAs.size() != GAsSup.size()) {
			cerr2("Supplementary abundance/sample count mismatch while building gene matrix.\n", 68);
		}
		int supcnt = static_cast<int>(GAsSup.size());
		smplSupIdx.resize(GAsSup.size(), -1);
		for (size_t i = 0; i < SmplNmsVec.size(); i++) {
			//if (!useSmpl[i]) { continue; }
			if (!GAsSup[i]->isEmpty()) {
				matrixHeader << "\t" << SmplNmsVec[i] + ".sup";
				smplSupIdx[i] = supcnt; supcnt++;
			}
		}
		smplNsup = supcnt;
		SmplSum.resize(smplNsup, 0.f);
	}


	matrixHeader << "\n";
	matO->write(matrixHeader.str());

	(*geneNames) << "#GID\tCluster\tRepSeq\n";
	//CCH = new ContigCrossHit((int)smplN);
	//CCH->setSmplNms(SmplNmsVec);
	//SparseMatrix * mat = new SparseMatrix();
	int64_t CLidx = 2; //start at number 2
	//const string sampleStrSep = "__";
	//sampleSep = "__";
	//vector<smat_fl> SmplSum(smplN, 0.f);
	
	int numthr = opts->threads-1;// -1 to include wrThr
	if (numthr <= 0) { numthr = 1; }
	ThreadPool pool(static_cast<size_t>(numthr));

	// Keep only a bounded window of jobs and consume it in submission order.
	// This preserves deterministic output without allowing the reorder buffer to grow with the input.
	deque<future<clusWrk*>> pendingClusters;
	string lastline = "";
	rtk_count tCnt(0);
	double total_Counts(0);
	bool abortWork(false);
	auto collectOldestCluster = [&]() {
		clusWrk* curClus = nullptr;
		try {
			curClus = pendingClusters.front().get();
		}
		catch (const std::exception& ex) {
			cerr << "geneMat worker exception: " << ex.what() << "\n";
			abortWork = true;
		}
		catch (...) {
			cerr << "geneMat worker exception: unknown\n";
			abortWork = true;
		}
		pendingClusters.pop_front();
		if (curClus != nullptr) {
			if (abortWork) { delete curClus; }
			else { manage_write(curClus); }
		}
	};

	while (true) {
		textBlock* clbl = getClusBlock(incl2, lastline);
		if (clbl->txt.size() == 0 || (clbl->txt.size() == 1 && clbl->txt[0].empty())) {
			delete clbl;
			break;
		}
		const bool hasMoreClusters = clbl->cont;
		tCnt++;
		try {
			pendingClusters.push_back(pool.enqueue(&ClStr2Mat::workClusBlock, this, clbl, CLidx, calcSupplCov));
		}
		catch (...) {
			delete clbl;
			throw;
		}
		CLidx++;
		if (pendingClusters.size() >= static_cast<size_t>(numthr)) { collectOldestCluster(); }
		if (abortWork || !hasMoreClusters) { break; }
	}
	while (!pendingClusters.empty()) { collectOldestCluster(); }
	if (abortWork) {
		cerr2("geneMat aborted due to worker failure; output may be incomplete\n", 70);
	}
	finish_write();
	if (!tmpSave.empty()) {
		cerr2("geneMat writer queue not empty after drain; missing cluster blocks or ordering gap detected\n", 73);
	}
	if (malformedCDHit > 0) {
		const rtk_count suppressed = malformedCDHit - malformedCDHitWarned;
		cerr << "Warning summary: geneMat encountered " << malformedCDHit << " malformed CD-HIT entries";
		if (suppressed > 0) {
			cerr << " (" << malformedCDHitWarned << " shown, " << suppressed << " suppressed)";
		}
		cerr << ".\n";
	}

	for (size_t i = 0; i < SmplNmsVec.size(); i++) {
		total_Counts += SmplSum[i];
	}

	std::cout << "Read/wrote " << tCnt << " total genes counting to " << total_Counts<< " counts in matrix\n";
	

	incl2.close();
	geneNames->flush();
	matO->flush();
	geneNames->close();
	matO->close();
	delete matO;
	delete geneNames;

	ofstream* matS = new ofstream((outF + ".mat.sum"), ofstream::out);
	if (!(*matS)) {
		delete matS;
		cerr2("Couldn't open output sum file " + outF + ".mat.sum\n", 71);
	}
	(*matS) << setprecision(numeric_limits<mat_fl>::max_digits10);
	for (size_t i = 0; i < SmplNmsVec.size(); i++) {
		(*matS) << SmplNmsVec[i] << "\t" << SmplSum[i] << endl;
	}
	if (smplN < smplNsup) {
		for (size_t i = 0; i < SmplNmsVec.size(); i++) {
			if (smplSupIdx[i] != -1) {
				(*matS) << SmplNmsVec[i] + ".sup" << "\t" << SmplSum[smplSupIdx[i]] << endl;
			}
		}
	}
	if (!(*matS)) {
		matS->close();
		delete matS;
		cerr2("I/O error while writing output sum file " + outF + ".mat.sum\n", 72);
	}
	matS->close();
	delete matS;
	std::cout << "Finished Gene matrix conversion\n";

}



ClStr2Mat::~ClStr2Mat() {
	for (size_t i = 0; i < GAs.size(); i++) { delete GAs[i]; }
	GAs.resize(0);
	for (size_t i = 0; i < GAsSup.size(); i++) { delete GAsSup[i]; }
	GAsSup.resize(0);
	
	//delete CCH;
}
void ClStr2Mat::finish_write() {
	//if (wrThr.inUse) {		wrThr.fut.get();		wrThr.inUse = false;	}
	if (tmpSave.size() > 0) {
		auto saveCl = tmpSave.begin();
		while (saveCl != tmpSave.end()) {
			if ((*saveCl)->Clnum == (lastClIdWr + 1)) {
				printVec((*saveCl), matO, geneNames); //useSmpl
				/*wrThr.fut.get();
				wrThr.fut = async(std::launch::async, printVec, (*saveCl), matO, geneNames,	useSmpl);
				wrThr.inUse = true;*/
				lastClIdWr++;
				tmpSave.erase(saveCl);
				saveCl = tmpSave.begin();
			}
			else {
				saveCl++;
			}
		}
	}

	//if (wrThr.inUse) {		wrThr.fut.get();		wrThr.inUse = false;	}
}
void ClStr2Mat::manage_write(clusWrk* curClus) {
	//push into wrThr.. (only one can run at a time, therefore wait)
	//bool asyncWr = false;
	//if (wrThr.inUse) {	wrThr.fut.get();	wrThr.inUse = false;	}
	finish_write();
	if (curClus->Clnum == (lastClIdWr + 1)) {
		//wrThr.fut = async(std::launch::async, printVec, curClus, matO, geneNames);
				//,useSmpl);
		//wrThr.inUse = true;
		printVec(curClus, matO, geneNames);
		lastClIdWr++;
	} else {
		tmpSave.push_back(curClus);
		//tmpSave[(lastClIdWr + 1)] = curClus;
	}

}

void ClStr2Mat::read_map(const string mapF, string& baseP, options* opts) {
	ifstream in;
	bool oldFolderStructure = opts->oldMapStyle;
	const bool calcSupplCov = opts->calcSupplCov;

	int map2folderIdx = 0; if (oldFolderStructure) { map2folderIdx = 1; }
	//curr++;//keep track of different maps and inPaths
	uint preMapSize((int)smplLoc.size());

	in.open(mapF.c_str());
	if (!in) {
		cerr << "Couldn't open mapping file " << mapF << endl;
		exit(56);
	}
	std::cout << "Reading map " << mapF << endl;// " on path " << baseP[curr] << endl;
	SmplOccurMult CntAssGrps;
	string line(""); int cnt(-1); int assGrpN(-1);
	//mapping group params
	int mapGrpN(-1); //bool fillMapGrp(false); 
	int artiCntAssGrps(0); int skSmplCol(-1);

	string baseP1 = ""; string baseP2 = "";//#OutPath #RunID
	//string baseP = "";


	//1st part: just read headers
	while (getline(in, line)) {
		cnt++; int sbcnt(-1);
		stringstream ss(line); string segments;
		if (line.substr(0, 1) == "#") { //read column position of smpl id, path & assGrps
			//if (cnt > 0) { continue; }
			if (cnt == 0) {
				while (getline(ss, segments, '\t')) {
					segments = trimWhitespace(segments);
					sbcnt++;
					if (sbcnt == 0 && segments != "#SmplID") {
						cerr << "Map has to start with tag \"#SmplID\"\n"; exit(83);
					}
					if (sbcnt == 1 && !(segments == "Path" || segments == "SmplPrefix")) {
						cerr << "Map has to have tag \"Path\" || \"SmplPrefix\" as second entry\n";
						exit(83);
					}
					if (segments == "AssmblGrps") {
						assGrpN = sbcnt;
						std::cout << "Found Assembly groups in map\n";
					}
					if (segments == "MapGrps") {
						mapGrpN = sbcnt;
						//fillMapGrp = true;
						std::cout << "Found Mapping groups in map\n";
					}
					if (segments == "ExcludeAssembly") {
						skSmplCol = sbcnt;
						std::cout << "Samples can be excluded from assembly\n";
					}
				}
			}
			while (getline(ss, segments, '\t')) {
				segments = trimWhitespace(segments);
				if (segments == "#OutPath") {
					getline(ss, segments, '\t'); baseP1 = trimWhitespace(segments);
				}
				if (segments == "#RunID") {
					getline(ss, segments, '\t'); baseP2 = trimWhitespace(segments);
				}
			}

			continue;
		}
	}
	if (baseP2 == "") { cerr << "No #RunID specified in map!?"; exit(123); }
	if (baseP1 == "") { cerr << "No #OutPath specified in map!?"; exit(123); }
	baseP = baseP1 + "/" + baseP2 + "/";

	//2nd part: read samples
	SmplOccurMult CntMapGrps ;
	

	in.clear();                 // clear fail and eof bits
	in.seekg(0, std::ios::beg); // back to the start!
	cnt = 0;
	while (getline(in, line)) {
		cnt++;
		if (line.substr(0, 1) == "#") {
			continue;
		}
		stringstream ss(line); string segments;
		vector<string> curLine(0);
		while (getline(ss, segments, '\t')) {
			curLine.push_back(trimWhitespace(segments));
		}
		if (curLine.size() < 2) {
			cerr << "Skipping malformed map line " << cnt << " in " << mapF << "\n";
			continue;
		}
		if (skSmplCol > -1 && (int)curLine.size() > skSmplCol && curLine[skSmplCol] == "1") { continue; }
		if (assGrpN >= 0 && (int)curLine.size() <= assGrpN) {
			curLine.resize(assGrpN + 1, "");
		}
		if (mapGrpN >= 0 && (int)curLine.size() <= mapGrpN) {
			curLine.resize(mapGrpN + 1, "");
		}
		if ((int)curLine.size() <= map2folderIdx) {
			cerr << "Skipping malformed map line " << cnt << ": missing Path/SmplPrefix in " << mapF << "\n";
			continue;
		}

		string smpID = curLine[0];
		if (smpID.empty()) {
			cerr << "Skipping malformed map line " << cnt << ": empty sample id in " << mapF << "\n";
			continue;
		}

		//idx 1 for old folder structure, 0 for new folder structure
		string subDir = curLine[map2folderIdx];


		//assembly groups
		string assGrp("");
		if (assGrpN != -1) {
			//handles assembly groups from here
			assGrp = curLine[assGrpN];
		}
		else {//simulate CntAssGrps
			assGrp = itos(artiCntAssGrps);
			artiCntAssGrps++;
		}
		if (CntAssGrps.find(assGrp) != CntAssGrps.end()) {
			CntAssGrps[assGrp].push_back((int)smplLoc.size());
		}
		else {
			CntAssGrps[assGrp] = vector<int>(1, (int)smplLoc.size());
		}

		if (assGrp != "" && CntAssGrps[assGrp].size() > 1) {
			string nsmpID = smpID + "M" + std::to_string(CntAssGrps[assGrp].size());
			if (smpls.find(nsmpID) != smpls.end()) {
				cerr << "Double sample ID: " << nsmpID << endl;
				exit(12);
			}
			smpls[nsmpID] = CntAssGrps[assGrp];//(int)smplLoc.size();
			smplRid[nsmpID] = smpID;
		}
		else {
			if (smpls.find(smpID) != smpls.end()) {
				cerr << "Double sample ID: " << smpID << endl;
				exit(12);
			}
			smpls[smpID] = vector<int>(1, (int)smplLoc.size());
			smplRid[smpID] = smpID;
		}

		//mapping groups
		string mapGrp("");
		if (mapGrpN != -1) {
			mapGrp = curLine[mapGrpN];
		}
		if (mapGrp != "" && CntMapGrps.find(mapGrp) != CntMapGrps.end()) {
			(CntMapGrps)[mapGrp].push_back((int)smplLoc.size());
		}
		else if (mapGrp != "") {
			(CntMapGrps)[mapGrp] = vector<int>(1, (int)smplLoc.size());
		}


		mapGr.push_back(mapGrp);
		//useSmpl.push_back(true);
		smplLoc.push_back(subDir);

	}
	in.close();

	smplN = smplLoc.size();

	//read the gene abundances sample-wise in, primary abundance
	std::cout << "Reading primary abundances from " << smplN << " samples\n";
	SmplOccur currCntMpGr;
	uint smplsIncl(0), smplsInclSup(0);
	mat_fl totalCnts(0);
	for (uint i = preMapSize; i < smplN; i++) {
		if (!mapGrpCnt(currCntMpGr,  CntMapGrps,i, false)) { continue; }
		string abPath = baseP + "/" + smplLoc[i];
#ifdef notRpackage
		cerr << abPath << endl;
#endif
		//read in abundance of sample
		GAs.push_back(new GeneAbundance(abPath, getPath2Ab(opts),false));
		totalCnts += GAs.back()->getTotalCnts();
		smplsIncl++;
	}
	std::cout << "Read primary abundances from " << smplsIncl << " samples, with "<< totalCnts<<" counts\n";// , using "<< geneN<<" genes\n";

	if (calcSupplCov) {
		std::cout << "\n\nScanning for supplementary abundances in  " << smplsIncl << " samples\n";// , using "<< geneN<<" genes\n";
		//read the gene abundances sample-wise in, supplementary abundance
		SmplOccur currCntMpGrSup;
		for (uint i = preMapSize; i < smplN; i++) {
			if (!mapGrpCnt(currCntMpGrSup,  CntMapGrps, i, true)) { continue; }

			string abPath = baseP + "/" + smplLoc[i];
			//read in abundance of sample
			GAsSup.push_back(new GeneAbundance(abPath, getPath2Ab(opts, true), true));
			if (!(GAsSup.back())->isEmpty()) {
#ifdef notRpackage
				cerr << abPath << endl;
#endif
				smplsInclSup++;
			}
		}
		if (smplsInclSup) {
			std::cout << "Read supplementary abundances from " << smplsInclSup << " samples\n";// , using "<< geneN<<" genes\n";
		}
	}



}

bool ClStr2Mat::mapGrpCnt(SmplOccur& currCntMpGr,  
	SmplOccurMult& CntMapGrps, int i, bool isSuppl) {
	if (mapGr[i] != "") {//only include last sample of mapping group..
		SmplOccurIT cMGcnts = currCntMpGr.find(mapGr[i]);
		if (cMGcnts == currCntMpGr.end()) {
			currCntMpGr[mapGr[i]] = 1;
		}
		else { (*cMGcnts).second++; }
		if ((CntMapGrps)[mapGr[i]].size() != (uint)currCntMpGr[mapGr[i]]) {
			if (isSuppl) {
				GAsSup.push_back(new GeneAbundance("", "", isSuppl));
			} else {
				GAs.push_back(new GeneAbundance("", "", isSuppl));
			}
			return false;
		}
	}
	return true;

}


string ClStr2Mat::getPath2Ab(options* opts, bool suppl) {
	string pa2ab(""); 
	//bool calcCoverage = opts->calcCoverage;
	//bool calcCovMedian = opts->calcCovMedian;
	if (suppl) {
		if (opts->calcCoverage) { pa2ab = path2abundanceSup; }
		else if (opts->calcCovMedian) { pa2ab = path2mediABSup; }
		else { pa2ab = path2countsSup; }
	}else {
		if (opts->calcCoverage) { pa2ab = path2abundance; }
		else if (opts->calcCovMedian) { pa2ab = path2mediAB; }
		else { pa2ab = path2counts; }
	}
	return pa2ab;
}
void ClStr2Mat::sealMap() {// only required for mapping group.. don't use
	/*SmplOccurMult nSmpls;
	for (auto smNum : smpls) {
		const vector<int>& smplLocs ( smNum.second);
		for (size_t jj = 0; jj < smplLocs.size(); jj++) { //the loop takes account for multiple samples being grouped together in map
			int idxM = smplLocs[jj];
			auto map2f = nSmpls.find(smNum.first);
			
			if (useSmpl[idxM]) {
				if (map2f == nSmpls.end()) {
					nSmpls[smNum.first] = vector<int>(1, idxM);
				}	else {
					nSmpls[smNum.first].push_back(idxM);
				}
			}
			//else {
			//	nSmpls[smNum.first] = vector<int>(0);
			//}
		}
	}
	smpls = nSmpls;
	*/
}

///////////////////////////////////////////////////
void ContigCrossHit::addHit(int Smpl, int Ctg) {
	(void)Smpl;
	(void)Ctg;
}


///////////////////////////////////////////////////

GeneAbundance::GeneAbundance(const string path, const string abunF, bool isSup):
	isPsAss(false), isSuppl(isSup), isEmp(false){
	if (path == "" && abunF == "") { isEmp = true; return; }//not required to read this (non-existant) file
	//FILE* in;
	//first test if this is a pseudoassembly
	//in = fopen((path + pseudoAssMarker).c_str(), "r");
	//if (in != NULL) {
//	in = fopen(newS.c_str(), "r");
	if (fileExists((path + pseudoAssMarker).c_str())){
		isPsAss = true; isEmp = true;
		return;
	}
	if (fileExists(path + path2skip)) {
		isEmp = true;
		return;
	}

	//not? then read abundances
	string newS = path + abunF;
	//test if input is gzipped
	string newSgz = newS + ".gz"; 
	if (fileExists(newSgz) && !fileExists(newS)) {
		newS = newSgz;
		//std::cout << "switch";
	}
	if (isSuppl && !fileExists(newS)) {
		isEmp = true;
		return;
	}

	//std::cout << "FILE" << newS << endl << newSgz << endl;

	//isGZfile

	istream* in = openInputStream(newS);
	if (!(*in)) {
#ifdef notRpackage
		cerr << "Can't open abundance infile " << newS << endl; std::exit(99);
#endif
	}

	if (in == NULL) {
		cerr << "Couldn't open gene abundance file " << newS << endl;
		exit(36);
	}
	// char buf[400];
	/*while (fgets(buf, sizeof buf, in) != NULL) {
		smat_fl abu; char gene[300];
		bool wrk = sscanf(buf, "%s\t%f",  gene, &abu);
		GeneAbu[gene] = abu;
		//string line(buf);
		*/
	string line;
	size_t lineNo = 0;
	bool hasPositiveAbundance = false;
	while (getline((*in), line, '\n')) {
		lineNo++;
		if (line.empty() || line[0] == '#') { continue; }
		size_t pos = line.find("\t");
		if (pos == string::npos || pos == 0 || pos + 1 >= line.size()) {
			closeInputStream(in);
			cerr2("Malformed gene abundance row at " + newS + ":" + to_string(lineNo) + "\n", 36);
		}
		string gene = line.substr(0, pos);
		string valueText = line.substr(pos + 1);
		char* end = nullptr;
		errno = 0;
		double value = strtod(valueText.c_str(), &end);
		while (end != nullptr && *end != '\0' && isspace(static_cast<unsigned char>(*end))) { end++; }
		if (end == valueText.c_str() || end == nullptr || *end != '\0' || errno == ERANGE || !isfinite(value) || value < 0 || value > numeric_limits<mat_fl>::max()) {
			closeInputStream(in);
			cerr2("Invalid gene abundance at " + newS + ":" + to_string(lineNo) + "\n", 36);
		}
		if (GeneAbu.find(gene) != GeneAbu.end()) {
			closeInputStream(in);
			cerr2("Duplicate gene abundance ID '" + gene + "' in " + newS + "\n", 36);
		}
		GeneAbu[gene] = static_cast<mat_fl>(value);
		if (value > 0.0) { hasPositiveAbundance = true; }
		
	}
	closeInputStream(in);
	if (isSuppl && !hasPositiveAbundance) {
		GeneAbu.clear();
		isEmp = true;
	}

	//fclose(in);
}
mat_fl GeneAbundance::getAbundance(const string x) {
	if (isPsAss) {
		return (mat_fl)1;//return one read count
	}
	SmplAbunIT fnd = GeneAbu.find(x);
	if (fnd == GeneAbu.end()) {
		return (mat_fl)0;

cerr << "Can't find " << x << endl;
exit(33);
	}
	return (*fnd).second;
}
