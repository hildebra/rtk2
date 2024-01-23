#include "ClStr2Mat.h"



void printVec(clusWrk * curClus, ostream* mO, ofstream* gN) {//, const vector<bool>& useSmpl) {
	long CLidx = curClus->Clnum;
	//string outStr;outStr.reserve(4000);outStr = to_string(CLidx);
	//important to have this after outStr creation!
	//(*gN) << outStr + curClus->geneNamesStr;
	(*gN) << CLidx<< curClus->geneNamesStr;
	(*mO) << CLidx << curClus->matStr;// CLidx;

	/*const vector<smat_fl>& pr = curClus->Vec;
	for (size_t i = 0; i < pr.size(); i++) {
		//if (!useSmpl[i]) { continue; }
		if (pr[i] == 0) {
			(*mO) << "\t0"; //outStr += 
		}
		else {
			(*mO) << "\t" << pr[i]; //outStr += 
		}
		//SmplSum[i] += pr[i];
	}
	//outStr += "\n"; (*mO)<< outStr
	(*mO) << endl;
	*/
	delete curClus;
}
clusWrk* ClStr2Mat::workClusBlock(textBlock* inVs, 
	 long CLidx) {
	clusWrk* ret = new clusWrk(smplN, CLidx);
	bool repFound = false;
	for (uint i = 0; i < inVs->txt.size(); i++) {
		string &line = inVs->txt[i];
		if (i == 0) {//new cluster, add headerInfo
			ret->geneNamesStr = "\t" + line;
			continue;
		}


		//1 get gene, sample (deparse)
		size_t pos = line.find("nt, >");
		size_t pos2 = line.find("...", pos + 4);
		string gene = line.substr(pos + 5, pos2 - pos - 5);

		if (!repFound && line.back() == '*') {//report representative gene
			ret->geneNamesStr += "\t" + gene + "\n";
			repFound = true;
		}

		//bool geneInAssembl(true);
		pos = gene.find(sampleSeq);
		string sample;
		//SmplOccurITmult smNum;
		pos2 = gene.find("_L", pos + 3);
		if (pos != string::npos && pos2 != string::npos) { //has the characteristic "__" sample separator
			sample = gene.substr(0, pos);
			auto smNum = smpls.find(sample);
			if (smNum == smpls.end()) {
#ifdef notRpackage
				cerr << "incorrect sample name: " << sample << "::" << gene << endl;
				continue;
				//exit(55);
#endif
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
		}
		else {
			//geneInAssembl = false;
		}

	}


	if (!repFound) {
		ret->geneNamesStr += "\t?\n";
		cerr << "No RepSeq found for " << ret->geneNamesStr << " -1" << endl;
	}
	delete inVs;

	this->addSums(ret);

	ret->generateMatStr();


	return (ret);

}







ClStr2Mat::ClStr2Mat(options* opts):
	lastClIdWr(1),GAs(0), CCH(NULL),smplLoc(0), //baseP(0), 
	mapGr(0), SmplSum(0), smplN(0), curr(-1),
	lastline(""), sampleSeq("__"){
	ifstream incl2;
	//FILE* incl;
	const string inF = opts->input;
	const string outF = opts->output;
	const string mapF = opts->map;
	//const string basePX = opts->referenceDir;
	const bool doGZ = opts->gzOut;
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
	if (doGZ) {
#ifdef _gzipread
		matO = new ogzstream(outF2.c_str(), ios::out);
		cout << "Writing gzip'd matrix\n";
#else
		cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
	}
	else {
		matO = new ofstream(outF2);
	}

	
	if (!(*matO)) {
		cerr2 ("Couldn't open matrix output " +outF + ".mat" + matCmpr +"\n",57);
	}
	geneNames = new ofstream(outF + ".genes2rows.txt");
	if (!(*geneNames)) {
		cerr2( "Couldn't open report file " + outF + ".genes2rows.txt" +"\n", 56);
	}
	


	//read map(s) and check that
	stringstream ss2(mapF);  string segments;
	while (getline(ss2, segments, ',')) {
		string baseP = ""; //SmplOccurMult CntMapGrps ;
		read_map(segments, opts->calcCoverage, opts->calcCovMedian, opts->oldMapStyle, baseP);
		//read_abundances(CntMapGrps, baseP, opts->calcCoverage, opts->calcCovMedian, opts->oldMapStyle);
		
	}


	//smplnames in out matrix
	vector<string> SmplNmsVec(smplN, "");
	for (auto it = smpls.begin(); it != smpls.end(); it++) {
		vector<int> XX = (*it).second;
		if (XX.size() == 0) { continue; }
		SmplNmsVec[ XX[XX.size()-1] ] = smplRid[(*it).first];
	}
	cout << "Calculating abundance matrix on "<< smplN <<" samples (this may take awhile)\n";
	(*matO) << "Genes";
	for (size_t i = 0; i < SmplNmsVec.size(); i++) {
		//if (!useSmpl[i]) { continue; }
		(*matO) << "\t" << SmplNmsVec[i];
	}
	(*matO) << endl;

	(*geneNames) << "#GID\tCluster\tRepSeq\n";
	CCH = new ContigCrossHit((int)smplN);
	CCH->setSmplNms(SmplNmsVec);
	//SparseMatrix * mat = new SparseMatrix();
	long CLidx = 2; //start at number 2
	//const string sampleStrSep = "__";
	sampleSeq = "__";
	//vector<smat_fl> SmplSum(smplN, 0.f);
	SmplSum.resize(smplN, 0.f);
	
	int numthr = opts->threads-1;// -1 to include wrThr
	if (numthr <= 0) { numthr = 1; }

	//thread pool
	vector < job > slots(numthr);
	


	//bool repFound = false;
	textBlock* clbl;
	string lastline = "";
	//std::future <textBlock*> readThr;
	//readThr = async(std::launch::async, getClusBlock, incl, lastline);
	//readThr = async(std::launch::async, test);

	int j = 0;
	int tCnt = 0;

	while (true) {
		clbl = getClusBlock(incl2, lastline);
		tCnt++;
		//clbl = readThr.get();
		//lastline = clbl->lastLine;
		//readThr = async(std::launch::async, getClusBlock, &incl, lastline);

		
		//clusWrk* curClus = workClusBlock(clbl, smplN, CLidx,sampleStrSep,GAs, &smpls);

		//search all threads for empty slot to push job into..
		while(true) {
			if (j >= numthr) { j = 0; }
			if (slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
				slots[j].inUse = false;
				clusWrk* curClus = slots[j].fut.get();
				//this->addSums(curClus);
				manage_write(curClus);
				//printVec( curClus, matO, geneNames);
			}
			if (slots[j].inUse == false) {
				slots[j].fut = async(std::launch::async, &ClStr2Mat::workClusBlock, this, clbl, CLidx);
				CLidx++;
				slots[j].inUse = true;
				break;//job submitted, break search for empty thread
			}
			j++;
		}
		if (!clbl->cont) { break; }
		
	}
	for (j = 0; j < numthr; j++) {//collect remaining jobs
		if (slots[j].inUse == true) {
			clusWrk* curClus = slots[j].fut.get();
			//this->addSums(curClus);
			manage_write(curClus);
			//printVec(curClus, matO, geneNames);
			slots[j].inUse = false;
			CLidx++;
		}
	}
	//wrThr.fut.get();
	finish_write();
	finish_write();

	cout << "Read/wrote " << tCnt << " total genes in matrix\n";
	

	incl2.close();
	//std::fclose(incl);
	//matO->close(); 
	geneNames->close();
	//matO->close();
	delete matO; delete geneNames;

	ofstream* matS = new ofstream((outF + ".mat.sum"), ofstream::out);
	
	//matS.open(outF + ".mat.sum", ofstream::out);
	//print sample sum
	for (size_t i = 0; i < SmplNmsVec.size(); i++) {
		//if (!useSmpl[i]) { continue; }
		//cout << SmplNmsVec[i] << "\t" << SmplSum[i] << endl;
		(*matS) << SmplNmsVec[i] << "\t" << SmplSum[i] << endl;
	}
	matS->close();  delete matS;
}
ClStr2Mat::~ClStr2Mat() {
	for (size_t i = 0; i < GAs.size(); i++) { delete GAs[i]; }
	delete CCH;
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

void ClStr2Mat::read_map(const string mapF, bool calcCoverage, bool calcCovMedian,  bool oldFolderStructure, string& baseP) {
	ifstream in;
	int map2folderIdx = 0; if (oldFolderStructure) { map2folderIdx = 1; }
	curr++;//keep track of different maps and inPaths
	uint preMapSize((int)smplLoc.size());

	in.open(mapF.c_str());
	if (!in) {
		cerr << "Couldn't open mapping file " << mapF << endl;
		exit(56);
	}
	cout << "Reading map " << mapF << endl;// " on path " << baseP[curr] << endl;
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
						cout << "Found Assembly groups in map\n";
					}
					if (segments == "MapGrps") {
						mapGrpN = sbcnt;
						//fillMapGrp = true;
						cout << "Found Mapping groups in map\n";
					}
					if (segments == "ExcludeAssembly") {
						skSmplCol = sbcnt;
						cout << "Samples can be excluded from assembly\n";
					}
				}
			}
			while (getline(ss, segments, '\t')) {
				if (segments == "#OutPath") {
					getline(ss, segments, '\t'); baseP1 = segments;
				}
				if (segments == "#RunID") {
					getline(ss, segments, '\t'); baseP2 = segments;
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
			curLine.push_back(segments);
		}
		if (skSmplCol > -1 && curLine[skSmplCol] == "1") { continue; }
		if (assGrpN >= 0 && int(curLine.size()) <= assGrpN) {
			curLine.resize(assGrpN + 1, "");
		}


		string smpID = curLine[0];


		//getline(ss, segments, '\t');
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

	//return CntMapGrps;

//}
//void ClStr2Mat::read_abundances(SmplOccurMult& CntMapGrps, string baseP, bool calcCoverage,
	//bool calcCovMedian, bool oldFolderStructure){
	smplN = smplLoc.size();
	//uint geneN = 0;
	//uint preMapSize((int)smplLoc.size());

	cout << "Reading abundances from " << smplN << " samples\n";
	//read the gene abundances sample-wise in
	SmplOccur currCntMpGr;
	uint smplsIncl = 0;
	for (uint i = preMapSize; i < smplN; i++) {
		string pa2ab = path2counts;
		
		if (calcCoverage) { pa2ab = path2abundance; }
		if (calcCovMedian) { pa2ab = path2mediAB; }

		//only include last sample of mapping group..
		if (mapGr[i] != "") {
			SmplOccurIT cMGcnts = currCntMpGr.find(mapGr[i]);
			if (cMGcnts == currCntMpGr.end()) {
				currCntMpGr[mapGr[i]] = 1;
			}
			else {
				(*cMGcnts).second++;
			}
			if ((CntMapGrps)[mapGr[i]].size() != (uint) currCntMpGr[mapGr[i]]) {
				GAs.push_back(new GeneAbundance("",""));
				//useSmpl[i] = false;
				continue;
			}
		}

		//DEBUG
		//continue;

		//read in abundance of sample
#ifdef notRpackage
		cerr << baseP + "/" + smplLoc[i] << endl;
#endif
		smplsIncl++;
		GAs.push_back(new GeneAbundance(baseP + "/" + smplLoc[i], pa2ab));
	}
	cout << "Including " << smplsIncl << " samples\n";// , using "<< geneN<<" genes\n";
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

}


///////////////////////////////////////////////////

GeneAbundance::GeneAbundance(const string path, const string abunF):
	isPsAss(false){
	if (path == "" && abunF == "") { return; }//not required to read this (non-existant) file
	FILE* in;
	//first test if this is a pseudoassembly
	in = fopen((path + pseudoAssMarker).c_str(), "r");
	if (in != NULL) {
		fclose(in);
		isPsAss = true;
		return;
	}
	//not? then read abundances
	string newS = path + abunF;
	in = fopen(newS.c_str(), "r");
	if (in == NULL) {
		cerr << "Couldn't open gene abundance file " << newS << endl;
		exit(36);
	}
	char buf[400];
	while (fgets(buf, sizeof buf, in) != NULL) {
		//buf[strcspn(buf, "\n")] = 0;
		smat_fl abu; char gene[300];
		bool wrk = sscanf(buf, "%s\t%f",  gene, &abu);
		GeneAbu[gene] = abu;
		//string line(buf);
		//size_t pos = line.find("\t");
		//string gene = line.substr(0, pos);
		//GeneAbu[gene] = (smat_fl)atof(line.substr(pos + 1).c_str());
		
	}
	fclose(in);
}
smat_fl GeneAbundance::getAbundance(const string x) {
	if (isPsAss) {
		return (smat_fl) 1.f;//return one read count
	}
	SmplAbunIT fnd = GeneAbu.find(x);
	if (fnd == GeneAbu.end()) {
		return (smat_fl) 0;

cerr << "Can't find " << x << endl;
exit(33);
	}
	return (*fnd).second;
}
