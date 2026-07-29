#include "Rare.h"


//2.01: 4.4.23: changed output for gene matrix calc
//2.07: 26.9.23: (-v) version print
//2.08: parallelized clstr2mat output, more error tolerant
//2.09: 7.12.23: added new delim flags for functions (funcHAnnoAND ..)
//2.10: 23.12.23: reworked "higher lvl" matrix calculations, adding ability to measure abundance by "mean/median" occurrence at lowest level but summing then for higher levels
//2.11: 24.4.24: `geneMat` subcommand can now accept .gz input coverage files
//2.12: 30.4.24: fix open .gz in geneMat mode
//2.13: 2.5.24: print suppl coverage of samples
//2.14: 13.1.26: just some additional counts for geneMatrix mode
//2.15: 25.4.26: extra safety for extractRows function
//2.16: 1.5.26: more tech hardening via AI, no real bug found, just stability
//2.17: 20.7.26: implement mergeMat and median aggregation; correct diversity estimators,
//                 deterministic ordering, wide-matrix handling, 64-bit global counters, parsing, and concurrency bugs
//2.18: 20.7.26: harden sumMat hierarchy/output semantics and geneMat parsing, precision,
//                 supplementary-column handling, and bounded deterministic worker output
//2.19: 24.7.26: centralize threaded workers via shared thread pool, parallel gzip output writes,
//                 and cap malformed CD-HIT warning spam with end-of-run summary
//2.20: 24.7.26: centralize and parallelize input reading
//2.21: 24.7.26: add program timers for key matrix stages (read/process/write)
//2.22: 24.7.26: add RTK2 gzip member index sidecar (.gzi) for indexed parallel gzip reads,
//                 and format timer display as HH:MM:SS
//2.23: 26.7.26: write standard binary GZI offsets and use 16 MiB parallel gzip members
//2.24: 26.7.26: harden parallel gzip finalization and normalization shutdown
//2.25: 26.7.26: normalize sparse matrices in memory by default; retain streaming mode by flag
//2.26: 27.7.26: scale parallel gzip compression and its bounded queue with -t
//2.27: 27.7.26: apply -t to indexed parallel gzip input decompression
//2.28: 29.7.26: add mosaic detection using reciprocal PAF relationships
//                 and MGS-normalized abundance complementation
const char* rar_ver="2.28";

static void reportRtkFinishTimer() {
	printProgramTimerStep("rtk2 finish");
}


rareStruct* calcDivRar(int i, Matrix* Mo, DivEsts* div, options* opts,
		global_count_cube* abundInRow, global_count_cube* occuencesInRow){

    smplVec* cur        = Mo->getSampleVec(i);
    string curS         = Mo->getSampleName(i);
    div->SampleName     = curS;
    std::vector<vector<uint>> cnts;
    vector<vector<rare_map>> cntsMap(opts->depth.size());
    vector<string> cntsName(opts->depth.size());
    string skippedNames;
    bool wrAtAll(opts->write > 0);

    cur->rarefy(opts->depth, opts->output, opts->repeats,
            div, cntsMap, cntsName, skippedNames, abundInRow, occuencesInRow,
            opts->write, false, wrAtAll);

    // put everything in our nice return container
    rareStruct* tmpRS       = new rareStruct();
    tmpRS->div              = div;
    tmpRS->cnts             = cntsMap;
    tmpRS->cntsName         = cntsName;
    tmpRS->skippedNames     = skippedNames;
    tmpRS->i                = i;

    delete cur;
    return tmpRS;
}


rareStruct* calcDivRarVec(int i, vector<string> fileNames, DivEsts* div, options* opts,
		global_count_cube* abundInRow, global_count_cube* occuencesInRow){
    smplVec* cur = new smplVec(fileNames[i],4);

    std::vector<vector<uint>> cnts;
    vector<vector<rare_map>> cntsMap(opts->depth.size());
    vector<string> cntsName(opts->depth.size());
    string skippedNames;
    bool wrAtAll(opts->write > 0);
    cur->rarefy(opts->depth, opts->output, opts->repeats,
            div, cntsMap, cntsName, skippedNames, abundInRow, occuencesInRow,
            opts->write, false, wrAtAll);

    rareStruct* tmpRS       = new rareStruct();
    tmpRS->div              = div;
    tmpRS->cnts             = cntsMap;
    tmpRS->cntsName         = cntsName;
    tmpRS->skippedNames     = skippedNames;
    tmpRS->IDs              = cur->getRowNames();
    tmpRS->i                = i;

    delete cur;
    if( remove( fileNames[i].c_str() ) != 0 ){
        cerr2("Swap Mode: Error deleting file: " + fileNames[i] +"\n");
    }
    return tmpRS;
}



void helpMsglegacy(){
    string  AvailableModes = "Available run modes:\nnormalize\nsplitMat\nlineExtr\nmergeMat\nsumMat\nrarefaction\nrare_inmat\nmodule\n";
    cerr << AvailableModes << "Provide two arguments\nexiting..\n";
    cerr << "------------------------------\nAuthor: falk.hildebrand@gmail.com\n";
    std::exit(2);
}

void stateVersion(){
    printf("rarefaction tool kit (rtk) %s\n", rar_ver);
	std::cout << "compiled on: " << __DATE__ << " at " << __TIME__ << endl;
	std::cout << "compiled with: " << __cplusplus << endl;
}

void versionMsg() {
	stateVersion();
	exit(0);
}

void helpMsg(){
	stateVersion();
	printf("\n");
	printf("USAGE\n");
	printf("    rtk <mode> -i <input.csv> -o <output> [options] \n");
	printf("\n");
	printf("MODE rarefaction\n");
	printf("\n");
	printf("OPTIONS\n");

	printf("<mode>      For rarefaction: mode can be either swap or memory.\n");
	printf("            Swap mode creates temporary files but uses less memory. \n");
	printf("            The speed of both modes is comparable.\n");
	printf("    -i      path to an .txt file (tab delimited) to rarefy\n");
	printf("    -o      path to a output directory\n");
	printf("    -d      Depth or multiple comma seperated depths to rarefy to. Default is 0.95 times the minimal column sum.\n");
	printf("    -r      Number of times to create diversity measures. Default is 10.\n");
	printf("    -w      Number of rarefied tables to write.\n");
	printf("    -t      Number of threads to use. Default: 1\n");
	printf("    -ns     If set, no temporary files will be used when writing rarefaction tables to disk.\n");
	printf("\n");
	printf("EXAMPLE\n");
	printf("    Rarefy a table to 1000 counts per sample with two threads. Create one table:\n");
	printf("        rtk swap -i table.csv -o outputdir/prefix. -d 1000 -r 10 -w 1 -t 2\n");
	printf("\n");
	printf("    Rarefy with most memory and least amount of IO:\n");
	printf("        rtk memory -i table.csv -o outputdir/prefix. -ns\n");
	printf("\n");
	printf("MODE: Colsums\n");
	printf("Reports the column sums of all columns in form of a sorted and an unsorted file.\n");
	printf("\n");
	printf("EXAMPLE\n");
	printf("    Repot column sums of file 'table.csv'\n");
	printf("        rtk colsums -i table.csv -o prefix\n");
	printf("\n");
	printf("MODE: mergeMat\n");
	printf("Outer-joins two or more matrices by the row ID in their first column and appends sample columns.\n");
	printf("Input rows missing from a matrix are filled with zero; row order follows first appearance.\n");
	printf("\n");
	printf("EXAMPLE\n");
	printf("        rtk mergeMat -i first.tsv,second.tsv,third.tsv -o merged.tsv\n");
	printf("Repeated -i options may be used instead of a comma-separated list.\n");
	printf("\nMODE: mosaic\n");
	printf("Reports sequence-related gene pairs in one MGS whose normalized abundances\n");
	printf("complement one another to approximately one MGS copy. No matrix is modified.\n");
	printf("    -i          gene abundance matrix (rows are genes; columns are samples)\n");
	printf("    -reference  MGS/core file, including six-column SB.clusters.core format\n");
	printf("    -paf        minimap2 PAF for the catalog (use minimap2 -c for CIGAR evidence)\n");
	printf("    -groups     optional gene-to-NOG/group mapping\n");
	printf("    -o          output prefix for .mosaic.tsv, .concat.list, and summary\n");
	printf("Thresholds: -mosaicMinIdentity 0.90, -mosaicMinQueryCoverage 0.80,\n");
	printf("    -mosaicMinTargetCoverage 0.80, -mosaicMinScoreMargin 0.02,\n");
	printf("    -mosaicMinLengthRatio 0.80, -mosaicMinSamples 5,\n");
	printf("    -mosaicMinVariantSamples 2, -mosaicMinBaselineGenes 5,\n");
	printf("    -mosaicDetectionFraction 0.10, -mosaicMinMergedCopy 0.70,\n");
	printf("    -mosaicMaxMergedCopy 1.30, -mosaicMaxMergedError 0.30,\n");
	printf("    -mosaicMinFitGain 0.15, -mosaicMaxSampleOverlapFraction 0.15,\n");
	printf("    -mosaicMaxAbundanceOverlapFraction 0.15; -mosaicRequireCigar requires cg:Z.\n");
	printf("    SB core markers are preferred for the MGS baseline; -mosaicAllCoreBaseline disables this.\n");
	printf("EXAMPLE\n");
	printf("    rtk2 mosaic -i genes.mat.gz -reference SB.clusters.core -paf catalog.self.paf.gz -groups gene2nog.tsv -o run\n");
	printf("\n");
	printf("MODE: normalize\n");
	printf("Loads the matrix into sparse memory, normalizes columns, and writes the result.\n");
	printf("    -normalizeReadPass  Use the lower-memory two-pass streaming implementation instead.\n");

	printf("\n");

	std::exit(0);
}




vector<double> parseDepths(const string& a){
    std::vector<double> vect;
    std::stringstream ss(a);
    string token;
    while (getline(ss, token, ',')) {
        if (token.empty()) {
            throw invalid_argument("empty rarefaction depth");
        }
        size_t used = 0;
        double value = stod(token, &used);
		if (used != token.size() || !isfinite(value) || value <= 0.0 ||
			value >= static_cast<double>(numeric_limits<rtk_count>::max())) {
            throw invalid_argument("invalid rarefaction depth: " + token);
        }
        vect.push_back(value);
    }
    if (vect.empty()) { throw invalid_argument("no rarefaction depth provided"); }
    return vect;
}

static unsigned int parseUnsignedOption(const char* text, const string& flag, bool allowZero) {
    if (text == nullptr || *text == '\0' || *text == '-') {
        throw invalid_argument("invalid value for " + flag);
    }
    size_t used = 0;
    unsigned long long value = stoull(text, &used);
    if (used != strlen(text) || value > numeric_limits<unsigned int>::max() || (!allowZero && value == 0)) {
        throw invalid_argument("invalid value for " + flag + ": " + text);
    }
    return static_cast<unsigned int>(value);
}

static double parseNumericOption(const char* text, const string& flag) {
    if (text == nullptr || *text == '\0') { throw invalid_argument("missing value for " + flag); }
    size_t used = 0;
    double value = stod(text, &used);
    if (used != strlen(text) || !isfinite(value)) {
        throw invalid_argument("invalid value for " + flag + ": " + text);
    }
    return value;
}

void cerr2(const std::string x, int ex ) {
#ifdef notRpackage
	cerr << x;
	if (ex != 0) { exit(ex); }
#else
	std::cout << x;
	if (ex != 0) { exit(ex); }
#endif

}

options::options(int argc, char** argv) :input(""), output(""), mode(""),
	referenceDir(""), referenceFile(""),
	depth(), depthMin(0), sepChar('\t'), repeats(10), write(0), threads(1), writeSwap(true), verbose(false),
	oldMapStyle(false), sparse(false), pval(1e-5f),
	modDB(""), modRedund(5), modEnzCompl(0.5f), modModCompl(0.5f), modWrXtraInfo(false),
	modCollapse(false), calcCoverage(false), calcCovMedian(false), mean(false),
	highFromLow(false), hieraSrtDown(false),
	median(false), check4idxMatch(false), gzOut(false), header(true),
	calcSupplCov(false), normalizeReadPass(false),
	modDescr(""), modHiera(""), xtra(""),
	funcHieraSep(";" ), funcAnnoAND(","), funcAnnoOR("|") {


        bool hasErr = false;
        if (argc == 0) { return; }//could be the pseudo mod opt object

        //bool newIDthrs = false; string newIDthStr("");
        mode = argv[1];

		if (mode == "help" || mode == "-help" || mode == "-h" || mode == "--help") {
			helpMsg();
		}else if (mode == "version" || mode == "-version" || mode == "-v" || mode == "--version") {
			versionMsg();
		}

		auto nextValue = [&](int& i, const string& flag) -> const char* {
			if (i + 1 >= argc) { throw invalid_argument("missing value for " + flag); }
			return argv[++i];
		};

		try {
		for (int i = 2; i < argc; i++)
        {
			if (!strcmp(argv[i], "-i")) {
				string value = nextValue(i, "-i");
				if (mode == "mergeMat" && !input.empty()) { input += "," + value; }
				else { input = value; }
			}
			else if (!strcmp(argv[i], "-o"))
				output = nextValue(i, "-o");
			///else if (!strcmp(argv[i], "-m"))
			else if (!strcmp(argv[i], "-d")) {
				// split by any given komma
				depth = parseDepths(nextValue(i, "-d"));
	depthMin = static_cast<rtk_count>(*std::min_element(depth.begin(), depth.end()));
			}
			else if (!strcmp(argv[i], "-r"))
				repeats = parseUnsignedOption(nextValue(i, "-r"), "-r", false);
			else if (!strcmp(argv[i], "-w"))
				write = parseUnsignedOption(nextValue(i, "-w"), "-w", true);
			else if (!strcmp(argv[i], "-t"))
				threads = parseUnsignedOption(nextValue(i, "-t"), "-t", false);
			else if (!strcmp(argv[i], "-ns"))   // no swap
				writeSwap = false;
			else if (!strcmp(argv[i], "-v"))
				verbose = true;
			else if (!strcmp(argv[i], "-calcSupplCov"))
				calcSupplCov = true;
			else if (!strcmp(argv[i], "-hieraSrtDown"))
				hieraSrtDown = true;		
			else if (!strcmp(argv[i], "-s")) {
				string tmp = nextValue(i, "-s");
				if (tmp.length() != 1) { throw invalid_argument("-s argument must be exactly one character"); }
				sepChar = char(tmp[0]);
			}
			else if (!strcmp(argv[i], "-h"))
				helpMsg();
			//geneMat specific args
			else if (!strcmp(argv[i], "-map"))
				map = nextValue(i, "-map");
			else if (!strcmp(argv[i], "-refD"))
				referenceDir = nextValue(i, "-refD");
			else if (!strcmp(argv[i], "-reference"))
				referenceFile = nextValue(i, "-reference");
			else if (!strcmp(argv[i], "-paf"))
				pafFile = nextValue(i, "-paf");
			else if (!strcmp(argv[i], "-groups"))
				groupFile = nextValue(i, "-groups");
			else if (!strcmp(argv[i], "-mosaicMinIdentity"))
				mosaicMinIdentity = parseNumericOption(nextValue(i, "-mosaicMinIdentity"), "-mosaicMinIdentity");
			else if (!strcmp(argv[i], "-mosaicMinQueryCoverage"))
				mosaicMinQueryCoverage = parseNumericOption(nextValue(i, "-mosaicMinQueryCoverage"), "-mosaicMinQueryCoverage");
			else if (!strcmp(argv[i], "-mosaicMinTargetCoverage"))
				mosaicMinTargetCoverage = parseNumericOption(nextValue(i, "-mosaicMinTargetCoverage"), "-mosaicMinTargetCoverage");
			else if (!strcmp(argv[i], "-mosaicMinScoreMargin"))
				mosaicMinScoreMargin = parseNumericOption(nextValue(i, "-mosaicMinScoreMargin"), "-mosaicMinScoreMargin");
			else if (!strcmp(argv[i], "-mosaicMinLengthRatio"))
				mosaicMinLengthRatio = parseNumericOption(nextValue(i, "-mosaicMinLengthRatio"), "-mosaicMinLengthRatio");
			else if (!strcmp(argv[i], "-mosaicMinBaseline"))
				mosaicMinBaseline = parseNumericOption(nextValue(i, "-mosaicMinBaseline"), "-mosaicMinBaseline");
			else if (!strcmp(argv[i], "-mosaicDetectionFraction"))
				mosaicDetectionFraction = parseNumericOption(nextValue(i, "-mosaicDetectionFraction"), "-mosaicDetectionFraction");
			else if (!strcmp(argv[i], "-mosaicMinMergedCopy"))
				mosaicMinMergedCopy = parseNumericOption(nextValue(i, "-mosaicMinMergedCopy"), "-mosaicMinMergedCopy");
			else if (!strcmp(argv[i], "-mosaicMaxMergedCopy"))
				mosaicMaxMergedCopy = parseNumericOption(nextValue(i, "-mosaicMaxMergedCopy"), "-mosaicMaxMergedCopy");
			else if (!strcmp(argv[i], "-mosaicMaxMergedError"))
				mosaicMaxMergedError = parseNumericOption(nextValue(i, "-mosaicMaxMergedError"), "-mosaicMaxMergedError");
			else if (!strcmp(argv[i], "-mosaicMinFitGain"))
				mosaicMinFitGain = parseNumericOption(nextValue(i, "-mosaicMinFitGain"), "-mosaicMinFitGain");
			else if (!strcmp(argv[i], "-mosaicMaxSampleOverlapFraction"))
				mosaicMaxSampleOverlapFraction = parseNumericOption(nextValue(i, "-mosaicMaxSampleOverlapFraction"), "-mosaicMaxSampleOverlapFraction");
			else if (!strcmp(argv[i], "-mosaicMaxAbundanceOverlapFraction"))
				mosaicMaxAbundanceOverlapFraction = parseNumericOption(nextValue(i, "-mosaicMaxAbundanceOverlapFraction"), "-mosaicMaxAbundanceOverlapFraction");
			else if (!strcmp(argv[i], "-mosaicMinSamples"))
				mosaicMinSamples = parseUnsignedOption(nextValue(i, "-mosaicMinSamples"), "-mosaicMinSamples", false);
			else if (!strcmp(argv[i], "-mosaicMinVariantSamples"))
				mosaicMinVariantSamples = parseUnsignedOption(nextValue(i, "-mosaicMinVariantSamples"), "-mosaicMinVariantSamples", false);
			else if (!strcmp(argv[i], "-mosaicMinBaselineGenes"))
				mosaicMinBaselineGenes = parseUnsignedOption(nextValue(i, "-mosaicMinBaselineGenes"), "-mosaicMinBaselineGenes", false);
			else if (!strcmp(argv[i], "-mosaicAllCoreBaseline"))
				mosaicPreferMarkerBaseline = false;
			else if (!strcmp(argv[i], "-mosaicRequireCigar"))
				mosaicRequireCigar = true;
			//module specific args
			else if (!strcmp(argv[i], "-refMods"))
				modDB = nextValue(i, "-refMods");
			else if (!strcmp(argv[i], "-redundancy"))
				modRedund = static_cast<int>(parseUnsignedOption(nextValue(i, "-redundancy"), "-redundancy", true));
			else if (!strcmp(argv[i], "-enzymeCompl"))
				modEnzCompl = static_cast<float>(parseNumericOption(nextValue(i, "-enzymeCompl"), "-enzymeCompl"));
			else if (!strcmp(argv[i], "-pval"))
				pval = static_cast<float>(parseNumericOption(nextValue(i, "-pval"), "-pval"));
			else if (!strcmp(argv[i], "-moduleCompl"))
				modModCompl = static_cast<float>(parseNumericOption(nextValue(i, "-moduleCompl"), "-moduleCompl"));
			else if (!strcmp(argv[i], "-oldMapStyle"))
				oldMapStyle = true;
			else if (!strcmp(argv[i], "-writeExtraModEstimates"))
				modWrXtraInfo = true;
			else if (!strcmp(argv[i], "-checkRowName2Idx"))
				check4idxMatch = true;
			else if (!strcmp(argv[i], "-gz"))
				gzOut = true;
			else if (!strcmp(argv[i], "-normalizeReadPass") || !strcmp(argv[i], "-normReadPass"))
				normalizeReadPass = true;
			else if (!strcmp(argv[i], "--no_header"))
				header = false;
			else if (!strcmp(argv[i], "-no_header"))
				header = false;
			else if (!strcmp(argv[i], "-filterOccPerSmpl"))
				occPerSmpl = parseUnsignedOption(nextValue(i, "-filterOccPerSmpl"), "-filterOccPerSmpl", true);
			else if (!strcmp(argv[i], "-filterMinOcc"))
				occMin = parseUnsignedOption(nextValue(i, "-filterMinOcc"), "-filterMinOcc", true);
			else if (!strcmp(argv[i], "-extHiera"))
				extendHierachy = true;
			else if (!strcmp(argv[i], "-mean"))
				mean = true;
			
			else if (!strcmp(argv[i], "-highFromLow"))
				highFromLow = true;
			else if (!strcmp(argv[i], "-median"))
				median = true;
			else if (!strcmp(argv[i], "-collapseDblModules"))
                modCollapse = true;
            else if (!strcmp(argv[i], "-useCoverage"))//for gene catalog, default is counts
                calcCoverage = true;
			else if (!strcmp(argv[i], "-useCovMedian"))//for gene catalog, default is counts
				calcCovMedian = true;
			else if (!strcmp(argv[i], "-description")) //description of single modules
				modDescr = nextValue(i, "-description");
			else if (!strcmp(argv[i], "-hiera")) // hierachy for modules
				modHiera = nextValue(i, "-hiera");
			else if (!strcmp(argv[i], "-xtra"))
				xtra = nextValue(i, "-xtra");
			else if (!strcmp(argv[i], "-funcHieraSep")) { // |
				funcHieraSep = nextValue(i, "-funcHieraSep");
				//funcHieraSep = (tmp.c_str());
			}else if (!strcmp(argv[i], "-funcHAnnoAND")) {// ,
				funcAnnoAND = nextValue(i, "-funcHAnnoAND");
				//funcHAnnoAND = char(tmp.c_str());
			}else if (!strcmp(argv[i], "-funcAnnoOR")) { // ;
				funcAnnoOR = nextValue(i, "-funcAnnoOR");
				//funcAnnoOR = char(tmp.c_str());
			} else {
				throw invalid_argument("unknown option: " + string(argv[i]));
			}
        }
		} catch (const exception& e) {
			cerr << "Option error: " << e.what() << "\nUse \"rtk -h\" for help.\n";
			exit(98);
		}

		if (calcCovMedian && calcCoverage) {
			std::cout << "Can not calc median and coverage at the same time, please use only one option\n(-useCoverage && -useCovMedian)\n";
			exit(734);
		}
		if (mean && median) {
			cerr << "Options -mean and -median are mutually exclusive.\n";
			exit(734);
		}
		if (modEnzCompl < 0.f || modEnzCompl > 1.f || modModCompl < 0.f || modModCompl > 1.f || pval < 0.f || pval > 1.f) {
			cerr << "Completeness and p-value options must be between 0 and 1.\n";
			exit(734);
		}
		if (mosaicMinIdentity < 0.0 || mosaicMinIdentity > 1.0 ||
			mosaicMinQueryCoverage < 0.0 || mosaicMinQueryCoverage > 1.0 ||
			mosaicMinTargetCoverage < 0.0 || mosaicMinTargetCoverage > 1.0 ||
			mosaicMinLengthRatio < 0.0 || mosaicMinLengthRatio > 1.0 ||
			mosaicMinScoreMargin < 0.0 || mosaicMinScoreMargin > 1.0 ||
			mosaicDetectionFraction < 0.0 || mosaicDetectionFraction > 1.0 ||
			mosaicMaxMergedError < 0.0 || mosaicMinFitGain < 0.0 ||
			mosaicMaxSampleOverlapFraction < 0.0 || mosaicMaxSampleOverlapFraction > 1.0 ||
			mosaicMaxAbundanceOverlapFraction < 0.0 || mosaicMaxAbundanceOverlapFraction > 1.0 ||
			mosaicMinBaseline < 0.0 || mosaicMinMergedCopy < 0.0 ||
			mosaicMaxMergedCopy < mosaicMinMergedCopy) {
			cerr << "Invalid mosaic threshold or merged-copy range.\n";
			exit(734);
		}
		if (calcCovMedian) {
			std::cout << "Using median coverage for abundance estimates\n";
		}
		if (calcCoverage) {
			std::cout << "Using mean coverage for abundance estimates\n";
		}

        // sanity checks
        // we need input
        if (input == "") {//just set some defaults
            cerr << "Input must be specified\n";
            hasErr = true;
        }
        if (output == "") {//just set some defaults
            cerr << "Output must be specified\n";
            hasErr = true;
        }
        // default to min*0.95
        if(depth.size() == 0){
            depth.push_back(0.95);
        }

        if (hasErr) {
            cerr << "Use \"rtk -h\" to get full help.\n";
            exit(98);
        }
        if (mode == "rarefaction") {
            if (writeSwap) { mode = "swap"; 
            } else { mode = "memory"; }
        }
        if( write > repeats){
            cerr << "Can not create more tables than repeats of rarefaction\n";
            cerr << "Writes is now equal to repeats ("<< repeats <<")\n";
            write = repeats;
        }
    }

void options::print_rare_details(){

    stateVersion();
    // print run mode:
    std::cout << "------------------------------------ "  << std::endl;
    std::cout << "Run information:" << std::endl;
    std::cout << "input file:     " << input  << std::endl;
    std::cout << "output file:    " << output  << std::endl;
    std::cout << "mode:           " << mode  << std::endl;
    //	if(depth != 0){
    //std::cout << "depth:          " << depth  << std::endl;
    //}else{
    //	std::cout << "depth:          0.95 x min. column sum"   << std::endl;
    //}
    std::cout << "repeats:        " << repeats  << std::endl;
    std::cout << "writes:         " << write  << std::endl;
    std::cout << "threads:        " << threads  << std::endl;
    std::cout << "use swap:       ";
    if(writeSwap == false){
        std::cout << "false" << std::endl;
    }else{
        std::cout << "true" << std::endl;
    }

    //std::cout << "mode:           " << mode  << std::endl;
    std::cout << std::endl;

}

void rareExtremLowMem(options * opts,  int writeFiles, string arg4, int repeats, int numThr = 1, bool storeBinary = false){
	(void)writeFiles;
	(void)arg4;
	(void)repeats;
	(void)numThr;
	(void)storeBinary;
    // this mode takes the file, reads it in memory
    // prints the columns to their own files
    // then it loads those files again and
    // rarefies each column
    // the measures are then combines again.

    //split mat code
    std::cout << "Splitting input matrix to disk" << std::endl;
    vector<string> fileNames;
	const string outF = opts->output;
    Matrix* Mo 	= new Matrix(opts, "", fileNames, false, true);
    vector < string > SampleNames 	= Mo->getSampleNames();
    vector < string > rowNames 		= Mo->getRowNames();
    std::cout << "Done loading matrix" << std::endl;

    // abundance vectors to hold the number of occurences of genes per row
    // this will be used for ICE, ACE and or Chao2 estimation
	global_count_cube occuencesInRow(opts->depth.size(), vector<vector<rtk_count>>(opts->repeats, vector<rtk_count>(Mo->rowNum(), 0)));
	global_count_cube abundInRow(opts->depth.size(), vector<vector<rtk_count>>(opts->repeats, vector<rtk_count>(Mo->rowNum(), 0)));

    for(uint i = 0; i < opts->depth.size(); i++){
		if (opts->depth[i] < 1.) {
			// rarefy to smallest colSum
			const double resolvedDepth = round(opts->depth[i] * Mo->getMinColSum());
			if (!isfinite(resolvedDepth) || resolvedDepth < 0.0 || resolvedDepth >= static_cast<double>(numeric_limits<rtk_count>::max())) { cerr2("Resolved rarefaction depth is too large.\n", 1); }
			opts->depth[i] = resolvedDepth;
            if (opts->depth[i] == 0.0) {
                cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
                exit(1);
            }

        } 
    }

    size_t smpls = Mo->smplNum();
    delete Mo;

    //int NoOfMatrices = writeFiles;
	vector< vector< vector< rare_map > >> MaRare(opts->depth.size(), vector< vector< rare_map> > (opts->write));
    std::vector<vector<string>> cntsNames(opts->depth.size(), vector<string>());
    vector < vector < vector < string >> > tmpMatFiles(opts->depth.size(), vector<vector <string>>(opts->write));
    vector<DivEsts*> divvs(fileNames.size(),NULL);
	const uint workerCount = std::max(1u, opts->threads);
	ThreadPool pool(static_cast<size_t>(workerCount));
	vector < job > slots(workerCount);

    // now start a async in each slot
    uint i          = 0; 
    bool breakpoint(true);

    while (breakpoint) {
        // check for any finished jobs
        for( uint j = 0; j < slots.size(); j++ ){
            if( i >= smpls){
                breakpoint = false;
                // break in case we have more slots than work
                break;
            }

            if(slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(20)) == std::future_status::ready){
                // move the information
                rareStruct* tmpRS;
                tmpRS               = slots[j].fut.get();
                divvs[tmpRS->i]     = tmpRS->div;
                string curS 	    = SampleNames[tmpRS->i];
                // add the matrices to the container
                if (opts->write > 0) {
                    if (opts->writeSwap) {
                        binaryStoreSample(opts, tmpMatFiles, tmpRS, rowNames, outF, cntsNames, true);
                    }
                    else {
                        memoryStoreSample(opts, tmpRS, MaRare, cntsNames, true);
                    }
                }

                delete tmpRS;
                // free slot
                slots[j].inUse = false;
            }

            // open new slots
            if( slots[j].inUse == false){
                slots[j].inUse = true;
                // launch an async task
                DivEsts * div   = new DivEsts();
                div->SampleName = SampleNames[i];
				slots[j].fut    = pool.enqueue(calcDivRarVec, i, fileNames, div, opts, &abundInRow, &occuencesInRow);
                // send user some feedback
                int thirds = 1;
                if(smpls > 5){
                    thirds = (int) ceil((double)(smpls - 3) / 3.0);
                }

                if(i < 3 || i % thirds == 0  ){
                    std::cout << "At Sample " << i+1 << " of " << smpls << " Samples" << std::endl;
                    if(i > 3 && i % thirds == 0 ){
                        std::cout << "..." << std::endl ;
                    }
                }else if( i == 3){
                    std::cout << "..." << std::endl ;
                }

                i++;

            }

        }


    }

    // wait for the threads to finish up.
    for(uint j = 0; j < slots.size(); j++){
        if(slots[j].inUse == false ){
            // only copy if there is work to be done
            continue;
        }
        slots[j].fut.wait();
        // move the information
        rareStruct* tmpRS;
        tmpRS               = slots[j].fut.get();
        divvs[tmpRS->i]     = tmpRS->div;
        string curS 	    = SampleNames[tmpRS->i];

        // add the matrices to the container
        if (opts->write > 0) {
            if (opts->writeSwap) {
                binaryStoreSample(opts, tmpMatFiles, tmpRS, rowNames, outF, cntsNames, true);
            }
            else {
                memoryStoreSample(opts, tmpRS, MaRare, cntsNames, true);
            }
        }

        delete tmpRS;
        // free slot
        slots[j].inUse = false;
    }



    // print the div estimates out into a file
    printDivMat(outF , divvs, true, opts);
    for (size_t i = 0; i < divvs.size(); i++){
        delete divvs[i];
    }
    // write rarefaction matrices to disk
    if(opts->write > 0){
        //vector< string > rowNames = Mo->getRowNames();
        if (opts->writeSwap) {
            printRarefactionMatrix(opts, tmpMatFiles, outF, cntsNames, rowNames);
        }
        else {
            printRarefactionMatrix(opts, MaRare, outF,  cntsNames, rowNames);
        }

    }

    // compute chao2, ACE, ICE and write to file
    vector<vector<mat_fl>> chao2(opts->depth.size());
    vector<vector<mat_fl>> ICE(opts->depth.size());
    vector<vector<mat_fl>> ACE(opts->depth.size());
	computeChao2(chao2, abundInRow, static_cast<rtk_count>(smpls));
	computeCE(ICE, abundInRow, true, static_cast<rtk_count>(smpls));
	computeCE(ACE, occuencesInRow, false, static_cast<rtk_count>(smpls));
    writeGlobalDiv(opts, ICE, ACE, chao2, outF + "global_diversity.tsv");

    std::cout << "Finished\n";
}

void binaryStoreSample(options* opts, vector<vector<vector< string >>>& tmpMatFiles, rareStruct* tmpRS, vector<string>& rowNames, string outF,  vector<vector<string>>& cntsNames, bool reshapeMap){
    // store vectors of rarefied matrix on hard disk for memory reduction
    if(reshapeMap){
        vector < string > rowIDs = tmpRS->IDs;
        vector < uint > nrowIDs(rowIDs.size());
        // convert ids into integer vector, for remapping the values
        for(uint i = 0; i < rowIDs.size(); i++){
            nrowIDs[i] = std::stoi(rowIDs[i]);
        }
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            // reshape each vector, as some are zero, and we need to rematch values and rows
            // we use the row Ids which we created correctly when splitting the vector from the input file
            rare_map tmpVec;
            for (auto const& x : tmpRS->cnts[i][ii]){
                tmpVec[nrowIDs[x.first]] = x.second;
            }
            string vecLocation = printSimpleMap(tmpVec,	outF + "tmp_"  + std::to_string(ii)+"_" + std::to_string(i) + tmpRS->cntsName[i] + ".binary",	tmpRS->cntsName[i], rowNames);

			if (tmpMatFiles[i][ii].size() <= static_cast<size_t>(tmpRS->i)) { tmpMatFiles[i][ii].resize(tmpRS->i + 1); }
			tmpMatFiles[i][ii][tmpRS->i] = vecLocation;
        }
        }
    }else{
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            string vecLocation = printSimpleMap(tmpRS->cnts[i][ii],	outF + "tmp_" + std::to_string(ii) + "_" + std::to_string(i) + tmpRS->cntsName[i] + ".binary",	tmpRS->cntsName[i], rowNames);
			if (tmpMatFiles[i][ii].size() <= static_cast<size_t>(tmpRS->i)) { tmpMatFiles[i][ii].resize(tmpRS->i + 1); }
			tmpMatFiles[i][ii][tmpRS->i] = vecLocation;
        } }
    }
    // save sample name
    for(uint di = 0; di < opts->depth.size(); di++){
        if(tmpRS->cntsName[di].size() != 0){
			if (cntsNames[di].size() <= static_cast<size_t>(tmpRS->i)) { cntsNames[di].resize(tmpRS->i + 1); }
			cntsNames[di][tmpRS->i] = tmpRS->cntsName[di];

        }
    }
}

void memoryStoreSample(options* opts, rareStruct* tmpRS, vector< vector< vector< rare_map >> >& MaRare,  vector<vector<string>>& cntsNames, bool reshapeMap){
    if(reshapeMap){
        vector < string > rowIDs = tmpRS->IDs;
        vector < uint > nrowIDs(rowIDs.size());
        // convert ids into integer vector, for remapping the values
        for(uint i = 0; i < rowIDs.size(); i++){
            nrowIDs[i] = std::stoi(rowIDs[i]);
        }
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
            // reshape each vector, as some are zero, and we need to rematch values and rows
            // we use the row Ids which we created correctly when splitting the vector from the input file
            rare_map tmpVec;
            for (auto const& x : tmpRS->cnts[i][ii]){
                tmpVec[nrowIDs[x.first]] = x.second;
            }
			if (MaRare[i][ii].size() <= static_cast<size_t>(tmpRS->i)) { MaRare[i][ii].resize(tmpRS->i + 1); }
			MaRare[i][ii][tmpRS->i] = tmpVec;
        }
    }
    }else{
        for(uint i = 0; i < tmpRS->cnts.size(); i++){
            for(uint ii = 0; ii < tmpRS->cnts[i].size(); ii++){
			if (MaRare[i][ii].size() <= static_cast<size_t>(tmpRS->i)) { MaRare[i][ii].resize(tmpRS->i + 1); }
			MaRare[i][ii][tmpRS->i] = tmpRS->cnts[i][ii];
            }
        }
    }
    // save sample name
    for(uint di = 0; di < opts->depth.size(); di++){
        if(tmpRS->cntsName[di].size() != 0){
			if (cntsNames[di].size() <= static_cast<size_t>(tmpRS->i)) { cntsNames[di].resize(tmpRS->i + 1); }
			cntsNames[di][tmpRS->i] = tmpRS->cntsName[di];
        }
    }
}

void printRarefactionMatrix(options* opts, vector<vector<vector< string >>>& tmpMatFiles, string outF,  vector<vector<string>>& cntsNames, vector<string>& rowNames){
    //vector < string > rowNames = Mo->getRowNames();
    // reassemble tmp fev files:
    for(uint i = 0; i < tmpMatFiles.size(); i++){
        for(uint ii = 0; ii < tmpMatFiles[i].size(); ii++){
        string matOut = outF + "rarefied_to_" + std::to_string(opts->depth[i]) + "_n_" +  std::to_string(ii) + ".tsv";
        reassembleTmpMat(tmpMatFiles[i][ii], rowNames, cntsNames[i], matOut);
        // delete tmp rarefaction files now
        for(uint j = 0; j < tmpMatFiles[i][ii].size(); j++){
            if( remove( tmpMatFiles[i][ii][j].c_str() ) != 0 ){
                cerr << "printRarefactionMatrix: Error deleting file: " << tmpMatFiles[i][ii][j] << std::endl;
            }
        }
    }
    }
}    

void printRarefactionMatrix(options* opts, const vector<vector<vector< rare_map>>>& MaRare, string outF, vector<vector<string>>& cntsNames, vector<string>& rowNames){
    for(uint i = 0; i < MaRare.size(); i++){
        for(uint ii = 0; ii < MaRare[i].size(); ii++){
            printRareMat(outF + "rarefied_to_" + std::to_string(opts->depth[i]) + "_n_" +  std::to_string(ii) + ".tsv", MaRare[i][ii], cntsNames[i], rowNames);
    }                                                                                                  }
}





int main(int argc, char* argv[])
{

	if (argc < 2) { cerr << "Not enough arguments. Use \"rtk -h\" for getting started.\n"; exit(3); }
	startProgramTimer();
	printProgramTimerStep("rtk2 startup");
	std::atexit(reportRtkFinishTimer);

	clock_t tStart = clock();

	options* opts = new options(argc, argv);

	string inF = opts->input;
	string outF = opts->output;
	string mode = opts->mode;
	uint numThr = opts->threads;
	string arg4 = "";//std::to_string(opts->depth[0]);
	string map = opts->map;
	string refD = opts->referenceDir;
	//bool verbose = opts->verbose;

	if (mode == "rarefaction" || mode == "rare_inmat") {
		//rareDep = atoi(arg4.c_str());
		mode = "memory";
	}
	printProgramTimerStep("mode " + mode + " start");


	//all modes that classify as rarefactions:
	if (mode == "swap" || mode == "memory") {
		opts->print_rare_details();
	}


	// start the modes
	//	if (argc>3){
	if (mode == "splitMat") {
		std::cout << mode << " mode" << std::endl;
		vector<string> empt;
		Matrix* Mo = new Matrix(opts, opts->xtra, empt, false);
		//Mo->splitOnHDD(opts->output);	//Mo->writeSums(opts->output);
		delete Mo;
		printProgramTimerStep("mode splitMat done");
		std::exit(0);
		//}else if (mode == "rare_lowMem") {
		//	rareLowMem(opts->input, opts->output, writeFiles,  arg4,  repeats, numThr);
	}
	else if (mode == "mosaic") {
		std::cout << "mosaic mode" << std::endl;
		runMosaicMode(opts);
		printProgramTimerStep("mode mosaic done");
		delete opts;
		std::exit(0);
	}
	else if (mode == "decluter" || mode == "decluter2") {
		std::cout << mode << " mode" << std::endl;
		//decluttering of matrix based on 90% gene clustering + co-exclusion of genes
		opts->sparse = true;
		Matrix* Mo = new Matrix(opts, "");
		if (mode == "decluter") {
			Mo->decluter(opts);
		}
		else {
			Mo->decluter2(opts);
		}
		delete Mo;
		printProgramTimerStep("mode " + mode + " done");
	}
	else if (mode == "module") {
		std::cout << mode << " mode" << std::endl;
		if (opts->modDB.empty()) { cerr2("module mode requires a module database via -refMods.\n", 3); }
		Matrix* Mo = new Matrix(opts, ""); //needs to be KO file
		cerr << "Estimate mod AB\n";
		Mo->estimateModuleAbund(opts);
		delete Mo;
		printProgramTimerStep("mode module done");
		std::exit(0);
	}
	else if (mode == "normalize") {
		std::cout << mode << " mode" << std::endl;
		if (argc < 4) {
			std::cout << "Usage: ./rare normalize [in matrix] [outputfile]\n";
			cerr << "Not enough arguments for \"normalize\" function\n";
			exit(3);
		}

		string matCmpr = "";
		if (opts->gzOut) { matCmpr = ".gz"; }

		if (!opts->normalizeReadPass) {
			std::cout << "[normalize] strategy=sparse in-memory\n";
			opts->sparse = true;
			Matrix* Mo = new Matrix(opts, "", false);
			const vector<double> colsums = Mo->getCSum();
			vector<double> sortedSums = colsums;
			sort(sortedSums.begin(), sortedSums.end(), greater<double>());
			if (sortedSums.size() > 5) {
				std::cout << "Top 5 colSums:" << sortedSums[0] << " " << sortedSums[1] << " "
					<< sortedSums[2] << " " << sortedSums[3] << " " << sortedSums[4] << endl;
				std::cout << "Bottom 5 colSums:" << sortedSums[sortedSums.size() - 5] << " "
					<< sortedSums[sortedSums.size() - 4] << " " << sortedSums[sortedSums.size() - 3]
					<< " " << sortedSums[sortedSums.size() - 2] << " "
					<< sortedSums[sortedSums.size() - 1] << endl;
			}
			const size_t zeroCnt = static_cast<size_t>(
				std::count(sortedSums.begin(), sortedSums.end(), 0.0));
			std::cout << "Total of " << zeroCnt << " zero colsum samples\n";
			Mo->normalize(static_cast<int>(opts->threads));
			printProgramTimerStep("matrix normalized in sparse memory");
			Mo->writeMatrix(opts->output + matCmpr, false,
				static_cast<int>(opts->threads), "Norm_rtk");
			delete Mo;
		}
		else {
			std::cout << "[normalize] strategy=two-pass streaming\n";
			// First calculate column sums, then reread and normalize while writing.
		vector<double> colsums;
		vector<string> colID;
		int totalRows(0);
		{
			vector<string> fileNames;
			Matrix* Mo = new Matrix(opts, "", fileNames, false, true, false);
			colsums  = Mo->getCSum();
			colID    = Mo->getSampleNames();
			totalRows = (int)Mo->getRowNames().size();
			delete Mo;
			std::cout << "Read matrix with " << totalRows << " rows and " << colID.size() << " columns\n";
		}
		vector<double> srtCS = colsums;
		sort(srtCS.begin(), srtCS.end(), greater<double>());
		if (srtCS.size() > 5) {
			std::cout << "Top 5 colSums:" << srtCS[0] << " " << srtCS[1] << " " << srtCS[2] << " " << srtCS[3] << " " << srtCS[4] << endl;
			std::cout << "Bottom 5 colSums:" << srtCS[srtCS.size() - 5] << " " << srtCS[srtCS.size() - 4] << " " << srtCS[srtCS.size() - 3] << " " << srtCS[srtCS.size() - 2] << " " << srtCS[srtCS.size() - 1] << endl;
		}
		int zeroCnt(0);
		for (size_t i = 0; i < srtCS.size(); i++) {
			if (srtCS[i] == 0) { zeroCnt++; }
		}
		std::cout << "Total of " << zeroCnt << " zero colsum samples\n";
		std::cout << "Starting with matrix normalization\n";
		const char sep = opts->sepChar;

		Matrix* Mo = new Matrix(opts->input, opts->output + matCmpr, colsums, colID, sep, opts);
		delete Mo;
		}
		printProgramTimerStep("mode normalize done");
		delete opts;
		return 0;
	}
	else if (mode == "help" || mode == "-help" || mode == "-h" || mode == "--help") {
		helpMsg();	
	}	else if (mode == "version" || mode == "-version" || mode == "-v" || mode == "--version") {
		versionMsg();
	}	else if (mode == "submatrices") {
		std::cout << mode << " mode" << std::endl;
		extractRowsMultiMat(opts);
		std::exit(0);

	} else if (mode == "lineExtr" || mode == "rowSubset") {
		//extracts specific lines from matrix, given by ref file 
		std::cout << mode << " mode" << std::endl;
		extractRows(opts);
		printProgramTimerStep("mode " + mode + " done");
		std::exit(0);
	}
	else if (mode == "mergeMat") {
		std::cout << mode << " mode" << std::endl;
		VecFiles* VFs = new VecFiles(opts->input, opts->output, arg4);
		delete VFs;
		printProgramTimerStep("mode mergeMat done");
		std::exit(0);
	}
	else if (mode == "sumMat") {
		//creates a matrix summed to reference system (NOG, Taxa, ..)
		std::cout << mode << " mode" << std::endl;
		if (opts->referenceDir == "") {
			cerr << "sumMat requires a hierarchy/reference file path via -refD\n";
			std::exit(3);
		}

		vector<string> empt;
		Matrix* Mo = new Matrix(opts, opts->referenceDir, empt, true, false, false);
		delete Mo;
		printProgramTimerStep("mode sumMat done");
		std::exit(0);
	}
	else if (mode == "colSums" || mode == "colsums" || mode == "colSum" || mode == "rowOcc") {
		// just load and discard the matrix and report back the colsums
		std::cout << mode<< " mode" << std::endl;
		vector<string> fileNames;
		Matrix* Mo = new Matrix(opts, "", fileNames, false, true, false);
		if (mode == "colSums" || mode == "colsums" || mode == "colSum") {
			column co = Mo->getMinColumn();
			Mo->writeColSums(opts->output);
			std::cout << "" << std::endl;
			std::cout << "------------------------------------" << std::endl;
			std::cout << "ColSums output" << std::endl;
			std::cout << "Smallest column: " << co.id << std::endl;
			std::cout << "With counts:     " << co.colsum << std::endl<<endl;
			std::cout << "Colsums were written into the files:" << std::endl;
			std::cout << "    " << opts->output << "colSums.txt" << std::endl;
			std::cout << "    " << opts->output << "colSums_sorted.txt" << std::endl;
		}
		else if (mode == "rowOcc") {
			//vector<int> RO = Mo->getRowOcc();
			Mo->writeRowOcc(opts);
		}
		delete Mo;
		printProgramTimerStep("mode " + mode + " done");
		std::exit(0);
	}
	else if (mode == "geneMat") {
		if (argc < 5) { cerr << "Needs at least 4 arguments\n"; std::exit(0); }
		ClStr2Mat* cl = new ClStr2Mat(opts);// opts->input, opts->output, opts->map, opts->referenceDir, opts->calcCoverage, opts->oldMapStyle);
		delete cl;
		printProgramTimerStep("mode geneMat done");
		std::exit(0);
	}
	else if (mode == "swap") {
		std::cout << "Swap mode" << std::endl;
		rareExtremLowMem(opts, opts->write, arg4, opts->repeats, numThr, opts->writeSwap);
		printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
		printProgramTimerStep("mode swap done");
		std::exit(0);
	}
	else if (mode == "memory") {
		std::cout << "Loading input matrix to memory" << std::endl;
		Matrix* Mo = new Matrix(opts, "");//no arg for outfile &  hierachy | gene subset
		vector<DivEsts*> divvs(Mo->smplNum(), NULL);
		vector< string > rowNames = Mo->getRowNames();
		std::cout << "Done loading matrix" << std::endl;

		// transform all percentage sizes into correct values
		for(uint i = 0; i < opts->depth.size(); i++){
			if (opts->depth[i] < 1.) {
				// rarefy to smallest colSum
				const double resolvedDepth = round(opts->depth[i] * Mo->getMinColSum());
			if (!isfinite(resolvedDepth) || resolvedDepth < 0.0 || resolvedDepth >= static_cast<double>(numeric_limits<rtk_count>::max())) { cerr2("Resolved rarefaction depth is too large.\n", 1); }
				opts->depth[i] = resolvedDepth;
				if (opts->depth[i] == 0.0) {
					cerr << "Minimal sample count is 0. This can not be the rarefaction depth. Please provide a rarefaction depth > 0." << std::endl;
					exit(1);
				}
			} 
		}

		// hold rarefied matrices
		// stores : repeats - sampels eg rows - vectors of columns
		//int NoOfMatrices = opts->write;
		vector< vector< vector< rare_map > >> MaRare(opts->depth.size(), vector< vector< rare_map> > (opts->write));
		std::vector<vector<string>> cntsNames(opts->depth.size(), vector<string>());


		// abundance vectors to hold the number of occurences of genes per row
		// this will be used for Chao2 estimation
		global_count_cube occuencesInRow(opts->depth.size(), vector<vector<rtk_count>>(opts->repeats, vector<rtk_count>(Mo->rowNum(), 0)));
		global_count_cube abundInRow(opts->depth.size(), vector<vector<rtk_count>>(opts->repeats, vector<rtk_count>(Mo->rowNum(), 0)));

		//object to keep matrices
		vector < vector < vector < string >> > tmpMatFiles(opts->depth.size(), vector<vector <string>>(opts->write));
		//cerr << "TH";
		// vector keeping all the slots
		vector < job > slots(opts->threads);

		// now start a async in each slot
		uint i          = 0; 

		size_t smpls = Mo->smplNum();
		bool breakpoint(true);
		while (breakpoint) {
			// check for any finished jobs
			for( uint j = 0; j < slots.size(); j++ ){
				if( i >= smpls){
					breakpoint = false;
					// break in case we have more slots than work
					break;
				}

				if(slots[j].inUse == true && slots[j].fut.wait_for(std::chrono::milliseconds(20)) == std::future_status::ready){

					// move the information
					rareStruct* tmpRS;
					tmpRS               = slots[j].fut.get();
					divvs[tmpRS->i]     = tmpRS->div;
					string curS         = Mo->getSampleName(tmpRS->i);

					// add the matrices to the container
					if (opts->write > 0) {
						if (opts->writeSwap) {
						   binaryStoreSample(opts, tmpMatFiles, tmpRS, rowNames, opts->output, cntsNames, false);
						}
						else {
							memoryStoreSample(opts, tmpRS, MaRare, cntsNames, false);
						}
					}

					delete tmpRS;
					// free slot
					slots[j].inUse = false;
				}

				// open new slots
				if( slots[j].inUse == false){

					slots[j].inUse = true;
					// launch an async task
					DivEsts * div   = new DivEsts();
					slots[j].fut    = async(std::launch::async, calcDivRar, i, Mo, div, opts, &abundInRow, &occuencesInRow);

					// send user some feedback
					int thirds = 1;
					if(smpls > 5){
						thirds = (int) ceil((double)(smpls - 3) / 3.0);
					}
					if(i < 3 || i % thirds == 0  ){
						std::cout << "At Sample " << i+1 << " of " << smpls << " Samples" << std::endl  ;
						if(i > 3 && i % thirds == 0 ){
							std::cout << "..." << std::endl ;
						}
					}else if( i == 3){
						std::cout << "..." << std::endl ;
					}

					i++;

				}

			}


		}

		// wait for the threads to finish up.
		for(uint j = 0; j < slots.size(); j++){
			if(slots[j].inUse == false ){
				// only copy if there is work to be done
				continue;
			}
			slots[j].fut.wait();
				// move the information
				rareStruct* tmpRS;
				tmpRS               = slots[j].fut.get();
				divvs[tmpRS->i]     = tmpRS->div;

				// add the matrices to the container
				if (opts->write > 0) {
					if (opts->writeSwap) {
					   binaryStoreSample(opts, tmpMatFiles, tmpRS, rowNames, opts->output, cntsNames, false);
					}
					else {
						memoryStoreSample(opts, tmpRS, MaRare, cntsNames, false);
					}
				}

				delete tmpRS;
				// free slot
				slots[j].inUse = false;
			}

			// output matrix
		printDivMat(opts->output, divvs, true, opts);
		for (size_t i = 0; i < divvs.size(); i++) {
			delete divvs[i];
		}

		// write rarefaction matrices to disk
		if (opts->write > 0) {
			if (opts->writeSwap) {
				printRarefactionMatrix(opts, tmpMatFiles, opts->output, cntsNames, rowNames);
			}
			else {
				printRarefactionMatrix(opts, MaRare, opts->output,  cntsNames, rowNames);
			}
		}

		delete Mo;

		// compute chao2, ACE, ICE and write to file
		vector<vector<mat_fl>> chao2(opts->depth.size());
		vector<vector<mat_fl>> ICE(opts->depth.size());
		vector<vector<mat_fl>> ACE(opts->depth.size());
		computeChao2(chao2, abundInRow, static_cast<rtk_count>(smpls));
		computeCE(ICE, abundInRow, true, static_cast<rtk_count>(smpls));
		computeCE(ACE, occuencesInRow, false, static_cast<rtk_count>(smpls));
		writeGlobalDiv(opts, ICE, ACE, chao2, opts->output + "global_diversity.tsv");

		printf("CPU time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
		printProgramTimerStep("mode memory done");
		//std::cout << "Finished\n";
		std::exit(0);
	}
	else {
		std::cout << "rtk run mode \"" << mode << "\" undefined.\nSee ./rtk -h for more information.\n";
		printProgramTimerStep("mode undefined");
		//helpMsg();
	}


	printProgramTimerStep("mode " + mode + " done");
return 0;
}
