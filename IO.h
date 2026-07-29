#pragma once

#include "options.h"
#include <condition_variable>
#include <deque>
#include <functional>
#include <memory>
#include <queue>
#include <thread>
#include <type_traits>


//#include <tchar.h>
//#include <string.h>

//multi threading
//#include <thread>
//#include <future>


const bool verbose=0;

#define LINELENGTH 20000;
#define OFBUFFER  800000;

std::string ReplaceAll(std::string& str, const std::string& from, const std::string& to);

std::ostream* openOutputStream(const std::string& filePath);
void closeOutputStream(std::ostream* stream);
std::istream* openInputStream(const std::string& filePath);
void closeInputStream(std::istream* stream);
bool forEachInputLine(const std::string& filePath,
	const std::function<bool(const std::string&)>& onLine,
	std::size_t gzipThreadCount = 0);
void startProgramTimer();
double getProgramElapsedSeconds();
void printProgramTimerStep(const std::string& stepName);

class ThreadPool;

class CentralOutputWriter {
public:
	explicit CentralOutputWriter(const std::string& filePath,
		std::size_t gzipThreadCount = 0);
	~CentralOutputWriter();
	void write(const std::string& chunk);
	void flush();
	void close();

private:
	struct GzipTask {
		std::future<std::vector<char>> compressed;
		std::size_t rawBytes;
		GzipTask(std::future<std::vector<char>>&& fut, std::size_t bytes)
			: compressed(std::move(fut)), rawBytes(bytes) {}
	};
	struct GzipIndexEntry {
		std::uint64_t offset;
		std::uint64_t compressedBytes;
		std::uint64_t rawBytes;
		GzipIndexEntry(std::uint64_t o, std::uint64_t c, std::uint64_t r)
			: offset(o), compressedBytes(c), rawBytes(r) {}
	};

	void workerLoop();
	void writeDirect(const std::string& chunk);
	void writeGzipSidecar();
	void enqueueGzipChunk(std::unique_lock<std::mutex>& lock, std::string&& chunk);
	void enqueueBufferedGzipData(std::unique_lock<std::mutex>& lock, bool flushPartial);

	std::string targetFile;
	std::ostream* out;
	bool asyncMode;
	bool parallelGzipMode;
	bool closed;
	bool stopWorker;
	std::size_t pendingBytes;
	std::size_t maxPendingBytes;
	std::unique_ptr<ThreadPool> gzipCompressionPool;
	std::uint64_t gzipOutputOffset;
	std::thread worker;
	std::mutex queueMutex;
	std::condition_variable queueCv;
	std::condition_variable drainCv;
	std::deque<std::string> queue;
	std::deque<GzipTask> gzipTasks;
	std::string gzipInputBuffer;
	std::vector<GzipIndexEntry> gzipMemberIndex;
};

class ThreadPool {
public:
	explicit ThreadPool(size_t threadCount);
	~ThreadPool();
	ThreadPool(const ThreadPool&) = delete;
	ThreadPool& operator=(const ThreadPool&) = delete;

	template <class F, class... Args>
	auto enqueue(F&& f, Args&&... args)
		-> std::future<typename std::result_of<F(Args...)>::type> {
		using return_type = typename std::result_of<F(Args...)>::type;
		auto task = std::make_shared<std::packaged_task<return_type()>>(
			std::bind(std::forward<F>(f), std::forward<Args>(args)...)
		);
		std::future<return_type> result = task->get_future();
		{
			std::lock_guard<std::mutex> lock(taskMutex);
			if (stopping) {
				throw std::runtime_error("enqueue on stopped ThreadPool");
			}
			tasks.emplace([task]() { (*task)(); });
		}
		taskCv.notify_one();
		return result;
	}

private:
	void workerLoop();
	std::vector<std::thread> workers;
	std::queue<std::function<void()>> tasks;
	std::mutex taskMutex;
	std::condition_variable taskCv;
	bool stopping;
};

class ofbufstream {
public:
	ofbufstream(void) { cerr << "Default ofbufstream, should not be called\n"; exit(023); }
	ofbufstream(const string IF, std::ios_base::openmode mif) :file(IF), modeIO(mif), used(0) {
		if (modeIO == ios::out) {
			remove(file.c_str());
		}
		keeper = new char[bufS];
		//test if file is writeable..
		ofstream of(file.c_str(), modeIO);
		if (!of) { cerr << "It appears outfile \"" << file << "\" is not writeable\n"; exit(763); }
		of.close();

	}
	~ofbufstream() {
		writeStream();
		delete[] keeper;
	}
	void operator<< (const string& X) {
		size_t lX(X.length());
		if (lX > bufS) {
			writeStream();
			ofstream of(file.c_str(), ios::app | ios::binary);
			if (!of) { cerr << "Could not append to output file \"" << file << "\"\n"; exit(763); }
			of.write(X.data(), static_cast<std::streamsize>(lX));
			if (!of) { cerr << "Could not write to output file \"" << file << "\"\n"; exit(763); }
			return;
		}
		if (lX + used > bufS) {
			writeStream();
		}
		memcpy(keeper + used, X.c_str(), lX);
		used += lX;
	}
private:
	void writeStream() {
		if (used == 0) { return; }
		ofstream of(file.c_str(), ios::app | ios::binary);
		of.write(keeper, used);
		if (!of) { cerr << "Could not write to output file \"" << file << "\"\n"; exit(763); }
		of.close();
		used = 0;
	}
	string file;
	char *keeper;
	std::ios_base::openmode modeIO;
	size_t used;
	static const size_t bufS = OFBUFFER;
};


typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
               // e.g. keep one global instance (per thread)


mat_fl median(std::vector<mat_fl> vec, bool ignoreZeros = false);


ulong thr_rng(unsigned long,MyRNG&);
std::istream& safeGetline(std::istream& is, std::string& t);
std::string safeGetline2(std::istream* is);

template<typename CharT, typename Traits, typename Alloc>
auto getline_n(std::basic_istream<CharT, Traits>& in, std::basic_string<CharT, Traits, Alloc>& str, std::streamsize n) -> decltype(in) {
	std::ios_base::iostate state = std::ios_base::goodbit;
	bool extracted = false;
	const typename std::basic_istream<CharT, Traits>::sentry s(in, true);
	if (s) {
		try {
			str.erase();
			typename Traits::int_type ch = in.rdbuf()->sgetc();
			for (; ; ch = in.rdbuf()->snextc()) {
				if (Traits::eq_int_type(ch, Traits::eof())) {
					// eof spotted, quit
					state |= std::ios_base::eofbit;
					break;
				}
				else if (str.size() == n) {
					// maximum number of characters met, quit
					extracted = true;
					in.rdbuf()->sbumpc();
					break;
				}
				else if (str.max_size() <= str.size()) {
					// string too big
					state |= std::ios_base::failbit;
					break;
				}
				else {
					// character valid
					str += Traits::to_char_type(ch);
					extracted = true;
				}
			}
		}
		catch (...) {
			in.setstate(std::ios_base::badbit);
		}
	}

	if (!extracted) {
		state |= std::ios_base::failbit;
	}

	in.setstate(state);
	return in;
}



template<typename T> T getMedian(vector<T> in){
	sort(in.begin(), in.end());
	size_t size = in.size();
	if (size == 0){ return (T)0; }
	if (size == 1){ return (in[0]) ; }
	if (size == 2){ return (in[0] + in[1]) / 2; }
	T median(in[size / 2]);
	if (size % 2 == 0)	{
		median = (in[size / 2 - 1] + in[size / 2]) / 2;
	}
	return median;
}
//extract rows from a matrix
void extractRows(options*);
void extractRowsMultiMat(options* opts);


inline std::string stringify(double x)
 {
   std::ostringstream o;
   o << x;
   return o.str();
 }
inline std::string itos(int number) {
	std::stringstream ss;
	ss << number;
	return ss.str();
}

class DivEsts{
public:
	DivEsts():richness(0),shannon(0),
		simpson(0),invsimpson(0),chao1(0),eve(0), depth(10){}
	~DivEsts(){}
	//void print2file(const string);
	//data vectors
	vector<vector<int64_t>> richness;
	vector<vector<double>> shannon,simpson,invsimpson,chao1,eve;
	string SampleName;
	int depth;
};
void printDivMat(const string outF, vector<DivEsts*>&, bool, options*);
void printRareMat(const string outF,const vector< rare_map>& rMat, vector< string >& sampleNames, vector < string >& rowId);
string printSimpleMap(const rare_map &vec, string outF, string id, vector<string> rowNames);
void reassembleTmpMat(vector<string> inF, vector< string > rowNames,vector< string > colNames, string outF);

class smplVec{
public:
	smplVec(const string, const int);
	smplVec(const vector<mat_fl>&, const int);
	~smplVec(){
		//delete[] arr;
	}
	void rarefy(vector<double> ,string o,int rep,DivEsts*, vector<vector<rare_map>>& RareSample,
		vector<string>& retCntsSampleName, string& skippedSample, global_count_cube*, global_count_cube*, int=0,bool=false, bool=false);
	int64_t getRichness(rare_map& cnts);
	int64_t getRichness(const vector<unsigned int>&);
	//int maxSiz(){return vector<unsigned short>::max_size();}
	vector < string > getRowNames(){ return(IDs); }

private:
	int binarySearch(vector<float>,const float x);
	//void shuffle();
	void shuffle_singl();

	//diversity indices
	//method: 1=shannon, 2=simpson, 3=invsimpson
	vector<double> calc_div(const vector<uint>& vec,int meth=1, float base=2.718282f);
	vector <double> calc_div(rare_map& , int meth=1, float base=2.718282f);
	double calc_chao1(const vector<uint> & vec,int corrBias=1);
	double calc_chao1(rare_map& , int corrBias=1); //corrBias: 0/1
	double calc_eveness(const vector<uint>& vec);
	double calc_eveness(rare_map& );

	void print2File(const vector<unsigned int>&,const string);
	//unsigned short * arr;
	vector<string> IDs;
	vector<unsigned int> arr;
	rtk_count totSum;
	vector<MyRNG> rng_P;
	MyRNG rng;
	int num_threads;
	int64_t richness;
	double Shannon;
	size_t numFeatures;

	//vector<float> vec;
};

void computeChao2(std::vector<vector<mat_fl>>& chao2, const global_count_cube& abundInRow, rtk_count sampleCount);
// compute ace or ice, depending on input data
void computeCE(vector<vector<mat_fl>>& CE, const global_count_cube& abundInRow, bool incidence, rtk_count sampleCount);


void writeGlobalDiv(options* opts, vector<vector<mat_fl>>& ICE, vector<vector<mat_fl>>& ACE, vector<vector<mat_fl>>& chao2, string outF);
//void writeGlobalDiv(vector<mat_fl>& ICE, vector<mat_fl>& ACE, vector<mat_fl>& chao2, string outF);


//utility functions for ClStr2Mat

class GeneAbundance
{
public:
	GeneAbundance(const string, const string, bool);
	mat_fl getAbundance(const string);
	bool isEmpty() { return isEmp; }
	mat_fl getTotalCnts() { mat_fl total(0); for (const auto& ix : GeneAbu) { total += ix.second; } return total; }
private:
	bool isPsAss; // is PsAssembly
	bool isSuppl;
	bool isEmp;
	SmplAbun GeneAbu;
	//vector<int> ContigID;

};

struct textBlock {
	textBlock() :txt(0), cont(true), lastLine("") { txt.reserve(300); }
	vector<string> txt;
	bool cont;
	string lastLine;
};

vector<string> getGenesInClstr(vector<string>& clbl,vector<string>&,vector<float>&);
textBlock* getClusBlock(istream&, string& lastline);
textBlock* getClusBlock(FILE*, string& lastline);
vector<string> getMMSeqsClus(FILE*, string &lastline,bool&);
