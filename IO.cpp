#include "IO.h"
std::mutex rarefyMutex;

// Roughly twenty times the historical 800 kB output buffer. Larger members
// substantially reduce both gzip framing overhead and the number of GZI entries.
static const size_t kParallelGzipChunkBytes = 16ull * 1024ull * 1024ull;

static std::chrono::steady_clock::time_point g_programStart = std::chrono::steady_clock::now();
static std::mutex g_programTimerMutex;

static bool readLittleEndian64(std::istream& in, std::uint64_t& value) {
    unsigned char bytes[8];
    in.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
    if (!in) { return false; }
    value = 0;
    for (size_t i = 0; i < sizeof(bytes); i++) {
        value |= static_cast<std::uint64_t>(bytes[i]) << (8 * i);
    }
    return true;
}

static bool writeLittleEndian64(std::ostream& out, std::uint64_t value) {
    unsigned char bytes[8];
    for (size_t i = 0; i < sizeof(bytes); i++) {
        bytes[i] = static_cast<unsigned char>((value >> (8 * i)) & 0xffu);
    }
    out.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
    return static_cast<bool>(out);
}

std::string ReplaceAll(std::string& str, const std::string& from, const std::string& to) {
	if (from.empty()) { return str; }
	size_t start_pos = 0;
	while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
	}
    return str;
}

void startProgramTimer() {
    std::lock_guard<std::mutex> lock(g_programTimerMutex);
    g_programStart = std::chrono::steady_clock::now();
}

double getProgramElapsedSeconds() {
    const auto now = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(now - g_programStart).count();
}

void printProgramTimerStep(const std::string& stepName) {
    std::lock_guard<std::mutex> lock(g_programTimerMutex);
    const auto now = std::chrono::steady_clock::now();
    const auto elapsedSeconds = static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::seconds>(now - g_programStart).count());
    const std::uint64_t hh = elapsedSeconds / 3600;
    const std::uint64_t mm = (elapsedSeconds % 3600) / 60;
    const std::uint64_t ss = elapsedSeconds % 60;
    std::ostringstream timeFmt;
    timeFmt << std::setfill('0')
        << std::setw(2) << hh << ":"
        << std::setw(2) << mm << ":"
        << std::setw(2) << ss;
    std::cout << "[timer] " << stepName << ": " << timeFmt.str() << std::endl;
}

bool isGZfile(const std::string fi) {
	if (fi.length() < 3) { return false; }
	std::string subst = fi.substr(fi.length() - 3);
	if (subst == ".gz") {
		return true;
	}
	return false;
}

std::ostream* openOutputStream(const std::string& filePath) {
    if (isGZfile(filePath)) {
#ifdef _gzipread
        std::cout << "Writing gzip output\n";
        return new ogzstream(filePath.c_str(), ios::out);
#else
        std::cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
    }
    return new ofstream(filePath, ios::out | ios::binary);
}

void closeOutputStream(std::ostream* stream) {
    if (stream == nullptr) { return; }
    stream->flush();
    if (auto fileStream = dynamic_cast<ofstream*>(stream)) {
        fileStream->close();
    }
    delete stream;
}

std::istream* openInputStream(const std::string& filePath) {
    if (isGZfile(filePath)) {
#ifdef _gzipread
        std::cout << "Reading gzip input\n";
        return new igzstream(filePath.c_str(), ios::in);
#else
        std::cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
    }
    return new ifstream(filePath.c_str());
}

void closeInputStream(std::istream* stream) {
    if (stream == nullptr) { return; }
    if (auto fileStream = dynamic_cast<ifstream*>(stream)) {
        fileStream->close();
    }
    delete stream;
}

bool forEachInputLine(
    const std::string& filePath,
    const std::function<bool(const std::string&)>& onLine,
    std::size_t gzipThreadCount) {
    if (isGZfile(filePath)) {
#ifdef _gzipread
        struct GzipMemberRange {
            std::uint64_t offset;
            std::uint64_t compressedBytes;
            std::uint64_t rawBytes;
        };

        auto processChunkText = [&](const std::string& data, std::string& carry) -> bool {
            carry += data;
            size_t start = 0;
            while (true) {
                const size_t nl = carry.find('\n', start);
                if (nl == std::string::npos) { break; }
                size_t len = nl - start;
                if (len > 0 && carry[start + len - 1] == '\r') { len--; }
                if (!onLine(carry.substr(start, len))) { return false; }
                start = nl + 1;
            }
            if (start > 0) { carry.erase(0, start); }
            return true;
        };

        auto scanGenericGzipMembers = [&](std::vector<GzipMemberRange>& members) -> bool {
            std::ifstream in(filePath.c_str(), ios::in | ios::binary);
            if (!in) { return false; }
            in.seekg(0, ios::end);
            const std::uint64_t fileSize = static_cast<std::uint64_t>(in.tellg());
            in.seekg(0, ios::beg);
            if (fileSize == 0) { return false; }

            const size_t inBufSize = 1024 * 1024;
            const size_t outBufSize = 64 * 1024;
            std::vector<unsigned char> inBuf(inBufSize);
            std::vector<unsigned char> outBuf(outBufSize);
            std::uint64_t memberStart = 0;

            while (memberStart < fileSize) {
                z_stream zs;
                memset(&zs, 0, sizeof(zs));
                if (inflateInit2(&zs, 15 + 16) != Z_OK) { return false; }
                bool gotEnd = false;
                std::uint64_t rawBytes = 0;

                while (true) {
                    if (zs.avail_in == 0) {
                        in.read(reinterpret_cast<char*>(inBuf.data()), static_cast<std::streamsize>(inBuf.size()));
                        const std::streamsize got = in.gcount();
                        if (got <= 0) {
                            inflateEnd(&zs);
                            return false;
                        }
                        zs.next_in = inBuf.data();
                        zs.avail_in = static_cast<uInt>(got);
                    }

                    zs.next_out = outBuf.data();
                    zs.avail_out = static_cast<uInt>(outBuf.size());
                    const int rc = inflate(&zs, Z_NO_FLUSH);
                    rawBytes += static_cast<std::uint64_t>(outBuf.size() - zs.avail_out);
                    if (rc == Z_STREAM_END) {
                        gotEnd = true;
                        break;
                    }
                    if (rc != Z_OK) {
                        inflateEnd(&zs);
                        return false;
                    }
                }

                const std::streamoff unread = static_cast<std::streamoff>(zs.avail_in);
                inflateEnd(&zs);
                if (!gotEnd) { return false; }
                if (unread > 0) {
                    in.seekg(-unread, ios::cur);
                    if (!in) { return false; }
                }
                const std::uint64_t memberEnd = static_cast<std::uint64_t>(in.tellg());
                if (memberEnd <= memberStart) { return false; }
                members.push_back({ memberStart, memberEnd - memberStart, rawBytes });
                memberStart = memberEnd;

                while (memberStart < fileSize) {
                    int c = in.peek();
                    if (c == EOF) { break; }
                    if (c == 0) {
                        in.get();
                        memberStart++;
                        continue;
                    }
                    break;
                }
            }
            return members.size() > 1;
        };

        auto inflateMemberBytes = [&](std::vector<char> compressedData) -> std::string {
            z_stream zs;
            memset(&zs, 0, sizeof(zs));
            const int initRc = inflateInit2(&zs, 15 + 16);
            if (initRc != Z_OK) {
                cerr2("Failed to initialize gzip decompressor for " + filePath + "\n", 11);
            }
            zs.next_in = reinterpret_cast<Bytef*>(compressedData.data());
            zs.avail_in = static_cast<uInt>(compressedData.size());

            std::string out;
            out.reserve(compressedData.size() * 3);
            std::vector<char> outBuf(1024 * 1024);
            while (true) {
                zs.next_out = reinterpret_cast<Bytef*>(outBuf.data());
                zs.avail_out = static_cast<uInt>(outBuf.size());
                const int rc = inflate(&zs, Z_NO_FLUSH);
                const size_t produced = outBuf.size() - zs.avail_out;
                if (produced > 0) { out.append(outBuf.data(), produced); }
                if (rc == Z_STREAM_END) { break; }
                if (rc != Z_OK) {
                    inflateEnd(&zs);
                    cerr2("Failed to decompress gzip member in " + filePath + "\n", 11);
                }
            }
            inflateEnd(&zs);
            return out;
        };

        auto decodeMembersInParallel = [&](const std::vector<GzipMemberRange>& members) -> bool {
            std::ifstream in(filePath.c_str(), ios::in | ios::binary);
            if (!in) {
                cerr2("Can't open gzip input file " + filePath + "\n", 11);
            }
            const size_t hardwareThreads = std::max<size_t>(
                1, std::thread::hardware_concurrency());
            const size_t workerCount = gzipThreadCount == 0
                ? hardwareThreads
                : std::max<size_t>(1, gzipThreadCount);
            ThreadPool pool(workerCount);
            // One full member per worker keeps every decompressor occupied while
            // bounding 16 MiB-member output memory to about 16 MiB per thread.
            const size_t window = workerCount;
            std::cout << "[gzip-reader] decompression threads=" << workerCount
                << ", in-flight members=" << window << std::endl;
            std::deque<std::future<std::string>> pending;
            std::string carry;
            carry.reserve(1024 * 1024);

            size_t memberIdx = 0;
            while (memberIdx < members.size()) {
                const auto& m = members[memberIdx];
                in.seekg(static_cast<std::streamoff>(m.offset), ios::beg);
                if (!in) { return false; }
                std::vector<char> compressed(static_cast<size_t>(m.compressedBytes));
                in.read(compressed.data(), static_cast<std::streamsize>(compressed.size()));
                if (!in) { return false; }
                pending.push_back(pool.enqueue(inflateMemberBytes, std::move(compressed)));
                memberIdx++;
                if (pending.size() >= window) {
                    if (!processChunkText(pending.front().get(), carry)) { return false; }
                    pending.pop_front();
                }
            }

            while (!pending.empty()) {
                if (!processChunkText(pending.front().get(), carry)) { return false; }
                pending.pop_front();
            }

            if (!carry.empty()) {
                if (!carry.empty() && carry.back() == '\r') { carry.pop_back(); }
                if (!onLine(carry)) { return false; }
            }
            return true;
        };

        const std::string sidecarPath = filePath + ".gzi";
        if (fileExists(sidecarPath)) {
            std::vector<GzipMemberRange> indexedMembers;
            std::ifstream gzipFile(filePath.c_str(), ios::in | ios::binary);
            gzipFile.seekg(0, ios::end);
            const std::streamoff gzipEnd = gzipFile.tellg();
            const std::uint64_t gzipFileSize = gzipEnd < 0 ? 0 : static_cast<std::uint64_t>(gzipEnd);
            gzipFile.close();

            // GZI stores an implicit first entry at compressed/raw offset (0,0),
            // followed by little-endian uint64 pairs for each later member.
            std::ifstream binaryIdx(sidecarPath.c_str(), ios::in | ios::binary);
            binaryIdx.seekg(0, ios::end);
            const std::streamoff indexEnd = binaryIdx.tellg();
            binaryIdx.seekg(0, ios::beg);
            std::uint64_t entryCount = 0;
            bool validBinaryIndex = indexEnd >= 8 && readLittleEndian64(binaryIdx, entryCount);
            const std::uint64_t indexSize = indexEnd < 0 ? 0 : static_cast<std::uint64_t>(indexEnd);
            if (entryCount > (std::numeric_limits<std::uint64_t>::max() - 8) / 16 ||
                indexSize != 8 + entryCount * 16) {
                validBinaryIndex = false;
            }

            std::vector<std::uint64_t> compressedOffsets(1, 0);
            std::vector<std::uint64_t> rawOffsets(1, 0);
            for (std::uint64_t i = 0; validBinaryIndex && i < entryCount; i++) {
                std::uint64_t compressedOffset = 0;
                std::uint64_t rawOffset = 0;
                if (!readLittleEndian64(binaryIdx, compressedOffset) ||
                    !readLittleEndian64(binaryIdx, rawOffset) ||
                    compressedOffset <= compressedOffsets.back() ||
                    rawOffset <= rawOffsets.back() ||
                    compressedOffset >= gzipFileSize) {
                    validBinaryIndex = false;
                    break;
                }
                compressedOffsets.push_back(compressedOffset);
                rawOffsets.push_back(rawOffset);
            }
            binaryIdx.close();

            if (validBinaryIndex && compressedOffsets.size() > 1 && gzipFileSize > 0) {
                indexedMembers.reserve(compressedOffsets.size());
                for (size_t i = 0; i < compressedOffsets.size(); i++) {
                    const std::uint64_t compressedEnd =
                        i + 1 < compressedOffsets.size() ? compressedOffsets[i + 1] : gzipFileSize;
                    const std::uint64_t rawBytes =
                        i + 1 < rawOffsets.size() ? rawOffsets[i + 1] - rawOffsets[i] : 0;
                    if (compressedEnd <= compressedOffsets[i]) {
                        indexedMembers.clear();
                        break;
                    }
                    indexedMembers.push_back({
                        compressedOffsets[i],
                        compressedEnd - compressedOffsets[i],
                        rawBytes
                    });
                }
            }

            // Backward compatibility for the text index emitted by RTK 2.22.
            if (indexedMembers.empty()) {
                std::ifstream idx(sidecarPath.c_str());
                std::string magic;
                if (idx && std::getline(idx, magic) && magic == "RTK2_PAR_GZIP_INDEX_V1") {
                    struct IndexEntry { std::uint64_t offset; std::uint64_t compressed; std::uint64_t raw; };
                    std::vector<IndexEntry> entries;
                    std::string line;
                    while (std::getline(idx, line)) {
                        if (line.empty()) { continue; }
                        std::stringstream ss(line);
                        IndexEntry e{};
                        if (!(ss >> e.offset >> e.compressed >> e.raw) || e.compressed == 0) {
                            entries.clear();
                            break;
                        }
                        entries.push_back(e);
                    }
                    if (!entries.empty()) {
                        indexedMembers.reserve(entries.size());
                        for (size_t i = 0; i < entries.size(); i++) {
                            indexedMembers.push_back({
                                entries[i].offset,
                                entries[i].compressed,
                                entries[i].raw
                            });
                        }
                    }
                }
            }

            if (!indexedMembers.empty()) {
                std::cout << "[gzip-reader] using GZI indexed parallel member decode: "
                    << indexedMembers.size() << " members\n";
                if (decodeMembersInParallel(indexedMembers)) {
                    return true;
                }
            }
        }

        auto mightHaveMultipleGzipMembers = [&]() -> bool {
            std::ifstream in(filePath.c_str(), ios::in | ios::binary);
            if (!in) { return false; }
            const size_t bufSize = 1024 * 1024;
            const std::uint64_t probeLimitBytes = 64ull * 1024ull * 1024ull;
            std::vector<unsigned char> buf(bufSize);
            std::uint64_t offset = 0;
            while (in && offset < probeLimitBytes) {
                in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
                const std::streamsize got = in.gcount();
                if (got <= 0) { break; }
                for (std::streamsize i = 0; i + 9 < got; i++) {
                    if (buf[static_cast<size_t>(i)] != 0x1f ||
                        buf[static_cast<size_t>(i + 1)] != 0x8b ||
                        buf[static_cast<size_t>(i + 2)] != 8) {
                        continue;
                    }
                    const unsigned char flg = buf[static_cast<size_t>(i + 3)];
                    if ((flg & 0xE0) != 0) { continue; }
                    const std::uint64_t pos = offset + static_cast<std::uint64_t>(i);
                    if (pos > 0) { return true; }
                }
                offset += static_cast<std::uint64_t>(got);
            }
            return false;
        };

        if (mightHaveMultipleGzipMembers()) {
            std::vector<GzipMemberRange> genericMembers;
            if (scanGenericGzipMembers(genericMembers) && decodeMembersInParallel(genericMembers)) {
                return true;
            }
        }

        gzFile file = gzopen(filePath.c_str(), "rb");
        if (file == nullptr) {
            cerr2("Can't open gzip input file " + filePath + "\n", 11);
        }
        gzbuffer(file, 1024 * 1024);
        const size_t chunkSize = 1024 * 1024;
        std::vector<char> chunk(chunkSize);
        std::string carry;
        carry.reserve(chunkSize * 2);
        while (true) {
            const int bytesRead = gzread(file, chunk.data(), static_cast<unsigned int>(chunk.size()));
            if (bytesRead < 0) {
                int gzErr = Z_OK;
                const char* gzMsg = gzerror(file, &gzErr);
                const std::string msg = gzMsg == nullptr ? "unknown zlib error" : std::string(gzMsg);
                gzclose(file);
                cerr2("Error while reading gzip input " + filePath + ": " + msg + "\n", 11);
            }
            if (bytesRead == 0) { break; }
            carry.append(chunk.data(), static_cast<size_t>(bytesRead));
            size_t start = 0;
            while (true) {
                const size_t nl = carry.find('\n', start);
                if (nl == std::string::npos) { break; }
                size_t len = nl - start;
                if (len > 0 && carry[start + len - 1] == '\r') { len--; }
                std::string line = carry.substr(start, len);
                if (!onLine(line)) {
                    gzclose(file);
                    return false;
                }
                start = nl + 1;
            }
            if (start > 0) {
                carry.erase(0, start);
            }
        }
        if (!carry.empty()) {
            if (!carry.empty() && carry.back() == '\r') { carry.pop_back(); }
            if (!onLine(carry)) {
                gzclose(file);
                return false;
            }
        }
        gzclose(file);
        return true;
#else
        std::cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
    }

    std::ifstream in(filePath.c_str());
    if (!in) {
        cerr2("Can't open input file " + filePath + "\n", 11);
    }
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') { line.pop_back(); }
        if (!onLine(line)) { return false; }
    }
    return true;
}

#ifdef _gzipread
static std::vector<char> compressGzipMemberChunk(const std::string& chunk, const std::string& targetFile) {
    std::vector<char> compressed;
    if (chunk.empty()) { return compressed; }
    if (chunk.size() > static_cast<size_t>(std::numeric_limits<uInt>::max())) {
        cerr2("Gzip chunk too large while writing " + targetFile + "\n", 69);
    }

    z_stream zs;
    memset(&zs, 0, sizeof(zs));
    const int initRc = deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY);
    if (initRc != Z_OK) {
        cerr2("Failed to initialize gzip compressor for " + targetFile + "\n", 69);
    }

    const uLong deflateCapacity = deflateBound(&zs, static_cast<uLong>(chunk.size()));
    if (deflateCapacity == 0 ||
        deflateCapacity > static_cast<uLong>(std::numeric_limits<uInt>::max())) {
        deflateEnd(&zs);
        cerr2("Invalid gzip compression bound while writing " + targetFile + "\n", 69);
    }
    compressed.resize(static_cast<size_t>(deflateCapacity));
    zs.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(chunk.data()));
    zs.avail_in = static_cast<uInt>(chunk.size());
    zs.next_out = reinterpret_cast<Bytef*>(compressed.data());
    zs.avail_out = static_cast<uInt>(compressed.size());

    const int defRc = deflate(&zs, Z_FINISH);
    if (defRc != Z_STREAM_END) {
        deflateEnd(&zs);
        cerr2("Failed to finalize gzip chunk for " + targetFile + "\n", 69);
    }
    compressed.resize(static_cast<size_t>(zs.total_out));
    deflateEnd(&zs);
    return compressed;
}
#endif

CentralOutputWriter::CentralOutputWriter(
    const std::string& filePath, std::size_t gzipThreadCount)
    : targetFile(filePath), out(nullptr), asyncMode(false), parallelGzipMode(false),
    closed(false), stopWorker(false), pendingBytes(0), maxPendingBytes(64ull * 1024ull * 1024ull), gzipCompressionPool(nullptr), gzipOutputOffset(0) {
    if (isGZfile(filePath)) {
#ifdef _gzipread
        std::cout << "Writing gzip output (parallel compression)\n";
        out = new ofstream(filePath.c_str(), ios::out | ios::binary);
        asyncMode = true;
        parallelGzipMode = true;
        const size_t hardwareThreads = std::max<size_t>(1, std::thread::hardware_concurrency());
        const size_t workerCount = gzipThreadCount == 0
            ? hardwareThreads
            : std::max<size_t>(1, gzipThreadCount);
        if (workerCount > std::numeric_limits<size_t>::max() / kParallelGzipChunkBytes) {
            cerr2("Requested gzip thread count is too large for the pending queue.\n", 69);
        }
        // Keep at least one member available per compressor. The original fixed
        // 64 MiB budget only allowed four 16 MiB tasks, regardless of -t.
        maxPendingBytes = std::max<size_t>(
            maxPendingBytes, workerCount * kParallelGzipChunkBytes);
        std::cout << "[gzip-writer] compression threads=" << workerCount
            << ", member bytes=" << kParallelGzipChunkBytes
            << ", pending MiB=" << (maxPendingBytes / (1024ull * 1024ull))
            << std::endl;
        gzipCompressionPool.reset(new ThreadPool(workerCount));
#else
        std::cout << "gzip not supported in your rtk build\n"; exit(50);
#endif
    }
    else {
        out = openOutputStream(filePath);
    }

    if (!(*out)) {
        closeOutputStream(out);
        out = nullptr;
        cerr2("Couldn't open output file " + filePath + "\n", 57);
    }
    if (asyncMode) {
        worker = std::thread(&CentralOutputWriter::workerLoop, this);
    }
}

CentralOutputWriter::~CentralOutputWriter() {
    try {
        close();
    }
    catch (...) {
    }
}

void CentralOutputWriter::writeDirect(const std::string& chunk) {
    (*out) << chunk;
    if (!(*out)) {
        cerr2("I/O error while writing output file " + targetFile + "\n", 69);
    }
}

void CentralOutputWriter::workerLoop() {
    while (true) {
        if (parallelGzipMode) {
            GzipTask task((std::future<std::vector<char>>()), 0);
            {
                std::unique_lock<std::mutex> lock(queueMutex);
                queueCv.wait(lock, [this]() { return stopWorker || !gzipTasks.empty(); });
                if (gzipTasks.empty()) {
                    if (stopWorker) { break; }
                    continue;
                }
                task = std::move(gzipTasks.front());
                gzipTasks.pop_front();
            }

            std::vector<char> compressed = task.compressed.get();
            if (!compressed.empty()) {
                out->write(compressed.data(), static_cast<std::streamsize>(compressed.size()));
                if (!(*out)) {
                    cerr2("I/O error while writing gzip output file " + targetFile + "\n", 69);
                }
                gzipMemberIndex.emplace_back(gzipOutputOffset,
                    static_cast<std::uint64_t>(compressed.size()),
                    static_cast<std::uint64_t>(task.rawBytes));
                gzipOutputOffset += static_cast<std::uint64_t>(compressed.size());
            }

            {
                std::lock_guard<std::mutex> lock(queueMutex);
                pendingBytes -= task.rawBytes;
                if (pendingBytes == 0 && gzipTasks.empty()) {
                    drainCv.notify_all();
                }
            }
        }
        else {
            std::string chunk;
            {
                std::unique_lock<std::mutex> lock(queueMutex);
                queueCv.wait(lock, [this]() { return stopWorker || !queue.empty(); });
                if (queue.empty()) {
                    if (stopWorker) { break; }
                    continue;
                }
                chunk = std::move(queue.front());
                queue.pop_front();
            }
            writeDirect(chunk);
            {
                std::lock_guard<std::mutex> lock(queueMutex);
                pendingBytes -= chunk.size();
                if (pendingBytes == 0 && queue.empty()) {
                    drainCv.notify_all();
                }
            }
        }
        queueCv.notify_all();
    }
}

void CentralOutputWriter::write(const std::string& chunk) {
    if (closed) {
        cerr2("Attempted write on closed output file " + targetFile + "\n", 69);
    }
    if (!asyncMode) {
        writeDirect(chunk);
        return;
    }

    if (parallelGzipMode) {
        std::unique_lock<std::mutex> lock(queueMutex);
        gzipInputBuffer.append(chunk);
        enqueueBufferedGzipData(lock, false);
        return;
    }

    std::unique_lock<std::mutex> lock(queueMutex);
    queueCv.wait(lock, [this, &chunk]() { return pendingBytes + chunk.size() <= maxPendingBytes; });
    pendingBytes += chunk.size();
    queue.push_back(chunk);
    lock.unlock();
    queueCv.notify_one();
}

void CentralOutputWriter::enqueueGzipChunk(std::unique_lock<std::mutex>& lock, std::string&& chunk) {
    const size_t rawBytes = chunk.size();
    queueCv.wait(lock, [this, rawBytes]() {
        return pendingBytes == 0 || rawBytes <= maxPendingBytes - pendingBytes;
    });
    pendingBytes += rawBytes;
    gzipTasks.emplace_back(
        gzipCompressionPool->enqueue(
#ifdef _gzipread
            compressGzipMemberChunk,
#else
            [](const std::string&, const std::string&) { return std::vector<char>(); },
#endif
            std::move(chunk),
            targetFile),
        rawBytes);
    queueCv.notify_one();
}

void CentralOutputWriter::enqueueBufferedGzipData(
    std::unique_lock<std::mutex>& lock, bool flushPartial) {
    while (gzipInputBuffer.size() >= kParallelGzipChunkBytes) {
        std::string chunk(gzipInputBuffer.data(), kParallelGzipChunkBytes);
        gzipInputBuffer.erase(0, kParallelGzipChunkBytes);
        enqueueGzipChunk(lock, std::move(chunk));
    }
    if (flushPartial && !gzipInputBuffer.empty()) {
        std::string chunk;
        chunk.swap(gzipInputBuffer);
        enqueueGzipChunk(lock, std::move(chunk));
    }
}

void CentralOutputWriter::flush() {
    if (closed) { return; }
    if (asyncMode) {
        std::unique_lock<std::mutex> lock(queueMutex);
        if (parallelGzipMode) {
            enqueueBufferedGzipData(lock, true);
        }
        drainCv.wait(lock, [this]() {
            if (parallelGzipMode) { return pendingBytes == 0 && gzipTasks.empty(); }
            return pendingBytes == 0 && queue.empty();
            });
    }
    out->flush();
    if (!(*out)) {
        cerr2("I/O error while flushing output file " + targetFile + "\n", 69);
    }
}

void CentralOutputWriter::close() {
    if (closed) { return; }
    flush();
    if (asyncMode) {
        {
            std::lock_guard<std::mutex> lock(queueMutex);
            stopWorker = true;
        }
        queueCv.notify_all();
        if (worker.joinable()) {
            worker.join();
        }
    }
    closeOutputStream(out);
    out = nullptr;
    // All futures have been consumed and the ordered writer has joined. Stop
    // the otherwise-idle compressor threads before touching the index sidecar,
    // keeping the finalization phases strictly separated.
    gzipCompressionPool.reset();
    if (parallelGzipMode) {
        writeGzipSidecar();
    }
    closed = true;
}

void CentralOutputWriter::writeGzipSidecar() {
    const std::string sidecarPath = targetFile + ".gzi";
    std::ofstream idx(sidecarPath.c_str(), ios::out | ios::binary | ios::trunc);
    if (!idx) {
        cerr2("Couldn't open gzip index sidecar " + sidecarPath + "\n", 57);
    }
    // Standard GZI layout: uint64 entry count, then compressed/uncompressed
    // offset pairs. The first member at (0,0) is implicit.
    const std::uint64_t entryCount = gzipMemberIndex.empty()
        ? 0
        : static_cast<std::uint64_t>(gzipMemberIndex.size() - 1);
    if (!writeLittleEndian64(idx, entryCount)) {
        cerr2("Couldn't write gzip index sidecar " + sidecarPath + "\n", 57);
    }
    std::uint64_t rawOffset = 0;
    for (size_t i = 0; i < gzipMemberIndex.size(); i++) {
        const auto& e = gzipMemberIndex[i];
        if (i > 0 &&
            (!writeLittleEndian64(idx, e.offset) ||
             !writeLittleEndian64(idx, rawOffset))) {
            cerr2("Couldn't write gzip index sidecar " + sidecarPath + "\n", 57);
        }
        rawOffset += e.rawBytes;
    }
    idx.close();
}

ThreadPool::ThreadPool(size_t threadCount) : workers(), tasks(), stopping(false) {
    if (threadCount < 1) { threadCount = 1; }
    workers.reserve(threadCount);
    for (size_t i = 0; i < threadCount; i++) {
        workers.emplace_back(&ThreadPool::workerLoop, this);
    }
}

ThreadPool::~ThreadPool() {
    {
        std::lock_guard<std::mutex> lock(taskMutex);
        stopping = true;
    }
    taskCv.notify_all();
    for (size_t i = 0; i < workers.size(); i++) {
        if (workers[i].joinable()) {
            workers[i].join();
        }
    }
}

void ThreadPool::workerLoop() {
    while (true) {
        std::function<void()> task;
        {
            std::unique_lock<std::mutex> lock(taskMutex);
            taskCv.wait(lock, [this]() { return stopping || !tasks.empty(); });
            if (stopping && tasks.empty()) { return; }
            task = std::move(tasks.front());
            tasks.pop();
        }
        task();
    }
}

void extractRowsMultiMat(options* opts) {
	string inF = opts->input;
	string outD = opts->output + "/";
	string arg4 = opts->referenceFile;
	const char sep = opts->sepChar;
	//bool check4idxMatch = opts->check4idxMatch;

    istream* in = openInputStream(inF);
    if (!(*in)) {
#ifdef notRpackage
        cerr << "Can't open infile " << inF << endl; std::exit(99);
#endif
        closeInputStream(in);
        return;
    }
	string headerStr = "";
	if (opts->header) {
        if (!getline((*in), headerStr, '\n')) {
            closeInputStream(in);
            return;
        }
	}

	string line;
	ifstream idxS(arg4.c_str());
	//guide to link gene to Bin (submatrix).. vector<string> to account for use of row in distinct submats
	robin_hood::unordered_map < string, vector<string> > srchH(0);
	//assign Bins to index to link to outS
	robin_hood::unordered_map <string, int> binH(0);
	vector<ofbufstream*> outS(0);
	if (!idxS) {
		cerr << "Can't open index file " << arg4 << endl; std::exit(99);
	}
    char sepRef = '\t';
    size_t guideLine(1);
	//read in guide file
	while (getline(idxS, line, '\n')) {
        if (line.empty()) {
            guideLine++;
            continue;
        }
		size_t pos = line.find(sepRef);
		if (pos == std::string::npos) {
            std::cout << "requires \"" << sepRef << "\" in index file: line " << guideLine << "\n" << line << "\n";
            for (size_t i = 0; i < outS.size(); i++) { delete outS[i]; }
            idxS.close();
            closeInputStream(in);
            return;
		}
		string Bin = line.substr(0, pos);
        if (Bin.empty()) {
            guideLine++;
            continue;
        }
		Bin = ReplaceAll(Bin, ":", "_");
		if (Bin == "." || Bin.find("..") != string::npos || Bin.find('/') != string::npos || Bin.find('\\') != string::npos) {
			cerr << "Unsafe output category in index file: " << Bin << "\n";
			for (size_t i = 0; i < outS.size(); i++) { delete outS[i]; }
            closeInputStream(in);
			exit(99);
		}
		size_t pos2 = line.find(sepRef,pos+1);
		string gene;
		if (pos2 == string::npos) {
			gene = line.substr(pos + 1);
		}
		else {
			gene = line.substr(pos + 1, pos2 - pos - 1);
		}
        if (gene.empty()) {
            guideLine++;
            continue;
        }
		auto idx = srchH.find(gene);
		if (idx == srchH.end()) {
			srchH[gene] = vector<string>(1,Bin);
		}	else {
			//check that  each bin is only once associated to Bin
			//bit complicated but computationally feasible.. 
			bool present = false;
			for (size_t i = 0; i < idx->second.size(); i++) {
				if (idx->second[i] == Bin) {
					present = true;
				}
			}
			if (!present) {
				idx->second.push_back(Bin);
			}
		}
		auto idx2 = binH.find(Bin);
		if (idx2 == binH.end()) {
            int idx = static_cast<int>(outS.size());
			binH[Bin] = idx;
			string outF = outD + Bin + ".mat";
			outS.push_back(new ofbufstream(outF, ios::out));
			if (headerStr != "") {
				(*outS[idx]) << headerStr +"\n";
			}
		}	
        guideLine++;
	}
	idxS.close();
	std::cout << "Found " << outS.size() << " categories, creating " << outS.size() << " submatrices in " << outD << endl;
	int geneCntFnd(0);
	//guide file finished, now separate big matrix on the fly by guide index
    if (srchH.empty() || outS.empty()) {
        for (size_t i = 0; i < outS.size(); i++) {
            delete outS[i];
        }
        closeInputStream(in);
        return;
    }
	/*
	vector < job4 > slots(num_threads);
	int j(0), cntNA(0);
	string line("");// = safeGetline2(in);
	//thread for ostream
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
	
	*/
	
    size_t matrixLine = opts->header ? 2 : 1;
    while (getline((*in), line, '\n')) {
		
		
		
		size_t pos = line.find(sep);
		if (pos == std::string::npos) {
            std::cout << "requires \"" << sep << "\" in matrix file: line " << matrixLine << "\n" << line << "\n";
            for (size_t i = 0; i < outS.size(); i++) {
                delete outS[i];
            }
            closeInputStream(in);
            return;
		}
		string gene = line.substr(0, pos);
		auto gid = srchH.find(gene);
        if (gid == srchH.end()) {
            matrixLine++;
            continue;
        }
		//sure hit after this point, something has to be written out
		geneCntFnd++;
		for (size_t j = 0; j < gid->second.size(); j++) {
			//find bins/submatrices that line should be written to
			auto bid = binH.find(gid->second[j]);
			if (bid == binH.end()) {
				cerr << "Error: could not find bin again:" << gid->second[j] << "(" << j << "), " << gene << endl;
                for (size_t i = 0; i < outS.size(); i++) {
                    delete outS[i];
                }
                closeInputStream(in);
                return;
			}
            if (bid->second < 0 || static_cast<size_t>(bid->second) >= outS.size()) {
                cerr << "Error: bin index out of range for " << gid->second[j] << endl;
                for (size_t i = 0; i < outS.size(); i++) {
                    delete outS[i];
                }
    closeInputStream(in);
                return;
            }
			(*outS[bid->second]) << line + "\n";
		}
        matrixLine++;

	}

	std::cout << "Done extracting " << geneCntFnd << " genes\n";
	//clean up
    for (size_t i = 0; i < outS.size(); i++) {
		delete outS[i];
	}
    delete in;

}
void extractRows(options* opts){

    string inF = opts->input;
    string outF = opts->output;
    string refFile = opts->referenceFile;
    const char sep = opts->sepChar;
    bool check4idxMatch = opts->check4idxMatch;

    istream* in = openInputStream(inF);
    if (!(*in)){
#ifdef notRpackage
        cerr << "Can't open infile " << inF << endl; std::exit(99);
#endif
        closeInputStream(in);
        return;
    }
    ofstream out(outF.c_str(), ios::out);
    if (!out){
#ifdef notRpackage
        cerr << "Can't open outfile " << outF << endl; std::exit(99);
#endif
        delete in;
        return;
    }
    //read file that contains nicely ordered all indexes of lines to be extracted
    string line;
    vector<size_t> srtTar;
    ifstream idxS(refFile.c_str());
    if (!idxS){
#ifdef notRpackage
        cerr << "Can't open outfile " << refFile << endl; std::exit(99);
#endif
        closeInputStream(in);
        out.close();
        return;
    }
    //char sepRef = '\t'; int cnt = 0; 
    while (getline(idxS, line, '\n')) {
        if (line.empty()){
            continue;
        }
        if (line[0] == '>'){
            line.erase(0,1);
        }
        if (line.empty()){
            continue;
        }
        try {
            srtTar.push_back((size_t)stoull(line));
        }
        catch (const std::exception&) {
            cerr << "Invalid reference index: " << line << "\n";
            idxS.close();
            closeInputStream(in);
            out.close();
            return;
        }
    }
    idxS.close();
    //sort ascending
    sort(srtTar.begin(), srtTar.end());

    if (srtTar.empty()) {
        cerr << "No valid target indexes found in reference file\n";
        closeInputStream(in);
        out.close();
        return;
    }

    //sort through new file
    cerr << "Searching for "<< srtTar.size()<<" targets in matrix\n";
    size_t cnt(1); size_t j(0); size_t extractedCnt(0);
    while (getline((*in), line, '\n')) {
        if (j < srtTar.size() && cnt == srtTar[j]){
            if (check4idxMatch && cnt!=1) { 
                //don't check in line 1, as this is the header..
                size_t pos = line.find(sep);
                if (pos == std::string::npos) {
                    std::cout << "requires \""<< sep<<"\" separated row name: line " << cnt << "\n"<<line<<"\n";
                    exit(956);
                }
                string rowN = line.substr(0,pos);
                size_t parsedRow;
                try {
                    parsedRow = (size_t)stoull(rowN);
                }
                catch (const std::exception&) {
                    cerr << "Invalid row name index at line " << cnt << ": " << rowN << "\n";
                    closeInputStream(in);
                    out.close();
                    return;
                }
                if (parsedRow != cnt) {
                    cerr << "mismatch "<<rowN<<" != "<<cnt<<"\n Consider using non sorted matrix mode\n";
                    exit(955);
                }
            }
            out << line + "\n";
            extractedCnt++;
            size_t cur = srtTar[j];
            while (j < srtTar.size() && srtTar[j] == cur){ j++; }
            if (j == srtTar.size()){ break; }
        }
        cnt++;
    }

    //(*in).close(); 
    closeInputStream(in);
    out.close();
    if (j != srtTar.size()){

#ifdef notRpackage
        cerr << "Missed " << (srtTar.size() - j) << " entries." << endl;
#endif
    }

    cerr << "Done extracting " << extractedCnt << " lines\n";
}
//****************************  smplVec::smplVec ***********
smplVec::smplVec(const vector<mat_fl>& vec, const int nt) :IDs(0),totSum(0),
    rng(std::random_device{}()), num_threads(nt), richness(-1), Shannon(-1.f){
		rtk_count cumSum(0);
		for (size_t i = 0; i < vec.size(); i++){
			if (!isfinite(vec[i]) || vec[i] < 0 || floor(vec[i]) != vec[i] || vec[i] > static_cast<mat_fl>(numeric_limits<uint>::max())) {
				cerr2("Rarefaction requires finite, non-negative integer counts; invalid value at feature " + to_string(i) + "\n", 64);
			}
			const rtk_count value = static_cast<rtk_count>(vec[i]);
			if (numeric_limits<rtk_count>::max() - cumSum < value) { cerr2("Sample count total overflow.\n", 64); }
			cumSum += value;
        }
        if (verbose){
#ifdef notRpackage
			cerr << cumSum << " allocating ";
#endif
        }
        //arr = (int*) malloc((int)cumSum * sizeof(int));
        //arr = new unsigned short[(int)cumSum];
		if (cumSum > arr.max_size()) { cerr2("Sample is too large to expand for rarefaction.\n", 64); }
        arr.resize(static_cast<size_t>(cumSum));
        if (verbose){
#ifdef notRpackage
            cerr << "memory";
#endif
        }
        totSum = cumSum;
        size_t k(0);
        //numFeatures = 0;
        for (size_t i = 0; i< vec.size(); i++){

            IDs.push_back( std::to_string(i));
			const size_t maxG = k + static_cast<size_t>(vec[i]);
            for (; k<maxG; k++){
                arr[k] = static_cast<uint>(i);
            }
        }
		numFeatures = vec.size();
        if (verbose){
#ifdef notRpackage
            cerr << "..\n";
#endif
        }
    }
smplVec::smplVec(const string inF, const int nt) :IDs(0),totSum(0), rng(std::random_device{}()),
    num_threads(nt), richness(-1),Shannon(-1.f) {

        vector<double> vec;
        ifstream in(inF.c_str());
		if (!in) { cerr2("Cannot open rarefaction sample file: " + inF + "\n", 64); }
		string line; rtk_count cumSum(0);
        while(getline(in,line,'\n')) {
        if (line.empty()) { continue; }
        string ID, extra; double num;
        stringstream ss(line);
        if (!(ss >> ID >> num) || (ss >> extra) || !isfinite(num) || num < 0 || floor(num) != num || num > static_cast<double>(numeric_limits<uint>::max())) {
            cerr2("Invalid count line in rarefaction sample file: " + line + "\n", 64);
        }
        if (numeric_limits<rtk_count>::max() - cumSum < static_cast<rtk_count>(num)) { cerr2("Sample count total overflow.\n", 64); }
        cumSum += static_cast<rtk_count>(num);
        vec.push_back(num); IDs.push_back(ID);
        }
        in.close();

        if (verbose){
#ifdef notRpackage
			cerr<<cumSum<<" allocating ";
#endif
        }
		if (cumSum > arr.max_size()) { cerr2("Sample is too large to expand for rarefaction.\n", 64); }
        arr.resize(static_cast<size_t>(cumSum));
        if (verbose){
#ifdef notRpackage
            cerr<<"memory";
#endif
        }
        totSum = cumSum;
		size_t k(0);
        for (size_t i = 0; i< vec.size(); i++){

			const size_t maxG = k + static_cast<size_t>(vec[i]);
            for (; k<maxG; k++){
				arr[k] = static_cast<uint>(i);
            }
        }
		numFeatures = vec.size();
        if (verbose){
#ifdef notRpackage
            cerr<<"..\n";
#endif
        }
    }


void smplVec::rarefy(vector<double> depts, string ofile, int rep,
        DivEsts* divs, std::vector<vector<rare_map>> & RareSample,
        vector<string>& retCntsSampleName, string& skippedSample,
		global_count_cube* abundInRow, global_count_cube* occuencesInRow,
        int writes,bool write, bool fillret){
	(void)write;
	rtk_count dep;
    // resize divvs
    divs->richness.resize(depts.size());
    divs->shannon.resize(depts.size());
    divs->simpson.resize(depts.size());
    divs->invsimpson.resize(depts.size());
    divs->chao1.resize(depts.size());
    divs->eve.resize(depts.size());
    
	for(uint i = 0; i < depts.size(); i++){
		const double requestedDepth = depts[i];
		// max() rounds up to 2^64 as a double, so equality must also be rejected.
		if (!isfinite(requestedDepth) || requestedDepth < 0.0 || floor(requestedDepth) != requestedDepth ||
			requestedDepth >= static_cast<double>(numeric_limits<rtk_count>::max())) {
			cerr2("Rarefaction depth must be a finite, non-negative integer within the 64-bit count range.\n", 64);
		}
		dep = static_cast<rtk_count>(requestedDepth);

		if (dep > totSum){
			if (skippedSample.empty()) { skippedSample = divs->SampleName; }
			if (verbose){std::cout<<"skipped sample, because rowSums < depth \n";}
			const double na = numeric_limits<double>::quiet_NaN();
			for (int curRep = 0; curRep < rep; curRep++) {
				divs->richness[i].push_back(-1);
				divs->shannon[i].push_back(na);
				divs->simpson[i].push_back(na);
				divs->invsimpson[i].push_back(na);
				divs->chao1[i].push_back(na);
				divs->eve[i].push_back(na);
				if (curRep < writes && fillret) { RareSample[i].push_back(rare_map()); }
			}
			if (writes > 0 && fillret) { retCntsSampleName[i] = divs->SampleName; }
			continue;
		}
		for (int curRep=0;curRep<rep;curRep++){
			shuffle_singl();


            //count up
            vector<uint> cnts(numFeatures, 0);
			for (rtk_count draw = 0; draw < dep; draw++){
				cnts[arr[static_cast<size_t>(draw)]]++;
			}
            string t_out = ofile;
            if (rep!=1){
                std::ostringstream oss;
                oss<<curRep;
                t_out += "_" +oss.str();
            }

            if (curRep < writes && fillret) {
                rare_map cntsMap;
                // fill map:
                for(uint i = 0; i < cnts.size(); i++){
                    if(cnts[i] != 0){
                        cntsMap.insert( std::make_pair(i, cnts[i]) );
                    }
                }
                RareSample[i].push_back(cntsMap);
                if(curRep == 0){
                    retCntsSampleName[i] = divs->SampleName; // safe the sample name as well
                }
            }
            richness = 0;
            divs->richness[i].push_back(this->getRichness(cnts));
            vector<double> three = this->calc_div(cnts, 4);
            divs->shannon[i].push_back(three[0]);
            divs->simpson[i].push_back(three[1]);
            divs->invsimpson[i].push_back(three[2]);
            divs->chao1[i].push_back(this->calc_chao1(cnts,1));
            divs->eve[i].push_back(this->calc_eveness(cnts));
            richness = 0;

            // save abundance for chao2 calculations later
			lock_guard<mutex> rarefyGuard(rarefyMutex);
			for(uint im = 0; im < IDs.size(); im++){
                //sparse convertions in swap mode
                int id = std::stoi(IDs[im]);
				if(im < cnts.size() && id >= 0 && static_cast<size_t>(id) < abundInRow->at(i)[curRep].size() && cnts[im] != 0){
					rtk_count& incidence = abundInRow->at(i)[curRep][id];
					rtk_count& abundance = occuencesInRow->at(i)[curRep][id];
					if (incidence == numeric_limits<rtk_count>::max() || abundance > numeric_limits<rtk_count>::max() - cnts[im]) {
						cerr2("Global diversity counter overflow.\n", 64);
					}
					incidence++;
					abundance += static_cast<rtk_count>(cnts[im]);
				}
			}
        }
    }
}


// vector version of diversity fuctions:
int64_t smplVec::getRichness(const vector<unsigned int>& cnts){
    richness = 0;
    for (size_t i = 0; i<cnts.size(); i++){
        //out<<IDs[i]<<"\t"<<cnts[i]<<endl;
        if (cnts[i]>0){
            richness++;
        }
    }
    return richness;
}


double smplVec::calc_chao1(const vector<uint> & vec,int corrBias){
    double Sobs((double)richness);
    double singl(0); double doubl(0);
    for (size_t i=0;i<vec.size();i++){
        if (vec[i]==1){singl++;}
        else if (vec[i]==2){doubl++;}
    }
    double est=0.0;
    if (corrBias==0){
		est = doubl > 0 ? Sobs + (singl*singl)/(2*doubl) : Sobs + singl * (singl - 1) / 2;
    } else {
		est = Sobs + (singl*(singl-1))/(2*(doubl+1));
    }
    /*if (conf.int){
      N = apply(M,2,sum)
      P = exp(-N/Sobs)
      P1 =Sobs/(1-P)
      P2 = 1.96*sqrt((Sobs*P)/(1-P))
      low =  P1 - P2
      idx = low<Sobs
      low[idx] = Sobs[idx]
      hi = P1 + P2
      }*/
    return (est);
}


vector<double> smplVec::calc_div(const vector<uint>& vec,int meth, float base){
    double sum = 0;
    for (size_t i=0; i<vec.size();i++){sum+=(double)vec[i];}
	if (sum <= 0) { return vector<double>(3, 0.0); }
    vector<double> x(vec.begin(),vec.end());
    for (size_t i=0; i<x.size();i++){x[i] /= sum;}
    bool doexp = false;
    if (base <= 2.718284f && base >= 2.718280f){ // account for machine imprecission
        doexp = true;
    }
    vector<double> H(3, 0.0); double H1(0.0), H2(0.0), H3(0.0);
    if (meth == 1 || meth == 4){
        if (doexp){
            for (size_t i=0; i<x.size();i++){if (x[i]>0){H1 += x[i] * -log(x[i])  ;}}
        } else {
            float div = -log10(base);
            for (size_t i = 0; i<x.size(); i++){ if (x[i]>0){ H1 += x[i] * log10(x[i]) / div; } }
        }
        Shannon = H1;
    }
    if (meth == 3 || meth == 4 || meth == 2) {
        for (size_t i = 0; i<x.size(); i++){ H2 += x[i] * x[i]; }
        H3 = H2;
        H2 = 1 - H2;//simpson
        H3 = 1 / H3; //invsimpson
    }
    //for (size_t i=0; i<x.size();i++){H += x[i];}
    //if (meth == (int)2) {		H = 1 - H;	}else if (meth == 3)		H = 1/H;	}
H[0] = H1; H[1] = H2; H[2] = H3;

return(H);
}
double smplVec::calc_eveness(const vector<uint>& vec){
    //double sha = calc_div(vec,1);
    if (Shannon == -1.f){ vector<double> tm = calc_div(vec, 1); }
	if (richness <= 1) { return 0.0; }
    return(Shannon / log((double)richness));
}

// map versions of diversity fuctions:
int64_t smplVec::getRichness(rare_map& cnts){
	if (cnts.size() > static_cast<size_t>(numeric_limits<int64_t>::max())) { cerr2("Richness counter overflow.\n", 64); }
	richness = static_cast<int64_t>(cnts.size()); // set richness for other functions here
    return richness;
}

double smplVec::calc_chao1(rare_map& cnts,int corrBias){
    double Sobs((double)richness);
    double singl(0); double doubl(0);
	double est(0);
    typedef rare_map::iterator it;
    for(it iterator = cnts.begin(); iterator != cnts.end(); iterator++) {
        if (iterator->second == 1){
            singl++;
        }else if(iterator->second == 2){
            doubl++;
        }
    }

    /*
       for (size_t i=0;i<vec.size();i++){
       if (vec[i]==1){singl++;}
       else if (vec[i]==2){doubl++;}
       }
       */
    if (corrBias==0){
		est = doubl > 0 ? Sobs + (singl*singl)/(2*doubl) : Sobs + singl * (singl - 1) / 2;
    } else {
		est = Sobs + (singl*(singl-1))/(2*(doubl+1));
    }

    return est;
}


vector<double> smplVec::calc_div(rare_map& cnts,int meth, float base){
    double sum = 0;
    vector<double> H(3, 0.0); double H1(0.0), H2(0.0), H3(0.0);
    //for (size_t i=0; i<vec.size();i++){sum+=(double)vec[i];}
	sum = std::accumulate(std::begin(cnts), std::end(cnts), 0.0, [] (double value, const rare_map::value_type& p){ return value + p.second; });
	if (sum <= 0) { return H; }

    // copy counts into x
    unordered_map <uint, double> x;
    x.insert(cnts.begin(), cnts.end());
    typedef unordered_map <uint, double>::iterator it;
    for(it iterator = x.begin(); iterator != x.end(); iterator++) {
        // iterator->first = key
        iterator->second /= sum;
    }
    /*
       vector<double> x(vec.begin(),vec.end());
       for (size_t i=0; i<x.size();i++){x[i] /= sum;}*/


    bool doexp = false;
    if (base <= 2.718284f && base >= 2.718280f){ // account for machine imprecission
        doexp = true;
    }

    if (meth == 1 || meth == 4){
        if (doexp){
            typedef unordered_map <uint, double>::iterator it;
            for(it iterator = x.begin(); iterator != x.end(); iterator++) {
                if (iterator->second>0){H1 += iterator->second * -log(iterator->second)  ;}
            }
            //for (size_t i=0; i<x.size();i++){if (x[i]>0){H1 += x[i] * -log(x[i])  ;}}
        } else {
            float div = -log10(base);
            typedef unordered_map <uint, double>::iterator it;
            for(it iterator = x.begin(); iterator != x.end(); iterator++) {
                if (iterator->second>0){H1 += iterator->second * -log10(iterator->second)/ div;}
            }
            //for (size_t i = 0; i<x.size(); i++){ if (x[i]>0){ H1 += x[i] * log10(x[i]) / div; } }
        }
        Shannon = H1;
    }

    if (meth == 3 || meth == 4 || meth == 2) {
        typedef unordered_map <uint, double>::iterator it;
        for(it iterator = x.begin(); iterator != x.end(); iterator++) {
            H2 += iterator->second * iterator->second;
        }
        //for (size_t i = 0; i<x.size(); i++){ H2 += x[i] * x[i]; }
        H3 = H2;
        H2 = 1 - H2;//simpson
		H3 = H3 > 0 ? 1 / H3 : 0; //invsimpson
    }
    H[0] = H1; H[1] = H2; H[2] = H3;

    return(H);
}
double smplVec::calc_eveness(rare_map& cnts){
    //double sha = calc_div(vec,1);
    if (Shannon == -1.f){ vector<double> tm = calc_div(cnts, 1); }
	if (richness <= 1) { return 0.0; }
    return(Shannon / log((double)richness));
}






void smplVec::print2File(const vector<unsigned int>& cnts,const string t_out){
    richness=0;
    ofstream out(t_out.c_str());
    for (size_t i=0;i<cnts.size();i++){
        //out<<IDs[i]<<"\t"<<cnts[i]<<endl;
        if (cnts[i]>0){
            richness++;
            out<<IDs[i]<<"\t"<<cnts[i]<<endl;
        }
    }
    out.close();
    //std::cout<<"Richness: "<<richness<<endl;
    //return richness;
}

ulong thr_rng(unsigned long pos,MyRNG& rng) {
    std::uniform_int_distribution<unsigned long> uint_distx(0,pos);
    return uint_distx(rng);
}

/*
   void smplVec::shuffle(){
   time_t seed_val=time(NULL);           // populate somehow
   rng_P.resize(num_threads);
   for (long t=0; t<num_threads; t++){ //initialize N random number generators
   rng_P[t].seed((long)seed_val-(t*13));
   }
   std::vector<std::future<ulong>> futures(num_threads);
//auto *thr = new async[num_threads];
//vector<unsigned long> rngs(num_threads,0);
ulong trng;

unsigned long pos = (unsigned long)totSum- 1;
while ( pos > 0) {
//std::cout << pos<<" ";
//Launch a group of threads
for (int t = 0; t < num_threads; ++t) {
futures[t] = async(thr_rng,pos-t,rng_P[t]);
}
for (int t = 0; t < num_threads; ++t) {
trng = futures[t].get();
unsigned int temp = arr[pos] ;
arr[pos] = arr[trng];
arr[trng] = temp;
pos--;
}
//std::uniform_int_distribution<unsigned long> uint_distx(0,pos);
//unsigned long j = uint_distx(rng);
//swap(arr[i],arr[j]);
}
//delete [] thr;
if (verbose){
#ifdef notRpackage
cerr<<"fini";
#endif
}
}
*/
/*void smplVec::shuffle_singl(){
  time_t seed_val=time(NULL);           // populate somehow
  rng.seed((long)seed_val);
  unsigned long j; unsigned int temp;
  for (unsigned long i = 0 ; i < (unsigned long)(totSum - 1); i++) {
  std::uniform_int_distribution<unsigned long> uint_distx(0,i);
  j = uint_distx(rng);
  temp = arr[i] ;
  arr[i] = arr[j];
  arr[j] = temp;
//swap(arr[i],arr[j]);
}
if (verbose){
#ifdef notRpackage
//cerr<<"fini";
#endif
}
}*/

void smplVec::shuffle_singl() {
	std::shuffle(arr.begin(), arr.end(), rng);
}


int smplVec::binarySearch( vector<float> vec, const float toFind)
{
    int len =(int) vec.size();
    // Returns index of toFind in sortedArray, or -1 if not found
    int low = 0;
    int high = len - 1;
    int mid;

    float l = vec[low];
    float h = vec[high];

    while (l < toFind && h > toFind) {
        mid = (low + high)>>1;

        float m = vec[mid];

        if (m < toFind) {
            if (vec[mid+1] > toFind){return mid;}
            l = vec[low = mid + 1];
        } else if (m > toFind) {
            if (vec[mid-1] < toFind){return mid-1;}
            h = vec[high = mid - 1];
        } else {
            return mid;
        }
        if (mid==len || mid==0){return mid;}
    }
    if (vec[low] == toFind)
        return low;
    else
        return -1; // Not found
}



/*
void DivEsts::print2file(const string file){
    if (richness.size()<1){return;}
    ofstream out(file.c_str());
    if (!out){
#ifdef notRpackage
        cerr << "Couldn't open diversity estimate file " << file << endl; std::exit(99);
#endif
    }
    out<<"Richness\t"<<richness[0];
    for (size_t i=1; i<richness.size();i++){
        out << "\t"<<richness[i];
    }
    out<<"\nShannon\t"<<shannon[0];
    for (size_t i=1; i<shannon.size();i++){
        out << "\t"<<shannon[i];
    }
    out<<"\nSimpson\t"<<simpson[0];
    for (size_t i=1; i<simpson.size();i++){
        out << "\t"<<simpson[i];
    }
    out<<"\nInv. Simpson\t"<<invsimpson[0];
    for (size_t i=1; i<invsimpson.size();i++){
        out << "\t"<<invsimpson[i];
    }
    out<<"\nChao1\t"<<chao1[0];
    for (size_t i=1; i<chao1.size();i++){
        out << "\t"<<chao1[i];
    }
    out<<"\nEveness\t"<<eve[0];
    for (size_t i=1; i<eve.size();i++){
        out << "\t"<<eve[i];
    }
    out.close();
}*/
void printDivMat(const string outF, vector<DivEsts*>& inD, bool printDIV, options* opts){

    string outFmedian = outF + "median_alpha_diversity.tsv";
    ofstream out(outFmedian.c_str());
    if (!out){
#ifdef notRpackage
        cerr << "Couldn't open diversity estimate matrix " << outF << endl; std::exit(99);
#endif
    }
    out << "depth";
    for( uint di = 0; di < opts->depth.size(); di++){
        for( uint j = 0; j < 6; j++){
                out << "\t" << opts->depth[di];   
        }
    }
    out << "\n";
    out << "Smpl";
    for( uint di = 0; di < opts->depth.size(); di++){
       out << "\tRichness\tShannon\tSimpson\tInv. Simpson\tChao1\tEveness";
    }
    out << "\n";
    for (size_t i = 0; i < inD.size(); i++){
        if (inD[i] == NULL){
#ifdef notRpackage
            cerr << "Empty vector at index " << i << " in div mat building.\n";
#endif
            out << "-1\t-1\t-1\t-1\t-1\t-1\n";
            continue;
        }
        out << inD[i]->SampleName;
        for( uint di = 0; di < opts->depth.size(); di++){
            out << "\t" << getMedian(inD[i]->richness[di]) ;
            out  << "\t"<< getMedian(inD[i]->shannon[di]);
            out  << "\t"<< getMedian(inD[i]->simpson[di]);
            out  << "\t"<< getMedian(inD[i]->invsimpson[di]) ;
            out  << "\t"<< getMedian(inD[i]->chao1[di]) << "";
            out  << "\t" << getMedian(inD[i]->eve[di]) << "";
        }
        out << "\n";
    }
    out.close();

    // print now each div estimate as well:
    if(printDIV){
        // open all files as streams
        vector<string> divNames;
        divNames.push_back("richness");
        divNames.push_back("shannon");
        divNames.push_back("simpson");
        divNames.push_back("invsimpson");
        divNames.push_back("chao1");
        divNames.push_back("eve");

        vector<ofstream> outFs(divNames.size());

        // open files
        for(uint i = 0; i < divNames.size(); i++){
            string outFdiv = outF + "alpha_" + divNames[i] + ".tsv";
            outFs[i].open(outFdiv.c_str(), ios_base::out);

            // write depth in first line of each file
            // usefull for when more than one depth was
            // requested
            outFs[i] << "depth" ;
            for(uint ii = 0; ii < opts->depth.size(); ii++){
                for(uint ij = 0; ij < opts->repeats ; ij++){
                    outFs[i] << "\t" << opts->depth[ii];
                }
            }
            outFs[i] << std::endl;
            
            outFs[i] << "repeat" ;
            for(uint ii = 0; ii < opts->depth.size(); ii++){
                for(uint ij = 0; ij < opts->repeats ; ij++){
                    outFs[i] << "\t" << ij+1;
                }
            }
            outFs[i] << std::endl;
        }

        // write the divvs to disk
        for (size_t i = 0; i < inD.size(); i++){
            // richness
            uint k = 0;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
                for( uint j = 0; j < inD[i]->richness[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->richness[di][j] ;
                }
            }
            outFs[k] << '\n';

            // shannon
            k = 1;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
                for( uint j = 0; j < inD[i]->shannon[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->shannon[di][j] ;
                }
            }
            outFs[k] << '\n';

            // simpson
            k = 2;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
                for( uint j = 0; j < inD[i]->simpson[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->simpson[di][j] ;
                }
            }
            outFs[k] << '\n';

            // invsimpson
            k = 3;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
                for( uint j = 0; j < inD[i]->invsimpson[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->invsimpson[di][j] ;
                }
            }
            outFs[k] << '\n';

            // chao1
            k = 4;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
				for( uint j = 0; j < inD[i]->eve[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->chao1[di][j] ;
                }
            }
            outFs[k] << '\n';

            // eve
            k = 5;
            outFs[k] << inD[i]->SampleName ;
            for( uint di = 0; di < opts->depth.size(); di++){
                for( uint j = 0; j < inD[i]->chao1[di].size(); j++){
                    outFs[k] << "\t" << inD[i]->eve[di][j] ;
                }
            }
            outFs[k] << '\n';
        }

        // close streams
        for(uint i = 0; i < divNames.size(); i++){
            outFs[i].close();
        }
    }

}
void printRareMat(const string outF, const vector< rare_map>& rMat, vector< string >& sampleNames, vector < string >& rowId){
    ofstream out(outF.c_str());
    if (!out){
#ifdef notRpackage
        cerr << "Couldn't open rarefy matrix file " << outF << endl; std::exit(99);
#endif
    }

    // write the header
    out << "Rarefied";
    for(uint i = 0; i < sampleNames.size(); i++){
        out << "\t"<<sampleNames[i];
    }
    out << "\n";

    // write the tsv body
    for(uint i = 0; i < rowId.size(); i++){
        out << rowId[i];
        for(uint j = 0; j < sampleNames.size(); j++){
            auto fnd = rMat[j].find(i);
            if(fnd != rMat[j].end()){
                out << "\t" << fnd->second;
            }else{
                out << "\t0";
            }
        }
        out << "\n";
    }

    out.close();
}



string printSimpleMap(const rare_map & vec, string outF, string id, vector<string> rowNames){
	(void)id;
    // takes a map from the rarefaction function and writes the vector
    // to the disk.
    // this way we dont need memory to do
    ofstream out(outF.c_str(),  ios::binary);
    if (!out){
#ifdef notRpackage
        cerr << "Couldn't open tmpvec file " << outF << endl; std::exit(99);
#endif
    }
    for(uint i = 0; i < rowNames.size(); i++){
        uint value = 0;
        auto fnd = vec.find(i);
        if(fnd != vec.end()){
            value = fnd->second;
        }
        out.write((char*) &value, sizeof(uint));

    }
    out.close();
    /*
    // write the header
    out << id;
    out << "\n";

    // write the vector
    for(uint i = 0; i < rowNames.size(); i++){
    auto fnd = vec.find(i);
    if(fnd != vec.end()){
    out << fnd->second;
    }else{
    out << "0";
    }
    out << "\n";
    }
    out.close();
    */
    return outF;
}

void reassembleTmpMat(vector<string> inF, vector< string > rowNames, vector< string > colNames, string outF){
    if(inF.empty() || inF.size() != colNames.size()){
        cerr2("Temporary rarefaction vectors do not match output columns.\n", 99);
    }

    vector<std::ifstream*> inFs(inF.size(), nullptr);
    for(uint i = 0; i < inFs.size(); i++){
        inFs[i] = new std::ifstream(inF[i].c_str(), std::ios::in | std::ios::binary);
        if (!(*inFs[i])) {
            for (uint j = 0; j <= i; j++) { delete inFs[j]; }
            cerr2("Couldn't open temporary rarefaction vector " + inF[i] + "\n", 99);
        }
    }

    ofstream out(outF.c_str());
    if (!out){
        for (size_t i = 0; i < inFs.size(); i++) { delete inFs[i]; }
        cerr2("Couldn't open rarefaction matrix output " + outF + "\n", 99);
    }
    out << "Rarefied";
    for(uint i = 0; i < colNames.size(); i++) { out << '\t' << colNames[i]; }
    out << '\n';

    for (size_t row = 0; row < rowNames.size(); row++) {
        out << rowNames[row];
        for (size_t i = 0; i < inFs.size(); i++) {
            uint value = 0;
            if (!inFs[i]->read(reinterpret_cast<char*>(&value), sizeof(value))) {
                for (size_t j = 0; j < inFs.size(); j++) { delete inFs[j]; }
                cerr2("Truncated temporary rarefaction vector " + inF[i] + "\n", 99);
            }
            out << '\t' << value;
        }
        out << '\n';
    }
    for (size_t i = 0; i < inFs.size(); i++) { delete inFs[i]; }
    out.close();

}




std::istream& safeGetline(std::istream& is, std::string& t)
{
	t.clear();
	//from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
	// The characters in the stream are read one-by-one using a std::streambuf.
	// That is faster than reading them one-by-one using the std::istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.

	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();


	for (;;) {
		int c = sb->sbumpc();
		switch (c) {
		case '\n':
			return is;
		case '\r':
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:
			// Also handle the case when the last line has no line ending
			if (t.empty())
				is.setstate(std::ios::eofbit);
			return is;
		default:
			t += (char)c;
		}
	}
}
std::string  safeGetline2(std::istream* is)
{
	string t("");
    t.reserve(256);
	//from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
	// The characters in the stream are read one-by-one using a std::streambuf.
	// That is faster than reading them one-by-one using the std::istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.

	std::istream::sentry se((*is), true);
	std::streambuf* sb = (*is).rdbuf();


	for (;;) {
		int c = sb->sbumpc();
		switch (c) {
		case '\n':
			return t;
		case '\r':
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return t;
		case EOF:
			// Also handle the case when the last line has no line ending
			if (t.empty())
				(*is).setstate(std::ios::eofbit);
			return t;
		default:
			t += (char)c;
		}
	}
	return t;
}



void computeChao2(std::vector<vector<mat_fl>>& chao2, const global_count_cube& abundInRow, rtk_count sampleCount){
    const size_t workerCount = std::max<size_t>(1, std::min<size_t>(abundInRow.size(), std::max(1u, std::thread::hardware_concurrency())));
    ThreadPool pool(workerCount);
    std::vector<std::future<void>> jobs;
    jobs.reserve(abundInRow.size());
    for (size_t i = 0; i < abundInRow.size(); i++) {
        chao2[i].resize(abundInRow[i].size(), 0);
        jobs.push_back(pool.enqueue([&chao2, &abundInRow, sampleCount, i]() {
            for (size_t ii = 0; ii < abundInRow[i].size(); ii++) {
                rtk_count NoOfSpec = 0;
                rtk_count singletons = 0;
                rtk_count doubletons = 0;
                for (size_t j = 0; j < abundInRow[i][ii].size(); j++) {
                    if (abundInRow[i][ii][j] != 0) {
                        NoOfSpec++;
                    }
                    if (abundInRow[i][ii][j] == 1) {
                        singletons++;
                    }
                    else if (abundInRow[i][ii][j] == 2) {
                        doubletons++;
                    }
                }
                const double correction = sampleCount > 0 ? static_cast<double>(sampleCount - 1) / sampleCount : 0.0;
                chao2[i][ii] = static_cast<double>(NoOfSpec) + correction * static_cast<double>(singletons) * (static_cast<double>(singletons) - 1.0) / (2.0 * (static_cast<double>(doubletons) + 1.0));
            }
        }));
    }
    for (size_t i = 0; i < jobs.size(); i++) { jobs[i].get(); }
}

void computeCE(vector<vector<mat_fl>>& CE, const global_count_cube& abundInRow, bool incidence, rtk_count sampleCount){
    // ACE uses abundance frequencies; ICE additionally uses the finite-sample incidence correction.
    const size_t workerCount = std::max<size_t>(1, std::min<size_t>(abundInRow.size(), std::max(1u, std::thread::hardware_concurrency())));
    ThreadPool pool(workerCount);
    std::vector<std::future<void>> jobs;
    jobs.reserve(abundInRow.size());
    for (size_t i = 0; i < abundInRow.size(); i++) {
        CE[i].resize(abundInRow[i].size(), 0);
        jobs.push_back(pool.enqueue([&CE, &abundInRow, incidence, sampleCount, i]() {
            for (size_t ii = 0; ii < abundInRow[i].size(); ii++) {
                vector<rtk_count> freq(11, 0);
                rtk_count nRare = 0, sAbundant = 0, sRare = 0, observed = 0;

                for (size_t j = 0; j < abundInRow[i][ii].size(); j++) {
                    const rtk_count val = abundInRow[i][ii][j];
                    if (val == 0) { continue; }
                    if (observed == numeric_limits<rtk_count>::max()) { cerr2("Diversity counter overflow while computing ACE/ICE.\n", 55); }
                    observed++;
                    if (val <= 10) {
                        if (freq[val] == numeric_limits<rtk_count>::max() || sRare == numeric_limits<rtk_count>::max() || nRare > numeric_limits<rtk_count>::max() - val) {
                            cerr2("Diversity counter overflow while computing ACE/ICE.\n", 55);
                        }
                        freq[val]++;
                        nRare += val;
                        sRare++;
                    }
                    else {
                        if (sAbundant == numeric_limits<rtk_count>::max()) { cerr2("Diversity counter overflow while computing ACE/ICE.\n", 55); }
                        sAbundant++;
                    }
                }
                const double f1 = static_cast<double>(freq[1]);
                const double nRareValue = static_cast<double>(nRare);
                const double sAbundantValue = static_cast<double>(sAbundant);
                const double sRareValue = static_cast<double>(sRare);
                const double observedValue = static_cast<double>(observed);
                if (sRare == 0) { CE[i][ii] = sAbundantValue; continue; }
                const double coverage = nRare > 0 ? 1.0 - f1 / nRareValue : 0.0;
                if (coverage <= 0.0 || nRare <= 1) {
                    const double f2 = static_cast<double>(freq[2]);
                    const double correction = incidence && sampleCount > 0 ? static_cast<double>(sampleCount - 1) / sampleCount : 1.0;
                    CE[i][ii] = observedValue + correction * f1 * (f1 - 1.0) / (2.0 * (f2 + 1.0));
                    continue;
                }
                double secondMoment = 0.0;
                for (size_t k = 1; k <= 10; k++) {
                    secondMoment += static_cast<double>(k) * (k - 1.0) * freq[k];
                }
                double gamma = (sRareValue / coverage) * secondMoment / (nRareValue * (nRareValue - 1.0)) - 1.0;
                if (incidence) {
                    if (sampleCount <= 1) { gamma = 0.0; }
                    else { gamma = (sRareValue / coverage) * (static_cast<double>(sampleCount) / (static_cast<double>(sampleCount) - 1.0)) * secondMoment / (nRareValue * nRareValue) - 1.0; }
                }
                gamma = max(0.0, gamma);
                CE[i][ii] = sAbundantValue + sRareValue / coverage + (f1 / coverage) * gamma;
            }
        }));
    }
    for (size_t i = 0; i < jobs.size(); i++) { jobs[i].get(); }
}


void writeGlobalDiv(options* opts, vector<vector<mat_fl>>& ICE, vector<vector<mat_fl>>& ACE, vector<vector<mat_fl>>& chao2, string outF){
    ofstream out(outF.c_str());
    out << "depth";

    for(uint j = 0; j < opts->depth.size(); j++){
        for(uint i = 0; i < opts->repeats; i++){
            out << "\t" << opts->depth[j];
        }
    }    
    out << '\n';
    
    
    out << "repeat";

    for(uint j = 0; j < opts->depth.size(); j++){
        for(uint i = 0; i < opts->repeats; i++){
            out << "\t" << i+1;
        }
    }    
    out << '\n';

    out << "Chao2";                      
    for(uint i = 0; i < chao2.size(); i++){
        for(uint j = 0; j < chao2[i].size(); j++){
            out << '\t' << chao2[i][j];
        }                            
    }
    out << '\n';

    out << "ICE";
    for(uint i = 0; i < ICE.size(); i++){
        for(uint j = 0; j < ICE[i].size(); j++){
            out << '\t' << ICE[i][j];
        }   }
    out << '\n';

    out << "ACE";
    for(uint i = 0; i < ACE.size(); i++){
        for(uint j = 0; j < ACE[i].size(); j++){
            out << '\t' << ACE[i][j];
        }    }
    out << '\n';
    out.close();
}


///ClStr2Mat functions
vector<string> getGenesInClstr(vector<string>& txt, vector<string>& ret, vector<float>& retF) {
	//string header;
	bool repFound = false;
	for (uint i = 0; i < txt.size(); i++) {
		string &line = txt[i];
		if (i == 0) {//new cluster, add headerInfo
			//header = line;
			continue;
		}
		//1 get gene, sample (deparse)
		size_t pos = line.find("nt, >");
		if (pos == string::npos) { cerr2("Malformed CD-HIT cluster line: " + line + "\n", 632); }
		size_t pos2 = line.find("...", pos + 5);
		if (pos2 == string::npos || pos2 <= pos + 5) { cerr2("Malformed CD-HIT cluster line: " + line + "\n", 632); }
		//string gene = line.substr(pos + 5, pos2 - pos - 5);
		ret[i-1] = line.substr(pos + 5, pos2 - pos - 5);

		if (!repFound && !line.empty() && line.back() == '*') {//report representative gene
			repFound = true;
			if (i != 1) {
				cerr2("Warning:: representative not the first gene\n");
			}
		}
		else {
			//find ID
			pos = line.find("at +/",pos2+3);
			if (pos == string::npos || line.length() <= pos + 6) { cerr2("Malformed CD-HIT identity field: " + line + "\n", 632); }
			string ID = line.substr(pos+5, line.length()-pos-6 );
			retF[i-1] = stof(ID);
			//cerr << line << " " << ID << endl;
		}
	}
	return ret;
}

static bool readFileLine(FILE* input, string& line) {
	line.clear();
	char buffer[4096];
	while (fgets(buffer, sizeof(buffer), input) != NULL) {
		line += buffer;
		if (!line.empty() && line.back() == '\n') {
			line.pop_back();
			if (!line.empty() && line.back() == '\r') { line.pop_back(); }
			return true;
		}
	}
	return !line.empty();
}

vector<string> getMMSeqsClus(FILE* incl, string &lastline,bool& cont) {
	vector<string> ret;
	string curLine;
	if (!lastline.empty()) { curLine.swap(lastline); }
	else if (!readFileLine(incl, curLine)) { cont = false; return ret; }
	size_t pos = curLine.find('\t');
	if (pos == string::npos || pos + 1 >= curLine.size()) { cerr2("Malformed MMseqs cluster line: " + curLine + "\n", 632); }
	const string curClus = curLine.substr(0, pos);
	ret.push_back(curLine.substr(pos + 1));
	while (readFileLine(incl, curLine)) {
		size_t pos = curLine.find_first_of("\t", 0);
		if (pos == std::string::npos) {
			cerr2("Malformed MMseqs cluster line: " + curLine + "\n", 632);
		}
		string tok1 =  curLine.substr(0, pos);
		if (curClus != tok1) {
			lastline = curLine;
			return ret;
		}

		ret.push_back(curLine.substr(pos+1));
		if (ret.back().empty()) { cerr2("Empty MMseqs cluster member: " + curLine + "\n", 232); }
	}
	cont = false;

	return ret;
}

textBlock* getClusBlock(istream& incl , string& lastline) {// ,  FILE* incl,
	string line;
	textBlock* ret = new textBlock;
	if (lastline.length() > 0) {
		ret->txt.push_back(lastline);
	}
	else {
		getline(incl, line);
		if (!line.empty() && line.back() == '\r') { line.pop_back(); }
		ret->txt.push_back(line);
	}
	while (getline(incl, line)) {
		if (!line.empty() && line.back() == '\r') { line.pop_back(); }
	    if (line.substr(0, 1) == ">") {//new cluster, add line to matrix
            lastline = line;
			return(ret);
		}
		ret->txt.push_back(line);
	}
	ret->cont = false;
	return ret;
}



textBlock* getClusBlock(FILE* incl, string& lastline) {// ,  FILE* incl,
    textBlock* ret = new textBlock;
    string line;

    if (lastline.length() > 0) {
        ret->txt.push_back(lastline);
    }
    else {
        if (!readFileLine(incl, line)) {
            ret->cont = false;
            return ret;
        }
		ret->txt.push_back(line);
    }
    while (readFileLine(incl, line)) {
        if (!line.empty() && line[0] == '>') {
            lastline = line;
            return(ret);
        }
		ret->txt.push_back(line);
    }
    ret->cont = false;
    return ret;
}
