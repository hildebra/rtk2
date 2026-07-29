#include "Mosaic.h"

#include "IO.h"

#include <cstdlib>
#include <set>

namespace {

struct GeneMeta {
	std::string mgs;
	std::string group;
	double copyFraction = 0.0;
	double assignmentScore = 1.0;
	bool ambiguousMgs = false;
	bool ambiguousGroup = false;
	bool marker = false;
	bool hasCoreMetadata = false;
	bool coreEligible = true;
};

struct PafHit {
	std::string query;
	std::string target;
	std::size_t queryLength = 0;
	std::size_t targetLength = 0;
	double identity = 0.0;
	double queryCoverage = 0.0;
	double targetCoverage = 0.0;
	double relationScore = 0.0;
	double lengthRatio = 0.0;
	bool hasCigar = false;
	char alignmentType = '?';
};

struct MosaicResult {
	std::string mgs;
	std::string group;
	std::string left;
	std::string right;
	std::string status = "rejected";
	std::string reason;
	double identity = 0.0;
	double queryCoverage = 0.0;
	double targetCoverage = 0.0;
	double relationScore = 0.0;
	double scoreMargin = 0.0;
	double lengthRatio = 0.0;
	bool reciprocalBest = false;
	bool baseAligned = false;
	std::size_t baselineGenes = 0;
	std::string baselineSource;
	std::size_t informativeSamples = 0;
	std::size_t leftSupport = 0;
	std::size_t rightSupport = 0;
	std::size_t jointSupport = 0;
	double jointFraction = 0.0;
	double abundanceOverlapFraction = 0.0;
	double normalizedCorrelation = 0.0;
	double mergedCopyMedian = 0.0;
	double mergedCopyMad = 0.0;
	double mergedMedianError = 0.0;
	double bestIndividualError = 0.0;
	double fitGain = 0.0;
};

struct PafStats {
	std::size_t rows = 0;
	std::size_t retained = 0;
	std::size_t malformed = 0;
	std::size_t filtered = 0;
};

typedef std::unordered_map<std::string, GeneMeta> GeneMap;
typedef std::unordered_map<std::string,
	std::unordered_map<std::string, PafHit> > HitMap;
typedef std::unordered_map<std::string, std::vector<double> > AbundanceMap;

std::vector<std::string> splitOnChar(const std::string& line, char delimiter) {
	std::vector<std::string> fields;
	std::size_t start = 0;
	for (;;) {
		const std::size_t pos = line.find(delimiter, start);
		fields.push_back(line.substr(start, pos == std::string::npos
			? std::string::npos : pos - start));
		if (pos == std::string::npos) { break; }
		start = pos + 1;
	}
	return fields;
}

std::vector<std::string> splitTabs(const std::string& line) {
	return splitOnChar(line, '\t');
}

std::vector<std::string> splitComma(const std::string& text) {
	return splitOnChar(text, ',');
}

std::string lowerCopy(std::string text) {
	std::transform(text.begin(), text.end(), text.begin(),
		[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	return text;
}

bool parseDoubleStrict(const std::string& text, double& value) {
	if (text.empty()) { return false; }
	char* end = nullptr;
	errno = 0;
	value = std::strtod(text.c_str(), &end);
	return errno == 0 && end == text.c_str() + text.size() && std::isfinite(value);
}

bool parseSizeStrict(const std::string& text, std::size_t& value) {
	if (text.empty() || text[0] == '-') { return false; }
	char* end = nullptr;
	errno = 0;
	const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
	if (errno != 0 || end != text.c_str() + text.size() ||
		parsed > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max())) {
		return false;
	}
	value = static_cast<std::size_t>(parsed);
	return true;
}

double medianOf(std::vector<double> values) {
	if (values.empty()) { return 0.0; }
	const std::size_t middle = values.size() / 2;
	std::nth_element(values.begin(), values.begin() + middle, values.end());
	const double upper = values[middle];
	if (values.size() % 2 != 0) { return upper; }
	std::nth_element(values.begin(), values.begin() + middle - 1, values.end());
	return (upper + values[middle - 1]) / 2.0;
}

double pearson(const std::vector<double>& left, const std::vector<double>& right) {
	if (left.size() != right.size() || left.size() < 2) { return 0.0; }
	const double leftMean = std::accumulate(left.begin(), left.end(), 0.0) / left.size();
	const double rightMean = std::accumulate(right.begin(), right.end(), 0.0) / right.size();
	double covariance = 0.0;
	double leftVariance = 0.0;
	double rightVariance = 0.0;
	for (std::size_t i = 0; i < left.size(); ++i) {
		const double l = left[i] - leftMean;
		const double r = right[i] - rightMean;
		covariance += l * r;
		leftVariance += l * l;
		rightVariance += r * r;
	}
	if (leftVariance <= 0.0 || rightVariance <= 0.0) { return 0.0; }
	return covariance / std::sqrt(leftVariance * rightVariance);
}

std::string pairKey(const std::string& first, const std::string& second) {
	return first < second ? first + "\t" + second : second + "\t" + first;
}

bool looksLikeMgsHeader(const std::vector<std::string>& fields) {
	if (fields.empty()) { return false; }
	const std::string first = lowerCopy(fields[0]);
	if (first == "mgs" || first == "genome" || first == "bin") { return true; }
	if (fields.size() > 1) {
		const std::string second = lowerCopy(fields[1]);
		return (first.find("mgs") != std::string::npos &&
			(second.find("gene") != std::string::npos || second == "id"));
	}
	return false;
}

void registerGene(GeneMap& genes,
	std::unordered_map<std::string, std::vector<std::string> >& genesByMgs,
	const std::string& mgs, const std::string& gene, bool coreEligible,
	bool marker, bool hasCoreMetadata, double copyFraction,
	double assignmentScore) {
	if (mgs.empty() || gene.empty()) { return; }
	GeneMeta& meta = genes[gene];
	if (!meta.mgs.empty() && meta.mgs != mgs) {
		meta.ambiguousMgs = true;
		return;
	}
	if (meta.mgs.empty()) {
		meta.mgs = mgs;
		genesByMgs[mgs].push_back(gene);
	}
	meta.coreEligible = meta.coreEligible && coreEligible;
	meta.marker = meta.marker || marker;
	meta.hasCoreMetadata = meta.hasCoreMetadata || hasCoreMetadata;
	meta.copyFraction = std::max(meta.copyFraction, copyFraction);
	meta.assignmentScore = assignmentScore;
}

void readMgsFile(const options* opts, GeneMap& genes,
	std::unordered_map<std::string, std::vector<std::string> >& genesByMgs) {
	std::size_t lines = 0;
	std::size_t skipped = 0;
	const bool ok = forEachInputLine(opts->referenceFile, [&](const std::string& raw) {
		std::string line = raw;
		if (!line.empty() && line.back() == '\r') { line.pop_back(); }
		if (line.empty() || line[0] == '#') { return true; }
		const std::vector<std::string> fields = splitTabs(line);
		if (looksLikeMgsHeader(fields)) { return true; }
		if (fields.size() < 2) { ++skipped; return true; }

		const std::string& mgs = fields[0];
		bool coreEligible = true;
		bool marker = false;
		bool hasCoreMetadata = false;
		double copyFraction = 0.0;
		double assignmentScore = 1.0;
		if (fields.size() > 2 && fields.size() < 6) {
			++skipped;
			return true;
		}
		if (fields.size() >= 6) {
			double observations = 0.0, multicopy = 0.0, markerValue = 0.0;
			if (!parseDoubleStrict(fields[2], observations) ||
				!parseDoubleStrict(fields[3], multicopy) ||
				!parseDoubleStrict(fields[4], assignmentScore) ||
				!parseDoubleStrict(fields[5], markerValue) ||
				observations <= 0.0 || multicopy < 0.0 || assignmentScore < 0.0) {
				++skipped;
				return true;
			}
			hasCoreMetadata = true;
			copyFraction = multicopy / observations;
			marker = markerValue != 0.0;
			coreEligible = assignmentScore > 0.9 && assignmentScore < 3.0 &&
				(copyFraction < 0.1 || (marker && copyFraction < 0.2));
		}
		for (const std::string& gene : splitComma(fields[1])) {
			registerGene(genes, genesByMgs, mgs, gene, coreEligible,
				marker, hasCoreMetadata, copyFraction, assignmentScore);
		}
		++lines;
		return true;
	}, opts->threads);
	if (!ok) { cerr2("Could not parse MGS/core file: " + opts->referenceFile + "\n", 3); }
	if (genes.empty()) { cerr2("MGS/core file contains no usable gene assignments.\n", 3); }
	if (opts->verbose) {
		std::cerr << "mosaic: loaded " << genes.size() << " genes from " << lines
			<< " MGS/core rows; skipped " << skipped << " malformed rows\n";
	}
}

void assignGroup(GeneMap& genes, const std::string& gene, const std::string& group) {
	GeneMap::iterator it = genes.find(gene);
	if (it == genes.end() || group.empty() || group == "-") { return; }
	if (it->second.ambiguousGroup) { return; }
	if (!it->second.group.empty() && it->second.group != group) {
		it->second.group.clear();
		it->second.ambiguousGroup = true;
		return;
	}
	it->second.group = group;
}

void readGroupFile(const options* opts, GeneMap& genes) {
	if (opts->groupFile.empty()) { return; }
	std::size_t assignments = 0;
	const bool ok = forEachInputLine(opts->groupFile, [&](const std::string& raw) {
		std::string line = raw;
		if (!line.empty() && line.back() == '\r') { line.pop_back(); }
		if (line.empty() || line[0] == '#') { return true; }
		const std::vector<std::string> fields = splitTabs(line);
		if (fields.size() < 2) { return true; }
		if (fields.size() >= 3 && genes.find(fields[1]) != genes.end()) {
			assignGroup(genes, fields[1], fields[2]);
			++assignments;
			return true;
		}
		if (fields.size() >= 3 && fields[1].find(',') == std::string::npos) {
			std::size_t count = 0;
			if (parseSizeStrict(fields[1], count)) {
				for (const std::string& gene : splitComma(fields[2])) {
					assignGroup(genes, gene, fields[0]);
					++assignments;
				}
				return true;
			}
		}
		if (genes.find(fields[0]) != genes.end()) {
			assignGroup(genes, fields[0], fields[1]);
			++assignments;
		} else if (genes.find(fields[1]) != genes.end()) {
			assignGroup(genes, fields[1], fields[0]);
			++assignments;
		}
		return true;
	}, opts->threads);
	if (!ok) { cerr2("Could not parse group file: " + opts->groupFile + "\n", 3); }
	if (opts->verbose) {
		std::cerr << "mosaic: read " << assignments << " group assignments\n";
	}
}

bool parsePafRow(const std::vector<std::string>& fields, PafHit& hit) {
	if (fields.size() < 12) { return false; }
	std::size_t qStart = 0, qEnd = 0, tStart = 0, tEnd = 0;
	std::size_t matches = 0, blockLength = 0;
	if (!parseSizeStrict(fields[1], hit.queryLength) ||
		!parseSizeStrict(fields[2], qStart) || !parseSizeStrict(fields[3], qEnd) ||
		!parseSizeStrict(fields[6], hit.targetLength) ||
		!parseSizeStrict(fields[7], tStart) || !parseSizeStrict(fields[8], tEnd) ||
		!parseSizeStrict(fields[9], matches) ||
		!parseSizeStrict(fields[10], blockLength)) {
		return false;
	}
	if (hit.queryLength == 0 || hit.targetLength == 0 || blockLength == 0 ||
		qStart > qEnd || qEnd > hit.queryLength || tStart > tEnd || tEnd > hit.targetLength ||
		matches > blockLength) {
		return false;
	}
	hit.query = fields[0];
	hit.target = fields[5];
	hit.identity = static_cast<double>(matches) / blockLength;
	hit.queryCoverage = static_cast<double>(qEnd - qStart) / hit.queryLength;
	hit.targetCoverage = static_cast<double>(tEnd - tStart) / hit.targetLength;
	hit.relationScore = static_cast<double>(matches) /
		static_cast<double>(std::max(hit.queryLength, hit.targetLength));
	hit.lengthRatio = static_cast<double>(std::min(hit.queryLength, hit.targetLength)) /
		static_cast<double>(std::max(hit.queryLength, hit.targetLength));
	for (std::size_t i = 12; i < fields.size(); ++i) {
		if (fields[i].compare(0, 5, "cg:Z:") == 0) { hit.hasCigar = true; }
		if (fields[i].compare(0, 5, "tp:A:") == 0 && fields[i].size() > 5) {
			hit.alignmentType = fields[i][5];
		}
	}
	return !hit.query.empty() && !hit.target.empty() && hit.query != hit.target;
}

bool strongerHit(const PafHit& left, const PafHit& right) {
	if (left.relationScore != right.relationScore) {
		return left.relationScore > right.relationScore;
	}
	if (left.identity != right.identity) { return left.identity > right.identity; }
	const double leftCoverage = std::min(left.queryCoverage, left.targetCoverage);
	const double rightCoverage = std::min(right.queryCoverage, right.targetCoverage);
	if (leftCoverage != rightCoverage) { return leftCoverage > rightCoverage; }
	return left.target < right.target;
}

HitMap readPaf(const options* opts, const GeneMap& genes, PafStats& stats) {
	HitMap hits;
	const bool ok = forEachInputLine(opts->pafFile, [&](const std::string& raw) {
		std::string line = raw;
		if (!line.empty() && line.back() == '\r') { line.pop_back(); }
		if (line.empty() || line[0] == '#') { return true; }
		++stats.rows;
		PafHit hit;
		if (!parsePafRow(splitTabs(line), hit)) { ++stats.malformed; return true; }
		if (genes.find(hit.query) == genes.end()) { return true; }
		if (hit.identity < opts->mosaicMinIdentity ||
			hit.queryCoverage < opts->mosaicMinQueryCoverage ||
			hit.targetCoverage < opts->mosaicMinTargetCoverage ||
			hit.lengthRatio < opts->mosaicMinLengthRatio ||
			(opts->mosaicRequireCigar && !hit.hasCigar)) {
			++stats.filtered;
			return true;
		}
		PafHit& stored = hits[hit.query][hit.target];
		if (stored.query.empty() || strongerHit(hit, stored)) { stored = hit; }
		++stats.retained;
		return true;
	}, opts->threads);
	if (!ok) { cerr2("Could not parse PAF file: " + opts->pafFile + "\n", 3); }
	return hits;
}

bool directionalPartner(const HitMap& hits, const std::string& query,
	std::string& target, double& margin) {
	HitMap::const_iterator q = hits.find(query);
	if (q == hits.end() || q->second.empty()) { return false; }
	std::vector<const PafHit*> ranked;
	for (const auto& item : q->second) { ranked.push_back(&item.second); }
	std::sort(ranked.begin(), ranked.end(), [](const PafHit* left, const PafHit* right) {
		return strongerHit(*left, *right);
	});
	target = ranked[0]->target;
	margin = ranked.size() == 1 || ranked[0]->relationScore <= 0.0 ? 1.0 :
		(ranked[0]->relationScore - ranked[1]->relationScore) /
		ranked[0]->relationScore;
	return true;
}

const PafHit* findHit(const HitMap& hits, const std::string& query,
	const std::string& target) {
	HitMap::const_iterator q = hits.find(query);
	if (q == hits.end()) { return nullptr; }
	std::unordered_map<std::string, PafHit>::const_iterator t = q->second.find(target);
	return t == q->second.end() ? nullptr : &t->second;
}

std::vector<MosaicResult> buildSequenceCandidates(const options* opts,
	const GeneMap& genes, const HitMap& hits) {
	std::map<std::string, MosaicResult> candidates;
	for (const auto& q : hits) {
		GeneMap::const_iterator queryMeta = genes.find(q.first);
		if (queryMeta == genes.end() || queryMeta->second.ambiguousMgs) { continue; }
		for (const auto& targetEntry : q.second) {
			GeneMap::const_iterator targetMeta = genes.find(targetEntry.first);
			if (targetMeta == genes.end() || targetMeta->second.ambiguousMgs ||
				queryMeta->second.mgs != targetMeta->second.mgs) {
				continue;
			}
			if (!opts->groupFile.empty() &&
				(queryMeta->second.ambiguousGroup || targetMeta->second.ambiguousGroup ||
				 queryMeta->second.group.empty() ||
				 queryMeta->second.group != targetMeta->second.group)) {
				continue;
			}
			const std::string key = pairKey(q.first, targetEntry.first);
			if (candidates.find(key) != candidates.end()) { continue; }
			MosaicResult result;
			result.mgs = queryMeta->second.mgs;
			result.group = queryMeta->second.group;
			if (q.first < targetEntry.first) {
				result.left = q.first;
				result.right = targetEntry.first;
			} else {
				result.left = targetEntry.first;
				result.right = q.first;
			}
			candidates[key] = result;
		}
	}

	for (auto& item : candidates) {
		MosaicResult& result = item.second;
		std::string leftBest, rightBest;
		double leftMargin = 0.0, rightMargin = 0.0;
		const bool hasLeftBest = directionalPartner(hits, result.left, leftBest, leftMargin);
		const bool hasRightBest = directionalPartner(hits, result.right, rightBest, rightMargin);
		const PafHit* forward = findHit(hits, result.left, result.right);
		const PafHit* reverse = findHit(hits, result.right, result.left);
		const PafHit* representative = forward != nullptr ? forward : reverse;
		if (representative != nullptr) {
			result.identity = representative->identity;
			result.queryCoverage = representative->queryCoverage;
			result.targetCoverage = representative->targetCoverage;
			result.relationScore = representative->relationScore;
			result.lengthRatio = representative->lengthRatio;
			result.baseAligned = representative->hasCigar;
		}
		if (forward != nullptr && reverse != nullptr) {
			result.identity = std::min(forward->identity, reverse->identity);
			result.queryCoverage = std::min(forward->queryCoverage, reverse->targetCoverage);
			result.targetCoverage = std::min(forward->targetCoverage, reverse->queryCoverage);
			result.relationScore = std::min(forward->relationScore, reverse->relationScore);
			result.lengthRatio = std::min(forward->lengthRatio, reverse->lengthRatio);
			result.baseAligned = forward->hasCigar && reverse->hasCigar;
		}
		result.scoreMargin = std::min(leftMargin, rightMargin);
		result.reciprocalBest = hasLeftBest && hasRightBest &&
			leftBest == result.right && rightBest == result.left;
		if (forward == nullptr || reverse == nullptr) {
			result.reason = "missing_reciprocal_alignment";
		} else if (!result.reciprocalBest) {
			result.reason = "not_reciprocal_best";
		} else if (result.scoreMargin < opts->mosaicMinScoreMargin) {
			result.reason = "ambiguous_sequence_partner";
		} else {
			result.status = "sequence_candidate";
			result.reason = "pending_abundance";
		}
	}

	std::vector<MosaicResult> output;
	for (auto& item : candidates) { output.push_back(item.second); }
	return output;
}

AbundanceMap readMatrixRows(const options* opts,
	const std::unordered_set<std::string>& wantedGenes, std::size_t& sampleCount) {
	AbundanceMap abundances;
	bool haveHeader = false;
	std::size_t lineNumber = 0;
	const bool ok = forEachInputLine(opts->input, [&](const std::string& raw) {
		++lineNumber;
		std::string line = raw;
		if (!line.empty() && line.back() == '\r') { line.pop_back(); }
		if (line.empty()) { return true; }
		if (line[0] == '#' && (haveHeader || line.find(opts->sepChar) == std::string::npos)) {
			return true;
		}
		const std::vector<std::string> fields = splitOnChar(line, opts->sepChar);
		if (!haveHeader) {
			if (fields.size() < 2) { cerr2("Mosaic matrix header needs sample columns.\n", 3); }
			sampleCount = fields.size() - 1;
			haveHeader = true;
			return true;
		}
		if (fields.size() != sampleCount + 1) {
			cerr2("Mosaic matrix row " + std::to_string(lineNumber) +
				" has the wrong number of columns.\n", 3);
		}
		if (wantedGenes.find(fields[0]) == wantedGenes.end()) { return true; }
		std::vector<double> values(sampleCount, 0.0);
		for (std::size_t i = 0; i < sampleCount; ++i) {
			if (!parseDoubleStrict(fields[i + 1], values[i]) || values[i] < 0.0) {
				cerr2("Invalid abundance at matrix row " + std::to_string(lineNumber) + ".\n", 3);
			}
		}
		if (!abundances.emplace(fields[0], std::move(values)).second) {
			cerr2("Duplicate matrix row for gene " + fields[0] + ".\n", 3);
		}
		return true;
	}, opts->threads);
	if (!ok || !haveHeader) { cerr2("Could not parse abundance matrix: " + opts->input + "\n", 3); }
	return abundances;
}

std::unordered_map<std::string, std::vector<double> > buildBaselines(
	const options* opts, const GeneMap& genes,
	const std::unordered_set<std::string>& excludedGenes,
	const std::unordered_map<std::string, std::vector<std::string> >& genesByMgs,
	const AbundanceMap& abundances, std::size_t sampleCount,
	std::unordered_map<std::string, std::size_t>& baselineGeneCounts,
	std::unordered_map<std::string, std::string>& baselineSources) {
	std::unordered_map<std::string, std::vector<double> > baselines;
	for (const auto& mgsEntry : genesByMgs) {
		std::vector<const std::vector<double>*> allRows;
		std::vector<const std::vector<double>*> markerRows;
		for (const std::string& gene : mgsEntry.second) {
			if (excludedGenes.find(gene) != excludedGenes.end()) { continue; }
			GeneMap::const_iterator meta = genes.find(gene);
			AbundanceMap::const_iterator abundance = abundances.find(gene);
			if (meta == genes.end() || meta->second.ambiguousMgs ||
				!meta->second.coreEligible || abundance == abundances.end()) {
				continue;
			}
			allRows.push_back(&abundance->second);
			if (meta->second.marker) { markerRows.push_back(&abundance->second); }
		}
		const bool useMarkers = opts->mosaicPreferMarkerBaseline &&
			markerRows.size() >= opts->mosaicMinBaselineGenes;
		const std::vector<const std::vector<double>*>& rows =
			useMarkers ? markerRows : allRows;
		baselineGeneCounts[mgsEntry.first] = rows.size();
		baselineSources[mgsEntry.first] = useMarkers ? "marker" : "core";
		if (rows.empty()) { continue; }
		std::vector<double>& baseline = baselines[mgsEntry.first];
		baseline.resize(sampleCount, 0.0);
		for (std::size_t sample = 0; sample < sampleCount; ++sample) {
			std::vector<double> column;
			column.reserve(rows.size());
			for (const std::vector<double>* row : rows) { column.push_back((*row)[sample]); }
			baseline[sample] = medianOf(std::move(column));
		}
	}
	return baselines;
}

void addAbundanceEvidence(const options* opts, const GeneMap& genes,
	const std::unordered_map<std::string, std::vector<std::string> >& genesByMgs,
	const AbundanceMap& abundances, std::size_t sampleCount,
	std::vector<MosaicResult>& results) {
	std::unordered_set<std::string> excludedGenes;
	for (const MosaicResult& result : results) {
		if (result.status != "sequence_candidate") { continue; }
		excludedGenes.insert(result.left);
		excludedGenes.insert(result.right);
	}
	std::unordered_map<std::string, std::size_t> baselineGeneCounts;
	std::unordered_map<std::string, std::string> baselineSources;
	const std::unordered_map<std::string, std::vector<double> > baselines =
		buildBaselines(opts, genes, excludedGenes, genesByMgs, abundances, sampleCount,
			baselineGeneCounts, baselineSources);

	for (MosaicResult& result : results) {
		if (result.status != "sequence_candidate") { continue; }
		result.baselineGenes = baselineGeneCounts[result.mgs];
		result.baselineSource = baselineSources[result.mgs];
		AbundanceMap::const_iterator left = abundances.find(result.left);
		AbundanceMap::const_iterator right = abundances.find(result.right);
		if (left == abundances.end() || right == abundances.end()) {
			result.status = "rejected";
			result.reason = "missing_matrix_row";
			continue;
		}
		if (baselineGeneCounts[result.mgs] < opts->mosaicMinBaselineGenes) {
			result.status = "rejected";
			result.reason = "insufficient_baseline_genes";
			continue;
		}
		auto baselineIt = baselines.find(result.mgs);
		if (baselineIt == baselines.end()) {
			result.status = "rejected";
			result.reason = "missing_mgs_baseline";
			continue;
		}

		std::vector<double> normalizedLeft;
		std::vector<double> normalizedRight;
		std::vector<double> mergedCopies;
		std::vector<double> mergedErrors;
		std::vector<double> leftErrors;
		std::vector<double> rightErrors;
		double totalLeftCopy = 0.0;
		double totalRightCopy = 0.0;
		double overlappingCopy = 0.0;
		for (std::size_t sample = 0; sample < sampleCount; ++sample) {
			const double baseline = baselineIt->second[sample];
			if (baseline <= opts->mosaicMinBaseline) { continue; }
			const double leftCopy = left->second[sample] / baseline;
			const double rightCopy = right->second[sample] / baseline;
			const double mergedCopy = leftCopy + rightCopy;
			normalizedLeft.push_back(leftCopy);
			normalizedRight.push_back(rightCopy);
			mergedCopies.push_back(mergedCopy);
			totalLeftCopy += leftCopy;
			totalRightCopy += rightCopy;
			overlappingCopy += std::min(leftCopy, rightCopy);
			mergedErrors.push_back(std::fabs(mergedCopy - 1.0));
			leftErrors.push_back(std::fabs(leftCopy - 1.0));
			rightErrors.push_back(std::fabs(rightCopy - 1.0));
			const double detection = opts->mosaicDetectionFraction;
			if (leftCopy >= detection) { ++result.leftSupport; }
			if (rightCopy >= detection) { ++result.rightSupport; }
			if (leftCopy >= detection && rightCopy >= detection) { ++result.jointSupport; }
		}
		result.informativeSamples = mergedCopies.size();
		if (result.informativeSamples == 0) {
			result.status = "rejected";
			result.reason = "no_informative_samples";
			continue;
		}
		const std::size_t smallerSupport = std::min(result.leftSupport, result.rightSupport);
		result.jointFraction = smallerSupport == 0 ? 0.0 :
			static_cast<double>(result.jointSupport) / smallerSupport;
		const double smallerTotal = std::min(totalLeftCopy, totalRightCopy);
		result.abundanceOverlapFraction = smallerTotal <= 0.0 ? 0.0 : overlappingCopy / smallerTotal;
		result.normalizedCorrelation = pearson(normalizedLeft, normalizedRight);
		result.mergedCopyMedian = medianOf(mergedCopies);
		std::vector<double> deviations;
		deviations.reserve(mergedCopies.size());
		for (double copy : mergedCopies) {
			deviations.push_back(std::fabs(copy - result.mergedCopyMedian));
		}
		result.mergedCopyMad = medianOf(std::move(deviations));
		result.mergedMedianError = medianOf(std::move(mergedErrors));
		result.bestIndividualError = std::min(medianOf(std::move(leftErrors)),
			medianOf(std::move(rightErrors)));
		result.fitGain = result.bestIndividualError - result.mergedMedianError;

		if (result.informativeSamples < opts->mosaicMinSamples) {
			result.status = "rejected";
			result.reason = "insufficient_informative_samples";
		} else if (result.leftSupport < opts->mosaicMinVariantSamples ||
			result.rightSupport < opts->mosaicMinVariantSamples) {
			result.status = "rejected";
			result.reason = "insufficient_variant_support";
		} else if (result.jointFraction > opts->mosaicMaxSampleOverlapFraction ||
			result.abundanceOverlapFraction > opts->mosaicMaxAbundanceOverlapFraction) {
			result.status = "rejected";
			result.reason = "excess_variant_overlap";
		} else if (result.mergedCopyMedian < opts->mosaicMinMergedCopy ||
			result.mergedCopyMedian > opts->mosaicMaxMergedCopy) {
			result.status = "rejected";
			result.reason = "merged_copy_out_of_range";
		} else if (result.mergedMedianError > opts->mosaicMaxMergedError) {
			result.status = "rejected";
			result.reason = "poor_merged_copy_fit";
		} else if (result.fitGain < opts->mosaicMinFitGain) {
			result.status = "rejected";
			result.reason = "insufficient_fit_gain";
		} else {
			result.status = "confirmed";
			result.reason = "sequence_and_abundance_support";
		}
	}
}

std::string formatDouble(double value) {
	std::ostringstream out;
	out << std::setprecision(8) << value;
	return out.str();
}

void writeResults(const options* opts, const PafStats& pafStats,
	const GeneMap& genes, const std::vector<MosaicResult>& results,
	std::size_t sampleCount) {
	const std::string detailPath = opts->output + ".mosaic.tsv";
	const std::string concatPath = opts->output + ".concat.list";
	const std::string summaryPath = opts->output + ".mosaic.summary.tsv";

	std::ostream* detail = openOutputStream(detailPath);
	if (detail == nullptr) { cerr2("Could not write " + detailPath + "\n", 3); }
	*detail << "mgs\tgroup\tleft\tright\tstatus\treason\tidentity\tquery_coverage"
		"\ttarget_coverage\trelation_score\tscore_margin\tlength_ratio"
		"\treciprocal_best\tbase_aligned\tbaseline_source\tbaseline_genes\tinformative_samples\tleft_support"
		"\tright_support\tjoint_support\tsample_overlap_fraction\tabundance_overlap_fraction\tnormalized_correlation"
		"\tmerged_copy_median\tmerged_copy_mad\tmerged_median_abs_error"
		"\tbest_individual_median_abs_error\tfit_gain\n";
	for (const MosaicResult& result : results) {
		*detail << result.mgs << '\t' << result.group << '\t' << result.left << '\t'
			<< result.right << '\t' << result.status << '\t' << result.reason << '\t'
			<< formatDouble(result.identity) << '\t' << formatDouble(result.queryCoverage) << '\t'
			<< formatDouble(result.targetCoverage) << '\t' << formatDouble(result.relationScore) << '\t'
			<< formatDouble(result.scoreMargin) << '\t' << formatDouble(result.lengthRatio) << '\t'
			<< (result.reciprocalBest ? 1 : 0) << '\t' << (result.baseAligned ? 1 : 0) << '\t'
			<< result.baselineSource << '\t' << result.baselineGenes << '\t'
			<< result.informativeSamples << '\t' << result.leftSupport << '\t'
			<< result.rightSupport << '\t' << result.jointSupport << '\t'
			<< formatDouble(result.jointFraction) << '\t'
			<< formatDouble(result.abundanceOverlapFraction) << '\t'
			<< formatDouble(result.normalizedCorrelation) << '\t'
			<< formatDouble(result.mergedCopyMedian) << '\t'
			<< formatDouble(result.mergedCopyMad) << '\t'
			<< formatDouble(result.mergedMedianError) << '\t'
			<< formatDouble(result.bestIndividualError) << '\t'
			<< formatDouble(result.fitGain) << '\n';
	}
	closeOutputStream(detail);

	std::ostream* concat = openOutputStream(concatPath);
	if (concat == nullptr) { cerr2("Could not write " + concatPath + "\n", 3); }
	for (const MosaicResult& result : results) {
		if (result.status == "confirmed") {
			*concat << result.left << ',' << result.right << '\n';
		}
	}
	closeOutputStream(concat);

	std::map<std::string, std::size_t> reasonCounts;
	std::size_t confirmed = 0;
	for (const MosaicResult& result : results) {
		++reasonCounts[result.reason];
		if (result.status == "confirmed") { ++confirmed; }
	}
	std::size_t ambiguousGenes = 0;
	for (const auto& gene : genes) {
		if (gene.second.ambiguousMgs) { ++ambiguousGenes; }
	}
	std::ostream* summary = openOutputStream(summaryPath);
	if (summary == nullptr) { cerr2("Could not write " + summaryPath + "\n", 3); }
	*summary << "metric\tvalue\n"
		<< "mgs_genes\t" << genes.size() << '\n'
		<< "ambiguous_mgs_genes\t" << ambiguousGenes << '\n'
		<< "matrix_samples\t" << sampleCount << '\n'
		<< "paf_rows\t" << pafStats.rows << '\n'
		<< "paf_rows_retained\t" << pafStats.retained << '\n'
		<< "paf_rows_filtered\t" << pafStats.filtered << '\n'
		<< "paf_rows_malformed\t" << pafStats.malformed << '\n'
		<< "candidate_pairs\t" << results.size() << '\n'
		<< "confirmed_pairs\t" << confirmed << '\n';
	for (const auto& reason : reasonCounts) {
		*summary << "reason." << reason.first << '\t' << reason.second << '\n';
	}
	closeOutputStream(summary);

	std::cout << "mosaic: wrote " << results.size() << " candidate pairs ("
		<< confirmed << " confirmed) to " << detailPath << '\n';
}

} // namespace

void runMosaicMode(options* opts) {
	if (opts == nullptr) { cerr2("mosaic mode received no options.\n", 3); }
	if (opts->input.empty() || opts->referenceFile.empty() || opts->pafFile.empty() ||
		opts->output.empty()) {
		cerr2("mosaic mode requires -i, -reference, -paf, and -o.\n", 3);
	}
	if (!fileExists(opts->input)) { cerr2("Mosaic matrix not found: " + opts->input + "\n", 3); }
	if (!fileExists(opts->referenceFile)) {
		cerr2("MGS/core file not found: " + opts->referenceFile + "\n", 3);
	}
	if (!fileExists(opts->pafFile)) { cerr2("PAF file not found: " + opts->pafFile + "\n", 3); }
	if (!opts->groupFile.empty() && !fileExists(opts->groupFile)) {
		cerr2("Group file not found: " + opts->groupFile + "\n", 3);
	}

	GeneMap genes;
	std::unordered_map<std::string, std::vector<std::string> > genesByMgs;
	readMgsFile(opts, genes, genesByMgs);
	readGroupFile(opts, genes);
	PafStats pafStats;
	const HitMap hits = readPaf(opts, genes, pafStats);
	std::vector<MosaicResult> results = buildSequenceCandidates(opts, genes, hits);

	std::unordered_set<std::string> wantedGenes;
	for (const MosaicResult& result : results) {
		wantedGenes.insert(result.left);
		wantedGenes.insert(result.right);
		std::unordered_map<std::string, std::vector<std::string> >::const_iterator mgs =
			genesByMgs.find(result.mgs);
		if (mgs != genesByMgs.end()) {
			wantedGenes.insert(mgs->second.begin(), mgs->second.end());
		}
	}
	std::size_t sampleCount = 0;
	const AbundanceMap abundances = readMatrixRows(opts, wantedGenes, sampleCount);
	addAbundanceEvidence(opts, genes, genesByMgs, abundances, sampleCount, results);
	writeResults(opts, pafStats, genes, results, sampleCount);
}
