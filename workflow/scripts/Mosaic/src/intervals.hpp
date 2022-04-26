/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#ifndef intervals_hpp
#define intervals_hpp

#include <iostream>
#include <fstream>

#include <algorithm>
#include <boost/tokenizer.hpp>
#include <htslib/sam.h>

namespace interval {


/**
 * @file
 * @defgroup interval Intervals along the genome
 * Summary of how the genome is binned.
 *
 * ### Overview of binning strategies
 *
 */



/**
 * BED interval (0-based, right-exclusive)
 * @ingroup interval
 *
 * Interval has a `less` operator and an `operator<<` for pringing.
 */
struct Interval {
    int32_t chr;  /**< chromosome id, which typically refers to `chrom_names` */
    int32_t start;
    int32_t end;
    Interval() : chr(0), start(0), end(0) {}
    Interval(int32_t tid, int32_t s, int32_t e): chr(tid), start(s), end(e) {}
    bool operator<(Interval const & other) const {
        return chr == other.chr ? start < other.start : chr < other.chr;
    }
};

bool invt_less(Interval const & a, Interval const & b) {
    if (a.chr == b.chr)
        // special condition: If intervals start at the same position, prefer larger one!
        if (a.start == b.start)
            return (a.end - a.start) < (b.end - b.start);
        else
            return a.start < b.start;
    else
        return a.chr < b.chr;
};

std::ostream &operator<<(std::ostream &out, Interval const & bin)
{
    return (out << bin.chr << ":" << bin.start << "-" << bin.end);
}



/** 
 * Calculate the chromosome index.
 * @ingroup interval
 *
 * Given a list of genomic intervals, check that they are non-overlapping
 * and calculate `chrom_map`. This contains indices of the first bin of each
 * chromosome. Input Intervals need to be sorted and non-overlapping
 *
 * @return false if Intervals overlap each other or are not sorted. Otherwise true.
 * @param intervals **sorted** list of bed intervals
 * @param chrom_map Vector to write `chrom_map`. Must be resized to the numnber 
 *                  of chromosomes (e.g. `hdr->n_targets` or 
 *                 `chrom_names.size()`) beforehands.
 */
bool make_chrom_map(std::vector<Interval> const & intervals,
                    std::vector<int32_t> & chrom_map) {

    // check that intervals don't overlap
    int32_t prev = -1, j = 0;
    for (unsigned i=0; i<intervals.size(); ++i) {
        if (intervals[i].chr == prev) {
            if (intervals[i-1].end > intervals[i].start) {
                std::cerr << "[Error] Intervals overlap. This is not supported!" << std::endl;
                return false;
            }
        } else {
            while (j<=intervals[i].chr)
                chrom_map[j++] = i;
            prev = intervals[i].chr;
        }
    }
    while (j < (int32_t)chrom_map.size())
        chrom_map[j++] = (int32_t)intervals.size();

    return true;
}




/** 
 * Read intervals from a BED file and create chrom_map.
 * @ingroup interval
 *
 * Intervals are read using `read_exclude_file`, sorted and checked for overlaps
 * by `make_chrom_map`.
 *
 * @return true if file could be read and intervals are correct.
 * @param intervals vector of intervals to be written (initially empty)
 * @param chrom_map vector of chromosome numbers. Must have correct size.
 * @param filename BED file name
 * @param hdr BAM header needed to know number of chromosomes and tids.
 */
template <typename TFilename>
bool read_dynamic_bins(std::vector<Interval> & intervals,
                       std::vector<int32_t> & chrom_map,
                       TFilename const & filename,
                       bam_hdr_t* hdr,
                       bool verbose = 0)
{
    // read intervals
    if (!read_exclude_file(filename, hdr, intervals, verbose))
        return false;

    if (intervals.empty()) {
        std::cerr << "[Error] No intervals in file" << std::endl;
        return false;
    }

    // sort intervals
    std::sort(intervals.begin(), intervals.end(), invt_less);

    // check that intervals don't overlap
    return make_chrom_map(intervals, chrom_map);
}



/** 
 * Create fixed-width intervals along the genome.
 * @ingroup interval
 *
 * When regions are excluded (via `exlc`), bins always start exactly at the 
 * right end of the exclude interval. The last bin of a chromosome (or before
 * another exclude interval) can be trucated on the right.
 *
 * @return true if all worked out
 * @param intervals Vector of intervals to be written (initially empty).
 * @param chrom_map Vector of chromosome numbers. Will be resized to `n_chrom`.
 * @param binwidth Size of the bins
 * @param excl List of intervals to be exlcuded.
 * @param n_chrom Number of chromosomes.
 * @param target_len List of chromosome lengths, e.g. `hdr->target_len`.
 */
template <typename TVec>
bool create_fixed_bins(std::vector<Interval> & intervals,
                       std::vector<int32_t> & chrom_map,
                       unsigned binwidth,
                       std::vector<Interval> const & excl,
                       int32_t n_chrom,
                       TVec const & target_len)
{
    auto excl_iter = excl.begin();

    for (int32_t chrom=0; chrom < n_chrom; ++chrom)
    {
        // store chrom-pointer in chrom_map
        chrom_map[chrom] = static_cast<int32_t>(intervals.size());

        // skip excl. chromosomes "left" of this one
        while(excl_iter != excl.end() && excl_iter->chr < chrom)
            ++excl_iter;

        unsigned pos = 0;
        while (pos < target_len[chrom]) {

            // skip excl. bins left of pos
            while(excl_iter != excl.end() && excl_iter->chr == chrom && excl_iter->end <= (int32_t)pos)
                ++excl_iter;

            Interval ivl;
            ivl.chr = chrom;

            // if pos is inside an excl. interval, go to its end
            if (excl_iter != excl.end() && excl_iter->chr == chrom && (int32_t)pos >= excl_iter->start) {
                pos = excl_iter->end;

            } // if pos is ok but next interval is closer than binwidth
            else if (excl_iter != excl.end() && excl_iter->chr == chrom && (int32_t)(pos+binwidth) >= excl_iter->start) {
                ivl.start = pos;
                ivl.end   = std::min((int32_t)(excl_iter->start), (int32_t)(target_len[chrom]));
                intervals.push_back(ivl);
                pos = excl_iter->end;

            } // normal interval
            else {
                Interval ivl;
                ivl.chr = chrom;
                ivl.start = pos;
                ivl.end   = std::min(pos+binwidth, (unsigned)target_len[chrom]);
                intervals.push_back(ivl);
                pos += binwidth;
            }
        }
    }
    return true;
}



/**
 * Read a BED file.
 * @ingroup interval
 *
 * @param filename Filename.
 * @param intervals List of `Intervals` to be written.
 */
bool read_exclude_file(std::string const & filename, bam_hdr_t* hdr, std::vector<Interval> & intervals, bool verbose = 0)
{
    std::ifstream interval_file(filename.c_str(), std::ifstream::in);
    if (interval_file.is_open()) {
        while (interval_file.good()) {
            std::string line;
            getline(interval_file, line);
            typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
            boost::char_separator<char> sep(" \t,;");
            Tokenizer tokens(line, sep);
            Tokenizer::iterator tokIter = tokens.begin();
            if (tokIter!=tokens.end()) {
                std::string chrName = *tokIter++;
                int32_t tid = bam_name2id(hdr, chrName.c_str());
                if (tid >= 0 && tid < hdr->n_targets) {
                    Interval ivl;
                    ivl.chr = tid;
                    if (tokIter == tokens.end()) {// exclude whole chrom
                        ivl.start = 0;
                        ivl.end   = hdr->target_len[tid];
                    } else {
                        ivl.start = boost::lexical_cast<int32_t>(*tokIter++);
                        if (tokIter == tokens.end()) {
                            std::cerr << "[Warning] Invalid line: " << line << std::endl;
                            continue;
                        }
                        ivl.end   = boost::lexical_cast<int32_t>(*tokIter++);
                        if (ivl.end <= ivl.start) {
                            std::cerr << "[Warning] Invalid line: " << line << std::endl;
                            continue;
                        }
                    }
                    intervals.push_back(ivl);
                } else {
                    if (verbose) std::cerr << "[Warning] Unknown chromosome " << chrName << " in file " << filename << std::endl;
                }
            }
        }
        interval_file.close();
    } else {
        std::cerr << "[Error] BED file cannot be read: " << filename << std::endl;
        return false;
    }
    return true;
}




/**
 * Find (inclusive) start and end bins of an interval.
 * @ingroup interval
 *
 * Such that `bin.start` *<=* `pos` *<=* `bin.end` is true for both start and 
 * end coordinate of `where`.
 *
 * @return Pair of indices referring to the `bins` vector.
 * @param where Interval of SV.
 * @param bins Chromosmal bins (sorted).
 * @param chrom_map Mapping to first bin of each chromosome.
 */
inline std::pair<int32_t, int32_t> locate_bins(Interval const & where,
                                               std::vector<Interval> const & bins,
                                               std::vector<int32_t> const & chrom_map)
{
    // Get bins
    int start = std::upper_bound(bins.begin(), bins.end(), Interval(where.chr, where.start, where.start)) - bins.begin() - 1;
    int end   = std::lower_bound(bins.begin(), bins.end(), Interval(where.chr, where.end,   where.end)) - bins.begin() - 1;

    assert(bins[start].chr == where.chr);
    assert(bins[start].start <= where.start && bins[start].end >= where.start);
    assert(bins[end].chr == where.chr);
    assert(bins[end].end >= where.end && bins[end].end >= where.start);

    return std::make_pair(start, end);
}

inline float _left_frac(Interval const & bin, int32_t pos) {
    assert(pos >= bin.start);
    assert(pos <= bin.end);
    assert(bin.end > bin.start);
    return (float)(pos - bin.start) / (bin.end - bin.start);
}

inline float _right_frac(Interval const & bin, int32_t pos) {
    assert(pos >= bin.start);
    assert(pos <= bin.end);
    assert(bin.end > bin.start);
    return (float)(bin.end - pos) / (bin.end - bin.start);
}


/**
 * Find (inclusive) start and end bins of an interval plus fractions.
 * @ingroup interval
 *
 * Like `locate_bins` but also report what fraction of start and end bins
 * are inside the SV. A fraction of close to 1 means that this bin is nearly 
 * completely part of the SV. The first fraction refers to the right part of 
 * the start bin, the second fraction refers to the left part of the end bin.
 *
 * **Note:** When start bin == end bin, then the SV is completely inside a 
 * single bin. The SV covers then a total fraction of \f$f = f_1 - (1-f_2)$ of
 * that bin.
 *
 * @return Pair of indices + pair of fractions referring to the `bins` vector.
 * @param where Interval of SV.
 * @param bins Chromosmal bins (sorted).
 * @param chrom_map Mapping to first bin of each chromosome.
 */
inline std::pair<std::pair<int32_t,int32_t>, std::pair<float,float>>
    locate_partial_bins(Interval const & where,
                        std::vector<Interval> const & bins,
                        std::vector<int32_t> const & chrom_map)
{
    // Get bins
    std::pair<std::pair<int32_t,int32_t>, std::pair<float,float>> tuple;
    tuple.first = locate_bins(where, bins, chrom_map);
    // determine the portion of the bins that are within the SV:
    // right_frac of start bin, and left_frac of end bin
    tuple.second = std::make_pair(_right_frac(bins[tuple.first.first], where.start),
                                  _left_frac(bins[tuple.first.second], where.end));
    return tuple;
}





}
#endif /* intervals_hpp */
