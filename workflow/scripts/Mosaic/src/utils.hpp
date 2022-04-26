/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#ifndef utils_hpp
#define utils_hpp


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>



inline uint32_t alignmentLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
      if (bam_cigar_op(cigar[i]) == BAM_CMATCH) alen+=bam_cigar_oplen(cigar[i]);
    return alen;
  }



struct CellInfo {
    unsigned median_bin_count;      /* raw counts */
    float mean_bin_count;           /* after removal of bad bins */
    std::string sample_name;
    std::string cell_name;
    std::string bam_file;
    int32_t id;                     /* position in conf.f_in */
    float nb_p, nb_r, nb_a;         /* NB parameters */
    bool pass_qc;                   /* set to false if no (good) SS library */
    
    /* Alignment statistics */
    unsigned n_mapped, n_pcr_dups, n_supplementary, n_low_mapq, n_read2s;
    unsigned n_counted, n_unmap;
    CellInfo() : median_bin_count(0), mean_bin_count(0), id(-1), nb_p(0),
                 nb_r(0), nb_a(0), pass_qc(true), n_mapped(0), n_pcr_dups(0),
                 n_supplementary(0), n_low_mapq(0), n_read2s(0), n_counted(0),
                 n_unmap(0)
    {}
};

struct SampleInfo {
    std::vector<float> means;
    std::vector<float> vars;
    float p;
    SampleInfo() : p(0.33) {}
};


/**
 * Write CellInfo to file
 */
bool write_cell_info(std::string const & f_out,
                     std::vector<CellInfo> const & cells)
{
    std::ofstream out(f_out);
    if (out.is_open()) {
        out << "# sample:  Sample (has multiple cells)" << std::endl;
        out << "# cell:    Name of the cell." << std::endl;
        out << "# mapped:  Total number of reads seen" << std::endl;
        out << "# suppl:   Supplementary, secondary or QC-failed reads (filtered out)" << std::endl;
        out << "# dupl:    Reads filtered out as PCR duplicates" << std::endl;
        out << "# mapq:    Reads filtered out due to low mapping quality" << std::endl;
        out << "# read2:   Reads filtered out as 2nd read of pair" << std::endl;
        out << "# good:    Reads used for counting." << std::endl;
        out << "# pass1:   Enough coverage? If false, ignore all columns from now" << std::endl;
        out << "# nb_p:    Negative Binomial parameter p. Constant for one sample." << std::endl;
        out << "# nb_r:    Negative Binomial parameter r. We use NB(p,r/2) * NB(p,r/2) in WC states, but NB(p,(1-a)*r)*NB(p,a*r) in WW or CC states." << std::endl;
        out << "# nb_a:    Negative Binomial parameter a (alpha) used for zero expectation (see above)." << std::endl;
        out << "# bam:     Bam file of this cell" << std::endl;
        out << "sample\tcell\tmedbin\tmapped\tsuppl\tdupl\tmapq\tread2\tgood\tpass1\tnb_p\tnb_r\tnb_a\tbam" << std::endl;

        for (CellInfo const & cell : cells) {
            out << cell.sample_name << "\t";
            out << cell.cell_name << "\t";
            out << cell.median_bin_count << "\t";
            out << cell.n_mapped << "\t";
            out << cell.n_supplementary << "\t";
            out << cell.n_pcr_dups << "\t";
            out << cell.n_low_mapq << "\t";
            out << cell.n_read2s << "\t";
            out << cell.n_counted << "\t";
            out << cell.pass_qc << "\t";
            out << cell.nb_p << "\t";
            out << cell.nb_r << "\t";
            out << cell.nb_a << "\t";
            out << cell.bam_file << std::endl;
        }
    } else {
        std::cerr << "[Warning] Cannot write to " << f_out << std::endl;
        return false;
    }
    return true;
}




template <typename TReturn>
using TMedianAccumulator = boost::accumulators::accumulator_set<TReturn, boost::accumulators::stats<boost::accumulators::tag::median> >;

template <typename TReturn>
using TMeanVarAccumulator = boost::accumulators::accumulator_set<TReturn, boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::variance> >;


double sum(std::vector<double> const & vec)
{
    double sum = 0;
    for (double d : vec)
        sum += d;
    return(sum);
}


// from Delly
inline bool get_RG_tag(std::string const & tag,
                       std::string const& header, 
                       std::string& output,
                       bool allow_multiple_matches = false)
{
    std::set<std::string> values;
    typedef std::vector<std::string> TStrParts;
    TStrParts lines;
    boost::split(lines, header, boost::is_any_of("\n"));
    TStrParts::const_iterator itH = lines.begin();
    TStrParts::const_iterator itHEnd = lines.end();
    for(;itH!=itHEnd; ++itH) {
        if (itH->find("@RG")==0) {
            TStrParts keyval;
            boost::split(keyval, *itH, boost::is_any_of("\t "));
            TStrParts::const_iterator itKV = keyval.begin();
            TStrParts::const_iterator itKVEnd = keyval.end();
            for(;itKV != itKVEnd; ++itKV) {
                size_t sp = itKV->find(":");
                if (sp != std::string::npos) {
                    std::string field = itKV->substr(0, sp);
                    if (field == tag) {
                        std::string value = itKV->substr(sp+1);
                        values.insert(value);
                    }
                }
            }
        }
    }
    if (values.size() == 1) {
        output = *(values.begin());
        return true;
    } else if (values.size() > 1) {
        if (allow_multiple_matches) {
            output = *(values.begin());
            return true;
        } else {
            output = "";
            return false;
        }
    } else {
        output = "";
        return false;
    }
}






/**
 * template <class InputIter, class ForwardIter, class BinaryPredicate, class BinaryFunction> ForwardIter reduce_adjacent(InputIter first, InputIter last,ForwardIter result, BinaryPredicate _bool_mergeable, BinaryFunction _merge_func)
 * STL-style algorithm to merge/reduce adjacent elements (if they are mergeable) and return a shortened list
 *
 * @param first input iterator.
 * @param end input iterator end.
 * @param result output iterator. Can point to *first to do inplace operations.
 * @param _bool_mergeable Function stating whether two consecutive elements are mergeable.
 * @param _merge_func Funciton to merge to consecutive elements. Must return same type.
 */
template <class InputIter, class ForwardIter, class BinaryPredicate, class BinaryFunction>
ForwardIter reduce_adjacent(InputIter first, InputIter last,
                            ForwardIter result,
                            BinaryPredicate _bool_mergeable,
                            BinaryFunction _merge_func)
{
    if (first == last) return result; // skip empty container

    auto elem = *first; // copy first element anyways.
    while (++first != last) { // loop skippes first elem
        if ( !_bool_mergeable(elem, *first) ) {
            *(result++) = elem;
            elem = *first;
        } else {
            elem = _merge_func(elem, *first);
        }
    }
    *(result++) = elem;
    return result;
}


#endif /* utils_hpp */
