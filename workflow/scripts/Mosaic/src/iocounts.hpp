/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#ifndef iocounts_hpp
#define iocounts_hpp

#include <iostream>
#include <vector>
#include <set>
#include <limits>       // MAX_UNSIGNED
#include <algorithm>    // std::all_of

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "counter.hpp"       // Counter
#include "intervals.hpp"    // Interval


namespace io {



const unsigned MAX_UNSIGNED = std::numeric_limits<unsigned int>::max();

using count::TGenomeCounts;
using count::Counter;
using interval::Interval;

/** write count table in long format to a gzip file.
  *
  * @param f_out      file name to write to (should end in .gz)
  * @param counts     count matrix incl. strand state labels
  * @param bins       Lookup table for intervals, which correspond to the columns of `counts`
  * @param chrom_name Lookup table for chromosome names, e.g. `hdr->target_name`
  * @param sample_cell_name Lookup table for pair of sample and cell name, correspond to the rows of `counts`
  */
template <typename TString, typename TVec, typename TPairVec>
bool write_counts_gzip(TString const & f_out,
                       std::vector<TGenomeCounts> const & counts,
                       std::vector<Interval> const & bins,
                       TVec const & chrom_name,
                       TPairVec const & sample_cell_name)
{
    boost::iostreams::filtering_ostream out;
    boost::iostreams::file_sink fout(f_out, std::ios_base::out | std::ios_base::binary);
    out.push(boost::iostreams::gzip_compressor());
    out.push(fout);

    if (fout.is_open()) {
        out << "chrom\tstart\tend\tsample\tcell\tc\tw\tclass" << std::endl;
        for(unsigned i = 0; i < counts.size(); ++i) {
            for (unsigned bin = 0; bin < counts[i].size(); ++bin) {
                Counter<unsigned> const & cc = counts[i][bin];
                out << chrom_name[bins[bin].chr];
                out << "\t" << bins[bin].start << "\t" << bins[bin].end;
                out << "\t" << sample_cell_name[i].first;
                out << "\t" << sample_cell_name[i].second;
                out << "\t" << cc.crick_count;
                out << "\t" << cc.watson_count;
                out << "\t" << cc.get_label();
                out << std::endl;
            }
        }
    } else {
        std::cerr << "[Error] Cannot open file: " << f_out << std::endl;
        return false;
    }
    return true;
}


/** read gzipped count table into data structures.
  *
  * @param f_in filename (.gz)
  * @param chromosomes empty vector. Write list of chromosomes found in file
  * @param sample_cell_names empty vector. Write sample and cell names found in file
  * @param bins empty vector. Write bins found in file
  */
template <typename TString, typename TPrec>
bool read_counts_gzip(TString const & f_in,
                      std::vector<std::vector<Counter<TPrec>>> & counts,
                      std::vector<std::string> & chromosomes,
                      std::vector<std::pair<std::string,std::string>> & sample_cell_names,
                      std::vector<Interval> & bins)
{
    struct Interval_ {
        std::string chr;
        int32_t start, end;
        Interval_(std::string const & chr, int32_t start, int32_t end) :
        chr(chr), start(start), end(end) {}
        bool operator==(Interval_ const & other) const {
            return chr == other.chr && start == other.start && end == other.end;
        }
        bool operator<(Interval_ const & other) const {
            return chr == other.chr ? start < other.start : chr < other.chr;
        }
    };

    // collect different factors
    std::set<std::pair<std::string,std::string>> sample_cell_set;
    std::set<Interval_> bin_set;
    std::set<std::string> chromosome_set;

    // Read Gzipped file a first time(!) to determine sample + cell names and bins
    {
        std::string last_sample;
        std::string last_cell;
        std::string last_chr;
        // Interval last_interval;     // don't need that, as I don't expect files
        // std::string last_chrom;     // to be sorted by intervals

        bool await_header = true;
        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(boost::iostreams::file(f_in));

        std::string line;
        while (in.good())
        {
            getline(in, line);
            if (line.size() < 1) continue; // ignore empty lines (e.g. last line)

            if (await_header && line == "chrom\tstart\tend\tsample\tcell\tc\tw\tclass") {
                await_header = false;
                continue;
            }
            if (await_header) {
                std::cerr << "[read_counts_gzip] No correct header in file" << std::endl;
                return false;
            }

            // split line
            std::vector<std::string> fields;
            boost::split(fields, line, boost::is_any_of("\t"));

            // if not 8 fields, error!
            if (fields.size() != 8) {
                std::cerr << "[read_counts_gzip] Line contains != 8 fields: \""  << line << "\"" << std::endl;
                return false;
            }

            std::string row_chr;
            int32_t     row_start;
            int32_t     row_end;
            std::string row_sample;
            std::string row_cell;
            //unsigned    row_w;
            //unsigned    row_c;
            //std::string row_type;

            try {
                row_chr    = fields[0];
                row_start  = std::stoi(fields[1]);
                row_end    = std::stoi(fields[2]);
                row_sample = fields[3];
                row_cell   = fields[4];
                //row_w      = std::stoi(fields[5]);
                //row_c      = std::stoi(fields[6]);
                //row_type   = fields[7];
            } catch (const std::exception&) {
                std::cerr << "[read_counts_gzip] Cannot read interger number in: " << line << std::endl;
                return false;
            }

            // Insert to sets
            if (last_sample != row_sample || last_cell != row_cell)
                sample_cell_set.insert(std::pair<std::string,std::string>(row_sample,row_cell));
            bin_set.insert(Interval_(row_chr, row_start, row_end));
            if (last_chr != row_chr)
                chromosome_set.insert(row_chr);
        }

        // All rows have been read. Now turn sample_cell_set into vectors.
        for (std::string const & c : chromosome_set)
            chromosomes.push_back(c);

        for (Interval_ const & bin : bin_set)
        {
            auto it = std::lower_bound(chromosomes.begin(), chromosomes.end(), bin.chr);
            int32_t chrom_id = (int32_t)(it - chromosomes.begin());
            assert(it != chromosomes.end());
            assert(*it == chromosomes[chrom_id]);
            bins.push_back(Interval(chrom_id, bin.start, bin.end));
        }

        for (std::pair<std::string,std::string> const & s_c : sample_cell_set)
            sample_cell_names.push_back(s_c);
    }

    if (sample_cell_names.size() > 10000) {
        std::cerr << "[read_counts_gzip] Over 10,000 samples recognized, this seems wrong (STOP)" << std::endl;
        return 9;
    }
    if (bins.size() > 3300000) {
        std::cerr << "[read_counts_gzip] Over 3.3M bins recognized, this seems wrong (STOP)" << std::endl;
        return 8;
    }

    // Create matrix of counts
    counts.resize(sample_cell_names.size());
    for (unsigned i=0; i<counts.size(); ++i)
        counts[i].resize(bins.size());

    // Read file again(!) and save all data into memory
    {
        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(boost::iostreams::file(f_in));

        std::pair<std::string,std::string> last_sample_cell;
        std::string last_cell, last_chr;

        bool await_header = true;
        while (in.good()) {
            std::string line;
            getline(in, line);
            if (line.size() < 1) continue; // ignore empty lines (e.g. last line)

            if (await_header && line == "chrom\tstart\tend\tsample\tcell\tc\tw\tclass") {
                await_header = false;
                continue;
            }
            if (await_header) {
                std::cerr << "[read_counts_gzip] No correct header in file" << std::endl;
                return 3;
            }
            // Read line
            std::vector<std::string> fields;
            boost::split(fields, line, boost::is_any_of("\t "));

            // if not 8 fields, error!
            if (fields.size() != 8) return 10;

            std::string row_chr;
            int32_t     row_start;
            int32_t     row_end;
            std::string row_sample;
            std::string row_cell;
            TPrec    row_c;
            TPrec    row_w;
            std::string row_type;

            try {
                row_chr    = fields[0];
                row_start  = std::stoi(fields[1]);
                row_end    = std::stoi(fields[2]);
                row_sample = fields[3];
                row_cell   = fields[4];
                row_c      = boost::lexical_cast<TPrec>(fields[5]);
                row_w      = boost::lexical_cast<TPrec>(fields[6]);
                row_type   = fields[7];
            } catch (const std::exception&) {
                std::cerr << "[read_counts_gzip] Cannot read interger number in: " << line << std::endl;
                return false;
            } catch (const boost::bad_lexical_cast&) {
                std::cerr << "[read_counts_gzip] Cannot read W or Crick number: " << line << std::endl;
                return false;
            }



            // Derive factor for each of the strings:
            unsigned sample_cell_id = MAX_UNSIGNED,
                     bin_id = MAX_UNSIGNED;
            int32_t chrom_id = -1;

            if (row_chr != last_chr)
                chrom_id = (int32_t)(std::lower_bound(chromosomes.begin(), chromosomes.end(), row_chr) - chromosomes.begin());

            bin_id = (unsigned)(std::lower_bound(bins.begin(), bins.end(), Interval(chrom_id, row_start, row_end)) - bins.begin());

            if (std::make_pair(row_sample, row_cell) != last_sample_cell)
                sample_cell_id = (unsigned) (std::lower_bound(sample_cell_names.begin(),
                                                              sample_cell_names.end(),
                                                              std::make_pair(row_sample, row_cell))
                                             - sample_cell_names.begin());

            Counter<TPrec> & cc = counts[sample_cell_id][bin_id];
            cc.crick_count  = row_c;
            cc.watson_count = row_w;
            cc.label        = row_type;
        }
    }


    // Validate entries, e.g. that all fields are filled
    for (unsigned i = 0; i<counts.size(); ++i) {
        for (unsigned bin = 0; bin < counts[i].size(); ++bin) {
            if (counts[i][bin].watson_count >= MAX_UNSIGNED || \
                counts[i][bin].crick_count  >= MAX_UNSIGNED || \
                (counts[i][bin].label != "WC" && counts[i][bin].label != "WW" && \
                 counts[i][bin].label != "CC" && counts[i][bin].label != "None")  )
            {
                std::cerr << "[read_counts_gzip] read matrix contains false/missing entries for ";
                std::cerr << sample_cell_names[i].first << ":" << sample_cell_names[i].second;
                std::cerr << " at bin " <<  chromosomes[bins[i].chr] << ":";
                std::cerr << bins[i].start << "-" << bins[i].end << std::endl;
                return false;
            }
        }
    }

    return true;
}

template <typename TPrec>
std::vector<unsigned> get_good_cells(std::vector<std::vector<Counter<TPrec>>> const & counts)
{
    std::vector<unsigned> good_cells;
    for (unsigned i = 0; i < counts.size(); ++i) {
        if (!std::all_of(counts[i].begin(),
                         counts[i].end(),
                         [](Counter<TPrec> const & x) { return (x.label == "None");}))
            good_cells.push_back(i);
    }
    // good_cells must be <= counts and sorted !!
    assert(std::is_sorted(good_cells.begin(), good_cells.end()));
    assert(good_cells.size() <= counts.size());
    return good_cells;
}


}
#endif /* iocounts_hpp */
