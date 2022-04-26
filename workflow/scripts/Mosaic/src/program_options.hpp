/*
 Copyright (C) 2017 Sascha Meiers
 Distributed under the MIT software license, see the accompanying
 file LICENSE.md or http://www.opensource.org/licenses/mit-license.php.
 */

#ifndef program_options_h
#define program_options_h

#include <fstream>
#include <boost/program_options/errors.hpp>

auto in_range = [](int min, int max, char const * const opt_name)
{
    using boost::program_options::validation_error;
    // returns another lambda function
    return [opt_name, min, max](unsigned short v)
    {
        if(v < min || v > max){
            throw validation_error(validation_error::invalid_option_value,
                                   opt_name,
                                   std::to_string(v));
        }
    };
};


// opt::value<unsigned short>()->default_value(5)->notifier(in(0, 10, "my_opt"));

template <typename TString>
bool file_exists(TString const & fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}




#endif
