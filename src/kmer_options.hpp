/**
 * \file kmer_options.hpp
 * \author Gokhan Tanisik
 * \brief Defines options and related functions for the kmer_counting problem.
 */


#ifndef kmer_options_hpp
#define kmer_options_hpp

#include "BloomFilter.hpp"
#include "FastQFileStreamer.hpp"
#include <string>

namespace kmers
{
    /**
      \brief Defines all of the command line argument options and parser for those options.
     */
    class Options
    {
    public:
        
        /// Checks the validity of the options
        bool checkOptions();
        
        /// Prints the usage for this program
        void printUsage() const;
        
        /// Parses the command line arguments
        void parse(int argc, char *argv[]);
        
    private:
        /// prints the given error and returns false.
        bool printErrorAndReturnFalse(const std::string&& error);
        
    public:
        std::string fileName; //!< FASTQ file to process
        int k = 0; //!< size of k-mer
        int t = 0; //!< most frequent "t" k-mers to list as output
        bool forceBloomFilter = false; //!< forces the program to use bloom filter kmer counter solution always.
        kmers::FastQFileStreamerOptions streamerOptions; //!< options for parsing fastq file
        kmers::BloomFilterOptions bfOptions; //!< options for Bloom Filter
    };
}


#endif /* kmer_options_hpp */
