/**
 * \file kmer_options.cpp
 * \author Gokhan Tanisik
 * \brief Defines options and related functions for the kmer_counting problem.
 */


#include "kmer_options.hpp"

#include <getopt.h>
#include <iostream>

using namespace kmers;

bool Options::printErrorAndReturnFalse(const std::string&& error)
{
    std::cerr << error << std::endl;
    return false;
}

bool Options::checkOptions()
{
    if(fileName.empty())
        return printErrorAndReturnFalse("A FASTQ file must be provided with -f option such as: -f <file_name>");
    if(k == 0)
        return printErrorAndReturnFalse("-kmersize (or -k) option is mandatory and it's value must be between [1,30]!");
    if(k < 1 || k > 30)
        return printErrorAndReturnFalse("kmersize must be between [1, 30]!");
    if(t == 0)
        return printErrorAndReturnFalse("-t option is mandatory and it's value must be greater than 0!");
    if(streamerOptions.maxReadSize == 0)
        return printErrorAndReturnFalse("Max read size (-R option) must be a valid integer number!");
    if(bfOptions.maxMemorySizeMB == 0)
        return printErrorAndReturnFalse("Max memory size (-M option) for number of bits in the bloom filter must be a valid integer number!");
    if(bfOptions.expectedFPRate == 0.0f)
        return printErrorAndReturnFalse("Expected false positive rate (-F option) for Bloom Filter must be a valid real number and greater than 0.0f!");
    if(streamerOptions.maxReadSize < k)
        return printErrorAndReturnFalse("Max read size (-R option) can not be less then kmersize!");
    if(bfOptions.numHashFunctions == 0)
        return printErrorAndReturnFalse("Number of hash functions in the bloom filter must be -1(auto calculated) or greater then 0!");
    if(bfOptions.layerCount < 1 || bfOptions.layerCount > 3)
        return printErrorAndReturnFalse("Number of bloom filter layers can be 1,2 or 3!.");
    if(bfOptions.cascadedKmerCountDecimation > 1.0f && bfOptions.cascadedKmerCountDecimation < 0.01f)
        return printErrorAndReturnFalse("Bloom filter size decimation rate must be in the interval [0.01, 1.0]");
    if(!(bfOptions.bitsPerEntry == 1 || bfOptions.bitsPerEntry == 8 || bfOptions.bitsPerEntry == 16 ||
         bfOptions.bitsPerEntry == 32 || bfOptions.bitsPerEntry == 64))
        return printErrorAndReturnFalse("Bits per entry (for single bloom filter version) must be between [1,5]");
    return true;
}

void Options::printUsage() const
{
    std::cout << "\n\n------------------------------------------------------------ \n" <<
    "This program counts the k-mers in a give FASTQ file. And prints the max frequent N k-mers in the file.\n" <<
    "  Usage: <program_name> [options]\n"
    "  -f <NAME>, --filename <NAMe>         FASTQ file to process (REQUIRED).\n"
    "  -k <SIZE>, --kmersize <SIZE>         Size of k-mers (value of k) to count (REQUIRED).\n"
    "  -t <COUNT>, --topcount <COUNT>       Number of most frequent kmers (REQUIRED).\n"
    "  -R <LENGTH>, --maxread <LENGTH>      Length of one read operation from file. Default value is 65536.\n"
    "  -M <SIZE>, --maxmemory <SIZE>        Max memory size(MB) allowed for the bloom filter. \n"
    "                                       Default value is 8192Mb.\n"
    "  -F <RATE>, --fprate <RATE>           Expected false positive rate for the bloom filter.\n"
    "                                       It may be increased depending on the -M option. Default value is 0.2\n"
    "  -K <COUNT>, --kmercount <COUNT>      Number of kmers in the file.\n"
    "                                       Decreases the execution time when it's provided.      \n"                      
    "  -h <COUNT>, --hashfunctions <COUNT>  Number of hash functions to be used in the bloom filter.\n"
    "                                       Default value is auto-calculated.\n"
    "  -L <COUNT>, --bflayers <COUNT>       Number of bloom filter layers. Valid values are: [1,2,3].\n"
    "                                       Default value is 1.\n"
    "  -D <RATE>, --decimation <RATE>       Decimation rate for cascaded bloom filters. \n"
    "                                       Used to decrease memory usage for each layer. Default value is 0.25\n"
    "  -B, --forcebloom                     Forces to use bloom filter solution in all cases. \n"
    "                                       Default value is false.\n"
    "  -b <COUNT>, --bitsperentry <COUNT>   Bits per entry to be used in single bloom filter version.\n"
    "                                       Valid values are [1,8,16,32,64]. Default value is 8." << std::endl;
}

void Options::parse(int argc, char *argv[])
{
    static struct option long_options[] =
    {
        {"filename",     required_argument, 0, 'f'},
        {"kmersize",     required_argument, 0, 'k'},
        {"topcount",     required_argument, 0, 't'},
        {"maxread",      required_argument, 0, 'R'},
        {"kmercount",    required_argument, 0, 'K'},
        {"fprate",       required_argument, 0, 'F'},
        {"maxmemory",    required_argument, 0, 'M'},
        {"hashcount",    required_argument, 0, 'h'},
        {"bflayers",     required_argument, 0, 'L'},
        {"decimation",   required_argument, 0, 'D'},
        {"forcebloom",   no_argument,       0, 'B'},
        {"bitsperentry", required_argument, 0, 'b'},
        {0, 0, 0, 0}
    };
    
    int opt;
    int option_index;
    
    while((opt = getopt_long(argc, argv, "f:k:t:R:M:F:K:h:L:D:Bb:", long_options, &option_index)) != -1)
    {
        switch (opt) {
            case 'f':
                fileName = std::string(optarg);
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 't':
                t = atoi(optarg);
                break;
            case 'R':
                streamerOptions.maxReadSize = atoi(optarg);
                break;
            case 'M':
                bfOptions.maxMemorySizeMB = atoi(optarg);
                break;
            case 'F':
                bfOptions.expectedFPRate = atof(optarg);
                break;
            case 'K':
                streamerOptions.numKmers = atol(optarg);
                break;
            case 'h':
                bfOptions.numHashFunctions = atoi(optarg);
                break;
            case 'L':
                bfOptions.layerCount = atoi(optarg);
                break;
            case 'D':
                bfOptions.cascadedKmerCountDecimation = atof(optarg);
                break;
            case 'B':
                forceBloomFilter = true;
                break;
            case 'b':
                bfOptions.bitsPerEntry = atoi(optarg);
                break;
            case '?':
                printUsage();
                exit(EXIT_SUCCESS);
            default: /* '?' */
                printUsage();
                exit(EXIT_FAILURE);
        }
    }
    
    if(!checkOptions())
    {
        printUsage();
        exit(EXIT_FAILURE);
    }
}
