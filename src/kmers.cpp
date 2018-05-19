/**
 * \file kmers.cpp
 * \author Gokhan Tanisik
 * \brief Main file for the k-mer counting problem.
 * 
 */
 
 /*! \mainpage k-mer Counter
 
  \section intro_sec Introduction
 
  This is my solution for the given k-mer counting problem.
  
  It implements three basic solutions to the problem:
	 - A simple lookup table counter which is used for small k values. It's also used as a baseline method to validate the second one.
	 - Bloom Filter counter which is a simple implementation of the <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-333">BFCounter</a> paper.
     - Cascaded Bloom Filter counter which resembles to <a href="https://www.hindawi.com/journals/mpe/2013/516298/">this</a> paper. 
        It's slower then the single layer Bloom Filter version (with 8 bits per entry) but it can use less memory and may be more accurate depending on the selected options. Only 3 layers is allowed for simplicity and each layer will use 1,8 and 16 bits per entry respectively.
     
  \section install_sec Installation
   A simple make file is provided with the source code. Just run the makefile to create the executable.
	\note Source code uses Google Sparse Hash library. To be able to compile the code, you must first install it on your system. 
	I can install it on my mac using <em>brew install google-sparsehash</em> command.
	 And also update the Include Search Directory in the make file(second line) to the correct location. 
  \section usage Usage
  
  \verbatim
  Usage: <program_name> [options]
	
      -f <NAME>, --filename <NAMe>         FASTQ file to process (REQUIRED).
      
      -k <SIZE>, --kmersize <SIZE>         Size of k-mers (value of k) to count (REQUIRED).
      
      -t <COUNT>, --topcount <COUNT>       Number of most frequent kmers (REQUIRED).
      
      -R <LENGTH>, --maxread <LENGTH>      Length of one read operation from file. Default value is 65536.
      
      -M <SIZE>, --maxmemory <SIZE>        Max memory size(MB) allowed for the bloom filter. 
                                           Default value is 8192Mb.
                                           
      -F <RATE>, --fprate <RATE>           Expected false positive rate for the bloom filter.
                                           It may be increased depending on the -M option. Default value is 0.2
                                           
      -K <COUNT>, --kmercount <COUNT>      Number of kmers in the file.
                                           Decreases the execution time when it's provided.
                                           
      -h <COUNT>, --hashfunctions <COUNT>  Number of hash functions to be used in the bloom filter.
                                           Default value is auto-calculated.
                                           
      -L <COUNT>, --bflayers <COUNT>       Number of bloom filter layers. Valid values are: [1,2,3].
                                           Default value is 1.
                                           
      -D <RATE>, --decimation <RATE>       Decimation rate for cascaded bloom filters. 
                                           Used to decrease memory usage for each layer. Default value is 0.25
                                           
      -B, --forcebloom                     Forces to use bloom filter solution in all cases. 
                                           Default value is false.
                                           
      -b <COUNT>, --bitsperentry <COUNT>   Bits per entry to be used in single bloom filter version.
                                           Valid values are [1,8,16,32,64]. Default value is 8.
 \endverbatim
  
  It uses lookup table solution for kmer size less then 15. Yet it's possible to choose a bloom filter solution with <em>-B</em> option.
  
  \subsection sample_usages Sample Usages
    \code
     ./kmer_counter -f myfile.fastq -k 32 -t 25
     ./kmer_counter --filename myfile.fastq --kmersize 32 --topcount 25
    \endcode
    Simple usages with short and long options. Single layer bloom filter solution is chosen because k > 14.

    \code
     ./kmer_counter -f myfile.fastq -k 5 -t 25
    \endcode
    Same as above but uses lookup table since k <= 14.

    \code
     ./kmer_counter -f myfile.fastq -k 5 -t 25 -B --bitsperentry 8
    \endcode
    Force to use bloom filter solution. Single layer filter is used as default. Also tells the filter to use 8 bits counter for each entry.
  
    \code
     ./kmer_counter -f myfile.fastq -k 5 -t 25 -B --bflayers 2 --kmercount 1000 --fprate 0.1 -h 3 --bitsperentry 8
    \endcode
    Forces to use a 2 layers (--bflayer 2) of bloom filters, and tells the filters to use 3 hash(-h 3) functions.
    Also tells that there are total 1000 kmers in the file. The last option \a bitsperentry is omitted because it is only used for single layered solution.
  
    code
    ./kmer_counter -f myfile.fastq -k 30 -t 25 --bflayers 3 --kmercount 1000 --decimation 0.2
    \endcode
    Uses a 3 layers bloom filter solution. Tells that there are 1000 kmers in the file. The first bloom filter will be initialized for 1000 kmers.
      Since the decimation value is 0.25, the second filter will be initialized for 1000*0.2 = 200 kmers, and the third filter with 250*0.2 = 40 kmers. The filters will use 1,8 and 16 bits per entry respectively.
  
  \section discussion Discussion
  
- k-mer size is limited to 31 to be able to use uin64_t values representing k-mers.
  It's not 32, because extra bits are needed for hash set to identify empty and deleted keys when selecting top frequent k-mers.
- The topcount parameter is expected to be a very small value with respect to the kmer count, such as 25, 256 or 1024.
 */

#include "kmer_options.hpp"
#include <iostream>

#include <sparsehash/sparse_hash_map>
#include <sparsehash/dense_hash_set>
#include <iostream>
#include <functional>
#include <queue>
#include <map>
#include <set>
#include <chrono>

// used to compress the nucleotides into 2 bits
static unsigned char lookup[128] = {0};

using pair = std::pair<uint64_t, uint64_t>;
auto cmp = [](const pair &left, const pair &right) { return left.second > right.second;};
using PriorityQueue = std::priority_queue<pair, std::vector<pair>, decltype(cmp)>;

/// This function returns the result of raising \a base to the power \a power.
uint64_t upow(uint64_t base, uint64_t power)
{
    uint64_t result = 1;
    while(power > 0)
    {
        result *= base;
        --power;
    }
    return result;
}

/**
 Used to store the counts of overflowed kmers. For example, for a 1 bit Bloom filter, if a kmer is seen twice, it is added into the overflow map.
 */
template <class CompressedKmerType, class CountType>
using SparseHashMap = google::sparse_hash_map<CompressedKmerType, CountType, std::hash<CompressedKmerType>, std::equal_to<CompressedKmerType>>;

/**
 \brief Defines a compressed kmer sequnceç
 */
template <class Type>
union compressed_kmer
{
    static void compress(compressed_kmer &kmer, const char *const rawKmer, const int &k)
    {
        kmer.compressed = lookup[static_cast<unsigned char>(rawKmer[0])];
        for(int kk = 1; kk < k; kk++)
        {
            kmer.compressed = kmer.compressed << 2;
            kmer.compressed = kmer.compressed | lookup[static_cast<unsigned char>(rawKmer[kk])];
        }
    }
    
    inline static int size(const int &k)
    {
        // 2bits per nucleotide=> k={1,2,3,4}->1byte, k={5,6,7,8}->2bytes, k={9,10,11,12}->3bytes...
        return ((k*2)+7) / 8;
    }
    
    inline static std::string toString(Type compressed, const int &k)
    {
        static char reverse_lookup[4] = {'A', 'C', 'G', 'T'};
        std::string str;
        for(int kk = 0; kk < k; kk++)
        {
            str = reverse_lookup[compressed & 0x3] + str;
            compressed = compressed >> 2;
        }
        return str;
    }
    
    Type compressed = 0;
    uint8_t bytes[sizeof(Type)];
};

/**
 \brief Returns the total count of given kmer in all of the bloom filters.
 */
template <class CompressedKmerType>
inline uint64_t countKmerFromBfCounters(std::vector<kmers::IBloomFilter*> &bfs, const compressed_kmer<CompressedKmerType>& kmer, const int& compressedKmerSize)
{
    uint64_t sum = 0;
    uint64_t count = 0;
    for(auto& bf: bfs)
    {
        if((count = bf->count(kmer.bytes, compressedKmerSize)) == 0)
            break;
        sum += count;
    }
    return sum;
}

/**
 \brief Selects top frequent t kmers from the container and pushes them into the priority queue. Since the queue does not provide find function, to know which kmers are in the queue, the topKmersSet is used.
 */
template <class AssociativeContainer, class CompressedKmerType>
void selectTopFrequentKmers(const AssociativeContainer& container, int T, PriorityQueue &queue, google::dense_hash_set<CompressedKmerType> topKmersSet)
{
    for(auto& it: container)
    {
        if(queue.size() < T)
        {
            queue.emplace(std::make_pair(it.first, it.second));
            topKmersSet.insert(it.first);
        }
        else if(queue.top().second < it.second)
        {
            topKmersSet.erase(queue.top().first);
            topKmersSet.insert(it.first);
            queue.pop();
            queue.emplace(std::make_pair(it.first, it.second));
        }
    }
}

/**
 \brief Selects top frequent kmers from the container
 */
template <class AssociativeContainer>
void selectTopFrequentKmers(const AssociativeContainer& container, int T, PriorityQueue &queue, uint64_t offset = 0)
{
    for(auto& it: container)
    {
        uint64_t value = offset + it.second;
        if(queue.size() < T)
        {
            queue.emplace(std::make_pair(it.first, value));
        }
        else if(queue.top().second < value)
        {
            queue.pop();
            queue.emplace(std::make_pair(it.first, value));
        }
    }
}

/**
 \brief Selects the top frequent kmers from the vector. Index corresponds to the kmer and value at the index is the count of that kmer.
 */
template <class Type>
void selectTopFrequentKmers(std::vector<Type>& vector, int T, PriorityQueue &queue)
{
    for(auto it = vector.begin(); it != vector.end(); ++it)
    {
        if(*it > 0)
        {
            auto idx = it - vector.begin();
            if(queue.size() < T)
            {
                queue.emplace(std::make_pair(idx, *it));
            }
            else if(queue.top().second < *it)
            {
                queue.pop();
                queue.emplace(std::make_pair(idx, *it));
            }
        }
    }
}

/**
 \brief Reads the whole file from the beginning and finds the most frequent kmers that are stored in the bloom filters.
 To avoid repeating kmers, topKmersSet is used.
 */
template <class CompressedKmerType>
void selectTopTFromBFCounters(kmers::FastQFileStreamer &fileStreamer, const kmers::Options &options ,std::vector<kmers::IBloomFilter*> &&bfs, PriorityQueue & queue, google::dense_hash_set<CompressedKmerType>& topKmersSet)
{
    using compressed_kmer = compressed_kmer<CompressedKmerType>;
    const unsigned int compressedKmerSize = compressed_kmer::size(options.k);
    compressed_kmer kmer;
    for(kmers::FastQFileStreamer::KmerIterator it = fileStreamer.begin(); it; ++it)
    {
        compressed_kmer::compress(kmer, *it, options.k);
        if(topKmersSet.find(kmer.compressed) == topKmersSet.end())
        {
            auto count = countKmerFromBfCounters(bfs, kmer, compressedKmerSize);
            if(queue.size() < options.t)
            {
                queue.push(std::make_pair(kmer.compressed, count));
                topKmersSet.insert(kmer.compressed);
            }
            else
            {
                if(queue.top().second < count)
                {
                    auto topIt = topKmersSet.find(queue.top().first);
                    if(topIt != topKmersSet.end())
                        topKmersSet.erase(topIt);
                    queue.pop();
                    topKmersSet.insert(kmer.compressed);
                    queue.push(std::make_pair(kmer.compressed, count));
                }
            }
        }
    }
}



/**
  \brief Prints the top frequent kmers which are stored in the queue.
 */
static void printTopFrequentKmers(PriorityQueue& queue, int k, const std::string& fileName)
{
    std::cout << "most frequent " << queue.size() << " kmers in the file:" << fileName<< ", for k=" << k << "(#  kmer  frequency):" << std::endl;
    int counter = 0;
    while(!queue.empty())
    {
        std::cout << "   " << ++counter << "  " << compressed_kmer<uint64_t>::toString(queue.top().first, k) << " --> " << queue.top().second << std::endl;
        queue.pop();
    }
}

/**
  \brief Kmer counter with single layer bloom filter.
  
  CompressedKmerType is 
*/
template <class CompressedKmerType, class BFCountType, class OverflowDataType>
void countKmersBloomFilter(const kmers::Options& options, kmers::FastQFileStreamer& fileStreamer)
{
    using compressed_kmer = compressed_kmer<CompressedKmerType>;
    const unsigned int compressedKmerSize = compressed_kmer::size(options.k);
    
    kmers::BloomFilter<BFCountType> bf(options.k, fileStreamer.kmerCount(), std::move(options.bfOptions));
    SparseHashMap<CompressedKmerType, OverflowDataType> hashMap;
    
    compressed_kmer kmer;
    
    for(kmers::FastQFileStreamer::KmerIterator it = fileStreamer.begin(); it; ++it)
    {
#ifdef DEBUG
        /*char *seq = *it;*/ static uint64_t counter; ++counter;
#endif
        
        compressed_kmer::compress(kmer, *it, options.k);
        if(!bf.add(kmer.bytes, compressedKmerSize))
        {// try to add the kmer into the bloom filter. if it fails, use the hash table.
            auto it = hashMap.find(kmer.compressed);
            if(it == hashMap.end())
            {// the kmer isn't added into the table before?
                hashMap[kmer.compressed] = static_cast<OverflowDataType>(bf.maxValuePerEntry()) + 1; // hash map stores the total count for each kmer entry
            }
            else
            {// it is already in the table, increase the counter for the kmer.
                it->second++;
            }
        }
    }
    
    std::cout << "Done counting the kmers! Selecting the top frequent ones to list.." << std::endl;
    {// select and print top frequent kmers
        PriorityQueue queue(cmp);
        google::dense_hash_set<CompressedKmerType> topKmersSet;
        topKmersSet.set_empty_key(std::numeric_limits<CompressedKmerType>::max());
        topKmersSet.set_deleted_key(std::numeric_limits<CompressedKmerType>::max()-1);
        selectTopFrequentKmers(hashMap, options.t, queue, topKmersSet);
        
        if(queue.size() < options.t)
        {
            selectTopTFromBFCounters(fileStreamer, options, {&bf}, queue, topKmersSet);
        }
        
        printTopFrequentKmers(queue, options.k, options.fileName);
    }
}

/**
  \brief Inits the bloom filters for cascaded kmer counter solution.
 */
void initCascadedBloomFilters(std::vector<kmers::IBloomFilter*>& bfs, uint64_t kmerCount, const kmers::Options& options)
{
    // init filters
    auto bfOptions = options.bfOptions;
    if(bfOptions.layerCount >= 1)
    {
        //bfOptions.maxMemorySizeMB *= options.bfOptions.cascadedKmerCountDecimation;
        //kmerCount *= options.bfOptions.cascadedKmerCountDecimation;
        bfs.push_back(new kmers::BloomFilter<bool>(options.k, kmerCount, std::move(bfOptions)));
        
        if(bfOptions.layerCount >= 2)
        {
            //bfOptions.maxMemorySizeMB *= options.bfOptions.cascadedKmerCountDecimation;
            kmerCount *= options.bfOptions.cascadedKmerCountDecimation;
            bfs.push_back(new kmers::BloomFilter<uint8_t>(options.k, kmerCount, std::move(bfOptions)));
            
            if(bfOptions.layerCount >= 3)
            {
                //bfOptions.maxMemorySizeMB *= options.bfOptions.cascadedKmerCountDecimation;
                kmerCount *= options.bfOptions.cascadedKmerCountDecimation;
                bfs.push_back(new kmers::BloomFilter<uint16_t>(options.k, kmerCount, std::move(bfOptions)));
            }
        }
    }
}

/**
  \brief Kmer counter with cascaded bloom filters.
 */
template <class CompressedKmerType, class OverflowDataType>
void countKmersCascadeBloomFilter(const kmers::Options& options, kmers::FastQFileStreamer& fileStreamer)
{
    using compressed_kmer = compressed_kmer<CompressedKmerType>;
    const unsigned int compressedKmerSize = compressed_kmer::size(options.k);
    
    std::vector<kmers::IBloomFilter*> bfs;
    initCascadedBloomFilters(bfs, fileStreamer.kmerCount(), options);
    
    SparseHashMap<CompressedKmerType, OverflowDataType> hashMap;
    
    OverflowDataType bfsMaxValue = 0;
    for(auto& bf: bfs)
        bfsMaxValue += bf->maxValuePerEntry();
    
    compressed_kmer kmer;
    
    for(kmers::FastQFileStreamer::KmerIterator it = fileStreamer.begin(); it; ++it)
    {
#ifdef DEBUG
        /*char *seq = *it; */static uint64_t counter; ++counter;
#endif
        compressed_kmer::compress(kmer, *it, options.k);
        bool success = false;
        
        for(auto& bf : bfs)
        {
            if((success = bf->add(kmer.bytes, compressedKmerSize)))
                break;
        }
        
        if(!success)
        {
            auto it = hashMap.find(kmer.compressed);
            if(it == hashMap.end())
            {// the kmer isn't added into the table before?
                hashMap[kmer.compressed] = bfsMaxValue + 1; // hash map stores the total count for each kmer entry
            }
            else
            {// it is already in the table, increase the counter for the kmer.
                it->second++;
            }
        }
    }
    
    std::cout << "Done counting the kmers! Selecting the top frequent ones to list.." << std::endl;
    compressed_kmer::compress(kmer, "TTTTTTTTTTTT", 24);
    {// select and print top frequent kmers
        PriorityQueue queue(cmp);
        google::dense_hash_set<CompressedKmerType> topKmersSet;
        topKmersSet.set_empty_key(std::numeric_limits<CompressedKmerType>::max());
        topKmersSet.set_deleted_key(std::numeric_limits<CompressedKmerType>::max()-1);
        selectTopFrequentKmers(hashMap, options.t, queue, topKmersSet);
        
        if(queue.size() < options.t)
        {
            selectTopTFromBFCounters(fileStreamer, options, std::move(bfs), queue, topKmersSet);
        }
        
        printTopFrequentKmers(queue, options.k, options.fileName);
    }
}

/*
  \brief Merges two priority queues into one.
 **/
void mergeQueues(PriorityQueue &queue1, PriorityQueue &queue2, PriorityQueue &merged, int maxSize);

/**
  \brief A simple lookup table counter for small k-mers.
 
   A lookup table is used to count the k-mers. Value of the compressed kmer is used as index to the table.
   To avoid high memory usage, small integral types (CountType template parameter) can be used. 
   In that case, if an overflow occurs for highly frequent k-mers, an extra map is used.
   
   This solution is very fast for small k-mers, and also used as a baseline method for the Bloom Filter solutions.
*/
template <class CompressedType, class CountType, class OverflowType>
void countShortKmers(const kmers::Options& options, kmers::FastQFileStreamer& fileStreamer)
{
    using compressed_kmer = compressed_kmer<CompressedType>;
    const CountType maxCount = static_cast<CountType>(upow(2, (8 * sizeof(CountType) - 1 ))) + static_cast<CountType>(upow(2, (8 * sizeof(CountType) - 1 )) - 1);
    
    std::vector<CountType> lookupTable;
    lookupTable.resize(upow(4, options.k), 0);
    
    SparseHashMap<CompressedType, OverflowType> overflows; //std::map<CompressedType, OverflowType> overflows;
    
    compressed_kmer kmer;
    
    for(kmers::FastQFileStreamer::KmerIterator it = fileStreamer.begin(); it; ++it)
    {
#ifdef DEBUG
        /*char *seq = *it;*/ static uint64_t counter; ++counter;
#endif
        compressed_kmer::compress(kmer, *it, options.k);
        
        if(lookupTable[kmer.compressed] < maxCount)
        {
            lookupTable[kmer.compressed]++;
        }
        else
        {
            overflows[kmer.compressed]++;
        }
    }
    
    std::cout << "Done counting the kmers! Selecting the top frequent ones to list.." << std::endl;
    {// select and print top freqeunt kmers
        PriorityQueue queue1(cmp);
        selectTopFrequentKmers(overflows, options.t, queue1, maxCount);
        
        if(queue1.size() < options.t)
        {
            PriorityQueue queue2(cmp);
            selectTopFrequentKmers(lookupTable, options.t, queue2);
            PriorityQueue mergedQueue(cmp);
            mergeQueues(queue1, queue2, mergedQueue, options.t);
            printTopFrequentKmers(mergedQueue, options.k, options.fileName);
        }
        else
        {
            printTopFrequentKmers(queue1, options.k, options.fileName);
        }
    }
}

/// \brief Main entry point of the program. Parses the arguments and runs the related kmer counter algorithm.
int main(int argc, char *argv[])
{
    lookup[static_cast<uint8_t>('a')] = 0; lookup[static_cast<uint8_t>('a')] = 0;
    lookup[static_cast<uint8_t>('c')] = 1; lookup[static_cast<uint8_t>('C')] = 1;
    lookup[static_cast<uint8_t>('g')] = 2; lookup[static_cast<uint8_t>('G')] = 2;
    lookup[static_cast<uint8_t>('t')] = 3; lookup[static_cast<uint8_t>('T')] = 3;
    
    kmers::Options options;
    options.parse(argc, argv);
    
    try
    {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        
        kmers::FastQFileStreamer fileStreamer(options.k, std::move(options.fileName), std::move(options.streamerOptions));
        
        std::cout << "Total kmer count in the file:" << options.fileName << " is " << fileStreamer.kmerCount() << std::endl;
        
        bool useBloomFilter = false;
        if(options.forceBloomFilter)
            useBloomFilter = true;
        else if(options.k <= 4)
            countShortKmers<uint8_t, uint64_t, uint64_t>(options, fileStreamer);
        else if(options.k <= 8)
            countShortKmers<uint16_t, uint64_t, uint64_t>(options, fileStreamer);
        else if(options.k <= 14)
            countShortKmers<uint32_t, uint32_t, uint64_t>(options, fileStreamer);
        else
            useBloomFilter = true;
        
        if(useBloomFilter)
        {
            if(options.bfOptions.layerCount > 1)
                countKmersCascadeBloomFilter<uint64_t, uint64_t>(options, fileStreamer);
            else
            {
                switch (options.bfOptions.bitsPerEntry) {
                    case 1: // single layer bf counter with 1 bits per entry
                        countKmersBloomFilter<uint64_t, bool, uint64_t>(options, fileStreamer);
                        break;
                    case 8:// single layer bf counter with 8 bits per entry
                        countKmersBloomFilter<uint64_t, uint8_t, uint64_t>(options, fileStreamer);
                        break;
                    case 16:// single layer bf counter with 16 bits per entry
                        countKmersBloomFilter<uint64_t, uint16_t, uint64_t>(options, fileStreamer);
                        break;
                    case 32:// single layer bf counter with 32 bits per entry
                        countKmersBloomFilter<uint64_t, uint32_t, uint64_t>(options, fileStreamer);
                        break;
                    case 64:// single layer bf counter with 60 bits per entry
                        countKmersBloomFilter<uint64_t, uint64_t, uint64_t>(options, fileStreamer);
                        break;
                    default:
                        break;
                }
            }
        }
        
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " seconds!"<<std::endl;
        
    }
    catch (std::runtime_error& err){
        std::cerr << err.what() << std::endl;
    }
    
    std::vector<bool> a;
    return 0;
}


void mergeQueues(PriorityQueue &queue1, PriorityQueue &queue2, PriorityQueue &merged, int maxSize)
{
    while(!queue1.empty() && queue2.empty())
    {
        if(queue1.top().first != queue2.top().first)
        {
            if(queue1.top().second > queue2.top().second)
            {
                merged.push(queue1.top());
                queue1.pop();
            }
            else
            {
                merged.push(queue2.top());
                queue2.pop();
            }
        }
        else
        {
            merged.push(std::make_pair(queue1.top().first, queue1.top().second + queue2.top().second));
            queue1.pop();
            queue2.pop();
        }
    }
    while(!queue1.empty())
    {
        merged.push(queue1.top());
        queue1.pop();
    }
    while (!queue2.empty()) {
        merged.push(queue2.top());
        queue2.pop();
    }
    while(merged.size() > maxSize)
        merged.pop();
}
