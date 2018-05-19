/**
 * \file BloomFilter.hpp
 * \author Gokhan Tanisik
 * \brief Bloom Filter implementation.
 *
 * Defines a Bloom Filter and structures to be used by the filter.
 */

#ifndef BloomFilter_hpp
#define BloomFilter_hpp

#include <vector>
#include <cmath>
#include <climits>
#include <iostream>
#include <limits>
#include "external/MurmurHash3.h"

/// Defines the namespace for kmer counting problem.
namespace kmers
{
	/// Defines the options to configure the bloom filter.
    struct BloomFilterOptions
    {
    	/// Expected false positive rate for the filter.
        float expectedFPRate = 0.2f;
        /// Max allowed memory size for the filter.
        unsigned int maxMemorySizeMB = 8192;
        /// Number of hash functions to be used in the filter. Auto-calculated if it's -1.
        int numHashFunctions = -1; 
        
        /// Number of bits for each entry. Allowed values are [1, 8, 16, 32, 64].
        /// \note Used for single layer version only.
        int bitsPerEntry = 8; 
        
        /// Number of filter layers to be used (Cascaded Bloom Filters are used).
        /// \note It's not used in the BloomFilter class directly.
        int layerCount = 1; 
        /// When using cascaded bf counter, the bit size is decimated for each layer to save memory.
        //  this leads to an increase in the amount collisions if it's not 1.0.
        /// \note It's not used in the BloomFilter class directly.
        float cascadedKmerCountDecimation = 0.25;
    };
    
    /**
    	\brief Interface for the Bloom Filter class.
    	
    	This class is used for multi-layer solution. It enables to hold all of the filters in the same container.
    */
    class IBloomFilter
    {
    public:
    
    	/// Destructor.
        virtual ~IBloomFilter()
        {
        }
        
        /**
         \brief Adds a new kmer into the filter.
         
         Returns true if the operation succeeds, false otherwise.
         If 1 bit is used for each kmer, failure means a collision or this kmer is already added to the filter.
         Else, failure means all of (hash) slots corresponding to this kmer has already reached it's max value.
         
         \param seq The raw bytes sequence of entry to be added into the filter.
         \pafam length Length of entry to be added.
         */
        virtual bool add(const void *seq, unsigned int length) = 0;
        
        /**
         \brief Checks if the given entry is (possibley) added into the filter.
         
         \param seq The raw bytes sequence of entry to be added into the filter.
         \pafam length Length of entry to be added.
        */
        virtual bool contains(const void *seq, unsigned int length) = 0;
        
        /**
         \brief Returns the count of given entry in the filter.
         
         \param seq The raw bytes sequence of entry to be added into the filter.
         \pafam length Length of entry to be added.
        */
        virtual uint64_t count(const void *seq, unsigned int length) = 0;
        
        /**
         \brief Returns the maximum counter value for an entry.
        */
        virtual uint64_t maxValuePerEntry() const = 0;
        
        /**
         \brief Returns the bit count for each entry.
        */
        virtual int bitsPerEntry() const = 0;
    };
    
    /**
    	\brief Implements the Bloom Filter.
    	
    	The template <em>class T </em> defines the entry type. Should be bool or an unsigned integral type.
    */
    template <class T>
    class BloomFilter : public IBloomFilter
    {
    public:
        BloomFilter(unsigned int k, unsigned long long numKmers, const BloomFilterOptions &&options = BloomFilterOptions())
        : mNumHashes(options.numHashFunctions), mK(k), mNumKmers(numKmers), mFPRate(options.expectedFPRate)
        {
            allocateBits(options.maxMemorySizeMB);
        }
        
        virtual ~BloomFilter()
        {
            
        }
        
        bool add(const void *seq, unsigned int length)
        {
            std::vector<typename std::vector<T>::size_type> idxes(mNumHashes);
            
            if(count(seq, length, idxes) == maxValuePerEntry())
            {// overflow detected
                return false;
            }
            
            for(auto idx : idxes)
            {
                mBits[idx] = saturatedAdd(mBits[idx], 1);
            }
            
            return true;
        }
        
        bool contains(const void *seq, unsigned int length)
        {
            return !!count(seq, length);
        }
        
        uint64_t count(const void *seq, unsigned int length)
        {
            std::vector<typename std::vector<T>::size_type> idxes(mNumHashes);
            
            return count(seq, length, idxes);
        }
        
        uint64_t maxValuePerEntry() const {
            return static_cast<uint64_t>(std::numeric_limits<T>::max());
        }
        
        inline int bitsPerEntry() const
        {
            return 8*sizeof(T);
        }
        
    private:
        constexpr T maxValuePerEntryFast() const
        {// to avoid the call to the maxValuePerEntry() function in the interface. it causes upcast to uint64
            return std::numeric_limits<T>::max();
        }
        
        inline T saturatedAdd(const T& a, const T& b) const
        {
            T c = a+b;
            if(c < a)
                c = maxValuePerEntryFast();
            return c;
        }
        
        T count(const void *seq, unsigned int length, std::vector<typename std::vector<T>::size_type> &idxes)
        {
            hash(seq, length, idxes);
            
            T count = maxValuePerEntryFast();
            for(auto idx : idxes)
            {
                if(!mBits[idx])
                    return 0;
                if(count > mBits[idx])
                    count = mBits[idx];
            }
            
            return count;
        }
        
        void allocateBits(unsigned int maxMemorySizeMB)
        {
            // allocate memory for bits
            while (true) {
                calcParams();
                
                float memoryRequiredMB = (mNumBits * bitsPerEntry()) / (8.0f * 1024.0f * 1024.0f);
                
                if(static_cast<unsigned int>(memoryRequiredMB) > maxMemorySizeMB)
                {
                    std::cout << "With given parameters (#k-mers:" << mNumKmers << " and expected FP rate: "
                    << mFPRate <<") " << mNumBits << " bits (" << memoryRequiredMB << " MB) required."
                    << " Yet, you wanted me to use maximum " << maxMemorySizeMB << " MB of memory. FP rate will be increased by 0.01 and I will try again..." << std::endl << std::endl;
                }
                else
                {
                    mBits.resize(mNumBits, 0);
                    if(mBits.size() < mNumBits)
                    {
                        std::cout << "With given parameters (#k-mers:" << mNumKmers << " and expected FP rate: "
                        << mFPRate <<") " << mNumBits << " bits (" << memoryRequiredMB << " MB) required."
                        << " Unfortunately this computer can not handle that much memory. FP rate will be increased by 0.01 and I will try again..." << std::endl << std::endl;
                    }
                    else
                    {
                        std::cout << "The Bloom Filter is created with " << mNumBits << " bits (" << memoryRequiredMB
                        << " MB), " << bitsPerEntry() << " bits per entry, " << mNumHashes << " hash functions, and " << mFPRate << " expected false positive rate.." << std::endl << std::endl << std::endl;
                        break;
                    }
                }
                mFPRate += 0.01;
            }
        }
        
        void calcParams()
        {
            mNumBits = std::ceil( (static_cast<double>(mNumKmers) * std::log(mFPRate)) /
                                 log(1.0 / std::pow(2.0, log(2.0))) );
            //mNumBits = (mNumBits + 7) - ((mNumBits + 7) % 8);
            
            if(mNumHashes <= 0 )
            {// it is not enforced by the user?
                mNumHashes = round(log(2.0) * (mNumBits / mNumKmers));
            }
        }
        
        void hash(const void *seq, unsigned int length, std::vector<typename std::vector<T>::size_type> &idxes)
        {
            unsigned int seed = 761914759;
            
            uint64_t h[2];
            MurmurHash3_x86_128(seq, length, seed, h);
            
            for(uint64_t i = 0; i < mNumHashes; i++)
            {
                idxes[i] = (h[0] + i*h[1]) % mNumBits;
            }
        }
        
    private:
        int mNumHashes;
        
        uint64_t mNumBits; // number of bits
        std::vector<T> mBits;
        
        unsigned int mK;
        
        uint64_t mNumKmers;
        float mFPRate;
    };
    
    template <>
    inline int BloomFilter<bool>::bitsPerEntry() const// specialize
    {
        return 1;
    }

    template<>
    inline bool BloomFilter<bool>::saturatedAdd(const bool& a, const bool& b) const
    {
        return a||b;
    }
}

#endif /* BloomFilter_hpp */
