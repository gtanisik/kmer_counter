 /**
 * \file FastQFileStreamer.hpp
 * \author Gokhan Tanisik
 * \brief Defines the class to parse a FASQ file.
 */

#ifndef FASTQFILESTREAMER_HPP_
#define FASTQFILESTREAMER_HPP_

#include <fstream>
#include <string>

namespace kmers {
    
    /// Options for the \a FastQFileStreamer class.
    struct FastQFileStreamerOptions
    {
    	/// Constructor.
        FastQFileStreamerOptions()
        {}
        
        /// Maximum read size from file for each read operation. Default value is 65536.
        int maxReadSize = 65536;
        
        // Number of kmers in the file. If the number of kmers is not provided by the user, it will be counted one by one.
        uint64_t numKmers = 0;
    };
    
    /// Class to read kmers from a FASQ file.
    class FastQFileStreamer {
    public:
    
    	/// Defines a simple iterator-like type to get k-mers from the streamer one-by-one.
        struct KmerIterator
        {
        private:
        	/// Makes the streamer to read a new sequence from the file.
        	/// \return false if end of file is reached.
            bool read()
            {
                sequenceLength = streamer.read(sequence);
                currentIdx = 0;
                return (k <= sequenceLength);
            }
            
        public:
        	/// Constructor for the iterator.
            KmerIterator(FastQFileStreamer& streamer, char *sequence, unsigned int sequenceLength, unsigned int k)
            : streamer(streamer), sequence(sequence), sequenceLength(sequenceLength), k(k)
            {
            }
            
            /// Returns true if there are k-mers left to process.
            operator bool() {
                if((currentIdx + k) <= sequenceLength)
                    return true;
                return read();
            }
            
            /// Skips the current k-mer
            /// Returns true if a new k-mer to process is available.
            bool operator++()
            {
                bool retVal = false;
                if((currentIdx + k) < sequenceLength)
                {
                    currentIdx++;
                    retVal = true;
                }
                else retVal = read();
                return retVal;
            }
            
            /// \return a pointer to the current k-mer.
            char *operator*()
            {
                if(operator bool()) return &sequence[currentIdx];
                return 0;
            }
            
            /// FASTQ file streamer class.
            FastQFileStreamer &streamer;
            /// Current sequence pointer.
            char *sequence;
            /// Length of the current sequence.
            unsigned int sequenceLength;
            /// Current idx at the sequence.
            unsigned int currentIdx = 0;
            /// k of k-mer.
            unsigned int k;
        };
        
        /**
        	\brief Constructor for the FastQFileStreamer class.
        	
        	\param k K value of k-mer.
        	\param fileName FASTQ file to process.
        	\param options Options to the streamer.
        */
        FastQFileStreamer(int k, const std::string &&fileName, const FastQFileStreamerOptions &&options);
        
        /// Destructor
        virtual ~FastQFileStreamer();
        
        /**
        	\brief Reads a sequence from the file.
        	
        	\param buffer Returns the pointer to the sequence read.
        	\return Length of the sequence.
        */
        uint32_t read(char*& buffer);

        /// Used for test purposes only
        uint32_t dummyRead(char*& buffer);
        
        /// Returns the total kmer count in the file. If it's not counted before, reads the whole file to count.
        uint64_t kmerCount();
        
        /// Resets the streamer to beginning of the file and return an iterator.
        KmerIterator begin();
        
    private:
        /**
         Loads new data from file if there are less data than required length.
         
         @return false if the remaining data to be processed is less then minRequiredDataLength
         */
        bool loadDataIfNeeded(const int &minRequiredDataLength);
        
        /**
         skips a line.
         */
        void skipLine(const int &skipSizeHint);
        
        void reset();
    private:
        //! k in k-mers
        int mK;
        
        //! max read size for each read operation from file.
        int mMaxReadSize;
        
        //! size of the read buffer
        const int mBufferSize;
        
        //! input file reader.
        std::ifstream mIfs;
        
        //! name of the input file. must be a fastq formatted file!
        std::string mFileName;
        
        //! buffer read from file.
        char *mReadBuffer;
        
        //! current position of buffer to be proccessed
        int mBufPosition = 0;
        
        //! number of bytes in buffer that is meaningful
        int mBufPayload = 0;
        
        //! indicates whether a new line will be read or not
        bool mIsNewLine = true;
        
        //! Holds the size of the current line. Used to skip the quality line.
        int mCurrentLineLength = 0;
        
        //! Number of total k-mers in the file.
        uint64_t mNumKmers;
    };
    
} /* namespace kmers */

#endif /* FASTQFILESTREAMER_HPP_ */
