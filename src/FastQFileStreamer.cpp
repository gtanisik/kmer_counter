 /**
 * \file FastQFileStreamer.cpp
 * \author Gokhan Tanisik
 * \brief Defines the class to parse a FASQ file.
 */

#include "FastQFileStreamer.hpp"

#include <iostream>
#include <exception>

inline int MIN(const int &x, const int &y)
{
    return x < y ? x : y;
}

using namespace kmers;

FastQFileStreamer::FastQFileStreamer(int k, const std::string &&fileName, const FastQFileStreamerOptions &&options) :
mK(k), mMaxReadSize(options.maxReadSize), mBufferSize(2*(options.maxReadSize + k)), mFileName(fileName), mNumKmers(options.numKmers)
{
    mReadBuffer = new char[mBufferSize];
    if(!mReadBuffer)
        throw std::runtime_error("not enough memory :(");
    
    mIfs.open(mFileName.c_str(), std::ios::in | std::ios::binary);
    if(!mIfs.is_open())
        throw std::runtime_error("input file \"" + mFileName + "\" can not be opened!");
    
    reset();
    kmerCount();
}

FastQFileStreamer::~FastQFileStreamer()
{
    delete[] mReadBuffer;
}

FastQFileStreamer::KmerIterator FastQFileStreamer::begin()
{
    reset();
    char *buffer;
    unsigned int length;
    length = read(buffer);
    return KmerIterator(*this, buffer, length, mK);
}

void FastQFileStreamer::reset()
{
    mIfs.clear();
    mIfs.seekg(0);
    mBufPayload = 0;
    mBufPosition = 0;
    mCurrentLineLength = 0;
    mIsNewLine = true;
    memset(mReadBuffer, 0, mBufferSize);
}

uint64_t FastQFileStreamer::kmerCount()
{
    // if it's already counted just return the count
    if(mNumKmers)
        return mNumKmers;
    
    uint64_t counter = 0;
    char *seq;
    uint32_t seqLen;
    
    while((seqLen = read(seq)))
    {
        counter += seqLen - mK + 1;
    }
    
    reset();
    
    mNumKmers = counter;
    return mNumKmers;
}

uint32_t FastQFileStreamer::dummyRead(char *&buffer)
{
    mIfs.read(mReadBuffer , mMaxReadSize);
    return static_cast<unsigned int>(mIfs.gcount());
}

uint32_t FastQFileStreamer::read(char *&buffer)
{
    if(!loadDataIfNeeded(mK))
    {// no available data left to process
        return 0;
    }
    
    if(mIsNewLine)
    {// skip "+", quality sequence and sequence identifier lines.
        
        if(mCurrentLineLength > 0) // it's 0 only at the beginning
        { // don't run for the first line of the file
            skipLine(mCurrentLineLength); // plus line
            skipLine(mCurrentLineLength); // quality sequence line
        }
        
        skipLine(mK); // skip the sequence identifier
        
        mIsNewLine = false;
        
        // mK-1 because of overlapping size between streams.
        mCurrentLineLength = mK-1;
        
        if(!loadDataIfNeeded(mK))
        {// no available data left to process
            return 0;
        }
    }
    
    // the buffer will start from the current position of the mReadBuf
    buffer = mReadBuffer + mBufPosition;
    unsigned int length = 0;
    
    {// determine the length of sequence to be processed
        char *newLinePtr = static_cast<char*>(memchr(buffer, '\n', mBufPayload));
        if(!newLinePtr)
        { // no new line encountered. all of the data is valid
            length = mBufPayload;
            mBufPayload = mK - 1;
            mBufPosition = (mBufPosition + length) - mBufPayload;
            mCurrentLineLength += length - mK + 1;
        }
        else
        {// data is valid to the new line
            mIsNewLine = true; // plus and quality sequence lines will be skipped on next call
            length = static_cast<int>(newLinePtr-buffer);
            mBufPosition += length + 1; // plus new line
            mBufPayload -= length + 1;
            mCurrentLineLength += length - mK + 1;
        }
    }
    
    return length;
}

bool FastQFileStreamer::loadDataIfNeeded(const int &minRequiredDataLength)
{
    // read new data if necessary. return if end of file reached.
    if(mBufPayload < minRequiredDataLength)
    {// need to read new data
        // the remaining payload may still be valid. shift it to the beginning of the buffer.
        memmove(mReadBuffer, mReadBuffer+mBufPosition, mBufPayload);
        // restart from the beginning
        mBufPosition = 0;
        mIfs.read(mReadBuffer + mBufPayload, mMaxReadSize);
        // update the payload with the size of data read.
        mBufPayload += static_cast<unsigned int>(mIfs.gcount());
        
        if(mBufPayload < minRequiredDataLength)
        {// end of file reached..
            return false;
        }
    }
    return true;
}

void FastQFileStreamer::skipLine(const int &skipSizeHint)
{
    do
    {
        if(!loadDataIfNeeded(MIN(mMaxReadSize+mK-mBufPayload, skipSizeHint)))
        {// no available data left to process
            mBufPayload = 0;
            break;
        }
        
        char *newLinePtr = static_cast<char*>(memchr(mReadBuffer+mBufPosition, '\n', mBufPayload));
        if(newLinePtr)
        { // new line found
            int length = static_cast<int>(newLinePtr-mReadBuffer) - mBufPosition;
            mBufPosition += length + 1; // plus new line
            mBufPayload -= length + 1;
            
            break;
        }
        else
        {// no new line encountered, invalidate current data
            mBufPayload = 0;
        }
    }while(true);
}
