CXX=clang++
INC=-I/usr/local/Cellar/google-sparsehash/2.0.3/include
#LIB=-L/usr/local/Cellar/google-sparsehash/2.0.3/lib/pkgconfig
#CXXFLAGS=-g -std=c++11  -Wall -pedantic
CXXFLAGS=-std=c++11  -Wall -O3

kmer_counter: kmers.o kmer_options.o FastQFileStreamer.o MurmurHash3.o
	$(CXX) -o kmer_counter kmers.o kmer_options.o FastQFileStreamer.o MurmurHash3.o

kmers.o : kmers.cpp BloomFilter.hpp
	$(CXX) $(CXXFLAGS) $(LIB) $(INC) -c kmers.cpp

kmer_options.o : kmer_options.cpp
	$(CXX) $(CXXFLAGS) -c kmer_options.cpp

MurmurHash3.o : external/MurmurHash3.cpp
	$(CXX) $(CXXFLAGS) -c external/MurmurHash3.cpp
	
FastQFileStreamer.o : FastQFileStreamer.cpp
	$(CXX) $(CXXFLAGS) -c FastQFileStreamer.cpp

clean : 
	rm *.o
	
.PHONY: clean