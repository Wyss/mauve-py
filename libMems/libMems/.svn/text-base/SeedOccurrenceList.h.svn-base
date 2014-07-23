#ifndef __SeedOccurrenceList_h__
#define __SeedOccurrenceList_h__

#include <vector>
#include "libMems/SortedMerList.h"
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include "libMems/Files.h"

namespace mems
{

class SeedOccurrenceList
{
public:
	typedef float32 frequency_type;

	SeedOccurrenceList(){}

	template< typename SMLType >
	void construct( SMLType& sml )
	{
//		std::vector<mems::bmer> mer_vec;
//		sml.Read( mer_vec, sml.SMLLength(), 0 );
		const size_t total_len = sml.Length();
		frequency_type* count = new frequency_type[total_len];
		size_t seed_start = 0;
		size_t cur_seed_count = 1;
		uint64 mer_mask = sml.GetSeedMask();
		size_t seedI = 1;
		bmer prevmer;
		bmer merI; 
		const size_t sml_length = sml.SMLLength();
		if( sml_length > 0 )
			merI = sml[0];
		for( seedI = 1; seedI < sml_length; seedI++ )
		{
			prevmer = merI;
			merI = sml[seedI];
			if( (merI.mer & mer_mask) == (prevmer.mer & mer_mask) )
			{
				++cur_seed_count;
				continue;
			}
			// set seed frequencies
			for( size_t i = seed_start; i < seedI; ++i )
				count[sml[i].position] = (frequency_type)cur_seed_count;
			seed_start = seedI;
			cur_seed_count = 1;
		}
		// set seed frequencies for the last few
		for( size_t i = seed_start; i < seedI && i < sml_length; ++i )
			count[sml[i].position] = (frequency_type)cur_seed_count;
		// hack: fudge the last few values on the end of the sequence, necessary when sequence isn't circular
		for( ; seedI < total_len; ++seedI )
			count[seedI]=1;

		smoothFrequencies( sml, count );

		// wipe out any stray zeros
		for( size_t i = 0; i < total_len; ++i )
			if( count[i]== 0 )
				count[i] = 1;

		// create a temporary memory-mapped file to store mer counts
		tmpfile = CreateTempFileName("sol");
		// resize the file to be big
		std::ofstream tfout;
		tfout.open(tmpfile.c_str(), std::ios::binary );
		tfout.write(((const char*)count), sml.Length()*sizeof(frequency_type));
		tfout.close();
		delete[] count;

		data.open( tmpfile );	// map the file
	}


	frequency_type getFrequency( gnSeqI position )
	{
		return ((frequency_type*)data.data())[position];
	}
	~SeedOccurrenceList()
	{
		if(data.is_open())
		{
			data.close();
			boost::filesystem::remove(tmpfile);
			tmpfile.clear();
		}
	}

//	SeedOccurrenceList( const SeedOccurrenceList& sol ){if(data.is_open()) throw "not copyable";};
//	SeedOccurrenceList& operator=( const SeedOccurrenceList& sol ){if(data.is_open()) throw "not copyable";}

protected:
	/**
	 * converts position freqs to the average freq of all k-mers containing that position
	 */
	template< typename SMLType >
	void smoothFrequencies( const SMLType& sml, frequency_type* count )
	{
		size_t seed_length = sml.SeedLength();
		// hack: for beginning (seed_length) positions assume that previous
		// containing seeds were unique
		double sum = seed_length - 1 + count[0];
		std::vector<frequency_type> buf(seed_length, 1);
		buf[0] = count[0];
		for( size_t i = 1; i < sml.Length(); i++ )
		{
			count[i-1] = sum / seed_length;
			sum += count[i];
			size_t bufI = i % seed_length;
			sum -= buf[bufI];
			buf[bufI] = count[i];
		}
	}
	
	boost::iostreams::mapped_file_source data;
	std::string tmpfile;

};

}

#endif	// __SeedOccurrenceList_h__

