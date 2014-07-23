#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnFASSource.h"
#include <iostream>


#include "libGenome/gnFilter.h"

int32 main( int32 argc, char* argv[])
{
	argc; argv;

	gnFASSource source;
	string sourceFile = "w:/developm/_data/cft073/v6/cft073v6.fas";
	cout << "Enter Source: ";
	cin >> sourceFile;
	cout << "Opening " << sourceFile << endl;
	if( !source.Open( sourceFile ) )
	{
		cout << "Failed to open " << sourceFile << endl;
		return 1;
	}

	cout << "Number of Contigs " << source.GetContigListLength() << endl;
	for( uint32 i=0; i < source.GetContigListLength(); ++i )
	{
		gnFileContig* fg = source.GetContig(i);
		cout << "\t" << fg->GetName() <<  endl;
		cout << "\t\t" << fg->GetSeqLength() << " : {" 
			<< fg->GetFileStartEnd().first << ", "
			<< fg->GetFileStartEnd().second << ": ";
		for( uint32 j=0; j<CONTIG_SECTION_SIZE; ++j )
			cout << "(" << fg->GetSectStartEnd( (gnContigSection)j ).first << ","
				<< fg->GetSectStartEnd( (gnContigSection)j ).second << ") ";
		cout << "}" << endl;
		cout << "\t\t" << (fg->HasRepeatSeqGap()?"TRUE":"FALSE")
			<< " (" << fg->GetRepeatSeqGapSize().first << ","
			<< fg->GetRepeatSeqGapSize().second << ")" << endl;
	}

	uint32 seqLen = 75;
	gnSeqC* seq = new gnSeqC[seqLen+1];
	for( uint32 i=0; i < seqLen; ++i )
		seq[i] = '@';
	seq[seqLen] = 0;

	cout << "SeqRead(0, seq, " << seqLen << ", 0)" << endl;
	if( source.SeqRead( 158, seq, seqLen, 0 ) )
	{
		seq[seqLen] = 0;
		cout << seqLen << "> " << seq << endl;
	}
	else
		cout << "Read failed " << endl;
	return 0;
}