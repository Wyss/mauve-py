#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSourceFactory.h"
#include "libGenome/gnFASSource.h"
#include "libGenome/gnSequence.h"

#include <iostream>

#include "libGenome/gnFilter.h"
#include "libGenome/gnStringTools.h"
#include "libGenome/gnFastTranslator.h"

int main( int32 argc, char* argv[])
{
	argc; argv;
try{
	string filename;
	cout << "Enter a filename to read bases from.\n";
	cin >> filename;
	cout << "Opening " + filename + "\n";
	gnSequence gpseq, gnsubseq;
	
//	string pathname = filename;
//	standarizePathString( pathname );
//	pathname = getPathString( pathname );
//	gnSourceFactory * msf = gnSourceFactory::GetSourceFactory();
//	msf->AddPath( pathname );
	gpseq.LoadSource(filename);
	cout << "Length " << gpseq.length() << " in " << gpseq.contigListSize() << " contigs\n";

	string pretrans, posttrans, posttrans2;
	const gnTranslator*  dna2pro_good = gnTranslator::DNAProteinTranslator();
	const gnFastTranslator*  dna2pro_fast = gnFastTranslator::DNAProteinTranslator();
	pretrans = gpseq.ToString();

	posttrans = pretrans;
	posttrans2 = pretrans;
	cout << "Doing correct translation:\n";
	dna2pro_good->Filter(posttrans);
	cout << "done\n";
	gnsubseq = posttrans;
	gnFASSource::Write(gnsubseq, "profile_good.fas");
	
//	posttrans = pretrans;
	cout << "Doing fast translation:\n";
	dna2pro_fast->Filter(posttrans2);
	cout << "done\n";
	gnsubseq = posttrans2;
	gnFASSource::Write(gnsubseq, "profile_fast.fas");


	cout << "base pairs 8 thru 68 are: \n";
	cout << gpseq.subseq(8, 60) << "\n";

	for(uint32 i=0; i < 15; i++){
		try{
			cout << "Contig " << i << " length " << gpseq.contig(i).length() << "\n";
		}catch(gnException& gne){
			cout << gne;
		}
	}
	cout << gpseq.subseq(1000000, 10);
	cin >> filename;

	gnSeqI midpoint = gpseq.length() / 2;
	gnSequence seqA = gpseq.subseq(1, midpoint);
	gnSequence seqB = gpseq.subseq(1 + midpoint, gpseq.length() - midpoint);
	cout << "Splitting " << gpseq.length() << " to " << seqA.length() << " and " << seqB.length() << "\n";
	cin >> filename;

	gnsubseq = gpseq.subseq(5, 10);
	cout << "subseq len: " << gnsubseq.length() << "\n";
	cout << "subseq: " << gnsubseq << "\n";
	gnSeqC* buf = new gnSeqC[gpseq.length()];
	gpseq.ToArray(buf, gpseq.length());

	cout << "Give a file name to output data: ";
	string outfilename;
	cin >> outfilename;
	gnFASSource::Write(gpseq, outfilename);
	
/*	gnSequence new_seq = gpseq.contig(3);
	cout << new_seq;
	new_seq += gpseq.contig(3);
	new_seq += gpseq.contig(2);
	new_seq += gpseq.contig(1);
	cout << new_seq;
	gnFASSource::Write(new_seq, outfilename);
	
*/	cout << "All done.  Contigs 3, 2, 1 are in " << outfilename << "\n";
	char bubba[50];
	cin >> bubba;
}catch(gnException& gne){
	cout << gne;
	string reality_bites;
	cin >> reality_bites;
}
}
