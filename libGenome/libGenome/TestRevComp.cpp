#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSourceFactory.h"
#include "libGenome/gnFASSource.h"
#include "libGenome/gnDNXSource.h"
#include "libGenome/gnSEQSource.h"
#include "libGenome/gnSequence.h"
#include <iostream>
#include <fstream>

#include "libGenome/gnFilter.h"

int main( int32 argc, char* argv[])
{

	argc; argv;

	string filename;
	cout << "Enter a filename to read bases from.\n";
	cin >> filename;
	cout << "Opening " + filename + "\n";
	gnSequence gnseq, smallseq;
	if(gnseq.LoadSource(filename))
		cout << "Sequence has " << gnseq.length() << " base pairs.\n";
	
	smallseq = gnseq.subseq(3836480, 10);
	cout << smallseq;
	string dump;
	cin >> dump;

	cout << "Give a file name to output reverse complement data: ";
	string outfilename;
	cin >> outfilename;
	cout << "Bases are:\n";
//	cout << gpseq;
	cout << "\nComplement Bases are:\n";
	gnBaseSpec* gpbs = gnseq.GetSpec();
	gpbs->SetReverseComplement(true);
//	cout << gpseq << "\n";
	
	gnFASSource::Write(gnseq, outfilename);
//	gnGBKSource::Write(gpbs, "testfile.seq");
//	gnDNXSource::Write(gpbs, "testfile.dnx");
	
	cout << "All done.  RevComp is in " << outfilename << "\n";
	char bubba[50];
	cin >> bubba;
}
