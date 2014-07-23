#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libGenome/gnFASSource.h"
#include "libGenome/gnGBKSource.h"
#include <iostream>
#include <string>
#include "libGenome/gnTranslator.h"
#include "libGenome/gnDebug.h"
#include "libGenome/gnDNASequence.h"
#include "libGenome/gnDebug.h"

void runTests(gnSequence& seq);
void runTests(gnSequence& seq){
	cout << "Testing with sequence:\n";
	cout << seq << endl;

	seq.setReverseComplement(true);
	cout << "Reverse Complement:\n";
	cout << seq << endl;

	seq.setReverseComplement(false);
	cout << "Forward sequence:\n";
	cout << seq << endl;

	seq.setFilter(gnTranslator::DNAProteinTranslator);
	cout << "Translated to protein sequence:\n";
	cout << seq << endl;

	seq.setReverseComplement(true);
	cout << "Reverse complement translated to protein sequence:\n";
	cout << seq << endl;

	list<const gnBaseFilter*> filter_list;
	filter_list.push_back(gnTranslator::DNAProteinTranslator);
	filter_list.push_back(gnTranslator::ProteinDNATranslator);
	seq.setFilterList(filter_list);
	cout << "Reverse Comp translated to protein sequence and back to dna:\n";
	cout << seq << endl;
};

int main(int argc, char* argv[]){

	gnDNASequence seq = "ACGTTAGCGTAT";
	gnDNASequence seq2 = "A";
	seq2 += "CG";
	seq2 += "TTA";
	seq2 += "GCGT";
	seq2 += "AT";
	string blah;
	DebugMsg("Hello World.\n");
	runTests(seq);
	cout << "And now for the same thing with the same sequence.\n";
	cin >> blah;
	runTests(seq2);
	cin >> blah;
	
	do{
		cout << "Finally, enter a sequence file to play with: ";
		cin >> blah;
	}while(!seq.LoadSource(blah));
	seq.setFilter(NULL);
	gnFASSource::Write(seq, "test_rc.fas");
	gnGBKSource::Write(seq, "test_rc.gbk");
	
	return 0;
}
