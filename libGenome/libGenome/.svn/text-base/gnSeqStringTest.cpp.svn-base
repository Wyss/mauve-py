#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "gnSequence.h"
#include "gnFASSource.h"
#include <string>
#include <iostream>

const uint32 SURROUND_BASES = 100;

int main(int argc, char* argv[]){

	gnSequence sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
	sequence.append("GGGGGGGGGGG");
	
	cout << sequence;
	char* dumbo;
	cin >> dumbo;


/*	gnSequence sequence, output_seq;
	string seq_filename;
	string restriction_filename;
	string output_filename;
	fstream restriction_file;
	fstream output_file;
	
	do{
		cout << "Enter the name of the restriction site list:\n";
		cin >> restriction_filename;
		restriction_file.open(restriction_filename.c_str(), ios::in);
		if(!restriction_file.is_open())
			cout << "Error opening file!  Try again.\n";
	}while(!restriction_file.is_open());

	do{
		cout << "Enter the name of the sequence data FastA file:\n";
		cin >> seq_filename;
	}while(!sequence.LoadSource(seq_filename));
	
	do{
		cout << "Enter the name of an output file:\n";
		cin >> output_filename;
		output_file.open(output_filename.c_str(), ios::out);
		if(!output_file.is_open())
			cout << "Error creating file!  Try again.\n";
	}while(!output_file.is_open());
	output_file.close();
	
	
	while(restriction_file.good()){
		string bubba;
		restriction_file >> bubba;	//first four are junk
		restriction_file >> bubba;
		restriction_file >> bubba;
		restriction_file >> bubba;
		restriction_file >> bubba;	//this contains contig name
		
		gnSequence contig_seq = sequence.contigByName(bubba);

		gnSeqI siteI = 0;
		restriction_file >> siteI;

		char dummy[1000];
		restriction_file.getline(dummy, 999);
		
		//if the restriction site is too close to the beginning, set the start to 0
		gnSeqI start_loc = siteI < SURROUND_BASES ? 0 : siteI - SURROUND_BASES;
		output_seq += contig_seq.subseq(start_loc, SURROUND_BASES*2);
	}

	if(gnFASSource::Write(output_seq, output_filename))
		cout << "All done, results are in " << output_filename << "\n";
	else
		cout << "An error occurred while writing " << output_filename << "\n";
	*/
}