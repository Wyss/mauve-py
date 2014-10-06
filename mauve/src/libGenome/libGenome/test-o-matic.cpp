#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnFASSource.h"
#include "libGenome/gnSequence.h"
#include "libGenome/gnTranslator.h"
#include "libGenome/gnStringSpec.h"
#include <fstream>
int main( int argc, char* argv[] )
{
	gnSequence seq_a, seq_b;
	string seq_a_name="styphi.fas";
	string seq_b_name="typhi.gbk";
	try{
		seq_a.LoadSource( seq_a_name );
		cout << seq_a_name << '\t' << seq_a.length() << " bp\n"; 
		seq_b.LoadSource( seq_b_name );
		cout << seq_b_name << '\t' << seq_b.length() << " bp\n"; 
	}catch(gnException& gne){
		cout << gne << "\n";
	}
	
	uint len = 50;
	uint increment = 1027;
	const gnCompare* dna_comp = gnCompare::DNASeqCompare();
	for( uint seqI = 1; seqI <= 200000; seqI+= increment ){
		string str_a, str_b;
		str_a = seq_a.ToString( len, seqI );
		str_b = seq_b.ToString( len, seqI );
		if( !dna_comp->Contains( str_a, str_b ) ){
			cout << str_a << endl << str_b << endl;
			seqI -= increment;
		}
	}

	gnSequence seq_in;
	string seq_name = "sequin.gbk";
	string out_seq = "seqout.gbk";
	try{
		seq_in.LoadSource( seq_name );
		cout << seq_name << '\t' << seq_in.length() << " bp\n"; 
		gnFASSource::Write( seq_in, out_seq );
	}catch(gnException& gne){
		cout << gne << "\n";
	}
/*	
	gnSequence erwinia_junk;
	gnSequence pruned_bush;
	string sequence_name = "mouse_small.fas";
	string output_sequence = "mouse_filt.fas";
//	cout << "Enter the name of the sequence file you want contig info on.\n";
//	cin >> sequence_name;
	try{
		erwinia_junk.LoadSource( sequence_name );
		for(uint32 contigI = 0; contigI < erwinia_junk.contigListSize(); contigI++){
			cout << erwinia_junk.contigName( contigI ) << '\t' << erwinia_junk.contigLength( contigI ) << '\n';
			if( erwinia_junk.contigLength( contigI ) > 1000 )
				pruned_bush += erwinia_junk.contig( contigI );
		}
		gnFASSource::Write(pruned_bush, output_sequence);
	}catch(gnException& gne){
		cout << gne << "\n";
	}
*/
/*
	ifstream gbk1("cpneuJ.gbk");
	ifstream gbk2("test_rc.gbk");

	uint32 badline = 0;
	char buf[30000];
	char buf2[30000];
	while(gbk1.good() && gbk2.good()){
		gbk1.getline(buf, 30000);
		gbk2.getline(buf2, 30000);
		if(strcmp(buf, buf2)){
			cout << "mismatched line: " << badline << "\n";
			cout << buf << endl << buf2 << endl;
			cout << "Length 1: " << strlen(buf) << "   Length 2: " << strlen(buf2) << endl;
			cin >> buf;
		}
		badline++;
	}
*/
/*	gnStringSpec s("");
	cout <<s.GetLength()<<endl;
	
	gnSequence t = "ACTATATA";
	if(t == string("TATATA")){
		cout << "ASDFASDFSADF";
	}*/
}
/*	// define a string to store the sequence file name
	
	string filename;
	string nseq,pseq;	
	// define a gnSequence to store the sequence
	gnSequence file_sequence;
	gnSequence prot_trans;
	
	//nseq = "BXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXXBXX";
	//nseq = "WSTYYYBBB";
    //nseq="ATGCUXGCNGGGTATGAATEQILPFP";
	//nseq="ATGGTXAGCGGNTAR";	
	cout<<nseq<<" (original)---->  \n";
	gnSequence aSeq(nseq);
	cout<<aSeq<<" (gnSequence)----> \n";

    gnTranslator::ProteinDNAConverter->ToProtein(nseq);
    
    cout << nseq << '\n';;
  
  
  
  
  //  prot_trans.append(nseq);
	// output the first ten bases of the sequence
//	cout << prot_trans.subseq(1,10);

	return 0;
	
}

/*#include "libGenome/gnSequence.h"
#include <iostream>
#include <string>
#include "libGenome/gnTranslator.h"
#include "libGenome/gnFASSource.h"

int main(int argc, char* argv[]){
	gnSequence gps;
	string seq_filename;
	string out_filename;
	
	/* Get the name of the dnx file to convert to FastA /
	cout << "Debug-o-matic!\n";
	cout << "Please give the name of the fof to test: ";
	cin >> seq_filename;
	
	gps.LoadSource(seq_filename);

	cout << "Output filename: ";
	cin >> out_filename;
	gnFASSource::Write(gps, out_filename);
	
	return 0;	

}*/