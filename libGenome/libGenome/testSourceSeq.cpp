#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include "libGenome/gnSequence.h"
#include <fstream>

#include "libGenome/gnFilter.h"
#include "libGenome/gnFASSource.h"

int32 main( int32 argc, char* argv[])
{
	argc; argv;
	gnSequence gps;
	
	string filename;
	cout << "Enter a filename to compute a base composition.\n";
	cin >> filename;
	cout << "Opening " + filename + "\n";
	gps.LoadSource(filename);
	cout << "Genome Length: " << gps.length() << " in " << gps.contigListLength() << " contigs\n";
	
	cout << "Give a file name to output composition data: ";
	string outfilename;
	cin >> outfilename;
	ofstream outfile(outfilename.c_str(), ios::out);
	outfile << "Contig_Name	A_count	C_count	G_count	T_count\n";

	ofstream ambig("ambigs.txt", ios::out);
	
	gnSeqI a_total = 0, c_total = 0, g_total = 0, t_total = 0, ambiguities = 0;
	for(int contigI=0; contigI < gps.contigListLength(); contigI++){
		gnSeqI seqLength = gps.contig(contigI).length();
		gnSeqC *bases = new gnSeqC[seqLength];
		gps.contig(contigI).ToArray(bases, seqLength);
		string contigName = gps.contigName(contigI);

		gnSeqI a = 0, c = 0, g = 0, t = 0;
		for(int i=0; i < seqLength; i++){
			if((bases[i] == 'a')||(bases[i] == 'A'))
				a++;
			else if((bases[i] == 'c')||(bases[i] == 'C'))
				c++;
			else if((bases[i] == 'g')||(bases[i] == 'G'))
				g++;
			else if((bases[i] == 't')||(bases[i] == 'T'))
				t++;
			else{
				ambiguities++;
				ambig << bases[i];
			}
		}
		delete[] bases;
		outfile << contigName << "	" << a << "	" << c << "	" << g << "	" << t << "\n";
		a_total += a; c_total += c; g_total += g; t_total += t;
	}
	outfile.close();
	ambig.close();
	cout << "Enter a file for me to write to: ";
	cin >> outfilename;
	gnFASSource::Write(gps, outfilename);
	cout << "A: " << a_total << "\nC: " << c_total << "\nG: " << g_total << "\nT: " << t_total << "\n";
	cout << "Ambiguities: " << ambiguities << "\n";
	cout << "For a total of: " << a_total + c_total + g_total + t_total + ambiguities << "\n";
	cout << "Base comp. completed.  Results are in " << outfilename << "\n";
	char bubba[50];
	cin >> bubba;
}