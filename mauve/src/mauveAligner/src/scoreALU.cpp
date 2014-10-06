/*******************************************************************************
 * $Id: scoreAlignment.cpp,v 1.14 2004/02/28 00:01:31 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MatchList.h"
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include "libMems/IntervalList.h"
#include "libGenome/gnFilter.h"
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;
using namespace genome;
using namespace mems;

class AluRecord
{

public:

	AluRecord();
	//Smith-Waterman score of the match
	int score;
	//% substitutions in matching region compared to the consensus
	float divergent;
	//% deleted bp
	float deleted;
	//%	inserted bp
	float inserted;
	//name of query sequence
	string queryId;
	//starting position of match in query sequence
	gnSeqI start;
	//ending position of match in query sequence
	gnSeqI end;
	//no. of bases in query sequence past the ending position of match
	gnSeqI remaining;
	//match is with the Complement of the consensus sequence in the database
	string strand;
	//name of the matching interspersed repeat
	string repeatId;
	//class of the matching repeat
	string repeatClass;
	//no. of bases in (complement of) the repeat consensus sequence 
    //prior to beginning of the match (so 0 means that the match extended 
    //all the way to the end of the repeat consensus sequence)
	gnSeqI prior;
	//starting position of match in database sequence
	gnSeqI startDB;
	//ending position of match in database sequence
	gnSeqI endDB;
	gnSeqI length(void);
};

gnSeqI AluRecord::length(void)
{
	gnSeqI len = 0;
	len = absolut((int64)end)-absolut((int64)start);
	return len;
}
AluRecord::AluRecord()
{
	score = 0;
	divergent = 0.0;
	deleted = 0.0;
	inserted = 0.0;
	queryId = "none";
	start = 0;
	end = 0;
	remaining = 0;
	strand = "+";
	repeatId = "none";
	repeatClass = "none";
	prior = 0;
	startDB = 0;
	endDB = 0;
}
void ReadAluFile( istream& in_stream, vector<AluRecord*>& alu_list, gnSeqI& lr ) 
{
	uint seq_count = 0;
	gnSeqI max_len = 0;
	string cur_line;
	//3 lines of header info
	getline( in_stream, cur_line );
	getline( in_stream, cur_line);
	getline( in_stream, cur_line);
	uint seqI = 0;
	vector< gnSeqI > lengths;
	//vector< AluRecord* > alu_list;
	
	string empty_line;
	vector< string > aln_mat;
	uint line_count = 1;

	
	while( getline( in_stream, cur_line) )
	{
		
		AluRecord* alu = new AluRecord();
		// read and parse first AluRecord line
		//stringstream line_str( cur_line );
		
		//first line of data
		
		// take off leading whitespace
		string::size_type loc = cur_line.find("(");
		if (loc != string::npos )
			cur_line.replace(loc,1," ");
		
		loc = cur_line.find(")");
		if (loc != string::npos )
			cur_line.replace(loc,1," ");
		stringstream parse_str( cur_line );
		
		parse_str >> alu->score;	
		parse_str >> alu->divergent;
		parse_str >> alu->deleted;
		parse_str >> alu->inserted;
		parse_str >> alu->queryId;
		parse_str >> alu->start;
		parse_str >> alu->end;
		parse_str >> alu->remaining;
		parse_str >> alu->strand;
		parse_str >> alu->repeatId;
		parse_str >> alu->repeatClass;
		//punt: rest of info not needed
		//parse_str >> alu->prior;
		//parse_str >> alu->startDB;
		//parse_str >> alu->endDB;
		
		//end of line
		alu_list.push_back(alu);
		lr+= alu->length();

	}
	cout << "number of ALU records in file: " << alu_list.size() << endl;
}

/**
 * program to score alignments
 * reads in a "correct" alignment and a calculated alignment
 * scores the calculated alignment based on the correct one
 */
int main( int argc, char* argv[] ){
	
	string alignment_fname;
	string alu_fname;
	
	
	if( argc < 2 ){
		cout << "scoreALU <procrastAligner alignment> <repeatmasker ALUs>\n";
		return -1;
	}
	// Declare the supported options.
	
	po::variables_map vm;
	try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "get help message")
            ("alignment", po::value<string>(&alignment_fname), "procrastAligner alignment")
			("alus", po::value<string>(&alu_fname), "repeatmasker ALUs")
        ;

                
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }

        
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
	
	

	ifstream align_in;
	align_in.open( alignment_fname.c_str() );
	if( !align_in.is_open() ){
		cerr << "Error opening " << alignment_fname << endl;
		return -1;
	}
	ifstream alu_in;
	alu_in.open( alu_fname.c_str() );
	if( !alu_in.is_open() ){
		cerr << "Error opening " << alu_fname << endl;
		return -1;
	}
try{
	cout << "Calclutating specificity and sensitivity of procrastAligner on dataset..." << endl;
	IntervalList alignment;
	vector<AluRecord*> alus;
	
	//total length of all aligned repeats found by procrastAligner
	gnSeqI lt = 0;
	//total length of all alignments found by procrastAligner
	gnSeqI ld = 0;
	//total length of repeats found by repeatmasker
	gnSeqI lr=0;
	gnSeqI ln=0;
	gnSeqI lp=0;
	//total length of all regions found only in procrastAligner
	gnSeqI lo=0;
	//total length of repeats masked by both programs
	gnSeqI lc = 0;
	ReadAluFile( alu_in, alus, lr );
	alu_in.close();
	string cur_line;
	uint seqI = 0;
    //this will suffice for now, but should plan on using
	//IntervalList::ReadStandardAlignment or equivalent
	//to read in XMFA formatted output from procrastAligner
	pair<int64,int64> pos;
	vector< vector< pair<int64,int64> > > align_list;
	vector< pair<int64,int64> > pos_list;
	map<int64,bool> alncoverage;
    map<int64,bool> coverage;
	//list of maps, one for each alignment
	vector< map<int64,bool> > totcoverage;
	int64 ccount = 0;
	while( getline( align_in, cur_line) )
	{
		vector< int64 > start_list;
		getline( align_in, cur_line);
		stringstream parse_str( cur_line );
		int64 start = 0;
		int64 end = 0;
		int64 length = 0;
		string aln_len_str;
		parse_str >> aln_len_str;
		while( parse_str >> start )
		{
			start_list.push_back(start);
		}
		getline( align_in, cur_line);
		stringstream parse_string(cur_line);
		//parse_str.( cur_line );
		string lens;
		parse_string >> lens;
		uint region_count = 0;	
		while( parse_string >> length )
		{
			//cout << length << endl;
			if ( region_count >= start_list.size() )
			{
				//something's wrong
				cout << "alu data failed!" << endl;
				break;
			}
			pos.first = start_list.at(region_count);
			if (start_list.at(region_count) < 0 )
			{
				pos.second = start_list.at(region_count)-length;
				//simply add up the alignment coverage in the map
				for(int i = 0; i < length; i++)
				{
					alncoverage[pos.first-i] = true;
					coverage[pos.first-i] = true;		
					ccount++;
				}
			}
			else
			{
				pos.second = start_list.at(region_count)+length;
				//for both strands
				for(int i = 0; i < length; i++)
				{	
					alncoverage[pos.first+i] = true;
					coverage[pos.first+i] = true;
					ccount++;
				}
			}
			pos_list.push_back(pos);
			region_count++;
		}
		totcoverage.push_back(alncoverage);
		alncoverage.clear();
		align_list.push_back(pos_list);
		pos_list.clear();
		
	}//end of read procrastAligner output hack
	//alignment.ReadStandardAlignment( align_in );
	align_in.close();
	cout << "alu data processed!" << endl;
	int aluhits = 0;
	int matches = 0;
	//a first attempt at generating the sensitivity & specificity of our method
	//for comparison with zhang&waterman's eulerian path method...
	//hopefully we pull out these ~290bp repeats in a nice chain in each case
	//FIXME: is this ok?

	map<int,bool> ignoreAlignment;
	map<int64,bool> mergedCoverage;
	map<int64,bool> aluCoverage;

	//Total length of unaligned repeats(false positives?) found by procrastAligner
	map<int64,bool> lpt;
	map<int64,bool> lpn;

	//Total length of regions found only in procrastAligner
	map<int64,bool> lpo;

	
	map<int,bool> hitlist;

	map<int64,bool> specificity;

	map< uint,pair<int,int> > best_borders;
	map< uint,pair<int,int> > worst_borders;
	int64 matchhits = 0;
    int64 matchhitmult = 0;
	cout << "checking which alus are aligned" << endl;
	for ( int j = 0; j < align_list.size(); j++)
	{
		//if alufound in any component of curr alignment, consider 'aligned'
		//if not, throw out to help our sr. specificity
		
		//then, for each ALU, see if it is 'covered' by our procrastAlignments.
		//if so, increase lc2 accordingly

		//for each alignment returned by procrastAligner, highest multiplicity first
		bool alufound = false;

		//cout << "checking alignment #" << j << " for ALUs..." << endl;
		for ( int i = 0; i < alus.size(); i++)
		{
			
			//lpt = 0;
			if (alus.at(i)->strand == "+" )
			{
				for ( int a = 0; a < alus.at(i)->length(); a++)
				{
					

					//column in alignment coincides with an alu
					if(totcoverage.at(j).find((alus.at(i)->start)+a) != totcoverage.at(j).end())
					{
						
						alufound = true;
						//this column in sequence not accounted for
						if(aluCoverage.find((alus.at(i)->start)+a) == aluCoverage.end())
						{	
							lc+=1;
							
							//now it is
							
						}
						hitlist[i] = true;
						aluCoverage[(alus.at(i)->start)+a] = true;

					}
					
				}
			}
			else
			{
				for ( int a = 0; a < alus.at(i)->length(); a++)
				{
					if(totcoverage.at(j).find(-1*((alus.at(i)->start)+a)) != totcoverage.at(j).end())
					{
						if(aluCoverage.find(-1*((alus.at(i)->start)+a)) == aluCoverage.end())
						{	
							lc+=1;
							
							
						}
						hitlist[i] = true;
						aluCoverage[-1*((alus.at(i)->start)+a)] = true;
						//lc+=1;
						alufound =true;
					}
					
				}
			}
		}
		if(!alufound)
		{
			ignoreAlignment[j] = true;
			cout << "ignoring alignment " << j << endl;
			
			//calculate regions only appearing in procrastAligner alignments
			for(int k = 0; k < align_list.at(j).size();k++)
			{			
				gnSeqI len = absolut((int64)align_list.at(j).at(k).second)-absolut((int64)align_list.at(j).at(k).first);
				for(int n = 0; n<len;n++)
				{
					if(align_list.at(j).at(k).first<0)
						lpo[align_list.at(j).at(k).first-n] = true;
					else
						lpo[align_list.at(j).at(k).first+n] = true;
				}
			}
		}
		else
		{
			//cout << "ALU was aligned!" << endl;
			bool hit = false;
			bool debug_pos = false;
			bool inall = true;
			uint rnum = 0;
			for(int k = 0; k < align_list.at(j).size();k++)
			{			
				gnSeqI len = absolut((int64)align_list.at(j).at(k).second)-absolut((int64)align_list.at(j).at(k).first);
				for(int n = 0; n<len;n++)
				{
					if(align_list.at(j).at(k).first<0)
					{
						// j = lma #
						// k = component #
						// first,second = start,end pos
						if(aluCoverage.find(align_list.at(j).at(k).first-n)!= aluCoverage.end())
						{
							//find which alu is hit
							for ( int i = 0; i < alus.size(); i++)
							{	
								//is this ok for reverse strand?
								if( (abs((int)align_list.at(j).at(k).first) >= alus.at(i)->start) && (abs((int)align_list.at(j).at(k).first) < alus.at(i)->end ) 
								||  (abs((int)align_list.at(j).at(k).second) > alus.at(i)->start) && (abs((int)align_list.at(j).at(k).second) < alus.at(i)->end ) )
								{
									//the repeat #
									if (rnum != i+1 && rnum != 0)
										inall = false;
									rnum = i+1;				
									break;
								}
							}
							//current component of alignment pertains to alu
							//spec.at(j).push_back(k)
							lpn[align_list.at(j).at(k).first-n] = true;
							hit = true;
						}
						//motif missed by procrastAligner
						else
						{
							lpt[align_list.at(j).at(k).first-n] = true;
							rnum = -1;
						}
						mergedCoverage[align_list.at(j).at(k).first-n] = true;
					}
					else
					{
						if(aluCoverage.find(align_list.at(j).at(k).first+n)!= aluCoverage.end())
						{
							//find out which alu is hit
							for ( int i = 0; i < alus.size(); i++)
							{
								
								if( (abs((int)align_list.at(j).at(k).first) >= alus.at(i)->start) && (abs((int)align_list.at(j).at(k).first) < alus.at(i)->end )
								||  (abs((int)align_list.at(j).at(k).second) > alus.at(i)->start) && (abs((int)align_list.at(j).at(k).second) <= alus.at(i)->end ) )
								{
									//the repeat #
									//cout << rnum << " " << i+1 <<  endl;
									if (rnum != i+1 && rnum != 0)
										inall = false;
									rnum = i+1;
									break;
								}
							}	
							//current component of alignment pertains to alu
							lpn[align_list.at(j).at(k).first+n] = true;
							hit = true;
						}
						//motif missed by procrastAligner
						else
						{
							lpt[align_list.at(j).at(k).first+n] = true;
							rnum = -1;
						}
						mergedCoverage[align_list.at(j).at(k).first+n] = true;
					}
				}
				if (rnum <= 0)
					inall = false;

				if(hit)
				{
					matchhits+=1;
					
				}
			}
			//punt: DONT need to first check if it hits all components!!
			if (inall || 1)
			{
				for(int k = 0; k < align_list.at(j).size();k++)
				{			
					gnSeqI len = absolut((int64)align_list.at(j).at(k).second)-absolut((int64)align_list.at(j).at(k).first);
					uint rnum = 0;
					
					if(align_list.at(j).at(k).first<0)
					{
						// j = lma #
						// k = component #
						// first,second = start,end pos
						
						//find which alu is hit
						for ( int i = 0; i < alus.size(); i++)
						{
							//is this ok for reverse strand?
							if( (abs((int)align_list.at(j).at(k).first) >= alus.at(i)->start) && (abs((int)align_list.at(j).at(k).first) < alus.at(i)->end ) 
							||  (abs((int)align_list.at(j).at(k).second) > alus.at(i)->start) && (abs((int)align_list.at(j).at(k).second) <= alus.at(i)->end ) )
							{
								//the repeat #
								rnum = i+1;
								//find overlap
								int leftend = 0;
								int rightend = 0;
								leftend = abs((int)alus.at(i)->start)-abs((int)align_list.at(j).at(k).first);
								rightend =   abs((int)alus.at(i)->end)-abs((int)align_list.at(j).at(k).second);
								if (debug_pos && (abs(leftend)>500 || abs(rightend)>500))
								{
									cout << "alu\talignment" << endl;
									cout << alus.at(i)->start << "\t" << align_list.at(j).at(k).first << endl;
									cout << alus.at(i)->end << "\t" << align_list.at(j).at(k).second << endl;

								}
								
								if ( worst_borders.find( rnum ) != worst_borders.end() )
								{
									// if component has worse boundaries for this alu, record them
									if ( abs((int)worst_borders[rnum].first) < abs((int)leftend) )
										worst_borders[rnum].first = leftend;
									if ( abs((int)worst_borders[rnum].second) < abs((int)rightend) )
										worst_borders[rnum].second = rightend;
									if ( abs((int)best_borders[rnum].first) > abs((int)leftend) )
										best_borders[rnum].first = leftend;
									if ( abs((int)best_borders[rnum].second) > abs((int)rightend) )
										best_borders[rnum].second = rightend;
								}
								else
								{
									worst_borders[rnum] = make_pair(leftend,rightend);
									best_borders[rnum] = make_pair(leftend,rightend);
								}
								
								break;
							}
						} 
					}
					else
					{
						//find out which alu is hit
						for ( int i = 0; i < alus.size(); i++)
						{
							//if( (abs((int)align_list.at(j).at(k).first) <= alus.at(i)->start) && (abs((int)align_list.at(j).at(k).second) >= alus.at(i)->end ) )
							//if( (abs((int)align_list.at(j).at(k).first) >= alus.at(i)->start) && (abs((int)align_list.at(j).at(k).second) <= alus.at(i)->end ) )
							//if( ((abs((int)align_list.at(j).at(k).first) >= alus.at(i)->start) && (abs((int)align_list.at(j).at(k).first) < alus.at(i)->end ) && (abs((int)align_list.at(j).at(k).second) > alus.at(i)->end) )
							//||  ((abs((int)align_list.at(j).at(k).second) > alus.at(i)->start) && (abs((int)align_list.at(j).at(k).second) <= alus.at(i)->end ) && (abs((int)align_list.at(j).at(k).first) < alus.at(i)->start) ) )
							if( (abs((int)align_list.at(j).at(k).first) >= alus.at(i)->start) && (abs((int)align_list.at(j).at(k).first) < alus.at(i)->end ) 
							||  (abs((int)align_list.at(j).at(k).second) > alus.at(i)->start) && (abs((int)align_list.at(j).at(k).second) <= alus.at(i)->end ) )
							{
								//the repeat #
								rnum = i+1;
								//find overlap
								int leftend = 0;
								int rightend = 0;
								
								leftend = abs((int)alus.at(i)->start) -abs((int)align_list.at(j).at(k).first);
								rightend =   abs((int)alus.at(i)->end)-abs((int)align_list.at(j).at(k).second);

								if (debug_pos && (abs(leftend)>500 || abs(rightend)>500))
								{
									cout << "alu\talignment" << endl;
									cout << alus.at(i)->start << "\t" << align_list.at(j).at(k).first << endl;
									cout << alus.at(i)->end << "\t" << align_list.at(j).at(k).second << endl;

								}
								
								if ( worst_borders.find( rnum ) != worst_borders.end() )
								{
									// if component has worse boundaries for this alu, record them		
									if ( abs((int)worst_borders[rnum].first) < abs((int)leftend) )
										worst_borders[rnum].first = leftend;
									if ( abs((int)worst_borders[rnum].second) < abs((int)rightend) )
										worst_borders[rnum].second = rightend;

									// if component has better boundaries for this alu, record them
									if ( abs((int)best_borders[rnum].first) > abs((int)leftend) )
										best_borders[rnum].first = leftend;
									if ( abs((int)best_borders[rnum].second) > abs((int)rightend) )
										best_borders[rnum].second = rightend;
									
								}
								else
								{
									worst_borders[rnum] = make_pair(leftend,rightend);
									best_borders[rnum] = make_pair(leftend,rightend);
								}
								
								break;
							}
						}
					
					}

				}
			}
		    
		}
		alufound = false;
	}
	gnSequence empty_seq;
	//this is the length of the repeats found by procrastAligner, 
	//with overlaps removed
	//remember the alignments to ignore!
	ofstream boundary_file;
	alignment_fname.append(".boundary");
	boundary_file.open(alignment_fname.c_str());
	map< uint,pair<int,int> >::iterator iter;
	uint avg_worst_left = 0;
	uint avg_worst_right = 0;
	uint avg_best_left = 0;
	uint avg_best_right = 0;
	for( iter = worst_borders.begin(); iter != worst_borders.end(); iter++ ) 
	{
		avg_worst_left += abs(iter->second.first);
		avg_worst_right += abs(iter->second.second);
		boundary_file << "worst boundaries for repeat copy #" << iter->first << "\t left: " << iter->second.first << "\t right: " << iter->second.second << endl;
	}
	for( iter = best_borders.begin(); iter != best_borders.end(); iter++ ) 
	{
		avg_best_left += abs( iter->second.first);
		avg_best_right += abs(iter->second.second);
		boundary_file << "best boundaries for repeat copy #" << iter->first << "\t left: " << iter->second.first << "\t right: " << iter->second.second << endl;
	}

	if (worst_borders.size() > 0 )
	{
		avg_worst_left /= worst_borders.size();
		avg_worst_right /= worst_borders.size();
	}
	else
	{
		avg_worst_left = -1;
		avg_worst_right = -1;

	}
	if ( best_borders.size() > 0)
	{
		avg_best_left /= best_borders.size();
		avg_best_right /= best_borders.size();
	}
	else
	{
		avg_best_left = -1;
		avg_best_right = -1;

	}
	boundary_file << "left best: \t" << avg_best_left << endl;
	boundary_file << "right best: \t" << avg_best_right << endl;
	boundary_file << "left worst: \t" << avg_worst_left << endl;
	boundary_file << "right worst: \t" << avg_worst_right << endl;
	boundary_file << "#" << endl;
	boundary_file.close();

	lt = mergedCoverage.size();
	//lt2 = coverage.size();
	ld = coverage.size();
	lp = lpt.size();
	ln = lpn.size();
	lo = lpo.size();

	//length of only ALUs hit by procrastAligner
	gnSeqI hitlength =0;
	for(int i =0; i< hitlist.size(); i++)
		hitlength+= alus.at(i)->length();

	cout << "\nprocrastAlignments processed: " << align_list.size() << endl;
	cout << "matches processed: " << matches << endl;
	cout << "Total ALUs found by repeatmasker: " << alus.size() << endl;
	cout << "Total ALUs hit by procrastAligner: " << hitlist.size() << endl;
	cout << "ALU hit percentage: " << (float)hitlist.size()/(float)alus.size() << endl;

	//cout << aluCoverage.size() << endl;
    cout << "\nTotal length of all repeats found by procrastAligner: " << ld << endl;
	 cout << "Total length of all regions found only in procrastAligner: " << lo << endl;
	cout << "Total length of all (partially) aligned repeats found by procrastAligner: lt = " << lt << endl;
	cout << "Total length of unaligned repeats(false positives?) found by procrastAligner: lp = " << lp << endl;
	//cout << "Total length of ???: ln = " << ln << endl;
	cout << "Total length of all repeats(ALU) found by repeatmasker: lr = " << lr << endl;
	cout << "Total length of repeats(ALU) found by repeatmasker hit by procrastAligner: lh =" << hitlength << endl;
	cout << "Total length of ALU repeats found by both methods: lc = " << lc << endl;
	
	//cout << "Sensitivity: lc / lr = " << (double)(lc) / (double)(lr) << endl;
	//cout << "Specificity: lc / lt = " << (double)(lc) / (double)(lt) << endl;
	
	//score changes per Sunday email, focus on filtration
	cout << "\nSensitivity-old: lc / lh = " << (double)(lc) / (double)(hitlength) << endl;
	cout << "Specificity-old: lc / lt = " << (double)(lc) / (double)(lt) << endl;

	
	cout << "\nSensitivity= " << (double)hitlist.size()/(double)alus.size() << endl;
	cout << "Specificity= " << (double)matchhits/(double)matchhitmult <<  endl;


	//TN = ltn
	//TP = lc
	//FN = lfn
	//FP = lp
}catch( gnException& gne ){
	cerr << gne << endl;
}catch( exception& e ){
	cerr << e.what() << endl;
}

}


