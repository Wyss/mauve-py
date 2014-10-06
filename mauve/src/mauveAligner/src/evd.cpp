#include "libMems/Islands.h"
#include "libMems/IntervalList.h"
#include "libMems/MatchList.h"
#include "libGenome/gnSequence.h"

#include <sstream>

using namespace std;
using namespace genome;
using namespace mems;

template< typename MatchVector >
void getLocalRecordHeights( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, vector< score_t >& lrh )
{
	typedef typename MatchVector::value_type MatchType;
	if( iv_list.size() == 0 )
		return;
	uint seq_count = seq_table.size();
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		const MatchType& iv = iv_list[ iv_listI ];
		std::vector< std::string > aln_table;
		GetAlignment( *iv, seq_table, aln_table );
		
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			uint seqJ;
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){

				std::vector< score_t > scores;
				PairwiseScoringScheme pss;
				computeMatchScores( aln_table[seqI], aln_table[seqJ], pss, scores );
				computeGapScores( aln_table[seqI], aln_table[seqJ], pss, scores );

				// Invert the scores since we're trying to detect rare bouts of non-homologous sequence
				for( size_t sI = 0; sI < scores.size(); ++sI )
					if( scores[sI] != INVALID_SCORE)
						scores[sI] = -scores[sI];

				score_t score_sum = 0;	// start in an hss
				score_t local_record_height = 0;
				for( size_t colI = 0; colI < scores.size(); ++colI )
				{
					if( scores[colI] == INVALID_SCORE )
						continue;

					if( score_sum > 0 && score_sum + scores[colI] < 0 )
					{
						// end of an excursion
						score_sum = 0;
						lrh.push_back( local_record_height );
						local_record_height = 0;
					}else if( score_sum == 0 && scores[colI] > 0 )
					{
						// start a new excursion
						score_sum += scores[colI];
						if( score_sum > local_record_height )
							local_record_height = score_sum;
					}else if( score_sum > 0 ){
						score_sum += scores[colI];
						if( score_sum > local_record_height )
							local_record_height = score_sum;
					}
				}
			}
		}
	}
}


// read each input file, write summary statistics about the EVD to stdout
int main( int argc, char* argv[] )
{
	vector< score_t > lrh_all;
	if( argc != 2 )
	{
		cerr << "Usage: evd <simulation run count>\n";
		cerr << "This program must be run from a directory which contains alignjob directories\n";
		return -1;
	}
	int run_count = atoi( argv[1] );
	int simu_count = 0;
	for( int runI = 0; runI < run_count; ++runI )
	{
		IntervalList iv_list;
		stringstream aln_fname;
		aln_fname << "alignjob." << runI << "/evolved.dat";
		ifstream in_file( aln_fname.str().c_str() );
		if( !in_file.is_open() )
		{
			cerr << "Error opening " << aln_fname.str() << endl;
			continue;
		}
		simu_count++;
		iv_list.ReadStandardAlignment(in_file);
		stringstream seq_fname;
		seq_fname << "alignjob." << runI << "/evolved_seqs.fas";
		MatchList ml;
		LoadMFASequences(ml, seq_fname.str(), &cout);
		iv_list.seq_table = ml.seq_table;

		vector< Interval* > iv_ptrs( iv_list.size() );
		for( size_t ivI = 0; ivI < iv_list.size(); ++ivI )
			iv_ptrs[ivI] = &iv_list[ivI];

		vector< score_t > lrh;
		getLocalRecordHeights( iv_ptrs, iv_list.seq_table, lrh );
		lrh_all.insert( lrh_all.end(), lrh.begin(), lrh.end() );
	}
	std::sort( lrh_all.begin(), lrh_all.end() );
	size_t index_95 = lrh_all.size() * .95;
	size_t index_99 = lrh_all.size() * .99;
	size_t index_999 = lrh_all.size() * .999;
	size_t index_9999 = lrh_all.size() * .9999;
	index_95 = std::min(index_95, lrh_all.size()-1);
	index_99 = std::min(index_99, lrh_all.size()-1);
	index_999 = std::min(index_999, lrh_all.size()-1);
	index_9999 = std::min(index_9999, lrh_all.size()-1);
	cout << "Total number of simulations: " << simu_count << endl;
	cout << "Total number of excursions: " << lrh_all.size() << endl;
	cout << "95% score threshold: " << lrh_all[index_95] << endl;
	cout << "Number excursions above 95%: " << lrh_all.size() - index_95 << endl;
	cout << "99% score threshold: " << lrh_all[index_99] << endl;
	cout << "Number excursions above 99%: " << lrh_all.size() - index_99 << endl;
	cout << "99.9% score threshold: " << lrh_all[index_999] << endl;
	cout << "Number excursions above 99.9%: " << lrh_all.size() - index_999 << endl;
	cout << "99.99% score threshold: " << lrh_all[index_9999] << endl;
	cout << "Number excursions above 99.99%: " << lrh_all.size() - index_9999 << endl;
}


