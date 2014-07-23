/**
 * multiEVD
 * (c)left 2007 aaron darling
 * A program to calculate the extreme value distribution of alignment drops in homologous sequence.
 * INPUT: a simulated multiple alignment as input
 * OUTPUT: the 95%ile, 99%ile, etc of scores in the extreme value distribution
 * THEORY:
 * computes inverse substitution and gap scores, never allowing
 * the inverse score to drop below 0.  Each time the score rises above 0, an "excursion" begins, and when the score
 * drops back to 0, the excursion has ended.  The highest score achieved by the excursion is the "extreme value".
 * Each extreme value is recorded, and the distribution of these extreme values is what gets output.
 */
#include "libMems/Islands.h"
#include "libMems/IntervalList.h"
#include "libMems/MatchList.h"
#include "libMems/MuscleInterface.h"
#include "libGenome/gnSequence.h"
#include "libMems/MatchProjectionAdapter.h"
#include "libMems/ProgressiveAligner.h"

#include <sstream>

#include <boost/multi_array.hpp>

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
		
		std::vector< score_t > scores;
		PairwiseScoringScheme pss;
		score_t total_score;

		stripGapColumns(aln_table);
		computeSPScore( aln_table, pss, scores, total_score );

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

//bad: copied from progressiveAligner.cpp
template< class BoostMatType >
void print2d_matrix( BoostMatType& mat, std::ostream& os )
{
	for( size_t i = 0; i < mat.shape()[0]; ++i )
	{
		for( size_t j = 0; j < mat.shape()[1]; ++j )
		{
			if( j > 0 )
				os << "\t";
			os << mat[i][j];
		}
		os << endl;
	}
}


// read each input file, write summary statistics about the EVD to stdout
int main( int argc, char* argv[] )
{
//	vector< score_t > lrh_all;
	if( argc != 2 )
	{
		cerr << "Usage: multiEVD <simulation run count>\n";
		cerr << "This program must be run from a directory which contains alignjob directories\n";
		return -1;
	}
	int run_count = atoi( argv[1] );
	int simu_count = 0;
	vector< vector< score_t > > lrh_all;
	size_t seq_count = 0;
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
		if( seq_count == 0 )
		{
			seq_count = iv_list.seq_table.size();
			lrh_all.resize(seq_count+1);
		}

		vector< Interval* > iv_ptrs( iv_list.size() );
		for( size_t ivI = 0; ivI < iv_list.size(); ++ivI )
			iv_ptrs[ivI] = &iv_list[ivI];

		vector< gnSequence* > seq_table = iv_list.seq_table;

		vector< uint > proj_seqs(seq_count);
		for( size_t sI = 0; sI < seq_count; ++sI )
			proj_seqs[sI] = sI;

		std::vector< std::vector< mems::MatchProjectionAdapter* > > LCB_list;
		std::vector< mems::LCB > projected_adjs;
		for( size_t mult = seq_count; mult > 1; mult-- )
		{
			vector< score_t > lrh;
			getLocalRecordHeights( iv_ptrs, seq_table, lrh );
			lrh_all[mult].insert( lrh_all[mult].end(), lrh.begin(), lrh.end() );
			// randomly pick a sequence to discard
			int disc = rand() % proj_seqs.size();
			proj_seqs.erase(proj_seqs.begin()+disc);
			seq_table.erase(seq_table.begin()+disc);
			// project the original alignment down to the remaining sequences
			projectIntervalList( iv_list, proj_seqs, LCB_list, projected_adjs );
			// free storage used by the previous set of projections
			if( mult != seq_count )
			{
				for( size_t ivI = 0; ivI < iv_ptrs.size(); ivI++ )
					iv_ptrs[ivI]->Free();	
			}
			// update iv_ptrs to contain the new projections
			iv_ptrs.resize(LCB_list.size());
			for( size_t lcbI = 0; lcbI < LCB_list.size(); lcbI++ )
			{
				Interval iv;
				iv_ptrs[lcbI] = iv.Copy();
				iv_ptrs[lcbI]->SetMatches(LCB_list[lcbI]);
			}
		}
	}

	boost::multi_array<score_t, 2> evd_table;
	evd_table.resize( boost::extents[4][seq_count-1] );
	boost::multi_array<size_t, 2> ss_table;
	ss_table.resize( boost::extents[4][seq_count-1] );
	for( size_t mult = 2; mult < seq_count + 1; mult++ )
	{
		std::sort( lrh_all[mult].begin(), lrh_all[mult].end() );
		size_t index_95 = lrh_all[mult].size() * .95;
		size_t index_99 = lrh_all[mult].size() * .99;
		size_t index_999 = lrh_all[mult].size() * .999;
		size_t index_9999 = lrh_all[mult].size() * .9999;
		index_95 = (std::min)(index_95, lrh_all[mult].size()-1);
		index_99 = (std::min)(index_99, lrh_all[mult].size()-1);
		index_999 = (std::min)(index_999, lrh_all[mult].size()-1);
		index_9999 = (std::min)(index_9999, lrh_all[mult].size()-1);
//		cout << "Total number of simulations: " << simu_count << endl;
//		cout << "Total number of excursions: " << lrh_all[mult].size() << endl;
//		cout << "95% score threshold: " << lrh_all[mult][index_95] << endl;
		evd_table[0][mult-2] = lrh_all[mult][index_95];
//		cout << "Number excursions above 95%: " << lrh_all[mult].size() - index_95 << endl;
		ss_table[0][mult-2] = lrh_all[mult].size() - index_95;
//		cout << "99% score threshold: " << lrh_all[mult][index_99] << endl;
		evd_table[1][mult-2] = lrh_all[mult][index_99];
//		cout << "Number excursions above 99%: " << lrh_all[mult].size() - index_99 << endl;
		ss_table[1][mult-2] = lrh_all[mult].size() - index_99;
//		cout << "99.9% score threshold: " << lrh_all[mult][index_999] << endl;
		evd_table[2][mult-2] = lrh_all[mult][index_999];
//		cout << "Number excursions above 99.9%: " << lrh_all[mult].size() - index_999 << endl;
		ss_table[2][mult-2] = lrh_all[mult].size() - index_999;
//		cout << "99.99% score threshold: " << lrh_all[mult][index_9999] << endl;
		evd_table[3][mult-2] = lrh_all[mult][index_9999];
//		cout << "Number excursions above 99.99%: " << lrh_all[mult].size() - index_9999 << endl;
		ss_table[3][mult-2] = lrh_all[mult].size() - index_9999;
	}
	cout << "Matrix of score thresholds:\n";
	print2d_matrix( evd_table, cout );
	cout << "\n\nMatrix of sample sizes:\n";
	print2d_matrix( ss_table, cout );
	cout << endl;
}


