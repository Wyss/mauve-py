#include "libMems/IntervalList.h"
#include "libMems/Islands.h"
#include "libMems/DistanceMatrix.h"
#include <sstream>
#include <fstream>
using namespace std;
using namespace genome;
using namespace mems;

int main( int argc, char* argv[] )
{
	if( argc != 2 )
	{
		cerr << "Usage: pairCompare <sequence count>\n";
		return -1;
	}
	int seq_count = atoi( argv[1] );
	cout << "SeqI\tSeqJ\tNTidentity\tAvgBBpct\tLCB count\n";
	for( size_t seqI = 10; seqI < seq_count; seqI++ )
	{
		for( size_t seqJ = 0; seqJ < seq_count; seqJ++ )
		{
			if( seqJ <= seqI )
				continue;
			cout << seqI << '\t' << seqJ << '\t';

			size_t lcb_count = 0;

			stringstream aln_in_fname;
			aln_in_fname << "all_pairs/pair_" << seqI << "." << seqJ << ".xmfa";
			ifstream alignment_in(aln_in_fname.str().c_str());
			IntervalList aligned_ivs;
			aligned_ivs.ReadStandardAlignment( alignment_in );


			LoadSequences(aligned_ivs, NULL);

			// add the sequence data to the interval list
			uint seq_count = aligned_ivs.seq_table.size();
			vector< GappedAlignment > backbone_data;
			simpleFindBackbone( aligned_ivs, 50, 50, backbone_data );

			IntervalList backbone_ivs;
			backbone_ivs.seq_table = aligned_ivs.seq_table;

			// count up the total length of backbone in each genome
			vector< gnSeqI > total_bb( seq_count, 0 );
			NumericMatrix< double > overall_identity;
			for( uint bbI = 0; bbI < backbone_data.size(); bbI++ ){
				vector<AbstractMatch*> tmp_iv(1, &backbone_data[ bbI ]);
				backbone_ivs.push_back( Interval( tmp_iv.begin(), tmp_iv.end() ) );
				for( uint seqI = 0; seqI < seq_count; seqI++ ){
					total_bb[ seqI ] += backbone_data[ bbI ].Length( seqI );
				}
			}

			vector< AbstractMatch* > bbivs;
			for( uint bbI = 0; bbI < backbone_ivs.size(); bbI++ )
				bbivs.push_back( &backbone_ivs[bbI] );
			BackboneIdentityMatrix( bbivs, aligned_ivs.seq_table, overall_identity );

			gnSeqI avg_bb = 0;
			double seq_len_average = 0;
			for( uint seqI = 0; seqI < aligned_ivs.seq_table.size(); seqI++ ){
				avg_bb += total_bb[ seqI ];
				seq_len_average += aligned_ivs.seq_table[seqI]->length();
			}
			avg_bb /= aligned_ivs.seq_table.size();
			seq_len_average /= (double)seq_count;


			for( size_t lcbI = 0; lcbI < aligned_ivs.size(); lcbI++ )
				if( aligned_ivs[lcbI].Multiplicity() > 1 )
					lcb_count++;


			cout << overall_identity(0,1) << '\t';
			cout << avg_bb / seq_len_average << '\t';
			cout << lcb_count << endl;


		}
	}

}
