#include "libMems/Backbone.h"
using namespace mems;
using namespace std;
using namespace genome;

typedef pair< bb_seqentry_t, size_t > labeled_bb_t;

class BbSorter
{
public:
	BbSorter( size_t seqI ){ m_seq = seqI; }
	bool operator()( const labeled_bb_t& a, const labeled_bb_t& b )
	{
		return genome::absolut(a.first[m_seq].first) < genome::absolut(b.first[m_seq].first);
	}
	size_t m_seq;
};



class ShorterThan {
public:
	bool operator()( const bb_seqentry_t& a )
	{
		size_t sc = 0;
		size_t tot = 0;
		for( size_t i = 0; i < a.size(); i++ )
			if( a[i].first != 0 )
			{
				tot += genome::absolut(a[i].second - a[i].first) + 1;
				sc++;
			}
		if( tot == 0 )
			return true;
		return (tot / sc) < 20;
	}
};

int main( int argc, char* argv[] )
{
#if	WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	if( argc < 4 )
	{
		cerr << "bbFilter <backbone file> <independent dist> <output file> <beast|gp> [<seq1> <seq2>...<seqN>]\n";
		cerr << "seq index starts at 0.\n";
		cerr << "\nExample:\n";
		cerr << "bbFilter my_alignment.backbone 50 my_feats.bin gp\n";
		cerr << "the above command extracts binary features from \"my_alignment.backbone\" which are separated by a minimum of 50nt sequence conserved among all taxa in the alignment.  The output is written to my_feats.bin in genoplast format\n";
		cerr << "\n\nExample 2:\nbbFilter aln.backbone 100 feats.xml beast 0 1 2 5 6\n";
		cerr << "the above command extracts binary features from \"aln.backbone\" which are separated by a minimum of 100nt sequence conserved among genomes 0,1,2,5, and 6 from the alignment.  The output is written to feats.xml in beast format\n";
		return -1;
	}
	string bbseq_fname( argv[1] );
	int indie_dist = atoi( argv[2] );
	string output_fname( argv[3] );
	string target_format( argv[4] );
	bool allow_alternalogs = true;
	bool check_independence = false;

	ifstream bbseq_input( bbseq_fname.c_str() );
	if( !bbseq_input.is_open() ){
		cerr << "Error opening \"" << bbseq_fname << "\"" << endl;
		return -4;
	}
	ofstream anal_output( output_fname.c_str() );
	if( !anal_output.is_open() ){
		cerr << "Error opening \"" << output_fname << "\" for writing" << endl;
		return -6;
	}
	
	// read the backbone column file	
	vector< bb_seqentry_t > bb_seq_list;
	readBackboneSeqFile( bbseq_input, bb_seq_list );

	// read the list of seqs of interest
	vector< int > seqs;
	for( int i = 5; i < argc; i++ )
		seqs.push_back(atoi(argv[i]));

	// assume all seqs are of interest
	if( seqs.size() == 0 && bb_seq_list.size() > 0 )
	{
		for( int i = 0; i < bb_seq_list[0].size(); i++ )
			seqs.push_back(i);
	}
	// add any genome-specific segments
	addUniqueSegments( bb_seq_list );

	// remove short segments
	ShorterThan st;
	vector< bb_seqentry_t >::iterator new_end = std::remove_if( bb_seq_list.begin(), bb_seq_list.end(), st );
	cout << "Removing " << bb_seq_list.end() - new_end << " features shorter than 20 nt\n";
	bb_seq_list.erase( new_end, bb_seq_list.end() );

	// now assign tracking IDs to the backbone segments
	vector< labeled_bb_t > bb_segs;
	for( size_t i = 0; i < bb_seq_list.size(); i++ )
	{
		bb_segs.push_back( make_pair( bb_seq_list[i], i ) );
	}

	// create a sorted list for each genome and a map to the segment ID
	vector< vector< labeled_bb_t > > sorted_segs( seqs.size(), bb_segs );
	vector< vector< size_t > > seg_id_maps( seqs.size(), vector< size_t >( bb_segs.size() ) );
	for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
	{
		BbSorter bbs(seqs[seqI]);
		std::sort( sorted_segs[seqI].begin(), sorted_segs[seqI].end(), bbs );
		for( size_t bbI = 0; bbI < sorted_segs[seqI].size(); bbI++ )
			seg_id_maps[ seqI ][ sorted_segs[seqI][bbI].second ] = bbI;
	}


	bitset_t good_bb( bb_seq_list.size() );
	bitset_t nway( bb_seq_list.size() );
	bitset_t nunya( bb_seq_list.size() );

	// mark anything that has all of the seqs or none of seqs as not useful
	for( size_t bbI = 0; bbI < bb_seq_list.size(); bbI++ )
	{
		bool all = true;
		bool none = true;
		for( size_t sI = 0; sI < seqs.size(); sI++ )
		{
			if( bb_seq_list[bbI][seqs[sI]].first == 0 )
				all = false;
			else
				none = false;
		}
		if(all)
			nway.set(bbI);
		if(none)
			nunya.set(bbI);
	}
	good_bb = nway | nunya;
	good_bb.flip();
	
	// now mark segs that are too close to each other to be considered independent
	for( size_t sI = 0; check_independence && sI < seqs.size(); sI++ )
	{
		BbSorter bbs(seqs[sI]);
		std::sort( bb_segs.begin(), bb_segs.end(), bbs );
		for( size_t bbI = 1; bbI < bb_segs.size()-1; bbI++ )
		{
			if( nway[bb_segs[bbI].second] )
				continue;
			if( bb_segs[bbI].first[seqs[sI]].first == 0 )
				continue;
			// ensure that it has n-way on both sides and that they are at least "indie_dist" long
			if( nway.test(bb_segs[bbI-1].second) && 
				nway.test(bb_segs[bbI+1].second) &&
				absolut(bb_segs[bbI-1].first[seqs[sI]].second - bb_segs[bbI-1].first[seqs[sI]].first) >= indie_dist &&
				absolut(bb_segs[bbI+1].first[seqs[sI]].second - bb_segs[bbI+1].first[seqs[sI]].first) >= indie_dist )
			{
				if( !allow_alternalogs ){
					// ensure that there is no other feature in the other genomes
					for( size_t k = 0; k < seqs.size(); k++ )
					{
						if( k == sI )
							continue;
						size_t oid = seg_id_maps[k][ bb_segs[bbI-1].second ];
						int parity = ((bb_segs[bbI-1].first[seqs[sI]].first > 0 && bb_segs[bbI-1].first[seqs[k]].first > 0) ||
							(bb_segs[bbI-1].first[seqs[sI]].first < 0 && bb_segs[bbI-1].first[seqs[k]].first < 0)) ? 1 : -1;
						size_t prev_in_sI = bb_segs[bbI-1].second;
						size_t cur_in_sI = bb_segs[bbI].second;
						size_t next_in_sI = bb_segs[bbI+1].second;
						size_t prev_in_k = sorted_segs[k][oid].second;
						size_t cur_in_k = sorted_segs[k][oid+parity].second;
						size_t next_in_k = sorted_segs[k][oid+parity*2].second;
						if( (cur_in_sI == cur_in_k && next_in_sI == next_in_k) ||
							(next_in_sI == cur_in_k))					
							; // it's good because no other segs intervene
						else
						{
							good_bb.set( bb_segs[bbI].second, false );
							break;	// it's an alternalog or overlapping, no sense in checking other seqs
						}
					}
				}
			}else
				good_bb.set(bb_segs[bbI].second, false);
		}
	}

	// create site patterns, then write out the good ones
	bitset_t empty( bb_seq_list.size() );
	vector< bitset_t > spa_seqs( seqs.size(), empty );
	for( size_t bbI = 0; bbI < bb_seq_list.size(); bbI++ )
		for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
			spa_seqs[seqI].set(bbI, bb_seq_list[bbI][seqs[seqI]].first != 0);

	vector< string > binseqs( seqs.size(), string( good_bb.count(), '0' ) );
	for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
	{
		size_t goodI = 0;
		for( size_t bbI = 0; bbI < good_bb.size(); bbI++ )
			if(good_bb.test(bbI))
			{
				if(spa_seqs[seqI].test(bbI))
					binseqs[seqI][goodI] = '1';
				goodI++;
			}
	}
	map< string, int > sitepattern_count;
	// count how many segments of each site pattern
	for( size_t bbI = 0; bbI < good_bb.size(); bbI++ )
	{
		if(!good_bb.test(bbI))	continue;
		size_t length=0;
		size_t sc=0;
		string sitepat( seqs.size(), '0' );
		for( int seqI = 0; seqI < seqs.size(); seqI++ )
		{
			sitepat[seqI] = spa_seqs[seqI][bbI] ? '1' : '0';
			if(spa_seqs[seqI][bbI]){
				length += genome::absolut(bb_seq_list[bbI][seqI].second - bb_seq_list[bbI][seqI].first);
				sc++;
			}
		}
		length /= sc;
		map< string, int >::iterator iter = sitepattern_count.find(sitepat);
		if(iter == sitepattern_count.end())
			sitepattern_count.insert( make_pair( sitepat, length ) );
		else
			iter->second+= length;
	}

	// write out the seqs!!
	if( target_format == "beast" )
	{
		anal_output << "\t<taxa id=\"taxa\">\n";
		for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
		{
			anal_output << "\t\t<taxon id=\"seq" << seqI << "\"/>\n";
	
		}
		anal_output << "\t</taxa>\n";
		anal_output << "\t<alignment id=\"alignment\" dataType=\"binary\">\n";
	
		for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
		{
			anal_output << "\t\t<sequence>\n";
			anal_output << "\t\t\t<taxon idref=\"seq" << seqI << "\"/>\n";
			anal_output << "\t\t\t" << binseqs[seqI] << endl;
			anal_output << "\t\t</sequence>\n";
//			anal_output << "> seq" << seqI << endl;
//			for( size_t i = 0; i < binseqs[seqI].size(); i+=80 )
//				anal_output << binseqs[seqI].substr(i, 80) << endl;
		}
		anal_output << "\t</alignment>\n";
	}else{
		// write out a header line with the number of times each site pattern is used.
		map<string,int>::iterator f = sitepattern_count.begin();
		for(; f!= sitepattern_count.end(); f++){
			if(f!=sitepattern_count.begin())	anal_output << ' ';
			anal_output << (f->second / 20);
		}
		anal_output << endl;
		// write genoplast format
		for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
		{
			f = sitepattern_count.begin();
			for(; f!= sitepattern_count.end(); f++){
				if(f!=sitepattern_count.begin())	anal_output << ' ';
				anal_output << f->first[seqI];
			}
			anal_output << endl;
		}
	}

	anal_output.close();

	string loc_fname = output_fname + ".locs";
	ofstream location_output( loc_fname.c_str() );
	for( size_t bbI = 0; bbI < good_bb.size(); bbI++ )
	{
		if( good_bb.test(bbI) )
		{
			for( size_t seqI = 0; seqI < seqs.size(); seqI++ )
			{
				if( seqI > 0 )
					location_output << '\t';
				location_output << bb_seq_list[bbI][seqI].first << '\t' << bb_seq_list[bbI][seqI].second;
			}
			location_output << std::endl;
		}
	}
}

