#include "libMems/Backbone.h"
#include "libMems/ProgressiveAligner.h"
#include <sstream>
using namespace mems;
using namespace std;
using namespace genome;


template< typename MatchVector >
void getBpList( MatchVector& mvect, uint seq, vector< gnSeqI >& bp_list )
{
	bp_list.clear();
	for( size_t ivI = 0; ivI < mvect.size(); ivI++ )
	{
		if( mvect[ivI]->LeftEnd(seq) == NO_MATCH )
			continue;
		bp_list.push_back( mvect[ivI]->LeftEnd(seq) );
		bp_list.push_back( mvect[ivI]->RightEnd(seq)+1 );
	}
	std::sort( bp_list.begin(), bp_list.end() );
}

template< typename MatchVector >
void createMap( const MatchVector& mv_from, const MatchVector& mv_to, vector< size_t >& map )
{
	typedef typename MatchVector::value_type MatchPtr;
	vector< pair< MatchPtr, size_t > > m1(mv_from.size());
	vector< pair< MatchPtr, size_t > > m2(mv_to.size());
	for( size_t i = 0; i < mv_from.size(); ++i )
		m1[i] = make_pair( mv_from[i], i );
	for( size_t i = 0; i < mv_to.size(); ++i )
		m2[i] = make_pair( mv_to[i], i );
	std::sort( m1.begin(), m1.end() );
	std::sort( m2.begin(), m2.end() );
	map.resize( m1.size() );
	for( size_t i = 0; i < m1.size(); ++i )
		map[m1[i].second] = m2[i].second;
}


void makeAllPairwiseGenomeHSSBreakOnGenes( IntervalList& iv_list, vector< CompactGappedAlignment<>* >& iv_ptrs, vector< CompactGappedAlignment<>* >& iv_orig_ptrs, pairwise_genome_hss_t& hss_cols, const HssDetector* detector, vector< vector< gnSeqI > >& gene_bounds )
{
	uint seq_count = iv_list.seq_table.size();
	// make pairwise projections of intervals and find LCBs...
	for( size_t seqI = 0; seqI < seq_count; ++seqI )
	{
		for( size_t seqJ = seqI+1; seqJ < seq_count; ++seqJ )
		{
			vector< uint > projection;
			projection.push_back( seqI );
			projection.push_back( seqJ );
			vector< vector< MatchProjectionAdapter* > > LCB_list;
			vector< LCB > projected_adjs;
			projectIntervalList( iv_list, projection, LCB_list, projected_adjs );
			// make intervals
			IntervalList pair_ivs;
			pair_ivs.seq_table.push_back( iv_list.seq_table[seqI] );
			pair_ivs.seq_table.push_back( iv_list.seq_table[seqJ] );
			pair_ivs.resize( LCB_list.size() );
			for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
				pair_ivs[lcbI].SetMatches( LCB_list[lcbI] );
			LCB_list.clear();

			vector< CompactGappedAlignment<>* > pair_cgas( pair_ivs.size() );
			for( size_t lcbI = 0; lcbI < pair_ivs.size(); ++lcbI )
			{
				CompactGappedAlignment<> tmp_cga;
				pair_cgas[lcbI] = tmp_cga.Copy();
				new (pair_cgas[lcbI])CompactGappedAlignment<>( pair_ivs[lcbI] );
			}

			vector< CompactGappedAlignment<>* > hss_list;
			// now find islands
			hss_array_t hss_array;
			(*detector)( pair_cgas, pair_ivs.seq_table, hss_array );
			HssArrayToCga(pair_cgas, pair_ivs.seq_table, hss_array, hss_list);

			for( size_t cgaI = 0; cgaI < pair_cgas.size(); ++cgaI )
				pair_cgas[cgaI]->Free();
			pair_cgas.clear();

			// now split up on iv boundaries
			vector< gnSeqI > bp_list;
			getBpList( iv_ptrs, seqI, bp_list );
			GenericMatchSeqManipulator< CompactGappedAlignment<> > gmsm(0);
			SingleStartComparator< CompactGappedAlignment<> > ssc(0);
			std::sort(hss_list.begin(), hss_list.end(), ssc );
			applyBreakpoints( bp_list, hss_list, gmsm );
			// break on gene bounds in seqI
			std::sort(hss_list.begin(), hss_list.end(), ssc );
//			if( !(seqI == 1 && seqJ == 15 ) )
			applyBreakpoints( gene_bounds[seqI], hss_list, gmsm );
			// and again on seqJ
			getBpList( iv_ptrs, seqJ, bp_list );
			GenericMatchSeqManipulator< CompactGappedAlignment<> > gmsm1(1);
			SingleStartComparator< CompactGappedAlignment<> > ssc1(1);
			std::sort(hss_list.begin(), hss_list.end(), ssc1 );
			applyBreakpoints( bp_list, hss_list, gmsm1 );
			// break on gene bounds in seqJ
			std::sort(hss_list.begin(), hss_list.end(), ssc1 );
//			if( !(seqI == 1 && seqJ == 15 ) )
			applyBreakpoints( gene_bounds[seqJ], hss_list, gmsm1 );

			// now transform into interval-specific columns
			std::sort(hss_list.begin(), hss_list.end(), ssc );

			SingleStartComparator< CompactGappedAlignment<> > ivcomp(seqI);
			std::sort( iv_ptrs.begin(), iv_ptrs.end(), ivcomp );
			vector< size_t > iv_map;
			createMap( iv_ptrs, iv_orig_ptrs, iv_map );
			size_t ivI = 0;
			while( ivI < iv_ptrs.size() && iv_ptrs[ivI]->LeftEnd(0) == NO_MATCH )
				++ivI;
			for( size_t hssI = 0; hssI < hss_list.size(); ++hssI )
			{
				if( hss_list[hssI]->LeftEnd(0) == NO_MATCH || hss_list[hssI]->Length(0) == 0 )
					continue;
				if( ivI == iv_ptrs.size() )
				{
					cerr << "huh?\n";
					cerr << hss_list[hssI]->LeftEnd(0) << endl;
					cerr << hss_list[hssI]->RightEnd(0) << endl;
					cerr << iv_ptrs.back()->LeftEnd(seqI) << endl;
					cerr << iv_ptrs.back()->RightEnd(seqI) << endl;
				}
				while( ivI < iv_ptrs.size() && 
					(iv_ptrs[ivI]->LeftEnd(seqI) == NO_MATCH ||
					hss_list[hssI]->LeftEnd(0) > iv_ptrs[ivI]->RightEnd(seqI) ) )
					++ivI;
				if( ivI == iv_ptrs.size() )
				{
					cerr << "hssI fit!!\n";
					genome::breakHere();
				}
				// check for containment in seqJ
				if( iv_ptrs[ivI]->LeftEnd(seqJ) == NO_MATCH ||
					iv_ptrs[ivI]->RightEnd(seqJ) < hss_list[hssI]->LeftEnd(1) ||
					hss_list[hssI]->RightEnd(1) < iv_ptrs[ivI]->LeftEnd(seqJ) )
					continue;	// this hss falls to an invalid range in seqJ

				if( hss_list[hssI]->RightEnd(0) < iv_ptrs[ivI]->LeftEnd(seqI) )
				{
					cerr << "huh 2?\n";
					cerr << hss_list[hssI]->LeftEnd(0) << endl;
					cerr << hss_list[hssI]->RightEnd(0) << endl;
					cerr << iv_ptrs[ivI]->LeftEnd(seqI) << endl;
					cerr << iv_ptrs[ivI]->RightEnd(seqI) << endl;
					hssI++;
					continue;
				}

				vector< pair< size_t, size_t > >& cur_hss_cols = hss_cols[seqI][seqJ][iv_map[ivI]];

				gnSeqI left_col = iv_ptrs[ivI]->SeqPosToColumn( seqI, hss_list[hssI]->LeftEnd(0) );
				gnSeqI right_col = iv_ptrs[ivI]->SeqPosToColumn( seqI, hss_list[hssI]->RightEnd(0) );
				if(left_col > right_col && iv_ptrs[ivI]->Orientation(seqI) == AbstractMatch::reverse )
				{
					swap(left_col, right_col);	// must have been a revcomp seq
				}
				else if(left_col > right_col)
				{
					cerr << "bad cols\n";
					cerr << hss_list[hssI]->LeftEnd(0) << endl;
					cerr << hss_list[hssI]->RightEnd(0) << endl;
					cerr << iv_ptrs[ivI]->LeftEnd(seqI) << endl;
					cerr << iv_ptrs[ivI]->RightEnd(seqI) << endl;
					genome::breakHere();
				}

				if( left_col > 2000000000 || right_col > 2000000000 )
				{
					cerr << "huh 2?\n";
					cerr << hss_list[hssI]->LeftEnd(0) << endl;
					cerr << hss_list[hssI]->RightEnd(0) << endl;
					cerr << iv_ptrs[ivI]->LeftEnd(seqI) << endl;
					cerr << iv_ptrs[ivI]->RightEnd(seqI) << endl;
					genome::breakHere();
				}
				cur_hss_cols.push_back( make_pair( left_col, right_col ) );
			}
			for( size_t hssI = 0; hssI < hss_list.size(); ++hssI )
				hss_list[hssI]->Free();
		}
	}
}


class IntervalSeqManipulator
{
public:
	IntervalSeqManipulator( uint seq ) : m_seq(seq) {}
	gnSeqI LeftEnd(Interval& m) const{ return m.LeftEnd(m_seq); }
	gnSeqI Length(Interval& m) const{ return m.Length(m_seq); }
	void CropLeft(Interval& m, gnSeqI amount ) const{ m.CropLeft(amount, m_seq); }
	void CropRight(Interval& m, gnSeqI amount ) const{ m.CropRight(amount, m_seq); }
	template< typename ContainerType >
	void AddCopy(ContainerType& c, Interval& m) const{ c.push_back( m ); }
private:
	uint m_seq;
};


void detectBackboneBreakOnGenes( IntervalList& iv_list, backbone_list_t& bb_list, const HssDetector* detector, vector< CompactGappedAlignment<>* >& iv_orig_ptrs, vector< vector< gnSeqI > >& gene_bounds )
{
	uint seq_count = iv_list.seq_table.size();

	// indexed by seqI, seqJ, ivI, hssI (left col, right col)
	pairwise_genome_hss_t hss_cols(boost::extents[seq_count][seq_count][iv_list.size()]);

	// ugg.  need CompactGappedAlignment for its SeqPosToColumn
	vector< CompactGappedAlignment<>* > iv_ptrs(iv_list.size());
	for( size_t i = 0; i < iv_list.size(); ++i )
	{
		CompactGappedAlignment<> tmp_cga;
		iv_ptrs[i] = tmp_cga.Copy();
		new (iv_ptrs[i])CompactGappedAlignment<>( iv_list[i] );
	}

	iv_orig_ptrs = iv_ptrs;
	makeAllPairwiseGenomeHSSBreakOnGenes( iv_list, iv_ptrs, iv_orig_ptrs, hss_cols, detector, gene_bounds );

	// merge overlapping pairwise homology predictions into n-way predictions
	mergePairwiseHomologyPredictions( iv_orig_ptrs, hss_cols, bb_list );
}

int main( int argc, char* argv[] )
{
#if	WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	if( argc < 4 )
	{
		cerr << "bbBreakOnGenes <xmfa file> <min bb gap size> <bb output>\n";
		return -1;
	}
	string xmfa_fname( argv[1] );
	int min_bb_gap = atoi( argv[2] );
	string output_fname( argv[3] );

	ifstream xmfa_input( xmfa_fname.c_str() );
	if( !xmfa_input.is_open() ){
		cerr << "Error opening \"" << xmfa_fname << "\"" << endl;
		return -4;
	}
	ofstream bb_output( output_fname.c_str() );
	if( !bb_output.is_open() ){
		cerr << "Error opening \"" << output_fname << "\" for writing" << endl;
		return -6;
	}


	// read the alignment
	IntervalList iv_list;
	iv_list.ReadStandardAlignment( xmfa_input );
	LoadSequences(iv_list, &cout);
	vector< vector< gnSeqI > > gene_bounds( iv_list.seq_table.size() );

	if( argc - 4 == iv_list.seq_filename.size() )
	{
		cerr << "Reading gene coordinates from .ptt files\n";
		// read ptt files instead
		for( size_t aI = 0; aI < iv_list.seq_filename.size(); aI++ )
		{
			ifstream ptt_in( argv[aI+4] );
			string bub;
			getline( ptt_in, bub );
			getline( ptt_in, bub );
			getline( ptt_in, bub );
			while( getline( ptt_in, bub ) )
			{
				stringstream line_str(bub);
				string buf;
				getline( line_str, buf, '.' );
				int64 lend = atoi(buf.c_str());
				getline( line_str, buf, '.' );
				getline( line_str, buf );
				int64 rend = atoi(buf.c_str());
				gene_bounds[aI].push_back( lend -1);
				gene_bounds[aI].push_back( lend );
				gene_bounds[aI].push_back( rend );
				gene_bounds[aI].push_back( rend+1 );
	
			}
		}
	}else{

		// get gene boundary coordinates, break bb segs on genes...
		for( size_t genomeI = 0; genomeI < iv_list.seq_table.size(); genomeI++ )
		{
			for( size_t featureI = 0; featureI < iv_list.seq_table[genomeI]->getFeatureListLength(); ++featureI )
			{
				gnBaseFeature* feat = iv_list.seq_table[genomeI]->getFeature( featureI );
				string feat_name = feat->GetName();
				if( feat_name != "CDS" )
					continue;	// don't deal with other feature types (source, misc_RNA, etc)
				gnLocation loc = feat->GetLocation(0);
				if( loc.GetFirst() > loc.GetLast() || loc.GetFirst() == 0 || loc.GetLast() == 0 )
					continue;	// a problem parsing annotation?
				gene_bounds[genomeI].push_back( loc.GetFirst() );
				gene_bounds[genomeI].push_back( loc.GetLast() +1 );
			}
//			IntervalSeqManipulator ism(genomeI);
			std::sort( gene_bounds[genomeI].begin(), gene_bounds[genomeI].end() );
			cerr << "Found " << gene_bounds[genomeI].size() / 2 << " genes for " << iv_list.seq_filename[genomeI] << endl;
		}
	}

	// detect big gaps
	backbone_list_t bb_list;
	vector< CompactGappedAlignment<>* > iv_orig_ptrs;
	BigGapsDetector bgd( min_bb_gap );
	detectBackboneBreakOnGenes( iv_list, bb_list, &bgd, iv_orig_ptrs, gene_bounds );

	writeBackboneSeqCoordinates( bb_list, iv_list, bb_output );
	std::vector< bb_seqentry_t > bb_seq_list;
	bb_output.close();
	std::ifstream bbseq_input( output_fname.c_str() );
	readBackboneSeqFile( bbseq_input, bb_seq_list );

	// testing:  check whether any gene boundaries are violated
	gene_bounds[0].push_back(31337);	// test the test:
	gene_bounds[0].push_back(31333);	// insert some bogus gene bounds to make sure
	gene_bounds[0].push_back(31341);	// they get found and reported
	gene_bounds[0].push_back(31345);
	for( uint seqI = 0; seqI < iv_list.seq_table.size(); seqI++ )
	{
		cerr << "Checking seq " << seqI << " for errors\n";
		std::sort( gene_bounds[seqI].begin(), gene_bounds[seqI].end() );
		BbSeqEntrySorter bs(seqI);
		std::sort( bb_seq_list.begin(), bb_seq_list.end(), bs );
		size_t gI = 0;
		size_t bI = 0;
		cerr << gene_bounds[seqI].size() << " gene boundaries and " << bb_seq_list.size() << " bb segs\n";
		for( ; gI < gene_bounds[seqI].size() && bI < bb_seq_list.size(); gI++ )
		{
			cout << "checking " << bb_seq_list[bI][seqI].first << ", " <<bb_seq_list[bI][seqI].second << endl;  
			while( bI < bb_seq_list.size() && gene_bounds[seqI][gI] > abs(bb_seq_list[bI][seqI].second) )
				bI++;
			if( bI == bb_seq_list.size() )
				break;
			if(abs(bb_seq_list[bI][seqI].first) + 1 < gene_bounds[seqI][gI] && gene_bounds[seqI][gI] < abs(bb_seq_list[bI][seqI].second) - 1)
			{
				cerr << "segment " <<bb_seq_list[bI][seqI].first << ", " <<bb_seq_list[bI][seqI].second << " violates gene boundary " << gene_bounds[seqI][gI] << " in seq " << seqI << endl;  
			}else
				cout << "segment " <<bb_seq_list[bI][seqI].first << ", " <<bb_seq_list[bI][seqI].second << " is okay for " << gene_bounds[seqI][gI] << " in seq " << seqI << endl;  
		}
	}

//	mergeAdjacentSegments( bb_seq_list );
//	addUniqueSegments( bb_seq_list );
	bbseq_input.close();
	bb_output.open(output_fname.c_str());
	writeBackboneSeqFile( bb_output, bb_seq_list );

	return 0;
}

