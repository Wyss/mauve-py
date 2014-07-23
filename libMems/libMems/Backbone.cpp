/*******************************************************************************
 * $Id: Backbone.cpp,v 1.12 2004/04/19 23:11:19 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/ProgressiveAligner.h"
#include "libMems/Backbone.h"
#include "libMems/Islands.h"
#include "libMems/CompactGappedAlignment.h"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/undirected_dfs.hpp>

using namespace std;
using namespace genome;
namespace mems {


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

typedef pair< size_t, Interval* > iv_tracker_t;
class IvTrackerComp
{
public:
	IvTrackerComp( uint seq ) : ssc( seq ) {}
	bool operator()( const iv_tracker_t& a, const iv_tracker_t& b )
	{
		return ssc(a.second, b.second);
	}
private:
	SingleStartComparator<Interval> ssc;
};

const int LEFT_NEIGHBOR = -1;
const int RIGHT_NEIGHBOR = 1;
typedef vector< size_t > neighbor_t;

neighbor_t& getNeighbor( pair< neighbor_t, neighbor_t >& entry, int direction )
{
	if( direction == RIGHT_NEIGHBOR )
		return entry.first;
	else
		return entry.second;
}


void collapseCollinear( IntervalList& iv_list )
{
	if( iv_list.size() == 0 )
		return;	// nothing to see here, move along...
	const size_t seq_count = iv_list.seq_table.size();
	std::vector< Interval* > iv_ptrs(iv_list.size());
	size_t lilI = 0;
	for( size_t i = 0; i < iv_list.size(); ++i )
	{
		// ignore unaligned regions
		if( iv_list[i].Multiplicity() < 2 )
			continue;
		iv_ptrs[lilI++] = &iv_list[i];
	}
	iv_ptrs.resize(lilI);
	const size_t NEIGHBOR_UNKNOWN = (std::numeric_limits<size_t>::max)();
	neighbor_t lefties_tmp( seq_count, NEIGHBOR_UNKNOWN );
	pair< neighbor_t, neighbor_t > neighbor_pair( lefties_tmp, lefties_tmp );
	vector< pair< neighbor_t, neighbor_t > > neighbor_list( iv_ptrs.size(), neighbor_pair );
	vector< iv_tracker_t > iv_tracker( iv_ptrs.size() );
	for( size_t i = 0; i < iv_ptrs.size(); ++i )
	{
		iv_tracker[i] = make_pair( i, iv_ptrs[i] );
	}
	for( size_t seqI = 0; seqI < seq_count; ++seqI )
	{
		IvTrackerComp ivc( seqI );
		sort( iv_tracker.begin(), iv_tracker.end(), ivc );
		size_t prev_i = NEIGHBOR_UNKNOWN;
		size_t cur_i = NEIGHBOR_UNKNOWN;
		for( size_t i = 0; i < iv_tracker.size(); ++i )
		{
			if( iv_tracker[i].second->LeftEnd(seqI) == NO_MATCH )
				continue;
			if( cur_i != NEIGHBOR_UNKNOWN )
			{
				neighbor_list[cur_i].first[seqI] = prev_i;
				neighbor_list[cur_i].second[seqI] = iv_tracker[i].first;
			}
			prev_i = cur_i;
			cur_i = iv_tracker[i].first;
		}
		// get the last one
		if( cur_i != NEIGHBOR_UNKNOWN )
		{
			neighbor_list[cur_i].first[seqI] = prev_i;
			neighbor_list[cur_i].second[seqI] = NEIGHBOR_UNKNOWN;
		}
	}

	// now look for neighbor pair entries which can be merged
	for( int d = -1; d < 2; d+= 2 )	// iterate over both directions
	{
		size_t unknown_count = 0;
		for( size_t nI = 0; nI < neighbor_list.size(); ++nI )
		{
			size_t nayb = NEIGHBOR_UNKNOWN;
			size_t seqI = 0;
			bool parity = false;
			size_t ct = 0;
			for( ; seqI < seq_count; ++seqI )
			{
				if( iv_ptrs[nI]->Orientation(seqI) == AbstractMatch::undefined )
					continue;
				int orient = iv_ptrs[nI]->Orientation(seqI) == AbstractMatch::forward ? 1 : -1;

				if( nayb == NEIGHBOR_UNKNOWN )
				{
					nayb = getNeighbor( neighbor_list[nI], d * orient * -1 )[seqI];
					if( nayb != NEIGHBOR_UNKNOWN )
						parity = iv_ptrs[nI]->Orientation(seqI) == iv_ptrs[nayb]->Orientation(seqI);
				}
				else if( nayb != getNeighbor( neighbor_list[nI], d * orient * -1 )[seqI] )
					break;
				else if( parity != (iv_ptrs[nI]->Orientation(seqI) == iv_ptrs[nayb]->Orientation(seqI)) )
					break;
				if( nayb != NEIGHBOR_UNKNOWN )
					ct++;
			}
			if( seqI < seq_count || ct < iv_ptrs[nI]->Multiplicity() )
				continue;	// not collinear
			if( nayb == NEIGHBOR_UNKNOWN )
				continue;

			// merge nI and nayb
			uint fs = iv_ptrs[nI]->FirstStart();
			gnSeqI nI_lend_fs = iv_ptrs[nI]->LeftEnd(fs);
			gnSeqI nayb_lend_fs = iv_ptrs[nayb]->LeftEnd(fs);
			AbstractMatch::orientation o = iv_ptrs[nI]->Orientation(fs);
			vector< AbstractMatch* > nI_matches;
			iv_ptrs[nI]->StealMatches( nI_matches );
			vector< AbstractMatch* > nayb_matches;
			iv_ptrs[nayb]->StealMatches( nayb_matches );
			if( !parity )
			{
				std::reverse( nI_matches.begin(), nI_matches.end() );
				for( size_t i = 0; i < nI_matches.size(); ++i )
					nI_matches[i]->Invert();
				o = o == AbstractMatch::forward ? AbstractMatch::reverse : AbstractMatch::forward;
			}
			if( (o == AbstractMatch::forward && nI_lend_fs > nayb_lend_fs) ||
				(o == AbstractMatch::reverse && nI_lend_fs < nayb_lend_fs))
				nayb_matches.insert( nayb_matches.end(), nI_matches.begin(), nI_matches.end() );
			else
				nayb_matches.insert( nayb_matches.begin(), nI_matches.begin(), nI_matches.end() );

			iv_ptrs[nayb]->SetMatches( nayb_matches );

			// update all pointers to point to nayb
			seqI = 0;
			for( ; seqI < seq_count; ++seqI )
			{
				if( getNeighbor( neighbor_list[nI], -1 )[seqI] == NEIGHBOR_UNKNOWN &&
					getNeighbor( neighbor_list[nI], 1 )[seqI] == NEIGHBOR_UNKNOWN )
					continue;
				int orient = iv_ptrs[nayb]->Orientation(seqI) == AbstractMatch::forward ? 1 : -1;
				size_t other_nayb = getNeighbor( neighbor_list[nI], d * orient * (parity ? 1 : -1) )[seqI];
				if( other_nayb != NEIGHBOR_UNKNOWN )
				{
					if( getNeighbor( neighbor_list[other_nayb], 1 )[seqI] == nI )
						getNeighbor( neighbor_list[other_nayb], 1 )[seqI] = nayb;
					else if( getNeighbor( neighbor_list[other_nayb], -1 )[seqI] == nI )
						getNeighbor( neighbor_list[other_nayb], -1 )[seqI] = nayb;
					else
					{
						cerr << "serious programmer error\n";
						genome::breakHere();
					}
				}
				if( getNeighbor( neighbor_list[nayb], 1 )[seqI] == nI )
					getNeighbor( neighbor_list[nayb], 1 )[seqI] = other_nayb;
				else if( getNeighbor( neighbor_list[nayb], -1 )[seqI] == nI )
					getNeighbor( neighbor_list[nayb], -1 )[seqI] = other_nayb;
				else
				{
					cerr << "inexcusable programmer error\n";
					genome::breakHere();
				}
				neighbor_list[nI].first[seqI] = NEIGHBOR_UNKNOWN;
				neighbor_list[nI].second[seqI] = NEIGHBOR_UNKNOWN;
			}
		}
	}

	IntervalList new_list;
	new_list.seq_filename = iv_list.seq_filename;
	new_list.seq_table = iv_list.seq_table;
	new_list.resize( iv_ptrs.size() );
	size_t newI = 0;
	for( size_t ivI = 0; ivI < iv_ptrs.size(); ++ivI )
	{
		vector< AbstractMatch* > matches;
		iv_ptrs[ivI]->StealMatches( matches );
		if( matches.size() > 0 )
			new_list[newI++].SetMatches( matches );
	}
	new_list.resize(newI);
	swap( iv_list, new_list );
	addUnalignedRegions(iv_list);
}


void checkForAllGapColumns( IntervalList& iv_list )
{
	// debug: sanity check whether there are all gap columns
	for( size_t ivI = 0; ivI < iv_list.size(); ivI++ )
	{
		vector< string > aln;
		mems::GetAlignment( iv_list[ivI], iv_list.seq_table, aln );
		for( size_t colI = 0; colI < aln[0].size(); ++colI )
		{
			size_t rowI = 0;
			for( ; rowI < aln.size(); ++rowI )
				if( aln[rowI][colI] != '-' )
					break;
			if( rowI == aln.size() )
			{
				cerr << "ERROR!  IV " << ivI << " COLUMN " << colI << " IS ALL GAPS!\n";
			}
		}
	}
}



void translateToPairwiseGenomeHSS( const hss_array_t& hss_array, pairwise_genome_hss_t& hss_cols )
{
	uint seq_count = hss_array.shape()[0];
	uint iv_count = hss_array.shape()[2];
	hss_cols.resize( boost::extents[seq_count][seq_count][iv_count] );

	// make pairwise projections of intervals and find LCBs...
	for( size_t seqI = 0; seqI < seq_count; ++seqI )
	{
		for( size_t seqJ = seqI+1; seqJ < seq_count; ++seqJ )
		{
			for( size_t ivI = 0; ivI < iv_count; ++ivI )
			{
				const hss_list_t& cur_list = hss_array[seqI][seqJ][ivI];
				hss_cols[seqI][seqJ][ivI].resize( cur_list.size() );
				for( size_t hssI = 0; hssI < cur_list.size(); hssI++ )
				{
					hss_cols[seqI][seqJ][ivI][hssI].first = cur_list[hssI].left_col;
					hss_cols[seqI][seqJ][ivI][hssI].second = cur_list[hssI].right_col;
				}
			}
		}
	}
}


double computeGC( std::vector< gnSequence* >& seq_table )
{
	const uint8* tab = SortedMerList::BasicDNATable();
	size_t counts[4];
	for( int i = 0; i < 4; i++ )
		counts[i] = 0;
	for( size_t seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		std::string seq;
		seq_table[seqI]->ToString( seq );
		for( size_t cI = 0; cI < seq.size(); cI++  )
			counts[ tab[ seq[cI] ] ]++;
	}
	return double(counts[1]+counts[2]) / double(counts[1]+counts[2] + counts[0]+counts[3]);
}


void makeAllPairwiseGenomeHSS( IntervalList& iv_list, vector< CompactGappedAlignment<>* >& iv_ptrs, vector< CompactGappedAlignment<>* >& iv_orig_ptrs, pairwise_genome_hss_t& hss_cols, const HssDetector* detector )
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

			// break up these alignments on contig and chromosome boundaries
			for(int ssI=0; ssI<2; ssI++){
				vector<gnSeqI> contig_bounds;
				for( size_t cI=0; cI < pair_ivs.seq_table[ssI]->contigListSize(); cI++ ){
					contig_bounds.push_back(pair_ivs.seq_table[ssI]->contigLength(cI));
					if( cI > 0 )
						contig_bounds[cI] += contig_bounds[cI-1];
				}
				GenericMatchSeqManipulator< CompactGappedAlignment<> > gmsm(ssI);
				applyBreakpoints(contig_bounds, pair_cgas, gmsm);
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
			// and again on seqJ
			getBpList( iv_ptrs, seqJ, bp_list );
			GenericMatchSeqManipulator< CompactGappedAlignment<> > gmsm1(1);
			SingleStartComparator< CompactGappedAlignment<> > ssc1(1);
			std::sort(hss_list.begin(), hss_list.end(), ssc1 );
			applyBreakpoints( bp_list, hss_list, gmsm1 );

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

void mergePairwiseHomologyPredictions( 	vector< CompactGappedAlignment<>* >& iv_orig_ptrs, pairwise_genome_hss_t& hss_cols, vector< vector< ULA* > >& ula_list )
{
	uint seq_count = hss_cols.shape()[0];
	uint iv_count = hss_cols.shape()[2];
	//
	// FINALLY!  ready to merge.  how to do it?
	// make an empty list of UngappedLocalAlignments
	// start with the first seq and create a ULA for every col
	// range.  Then continue to the second seq, and when
	// a col range overlaps a pre-existing ULA, create a new ULA
	// for the intersected region and a smaller ULA for the non-intersected region
	ula_list.resize( iv_count );
	for( size_t ivI = 0; ivI < iv_count; ++ivI )
	{
		vector< ULA* >& iv_ulas = ula_list[ivI];
		for( size_t seqI = 0; seqI < seq_count; ++seqI )
		{
			for( size_t seqJ = seqI+1; seqJ < seq_count; ++seqJ )
			{
				vector< pair< size_t, size_t > >& cur_hss_cols = hss_cols[seqI][seqJ][ivI];
				vector< ULA* > cur_ulas( cur_hss_cols.size() );
				ULA tmp_ula(seq_count);
				for( size_t hssI = 0; hssI < cur_hss_cols.size(); ++hssI )
				{
					cur_ulas[hssI] = tmp_ula.Copy();
					cur_ulas[hssI]->SetStart(seqI, cur_hss_cols[hssI].first+1);
					cur_ulas[hssI]->SetStart(seqJ, cur_hss_cols[hssI].first+1);
					cur_ulas[hssI]->SetLength( cur_hss_cols[hssI].second - cur_hss_cols[hssI].first + 1 );
				}

				vector< gnSeqI > iv_bp_list;
				vector< gnSeqI > cur_bp_list;
				SingleStartComparator<ULA> ulacompI(seqI);
				std::sort( iv_ulas.begin(), iv_ulas.end(), ulacompI );
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompI );
				getBpList( iv_ulas, seqI, iv_bp_list );
				getBpList( cur_ulas, seqI, cur_bp_list );
				GenericMatchSeqManipulator< ULA > gmsm(seqI);
				applyBreakpoints( iv_bp_list, cur_ulas, gmsm );
				applyBreakpoints( cur_bp_list, iv_ulas, gmsm );

				SingleStartComparator<ULA> ulacompJ(seqJ);
				std::sort( iv_ulas.begin(), iv_ulas.end(), ulacompJ );
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompJ );
				getBpList( iv_ulas, seqJ, iv_bp_list );
				getBpList( cur_ulas, seqJ, cur_bp_list );
				GenericMatchSeqManipulator< ULA > gmsmJ(seqJ);
				applyBreakpoints( iv_bp_list, cur_ulas, gmsmJ );
				applyBreakpoints( cur_bp_list, iv_ulas, gmsmJ );

				// do seqI a second time to propagate any breakpoints introduced by seqJ
				std::sort( iv_ulas.begin(), iv_ulas.end(), ulacompI );
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompI );
				getBpList( iv_ulas, seqI, iv_bp_list );
				getBpList( cur_ulas, seqI, cur_bp_list );
				applyBreakpoints( iv_bp_list, cur_ulas, gmsm );
				applyBreakpoints( cur_bp_list, iv_ulas, gmsm );

				std::sort( iv_ulas.begin(), iv_ulas.end(), ulacompI );
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompI );
				// now that cur_ulas and iv_ulas are all broken up according to each other's boundaries
				// we can simply scan along and add
				size_t iv_ulas_size = iv_ulas.size();
				size_t ivuI = 0;
				size_t curuI = 0;
				vector< ULA* > added_to( cur_ulas.size(), NULL );	// this tracks which of iv_ulas a cur_ula was added to
				vector< ULA* > to_delete;
				while( ivuI < iv_ulas_size && curuI < cur_ulas.size() )
				{
					if( iv_ulas[ivuI]->LeftEnd(seqI) == cur_ulas[curuI]->LeftEnd(seqI) )
					{
						if( added_to[curuI] == iv_ulas[ivuI] )
						{
							// do nothing
						}else if( added_to[curuI] == NULL )
						{
							iv_ulas[ivuI]->SetLeftEnd(seqJ, cur_ulas[curuI]->LeftEnd(seqJ));
							added_to[curuI] = iv_ulas[ivuI];
						}else{
							ULA* merge = added_to[curuI];
							for( size_t seqK = 0; seqK < seq_count; ++seqK )
							{
								if( merge->Start(seqK) == NO_MATCH )
									continue;
								iv_ulas[ivuI]->SetStart( seqK, merge->Start(seqK) );
							}
							to_delete.push_back( merge );
						}
						ivuI++;
					}else if( iv_ulas[ivuI]->LeftEnd(seqI) < cur_ulas[curuI]->LeftEnd(seqI) )
					{
						ivuI++;
					}else
						curuI++;
				}

				// delete to_delete...
				std::sort( to_delete.begin(), to_delete.end() );
				vector< ULA* >::iterator last = std::unique( to_delete.begin(), to_delete.end() );
				to_delete.erase( last, to_delete.end() );
				vector< ULA* > new_iv_ulas( iv_ulas.size() - to_delete.size() );
				std::sort( iv_ulas.begin(), iv_ulas.end() );
				std::set_difference( iv_ulas.begin(), iv_ulas.end(), to_delete.begin(), to_delete.end(), new_iv_ulas.begin() );
				swap( iv_ulas, new_iv_ulas );
				for( size_t delI = 0; delI < to_delete.size(); ++delI )
					to_delete[delI]->Free();

				vector< ULA* > orig_ula_order = cur_ulas;
				// now do something similar for seqJ
				std::sort( iv_ulas.begin(), iv_ulas.end(), ulacompJ );
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompJ );

				vector< size_t > added_map;
				createMap( cur_ulas, orig_ula_order, added_map );

				ivuI = 0;
				curuI = 0;
				to_delete.clear();
				while( ivuI < iv_ulas_size && curuI < cur_ulas.size() )
				{
					if( iv_ulas[ivuI]->LeftEnd(seqJ) == cur_ulas[curuI]->LeftEnd(seqJ) )
					{
						if( added_to[added_map[curuI]] == iv_ulas[ivuI] )
						{
							// do nothing
						}else if( added_to[added_map[curuI]] == NULL )
						{
							iv_ulas[ivuI]->SetLeftEnd(seqI, cur_ulas[curuI]->LeftEnd(seqI));
							added_to[added_map[curuI]] = iv_ulas[ivuI];
						}else{
							ULA* merge = added_to[added_map[curuI]];
							for( size_t seqK = 0; seqK < seq_count; ++seqK )
							{
								if( merge->Start(seqK) == NO_MATCH )
									continue;
								iv_ulas[ivuI]->SetStart( seqK, merge->Start(seqK) );
							}
							to_delete.push_back( merge );
						}
						ivuI++;
					}else if( iv_ulas[ivuI]->LeftEnd(seqJ) < cur_ulas[curuI]->LeftEnd(seqJ) )
					{
						ivuI++;
					}else
					{
						curuI++;
					}
				}

				// anything with a null added_to entry needs to be added to iv_ulas
				// everything else needs to get freed
				std::sort( cur_ulas.begin(), cur_ulas.end(), ulacompI );
				for( curuI = 0; curuI < cur_ulas.size(); ++curuI )
				{
					if( added_to[curuI] == NULL )
						iv_ulas.push_back( cur_ulas[curuI] );
					else
						cur_ulas[curuI]->Free();
				}
				// delete to_delete...
				std::sort( to_delete.begin(), to_delete.end() );
				last = std::unique( to_delete.begin(), to_delete.end() );
				to_delete.erase( last, to_delete.end() );
				new_iv_ulas = vector< ULA* >( iv_ulas.size() - to_delete.size() );
				std::sort( iv_ulas.begin(), iv_ulas.end() );
				std::set_difference( iv_ulas.begin(), iv_ulas.end(), to_delete.begin(), to_delete.end(), new_iv_ulas.begin() );
				swap( iv_ulas, new_iv_ulas );
				for( size_t delI = 0; delI < to_delete.size(); ++delI )
					to_delete[delI]->Free();
			}
		}
	}

	// Eliminate segments that have no representation in a genome
	for( size_t ivI = 0; ivI < ula_list.size(); ++ivI )
	{
		for( size_t mI = 0; mI < ula_list[ivI].size(); ++mI )
		{
			size_t seqI = ula_list[ivI][mI]->FirstStart();
			std::vector<gnSeqI> l_pos;
			std::vector<bool> l_column;
			std::vector<gnSeqI> r_pos;
			std::vector<bool> r_column;
			gnSeqI left_col = ula_list[ivI][mI]->LeftEnd(seqI)-1;
			gnSeqI right_col = ula_list[ivI][mI]->RightEnd(seqI)-1;
			iv_orig_ptrs[ivI]->GetColumn(left_col, l_pos, l_column);
			iv_orig_ptrs[ivI]->GetColumn(right_col, r_pos, r_column);
			for( ; seqI < ula_list[ivI][mI]->SeqCount(); ++seqI )
			{
				if( ula_list[ivI][mI]->LeftEnd(seqI) == NO_MATCH )
					continue;
				if( l_pos[seqI] == r_pos[seqI] && !l_column[seqI] && !r_column[seqI] )
					ula_list[ivI][mI]->SetStart(seqI, NO_MATCH);	// no match in this col
			}
			if( ula_list[ivI][mI]->Multiplicity() < 2 )
			{
				ula_list[ivI][mI]->Free();
				ula_list[ivI][mI] = NULL;
			}
		}
		// clean out any NULL ptrs
		std::vector< ULA* >::iterator last = std::remove( ula_list[ivI].begin(), ula_list[ivI].end(), (ULA*)NULL );
		ula_list[ivI].erase( last, ula_list[ivI].end() );
	}
}

void unalignIslands( IntervalList& iv_list, vector< CompactGappedAlignment<>* >& iv_orig_ptrs, vector< vector< ULA* > >& ula_list )
{
	uint seq_count = iv_list.seq_table.size();
	// unalign regions in the iv list that aren't contained in backbone
	for( size_t ivI = 0; ivI < ula_list.size(); ++ivI )
	{
		vector< AbstractMatch* > new_matches(ula_list[ivI].size());
		for( size_t mI = 0; mI < ula_list[ivI].size(); ++mI )
		{
			size_t seqI = ula_list[ivI][mI]->FirstStart();
			gnSeqI left_col = ula_list[ivI][mI]->LeftEnd(seqI)-1;
			CompactGappedAlignment<> tmp_cga;
			CompactGappedAlignment<>* new_cga = tmp_cga.Copy();
			iv_orig_ptrs[ivI]->copyRange(*new_cga, left_col, ula_list[ivI][mI]->Length(seqI));
			for( seqI = 0; seqI < ula_list[ivI][mI]->SeqCount(); ++seqI )
			{
				if( ula_list[ivI][mI]->LeftEnd(seqI) == NO_MATCH )
					new_cga->SetLeftEnd(seqI, NO_MATCH);
			}
			new_cga->CondenseGapColumns();
			new_matches[mI] = new_cga;
		}
		if( new_matches.size() > 0 )
		{

			vector< vector< AbstractMatch* > > disjoint_subsets;
			{
				// split into multiple intervals if some sequences are completely unaligned
				// use a union-find structure to quickly figure out how many subgroups there are
				vector< uint > seq_map( seq_count );
				for( size_t sI = 0; sI < seq_map.size(); ++sI )
					seq_map[sI] = sI;
				for( size_t mI = 0; mI < new_matches.size(); ++mI )
				{
					uint sI = new_matches[mI]->FirstStart();
					uint map_to = seq_map[sI];
					while( map_to != seq_map[map_to] )
						map_to = seq_map[map_to];
					seq_map[sI] = map_to;
					for( ++sI; sI < seq_count; ++sI )
					{
						if( new_matches[mI]->LeftEnd(sI) == NO_MATCH )
							continue;
						uint map_from = seq_map[sI];
						while( map_from != seq_map[map_from] )
							map_from = seq_map[map_from];
						seq_map[map_from] = map_to;
					}
				}
				vector< vector< AbstractMatch* > > mapped_lists( seq_count );
				for( size_t mI = 0; mI < new_matches.size(); ++mI )
				{
					uint sI = new_matches[mI]->FirstStart();
					uint map_to = seq_map[sI];
					while( map_to != seq_map[map_to] )
						map_to = seq_map[map_to];
					mapped_lists[map_to].push_back( new_matches[mI] );
				}
				for( uint sI = 0; sI < seq_count; ++sI )
				{
					if( mapped_lists[sI].size() > 0 )
						disjoint_subsets.push_back( mapped_lists[sI] );
				}
			}

			for( size_t dI = 0; dI < disjoint_subsets.size(); ++dI )
			{
				vector< AbstractMatch* >& cur_d_matches = disjoint_subsets[dI];
				vector< AbstractMatch* > orig_order = cur_d_matches;
				// need to sort these.  use boost's topological sort.
				vector< size_t > id_map;
				typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, boost::property<boost::vertex_color_t, boost::default_color_type> > Graph;
				typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
				typedef std::pair< int, int > Pair;
				vector< Pair > edges;
				for( size_t seqI = 0; seqI < seq_count; ++seqI )
				{
					SingleStartComparator<AbstractMatch> ssc(seqI);
					std::sort( cur_d_matches.begin(), cur_d_matches.end(), ssc );
					createMap( cur_d_matches, orig_order, id_map );
					int prev = -1;
					int first = -1;
					bool reverse = false;
					for( int mI = 0; mI < cur_d_matches.size(); ++mI )
					{
						if( cur_d_matches[mI]->LeftEnd(seqI) == NO_MATCH )
							continue;
						if( prev != -1 )
						{
							Pair edge( id_map[prev], id_map[mI] );
							if( reverse )
								swap( edge.first, edge.second );
							edges.push_back(edge);
						}else
						{
							reverse = cur_d_matches[mI]->Start(seqI) < 0;
							first = mI;
						}
						prev = mI;
					}
					if( prev != -1 && !reverse )
						edges.push_back( Pair( id_map[prev], cur_d_matches.size() ) );
					else if( prev != -1 && reverse )
						edges.push_back( Pair( id_map[first], cur_d_matches.size() ) );
				}
				std::sort( edges.begin(), edges.end() );
				vector< Pair >::iterator ee_iter = std::unique( edges.begin(), edges.end() );
				edges.erase( ee_iter, edges.end() );
				Pair* edge_array = new Pair[edges.size()];
				for( size_t eI = 0; eI < edges.size(); ++eI )
					edge_array[eI] = edges[eI];
				typedef boost::graph_traits<Graph>::vertices_size_type v_size_t;
				Graph G(edge_array, edge_array + edges.size(), v_size_t(edges.size()));
				typedef std::vector< Vertex > container;
				container c;
				topological_sort(G, std::back_inserter(c));
				cur_d_matches.clear();
				for ( container::reverse_iterator ii=c.rbegin(); ii!=c.rend(); ++ii)
				{
					if( *ii < orig_order.size() )
						cur_d_matches.push_back( orig_order[ *ii ] );
				}
				if( dI == 0 )
					iv_list[ivI].SetMatches(cur_d_matches);
				else
				{
					Interval new_iv( cur_d_matches.begin(), cur_d_matches.end() );
					iv_list.push_back(new_iv);
				}
				delete[] edge_array;
			}
		}
		else
		{
			iv_orig_ptrs[ivI]->Free();
			iv_orig_ptrs[ivI] = NULL;
		}
	}


	// update iv_list to match the filtered iv_orig_ptrs
	size_t givI = 0;
	for( size_t iI = 0; iI < iv_orig_ptrs.size(); ++iI )
	{
		if( iv_orig_ptrs[iI] != NULL )
		{
			swap( iv_list[givI], iv_list[iI] );
			iv_orig_ptrs[iI]->Free();	// done with the CompactGappedAlignments
			iv_orig_ptrs[iI] = NULL;
			givI++;
		}
	}
	// pick up any intervals that were split in half
	for( size_t iI = iv_orig_ptrs.size(); iI < iv_list.size(); ++iI )
		swap( iv_list[givI++], iv_list[iI] );
	iv_list.erase( iv_list.begin()+givI, iv_list.end() );

	// collapse any intervals that are trivially collinear
	collapseCollinear( iv_list );
}

void createBackboneList( const IntervalList& iv_list, backbone_list_t& ula_list ) 
{
	ula_list.resize( iv_list.size() );
	for( size_t ivI = 0; ivI < iv_list.size(); ++ivI )
	{
		if( iv_list[ivI].Multiplicity() < 2 )
			continue;
		const vector< AbstractMatch* >& matches = iv_list[ivI].GetMatches();
		int64 right_col = 0;
		int64 left_col = 0;
		for( size_t mI = 0; mI < matches.size(); ++mI )
		{
			left_col = right_col;
			right_col += matches[mI]->AlignmentLength();
			if( matches[mI]->Multiplicity() < 2 )
				continue;
			ULA tmp_ula(matches[mI]->SeqCount());
			ULA* mula = tmp_ula.Copy();
			for( size_t seqI = 0; seqI < matches[mI]->SeqCount(); ++seqI )
				if( matches[mI]->LeftEnd(seqI) != NO_MATCH )
					mula->SetLeftEnd( seqI, left_col+1 );
			mula->SetLength( right_col - left_col );
			ula_list[ivI].push_back(mula);
		}
		// merge neighbors that cover identical match components
		for( size_t ulaI = 1; ulaI < ula_list[ivI].size(); ulaI++ )
		{
			size_t seqI = 0;
			for( ; seqI < ula_list[ivI][ulaI]->SeqCount(); ++seqI )
			{
				int64 s1 = ula_list[ivI][ulaI-1]->Start(seqI);
				int64 s2 = ula_list[ivI][ulaI]->Start(seqI);
				if( s1 == mems::NO_MATCH && s2 == mems::NO_MATCH )
					continue;
				if( s1 == mems::NO_MATCH && s2 != mems::NO_MATCH )
					break;
				if( s1 != mems::NO_MATCH && s2 == mems::NO_MATCH )
					break;
				int64 r1 = ula_list[ivI][ulaI-1]->RightEnd(seqI);
				if( r1 + 1 != s2 )
					break;	// must be adjacent to each other
			}
			if( seqI == ula_list[ivI][ulaI]->SeqCount() )
			{
				// ulaI-1 needs to be swallowed up by ulaI
				ula_list[ivI][ulaI]->ExtendStart( ula_list[ivI][ulaI-1]->Length() );
				ula_list[ivI][ulaI-1]->SetLength(0);
			}
		}
		// get rid of matches that were swallowed up
		vector< ULA* > condensed_list;
		for( size_t ulaI = 0; ulaI < ula_list[ivI].size(); ulaI++ )
		{
			if( ula_list[ivI][ulaI]->Length() > 0 )
				condensed_list.push_back(ula_list[ivI][ulaI]);
			else
				ula_list[ivI][ulaI]->Free();
		}
		swap( ula_list[ivI], condensed_list );
	}
}

void detectAndApplyBackbone( AbstractMatch* m, vector< gnSequence* >& seq_table, CompactGappedAlignment<>*& result, backbone_list_t& bb_list, const Params& hmm_params, boolean left_homologous, boolean right_homologous )
{
	vector< AbstractMatch* > mlist( 1, m );
	uint seq_count = seq_table.size();

	// indexed by seqI, seqJ, ivI, hssI (left col, right col)
	pairwise_genome_hss_t hss_cols(boost::extents[seq_count][seq_count][1]);

	// ugg.  need CompactGappedAlignment for its SeqPosToColumn
	vector< CompactGappedAlignment<>* > iv_ptrs(1);
	CompactGappedAlignment<> tmp_cga;
	iv_ptrs[0] = tmp_cga.Copy();
	new (iv_ptrs[0])CompactGappedAlignment<>( *m );	// this will be freed when unalignIslands() gets called

	vector< CompactGappedAlignment<>* > iv_orig_ptrs(iv_ptrs);
	hss_array_t island_array, hss_array;

	findHssHomologyHMM( mlist, seq_table, island_array, hmm_params, left_homologous, right_homologous );
	translateToPairwiseGenomeHSS( island_array, hss_cols );

	// merge overlapping pairwise homology predictions into n-way predictions
	backbone_list_t ula_list;
	mergePairwiseHomologyPredictions( iv_orig_ptrs, hss_cols, ula_list );

	// unalignIslands wants an IntervalList
	IntervalList iv_list;
	iv_list.seq_table = seq_table;
	iv_list.resize(1);
	vector<AbstractMatch*> asdf(1, iv_orig_ptrs.front()->Copy() );
	iv_list[0].SetMatches( asdf );
	// unalign regions found to be non-homologous
	unalignIslands( iv_list, iv_orig_ptrs, ula_list );

	// free all ULAs and reconstruct them from the new alignment column coordinates
	for( size_t ulaI = 0; ulaI < ula_list.size(); ++ulaI )
		for( size_t i = 0; i < ula_list[ulaI].size(); ++i )
			ula_list[ulaI][i]->Free();
	ula_list.clear();


	createBackboneList( iv_list, ula_list );

	iv_orig_ptrs.clear();

	bb_list.clear();
	bb_list = ula_list;

	result = tmp_cga.Copy();
	if( iv_list.size() > 0 )
		new (result)CompactGappedAlignment<>( iv_list[0] );
}



void applyBackbone( IntervalList& iv_list, vector< CompactGappedAlignment<>* >& iv_orig_ptrs, backbone_list_t& bb_list )
{
	// unalign regions found to be non-homologous
	unalignIslands( iv_list, iv_orig_ptrs, bb_list );

	// need to add in all the unaligned regions so the viewer doesn't throw a fit
	addUnalignedRegions( iv_list );

	// free all ULAs and reconstruct them from the new alignment column coordinates
	for( size_t ulaI = 0; ulaI < bb_list.size(); ++ulaI )
		for( size_t i = 0; i < bb_list[ulaI].size(); ++i )
			bb_list[ulaI][i]->Free();
	bb_list.clear();

	createBackboneList( iv_list, bb_list );
}

void detectBackbone( IntervalList& iv_list, backbone_list_t& bb_list, const HssDetector* detector, vector< CompactGappedAlignment<>* >& iv_orig_ptrs )
{
	// collapse any intervals that are trivially collinear
	collapseCollinear( iv_list );

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
	makeAllPairwiseGenomeHSS( iv_list, iv_ptrs, iv_orig_ptrs, hss_cols, detector );

	// merge overlapping pairwise homology predictions into n-way predictions
	mergePairwiseHomologyPredictions( iv_orig_ptrs, hss_cols, bb_list );
}


// add unique segments of some minimum length
// FIXME: does not add begin and end segments!
void addUniqueSegments( std::vector< bb_seqentry_t >& bb_seq_list, size_t min_length )
{
	if( bb_seq_list.size() == 0 )
		return;
	vector< bb_seqentry_t > new_segs;
	uint seq_count = bb_seq_list[0].size();
	// now mark segs that are too close to each other to be considered independent
	for( size_t sI = 0; sI < seq_count; sI++ )
	{
		BbSeqEntrySorter bbs(sI);
		std::sort( bb_seq_list.begin(), bb_seq_list.end(), bbs );
		for( size_t bbI = 1; bbI < bb_seq_list.size(); bbI++ )
		{
			if( bb_seq_list[bbI][sI].first == 0 )
				continue;
			int64 diff = genome::absolut(bb_seq_list[bbI][sI].first) - genome::absolut(bb_seq_list[bbI-1][sI].second); 
			if( genome::absolut(diff) > min_length )
			{
				bb_seqentry_t newb( seq_count, make_pair( 0,0 ) );
				newb[sI].first = genome::absolut(bb_seq_list[bbI-1][sI].second) + 1;
				newb[sI].second = genome::absolut(bb_seq_list[bbI][sI].first) - 1;
				new_segs.push_back( newb );
			}
		}
	}
	bb_seq_list.insert( bb_seq_list.end(), new_segs.begin(), new_segs.end() );
}


void mergeAdjacentSegments( std::vector< bb_seqentry_t >& bb_seq_list )
{
	if( bb_seq_list.size() == 0 )
		return;
	uint seq_count = bb_seq_list[0].size();
	// now mark segs that are too close to each other to be considered independent
	for( size_t sI = 0; sI < seq_count; sI++ )
	{
		BbSeqEntrySorter bbs(sI);
		std::sort( bb_seq_list.begin(), bb_seq_list.end(), bbs );
		bitset_t merged;
		merged.resize( bb_seq_list.size() );
		for( size_t bbI = 1; bbI < bb_seq_list.size(); bbI++ )
		{
			if( bb_seq_list[bbI][sI].first == 0 )
				continue;
			size_t j = 0;
			for( ; j < seq_count; j++ )
			{
				if( bb_seq_list[bbI][j].first == 0 ^ bb_seq_list[bbI-1][j].first == 0)
					break;
				if( bb_seq_list[bbI][j].first == 0)
					continue;
				int64 diff = 0;
				if( bb_seq_list[bbI][j].first > 0 )
					diff = bb_seq_list[bbI][j].first - bb_seq_list[bbI-1][j].second; 
				else
					diff = bb_seq_list[bbI][j].second - bb_seq_list[bbI-1][j].first;
				if( diff != 1 )
					break;
			}
			if(j == seq_count)
			{	// they can be merged!
				merged.set(bbI-1);
				for( j = 0; j < seq_count; j++ )
					if( bb_seq_list[bbI][j].first > 0 )
						bb_seq_list[bbI][j].first = bb_seq_list[bbI-1][j].first;
					else
						bb_seq_list[bbI][j].second = bb_seq_list[bbI-1][j].second;
			}
		}
		// remove merged entries
		size_t cur = 0;
		for( size_t bbI = 0; bbI < bb_seq_list.size(); bbI++ )
			if( !merged.test( bbI ) )
				swap( bb_seq_list[cur++], bb_seq_list[bbI] );
		bb_seq_list.erase( bb_seq_list.begin() + cur, bb_seq_list.end() );
	}
}


void detectBackbone( IntervalList& iv_list, backbone_list_t& bb_list, const HssDetector* detector )
{
	vector< CompactGappedAlignment<>* > iv_orig_ptrs;
	detectBackbone( iv_list, bb_list, detector, iv_orig_ptrs );
	// FIXME: clean up iv_orig_ptrs
}

void detectAndApplyBackbone( IntervalList& iv_list, backbone_list_t& bb_list, const Params& hmm_params )
{
	HomologyHmmDetector* hmm_detector = new HomologyHmmDetector( hmm_params, true, true );
	vector< CompactGappedAlignment<>* > iv_orig_ptrs;
	detectBackbone( iv_list, bb_list, hmm_detector, iv_orig_ptrs );
	applyBackbone( iv_list, iv_orig_ptrs, bb_list );
	delete hmm_detector;
}


void writeBackboneColumns( ostream& bb_out, backbone_list_t& bb_list )
{
	//
	// At last! write out the backbone list
	//
	for( size_t ivI = 0; ivI < bb_list.size(); ++ivI )
	{
		for( size_t mI = 0; mI < bb_list[ivI].size(); ++mI )
		{
			size_t seqI = bb_list[ivI][mI]->FirstStart();
			bb_out << ivI << '\t' << bb_list[ivI][mI]->LeftEnd(seqI) << '\t' << bb_list[ivI][mI]->Length();
			for( ; seqI < bb_list[ivI][mI]->SeqCount(); ++seqI )
			{
				if( bb_list[ivI][mI]->LeftEnd(seqI) == NO_MATCH )
					continue;
				bb_out << '\t' << seqI;
			}
			bb_out << endl;
		}
	}
}

void writeBackboneSeqCoordinates( backbone_list_t& bb_list, IntervalList& iv_list, ostream& bb_out )
{
	if( bb_list.size() == 0 )
		return;
	// find seq_count
	uint seq_count = 0;
	for( size_t bbI = 0; bbI < bb_list.size(); ++bbI )
		if( bb_list[bbI].size() > 0 )
		{
			seq_count = bb_list[bbI].front()->SeqCount();
			break;
		}

	// different format -- use real sequence coordinates...
	// print a header line first
	for( size_t seqI = 0; seqI < seq_count; ++seqI )
	{
		if( seqI > 0 )
			bb_out << '\t';
		bb_out << "seq_" << seqI << "_leftend\t";
		bb_out << "seq_" << seqI << "_rightend";
	}
	bb_out << endl;
	for( size_t ivI = 0; ivI < bb_list.size(); ++ivI )
	{
		// there seems to be a bug in the backbone creation code that causes the CGA that gets
		// stuffed into the interval to have the wrong coordinates internally, while the interval
		// maintains the correct coordinates.  work around it by converting the whole interval to a cga
		CompactGappedAlignment<> iv_cga( iv_list[ivI] );
		for( size_t mI = 0; mI < bb_list[ivI].size(); ++mI )
		{
			uint fs = bb_list[ivI][mI]->FirstStart();
			// get the sequence positions out of the alignment
			vector< gnSeqI > left_pos;
			vector< bool > left_cols;
			iv_cga.GetColumn( bb_list[ivI][mI]->LeftEnd(fs)-1, left_pos, left_cols );
			vector< gnSeqI > right_pos;
			vector< bool > right_cols;
			iv_cga.GetColumn( bb_list[ivI][mI]->RightEnd(fs)-1, right_pos, right_cols );
			for( size_t seqI = 0; seqI < bb_list[ivI][mI]->SeqCount(); ++seqI )
			{
				if( seqI > 0 )
					bb_out << '\t';
				if( bb_list[ivI][mI]->LeftEnd(seqI) == NO_MATCH )
				{
					bb_out << "0\t0";
					continue;
				}else{
					int64 leftI = left_pos[seqI];
					int64 rightI = right_pos[seqI];
					if( iv_cga.Orientation(seqI) == AbstractMatch::forward && leftI != 0 && !left_cols[seqI] )
						leftI++;
					if( iv_cga.Orientation(seqI) == AbstractMatch::reverse && rightI != 0 && !right_cols[seqI] )
						rightI++;
					if( iv_cga.Orientation(seqI) == AbstractMatch::reverse )
					{
						swap( leftI, rightI );	// must be reverse complement
					}
					if( rightI + 1 == leftI )
					{
						bb_out << "0\t0";
						continue;
					}
					if( leftI > rightI )
					{
						cerr << "oh crahpey!\n";
						cerr << "leftI: " << leftI << endl;
						cerr << "rightI: " << rightI << endl;
						cerr << "seqI: " << seqI << endl;
						cerr << "ivI: " << ivI << endl;
					}
					if( leftI == 0 )
						leftI = iv_cga.LeftEnd(seqI);
					if( rightI == iv_cga.RightEnd(seqI)+1 )
						rightI--;
					if( iv_cga.Orientation(seqI) == AbstractMatch::reverse )
					{
						leftI *= -1;
						rightI *= -1;
					}
					bb_out << leftI << '\t' << rightI;
				}
			}
			bb_out << endl;
		}
	}
}


}  // namespace mems

