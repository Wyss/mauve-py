#include <fstream>
#include <sstream>
#include <iomanip>
#include "libMems/Backbone.h"
#include "libGenome/gnFeature.h"
#include "libGenome/gnBaseQualifier.h"
#include "libMems/IntervalList.h"
#include "libMems/AbstractMatch.h"
#include "libMems/MatchList.h"
#include "libMems/PhyloTree.h"
#include "libMems/ProgressiveAligner.h"
#include <boost/algorithm/string/erase.hpp>
#include <boost/tuple/tuple.hpp>

using namespace std;
using namespace genome;
using namespace mems;

// important constants that affect inference
const uint SHORT_SEGMENT = 5;	// when considering overlaps to genes, ignore overlaps less than this amount
const uint DISCARD_SEGMENT = 20;	// do not consider segments shorter than this amount
const double ALTERNALOG_MIN_SIZE = 15.0;


class BbSeqComp
{
public:
	BbSeqComp( uint seq ) : m_seq( seq ) {}
	bool operator()( const bb_seqentry_t* a, const bb_seqentry_t* b )
	{
		return genome::absolut( (*a)[m_seq].first ) < genome::absolut( (*b)[m_seq].first );
	}
private:
	uint m_seq;
};

template< typename PtrVector >
void createMap( const PtrVector& mv_from, const PtrVector& mv_to, vector< size_t >& map )
{
	typedef typename PtrVector::value_type PtrType;
	vector< pair< PtrType, size_t > > m1(mv_from.size());
	vector< pair< PtrType, size_t > > m2(mv_to.size());
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

size_t getCDScount( gnSequence* anno_seq )
{
	size_t count = 0;
	for( size_t featureI = 0; featureI < anno_seq->getFeatureListLength(); ++featureI )
	{
		gnBaseFeature* feat = anno_seq->getFeature( featureI );
		string feat_name = feat->GetName();
		if( feat_name == "CDS"  )
			count++;
		delete feat;
	}
	return count;
}

void featureIntersect( vector< bb_seqentry_t >& bb_list, uint seqI, vector< vector< size_t > >& intersecting, gnSequence* anno_seq )
{
	// stores the bb segs that overlap each feature
	intersecting.resize( anno_seq->getFeatureListLength() );

	uint seq_count = bb_list.front().size();

	vector< bb_seqentry_t* > bb_ptrs( bb_list.size() );
	for( size_t i = 0; i < bb_list.size(); ++i )
		bb_ptrs[i] = &bb_list[i];
	vector< bb_seqentry_t* > orig_ptrs( bb_ptrs );
	BbSeqComp bsc( seqI );
	std::sort( bb_ptrs.begin(), bb_ptrs.end(), bsc );
	vector< size_t > ptr_map;
	createMap( bb_ptrs, orig_ptrs, ptr_map );

	for( size_t featureI = 0; featureI < anno_seq->getFeatureListLength(); ++featureI )
	{
		gnBaseFeature* feat = anno_seq->getFeature( featureI );
		string feat_name = feat->GetName();
		if( feat_name != "CDS" && 
			feat_name != "tRNA" &&
			feat_name != "rRNA" &&
			feat_name != "misc_rna" )
			continue;	// don't deal with other feature types (source, misc_RNA, etc)
		gnLocation loc = feat->GetLocation(0);
		if( loc.GetFirst() > loc.GetLast() || loc.GetFirst() == 0 || loc.GetLast() == 0 )
			continue;	// a problem parsing annotation?
		// find where feature lands in our list
		bb_seqentry_t tmp_bb( seq_count );
		tmp_bb[seqI].first = loc.GetFirst();
		tmp_bb[seqI].second = loc.GetFirst();
		vector< bb_seqentry_t* >::iterator liter = std::lower_bound( bb_ptrs.begin(), bb_ptrs.end(), &tmp_bb, bsc );
		tmp_bb[seqI].first = loc.GetLast();
		tmp_bb[seqI].second = loc.GetLast();
		vector< bb_seqentry_t* >::iterator uiter = std::lower_bound( bb_ptrs.begin(), bb_ptrs.end(), &tmp_bb, bsc );
		if( liter == bb_ptrs.end() &&
			bb_ptrs.size() > 0 &&
			genome::absolut( (*bb_ptrs.back())[seqI].second ) >= loc.GetFirst() )
			liter--;
		while( liter != bb_ptrs.end() &&
			liter != bb_ptrs.begin() &&
			genome::absolut( (**liter)[seqI].second ) >= loc.GetFirst() )
			--liter;
		if( liter != bb_ptrs.end() &&
			genome::absolut( (**liter)[seqI].second ) < loc.GetFirst() )
			++liter;
		for( ; liter != uiter; ++liter )
		{
			if( (**liter)[seqI].first == 0 )
				continue;
			// only add the bbseg if the intersection is larger than SHORT_SEGMENT
			gnLocation bb_loc;
			if( (**liter)[seqI].first > 0 )
				bb_loc = gnLocation((**liter)[seqI].first, (**liter)[seqI].second);
			else
				bb_loc = gnLocation(-(**liter)[seqI].first, -(**liter)[seqI].second);

			gnLocation intersect = loc.GetIntersection( bb_loc, gnLocation::determinedRegions );
			if( intersect.GetLast() - intersect.GetFirst() <= SHORT_SEGMENT )
				continue;

			intersecting[ featureI ].push_back( ptr_map[ liter - bb_ptrs.begin() ] );
		}
		delete feat;
	}
}

void getFeatureHits( const vector< vector< size_t > >& intersecting, const bitset_t& segs, bitset_t& features_hit )
{
	features_hit.resize(intersecting.size());
	features_hit.reset();
	for( size_t featI = 0; featI < intersecting.size(); featI++ )
	{
		for( size_t i = 0; i < intersecting[featI].size(); ++i )
		{
			if( segs.test( intersecting[featI][i] ) )
				features_hit.set( featI );
		}		
	}
}

typedef map< string, map< string, double > > multifun_map_t;
typedef map< string, map< string, string > > multifun_names_t;

void makeMultiFunCount( gnSequence* anno_seq, multifun_map_t& mf_count, multifun_names_t& mf_names, bitset_t& feature_mask )
{
	for( size_t featureI = 0; featureI < anno_seq->getFeatureListLength(); ++featureI )
	{
		if( !feature_mask.test( featureI ) )
			continue;	// skip this feature if we're not supposed to include it
		gnBaseFeature* feat = anno_seq->getFeature( featureI );
		string feat_name = feat->GetName();
		if( feat_name != "CDS"  )
		{
			delete feat;
			continue;	
		}
		bool found_multifun = false;
		for( size_t qualI = 0; qualI < feat->GetQualifierListLength(); ++qualI )
		{
			gnBaseQualifier* gnq = feat->GetQualifier(qualI);
			string qual_name = gnq->GetName();
			if( qual_name != "function" )
			{
				delete gnq;
				continue;
			}
			string qual_value = gnq->GetValue();
			if( qual_value[0] == '"' )
				qual_value = qual_value.substr(1);
			stringstream qv_str( qual_value );
			string mf_level1;
			getline( qv_str, mf_level1, '.' );
			if( mf_level1.size() > 1 )
			{
				// not a multifun tag
				delete gnq;
				continue;
			}
			string mf_level2;
			mf_level2 += qv_str.get();
			mf_count[mf_level1][mf_level2]++;

			// get the name
			string l1_name;
			getline( qv_str, l1_name, ' ' );
			getline( qv_str, l1_name, ';' );
			string l2_name;
			getline( qv_str, l2_name, ';' );
			string cur_name = l1_name + ';' + l2_name;
			std::remove( cur_name.begin(), cur_name.end(), '\r' );
			std::remove( cur_name.begin(), cur_name.end(), '\n' );
			string space_str = "  ";
			boost::algorithm::erase_all( cur_name, space_str );
			mf_names[mf_level1][mf_level2] = cur_name;
			delete gnq;

			found_multifun = true;
		}
		// if we didn't find multifun, call it an "Unknown"
		if( !found_multifun )
		{
			string q = "?";
			mf_names[q][q] = "Unknown; No MultiFun Tag";
			mf_count[q][q]++;
		}

		delete feat;
	}
}

typedef boost::tuple< size_t, size_t, double, double, string > anal_row_t;
class AnalRowComp
{
public:
	bool operator()( const anal_row_t& a, const anal_row_t& b )
	{
		return a.get<2>() < b.get<2>();
	}
};

double chi_square_threshold = 5;
double min_expected_threshold = 5;
void mfAnalyze( ofstream& anal_output, multifun_map_t& all_mf, multifun_map_t& subset_mf, multifun_names_t& mf_names, double expect_freq )
{
	vector< anal_row_t > rows;
	multifun_map_t::iterator l1_iter = subset_mf.begin();
	for( ; l1_iter != subset_mf.end(); ++l1_iter )
	{
		multifun_map_t::iterator all_l1_iter = all_mf.find(l1_iter->first);
		multifun_names_t::iterator names_l1_iter = mf_names.find(l1_iter->first);
		map< string, double >::iterator l2_iter = l1_iter->second.begin();
		for( ; l2_iter != l1_iter->second.end(); ++l2_iter )
		{
			map< string, double >::iterator all_l2_iter = all_l1_iter->second.find(l2_iter->first);
			map< string, string >::iterator names_l2_iter = names_l1_iter->second.find(l2_iter->first);

			// percent in this category:
			double pct = (l2_iter->second / all_l2_iter->second) * 100;
			// category number:
			string cat_num = l1_iter->first + "." + l2_iter->first;
			// chi-square
			double chi_square = (l2_iter->second - (all_l2_iter->second*expect_freq));
			chi_square *= chi_square;
			chi_square /= (all_l2_iter->second*expect_freq);
			// total category gene count
			// category name
			if( chi_square < chi_square_threshold )
				continue;	// not significantly different
			if( (all_l2_iter->second*expect_freq) < min_expected_threshold )
				continue;	// don't have enough elements to make reliable estimation
			rows.push_back( boost::make_tuple( l2_iter->second, all_l2_iter->second, pct, chi_square, names_l2_iter->second ) );
		}
	}
	AnalRowComp arc;
	string col_delim = " & ";
	string new_row = "\\\\\n\\hline\n";
	std::sort( rows.begin(), rows.end(), arc );
	anal_output << "NumGenes" << col_delim << "GenesInCat" << col_delim << "Percent" << col_delim;
	anal_output << "Chi_square" << col_delim << "Mf_Level_2_name" << new_row;
	for( size_t rI = 0; rI < rows.size(); ++rI )
	{
		// if we transition from under to over-represented, output an empty row
		if( rI > 0 && rows[rI-1].get<2>() < expect_freq * 100 && rows[rI].get<2>() > expect_freq * 100 )
			anal_output << new_row;
		anal_output << rows[rI].get<0>() << col_delim;
		anal_output << rows[rI].get<1>() << col_delim;
		anal_output << setprecision(3) << rows[rI].get<2>() << col_delim;
		anal_output << setprecision(3) << rows[rI].get<3>() << col_delim;
		anal_output << rows[rI].get<4>() << new_row;
	}
}


void featureNearestNeighbors( const vector< bb_seqentry_t >& bb_list, const bitset_t& filter, uint seqI, vector< pair< size_t, size_t > >& neighbors, gnSequence* anno_seq, const vector< string >& feature_types )
{
	// stores the bb segs that overlap each feature
	neighbors.resize( bb_list.size() );

	uint seq_count = bb_list.front().size();

	vector< gnBaseFeature* > feats( anno_seq->getFeatureListLength() );
	vector< gnLocation > locs( anno_seq->getFeatureListLength() );
	vector< string > names( anno_seq->getFeatureListLength() );
	for( size_t featureI = 0; featureI < anno_seq->getFeatureListLength(); ++featureI )
	{
		feats[featureI] = anno_seq->getFeature( featureI );
		locs[featureI] = feats[featureI]->GetLocation(0);
		names[featureI] = feats[featureI]->GetName();
	}
	for( size_t bbI = 0; bbI < bb_list.size(); ++bbI )
	{
		// find the nearest feature
		size_t best_left = (std::numeric_limits<size_t>::max)();
		size_t best_right = (std::numeric_limits<size_t>::max)();
		size_t best_left_dist = (std::numeric_limits<size_t>::max)();
		size_t best_right_dist = (std::numeric_limits<size_t>::max)();
		if( !filter.test(bbI) )
		{
			neighbors[bbI].first = best_left;
			neighbors[bbI].second = best_right;
			continue;
		}
		for( size_t featI = 0; featI < feats.size(); ++featI )
		{
			size_t ntype = 0;
			for( ; ntype < feature_types.size(); ++ntype )
				if( names[featI] == feature_types[ntype] )
					break;
			if( ntype == feature_types.size() )
				continue;	// apparently this type of feature isn't interesting...
			if( locs[featI].GetFirst() > locs[featI].GetLast() || locs[featI].GetFirst() == 0 || locs[featI].GetLast() == 0 )
				continue;	// a problem parsing annotation?
			if( genome::absolut(bb_list[bbI][seqI].first) > locs[featI].GetLast() - ALTERNALOG_MIN_SIZE &&
				(int64)genome::absolut(bb_list[bbI][seqI].first) - (int64)locs[featI].GetLast() < best_left_dist )
			{
				best_left_dist = (int64)genome::absolut(bb_list[bbI][seqI].first) - (int64)locs[featI].GetLast();
				best_left = featI;
			}
			if( genome::absolut(bb_list[bbI][seqI].second) < locs[featI].GetFirst() + ALTERNALOG_MIN_SIZE &&
				(int64)locs[featI].GetFirst() - (int64)genome::absolut(bb_list[bbI][seqI].second) < best_right_dist )
			{
				best_right_dist = (int64)locs[featI].GetFirst() - (int64)genome::absolut(bb_list[bbI][seqI].second);
				best_right = featI;
			}
		}
		neighbors[bbI].first = best_left;
		neighbors[bbI].second = best_right;
	}
	// clean up
	for( size_t featureI = 0; featureI < feats.size(); ++featureI )
		delete feats[featureI];
}

void printFilteredBbSeqList( ostream& os, const vector< bb_seqentry_t >& bb_seq_list, const bitset_t& filter )
{
	for( size_t aI = 0; aI < bb_seq_list.size(); ++aI )
	{
		if( filter.test(aI) )
		{
			printBbSeq( os, bb_seq_list[aI] );
			os << endl;
		}
	}
}

void classifyIntergenic( ostream& os, const vector< bb_seqentry_t >& bbseq_list, const bitset_t& intergenic, 
						uint anno_seqI, gnSequence* anno_seq, bitset_t& trna_neighbor, bitset_t& miscrna_neighbor, 
						bitset_t& converging_cds, bitset_t& diverging_cds, bitset_t& inline_cds, 
						bitset_t& variable_miscrna, bitset_t& variable_trna )
{
	vector< pair< size_t, size_t > > all_neighbors;
	vector< string > all_types;
	all_types.push_back( "CDS" );
	all_types.push_back( "rRNA" );
	all_types.push_back( "tRNA" );
	all_types.push_back( "misc_RNA" );
	trna_neighbor.resize( bbseq_list.size() );
	miscrna_neighbor.resize( bbseq_list.size() );
	variable_miscrna.resize(anno_seq->getFeatureListLength());
	variable_trna.resize(anno_seq->getFeatureListLength());
	featureNearestNeighbors( bbseq_list, intergenic, anno_seqI, all_neighbors, anno_seq, all_types );
	for( size_t bbI = 0; bbI < bbseq_list.size(); ++bbI )
	{
		if( !intergenic.test(bbI) )
			continue;
		if( all_neighbors[bbI].first >= anno_seq->getFeatureListLength() ||
			all_neighbors[bbI].second >= anno_seq->getFeatureListLength() )
			continue;
		gnBaseFeature* lfeat = anno_seq->getFeature(all_neighbors[bbI].first);
		gnBaseFeature* rfeat = anno_seq->getFeature(all_neighbors[bbI].second);
		if( lfeat->GetName() == "tRNA" || rfeat->GetName() == "tRNA" )
			trna_neighbor.set(bbI);
		if( lfeat->GetName() == "tRNA" )
			variable_trna.set(all_neighbors[bbI].first);
		if( rfeat->GetName() == "tRNA" )
			variable_trna.set(all_neighbors[bbI].second);
		if( lfeat->GetName() == "misc_RNA" || rfeat->GetName() == "misc_RNA" )
			miscrna_neighbor.set(bbI);
		if( lfeat->GetName() == "misc_RNA" )
			variable_miscrna.set(all_neighbors[bbI].first);
		if( rfeat->GetName() == "misc_RNA" )
			variable_miscrna.set(all_neighbors[bbI].second);
		delete lfeat;
		delete rfeat;
	}
	
	vector< pair< size_t, size_t > > cds_neighbors;
	vector< string > cds_types;
	cds_types.push_back( "CDS" );
	featureNearestNeighbors( bbseq_list, intergenic, anno_seqI, cds_neighbors, anno_seq, cds_types );
	
	converging_cds.resize( bbseq_list.size() );
	diverging_cds.resize( bbseq_list.size() );
	inline_cds.resize( bbseq_list.size() );
	for( size_t bbI = 0; bbI < bbseq_list.size(); ++bbI )
	{
		if( !intergenic.test(bbI) )
			continue;
		if( cds_neighbors[bbI].first >= anno_seq->getFeatureListLength() ||
			cds_neighbors[bbI].second >= anno_seq->getFeatureListLength() )
			continue;
		gnBaseFeature* lfeat = anno_seq->getFeature(cds_neighbors[bbI].first);
		gnBaseFeature* rfeat = anno_seq->getFeature(cds_neighbors[bbI].second);
		if( lfeat->GetLocationType() == gnLocation::LT_Complement &&
			rfeat->GetLocationType() != gnLocation::LT_Complement )
			diverging_cds.set(bbI);
		else if( lfeat->GetLocationType() != gnLocation::LT_Complement &&
			rfeat->GetLocationType() == gnLocation::LT_Complement )
			converging_cds.set(bbI);
		else
			inline_cds.set(bbI);
		delete lfeat;
		delete rfeat;
	}
}

void findVariableSegmentsWithFlankingBB( const vector< bb_entry_t >& bb_list, const vector< double >& avg_lens, vector< pair< size_t, size_t > >& variable_segs, size_t min_bb_size = ALTERNALOG_MIN_SIZE, size_t min_variable_size = ALTERNALOG_MIN_SIZE, bool alternalogs = false )
{
	// find alternalogs (only at root node)
	const size_t NO_PREVIOUS = (std::numeric_limits<size_t>::max)();
	size_t prev_bb_seg = NO_PREVIOUS;
	uint seq_count = bb_list.front().bb_seq.size();
	for( size_t bbI = 0; bbI < bb_list.size(); ++bbI )
	{
		if( bb_list[bbI].bb_cols.Multiplicity() != seq_count ||
			avg_lens[bbI] < min_bb_size )
			continue;	// don't count this as n-way backbone
		if( prev_bb_seg == NO_PREVIOUS ||
			(bb_list[prev_bb_seg].iv != bb_list[bbI].iv)
			)
		{
			// no intervening alternalog...
			prev_bb_seg = bbI;
			continue;
		}
		// was there an alternalog?
		uint subset_count = 0;	// count the subset backbone of substantial size 
		bitset_t in_bb( seq_count );
		for( size_t segI = prev_bb_seg + 1; segI < bbI; ++segI )
		{
			if( avg_lens[segI] < min_variable_size )
				continue;
			bool found_new = false;
			for( size_t seqI = 0; seqI < seq_count; ++seqI )
			{
				if( bb_list[segI].bb_seq[seqI].first != 0 )
				{
					if( !in_bb.test(seqI) )
						found_new = true;
					in_bb.set(seqI);
				}
			}
			if( found_new )
				subset_count++;
		}
		for( size_t seqI = 0; seqI < seq_count; ++seqI )
		{
			// debug:
			if( (bb_list[bbI].bb_seq[seqI].first < 0 && bb_list[bbI].bb_seq[seqI].second > 0) ||
				(bb_list[bbI].bb_seq[seqI].first > 0 && bb_list[bbI].bb_seq[seqI].second < 0) ||
				(bb_list[bbI].bb_seq[seqI].first < 0 && bb_list[prev_bb_seg].bb_seq[seqI].first > 0) ||
				(bb_list[bbI].bb_seq[seqI].first > 0 && bb_list[prev_bb_seg].bb_seq[seqI].first < 0) ||
				(bb_list[bbI].bb_seq[seqI].first < 0 && bb_list[prev_bb_seg].bb_seq[seqI].second > 0) ||
				(bb_list[bbI].bb_seq[seqI].first > 0 && bb_list[prev_bb_seg].bb_seq[seqI].second < 0) )
			{
				cerr << "mismatch parity!!\n";
				genome::breakHere();
			}
			// normal:
			if( in_bb.test(seqI) )
				continue;
			int64 diff = 0;
			if( bb_list[bbI].bb_seq[seqI].first < 0 )
				diff = genome::absolut( bb_list[prev_bb_seg].bb_seq[seqI].first - bb_list[bbI].bb_seq[seqI].second );
			else
				diff = bb_list[bbI].bb_seq[seqI].first - bb_list[prev_bb_seg].bb_seq[seqI].second;
			if( diff >= min_variable_size )
				subset_count++;
		}
		if( alternalogs && subset_count > 1 )
			variable_segs.push_back( make_pair( prev_bb_seg, bbI ) );
		else if( !alternalogs && subset_count > 0 )
			variable_segs.push_back( make_pair( prev_bb_seg, bbI ) );
		prev_bb_seg = bbI;
	}
}

void makeVariableSegmentsCoordinateList( const vector< bb_entry_t >& bb_list, const vector< pair< size_t, size_t > >& alternalogs, vector< bb_seqentry_t >& alternabb_list )
{
	uint seq_count = bb_list.front().bb_seq.size();
	alternabb_list.resize( alternalogs.size() );
	for( size_t aI = 0; aI < alternalogs.size(); ++aI )
	{
		const bb_seqentry_t& a = bb_list[ alternalogs[aI].first ].bb_seq;
		const bb_seqentry_t& b = bb_list[ alternalogs[aI].second ].bb_seq;
		bb_seqentry_t alternabb = a;
		for( size_t seqI = 0; seqI < seq_count; ++seqI )
		{
			if( alternabb[seqI].first < 0 )
			{
				alternabb[seqI].first = b[seqI].second;
				alternabb[seqI].second = a[seqI].first;
			}
			else
			{
				alternabb[seqI].first = a[seqI].second;
				alternabb[seqI].second = b[seqI].first;
			}
		}
		alternabb_list[aI] = alternabb;
	}
}

class LocComp {
public:
	bool operator()( const gnLocation& a, const gnLocation& b ) const
	{
		return a.GetFirst() < b.GetFirst();
	}
};

void identifyIntergenicRanges( 	vector< gnSequence* >& seq_table, vector< vector< pair< size_t, size_t > > >& ranges )
{
	ranges.resize(seq_table.size());
	for( size_t seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		vector< gnLocation > loc_list;
		for( size_t featI = 0; featI < seq_table[seqI]->getFeatureListLength(); featI++ )
		{
			gnBaseFeature* feat = seq_table[seqI]->getFeature(featI);
			string feat_name = feat->GetName();
			if( feat_name != "CDS" )
				continue;	// don't deal with other feature types (source, etc)
			loc_list.push_back( feat->GetLocation(0) );
			delete feat;
		}

		size_t sum = 0;
		LocComp lc;
		std::sort( loc_list.begin(), loc_list.end(), lc );
		size_t fI = 0; 
		size_t lI = 1; 
		while( fI < loc_list.size() && lI < loc_list.size() )
		{
			if( loc_list[fI].GetLast() < loc_list[lI].GetFirst() )
			{
				ranges[seqI].push_back( make_pair( loc_list[fI].GetLast(), loc_list[lI].GetFirst() ) );
				sum += loc_list[lI].GetFirst() - loc_list[fI].GetLast() - 1;
			}
			fI++; lI++;
			while( fI < loc_list.size() && lI < loc_list.size() &&
				loc_list[fI].GetLast() >= loc_list[lI].GetFirst() )
			{
				if( loc_list[fI].GetLast() >= loc_list[lI].GetLast() )
				{
					fI++; lI++;
					cerr << "danger, complete containment in seq " << seqI << endl;
				}
				fI++; lI++;
			}
		}
	}
}

//big_coli_sam_fixed_goh0001_gou000001.xmfa guide.tre big_coli_sam_fixed_goh0001_gou000001.xmfa.backbone big_coli_sam_fixed_goh0001_gou000001.xmfa.bbcols 5 bb.out

void classifyCoordinateRanges( 
			const vector< bb_seqentry_t >& alternabb_list,			
			gnSequence* annotated_seq,
			vector< gnSequence* >& seq_table,
			vector< bitset_t >& genic, 
			vector< bitset_t >& genic_fudge, 
			vector< bitset_t >& overlaps_cds_upstream, 
			vector< bitset_t >& overlaps_cds_upstream_fudge, 
			vector< bitset_t >& overlaps_cds_downstream, 
			vector< bitset_t >& overlaps_cds_downstream_fudge, 
			vector< bitset_t >& intergenic, 
			vector< bitset_t >& spanner,
			vector< bitset_t >& trna, 
			vector< bitset_t >& rrna,
			vector< bitset_t >& miscrna,
			vector< bitset_t >& pseudogenized,
			vector< bitset_t >& variable_miscrna,
			vector< bitset_t >& variable_trna,
			vector< bitset_t >& intergenic_segs
			)
{
	if( alternabb_list.size() == 0 )
		return;
	uint seq_count = seq_table.size();
	// count genic vs. intergenic alternalogs
	// classify alternalogs as genic, intergenic, multigenic
	// and pseudogenizing
	bitset_t bbclass_tmp( alternabb_list.size() );
	// all of these classifications should be mutually exclusive
	genic.resize( seq_count, bbclass_tmp );
	genic_fudge.resize( seq_count, bbclass_tmp );

	// set to true if a variable segment ends in a CDS, but isn't contained by the CDS
	overlaps_cds_upstream.resize( seq_count, bbclass_tmp );
	overlaps_cds_upstream_fudge.resize( seq_count, bbclass_tmp );
	overlaps_cds_downstream.resize( seq_count, bbclass_tmp );
	overlaps_cds_downstream_fudge.resize( seq_count, bbclass_tmp );
	intergenic.resize( seq_count, bbclass_tmp );
	spanner.resize( seq_count, bbclass_tmp );
//	vector< bitset_t > multigenic( seq_count, bbclass_tmp );
	// these are true if trna or rrna are intersected
	trna.resize( seq_count, bbclass_tmp );
	rrna.resize( seq_count, bbclass_tmp );
	miscrna.resize( seq_count, bbclass_tmp );
	variable_miscrna.resize( seq_count );
	variable_trna.resize( seq_count );
	// an alternalog is pseudogenizing if it's genic in other sequences but not in the subject
	pseudogenized.resize( seq_count, bbclass_tmp );

	vector< vector< pair< size_t, size_t > > > ranges;
	identifyIntergenicRanges( seq_table, ranges );
	intergenic_segs.resize(seq_table.size());

	vector< const bb_seqentry_t* > alterna_ptrs( alternabb_list.size() );
	for( size_t i = 0; i < alternabb_list.size(); ++i )
		alterna_ptrs[i] = &alternabb_list[i];
	vector< const bb_seqentry_t* > orig_ptrs( alterna_ptrs );
	for( size_t seqI = 0; seqI < seq_count; ++seqI )
	{
		BbSeqComp bsc( seqI );
		std::sort( alterna_ptrs.begin(), alterna_ptrs.end(), bsc );
		vector< size_t > ptr_map;
		createMap( alterna_ptrs, orig_ptrs, ptr_map );

		vector< vector< unsigned > > bb_features( alternabb_list.size() );	// stores feature IDs of overlapping features
		variable_miscrna[seqI].resize(seq_table[seqI]->getFeatureListLength());
		variable_trna[seqI].resize(seq_table[seqI]->getFeatureListLength());
		for( size_t featureI = 0; featureI < seq_table[seqI]->getFeatureListLength(); ++featureI )
		{
			gnBaseFeature* feat = seq_table[seqI]->getFeature( featureI );
			string feat_name = feat->GetName();
			if( feat_name != "CDS" && 
				feat_name != "tRNA" &&
				feat_name != "rRNA" &&
				feat_name != "misc_RNA" )
				continue;	// don't deal with other feature types (source, etc)
			if( feat->GetLocationListLength() > 1 )
				continue;	// any multi-part CDS features are likely to be pseudogene annotations
							// which we don't want to bias our results.  There are only a couple true multi-part
							// CDS in enteric bacteria

			gnLocation loc = feat->GetLocation(0);
			if( loc.GetFirst() > loc.GetLast() || loc.GetFirst() == 0 || loc.GetLast() == 0 )
				continue;	// a problem parsing annotation?
			// find where feature lands in our list
			bb_seqentry_t tmp_bb( seq_count );
			tmp_bb[seqI].first = loc.GetFirst();
			tmp_bb[seqI].second = loc.GetFirst();
			vector< const bb_seqentry_t* >::iterator liter = std::lower_bound( alterna_ptrs.begin(), alterna_ptrs.end(), &tmp_bb, bsc );
			tmp_bb[seqI].first = loc.GetLast();
			tmp_bb[seqI].second = loc.GetLast();
			vector< const bb_seqentry_t* >::iterator uiter = std::lower_bound( alterna_ptrs.begin(), alterna_ptrs.end(), &tmp_bb, bsc );
			if( liter == alterna_ptrs.end() &&
				alterna_ptrs.size() > 0 &&
				genome::absolut( (*alterna_ptrs.back())[seqI].second ) >= loc.GetFirst() )
				liter--;
			while( liter != alterna_ptrs.end() &&
				liter != alterna_ptrs.begin() &&
				genome::absolut( (**liter)[seqI].second ) >= loc.GetFirst() )
				--liter;
			if( liter != alterna_ptrs.end() &&
				genome::absolut( (**liter)[seqI].second ) < loc.GetFirst() )
				++liter;
			for( ; liter != uiter; ++liter )
			{
				bb_features[ liter - alterna_ptrs.begin() ].push_back( featureI );
			}
			delete feat;
		}

		intergenic_segs[seqI].resize(ranges[seqI].size());
		for( size_t bbI = 0; bbI < alterna_ptrs.size(); ++bbI )
		{
			size_t l = (*alterna_ptrs[bbI])[seqI].first;
			size_t r = (*alterna_ptrs[bbI])[seqI].second;
			for( size_t rI = 0; rI < ranges[seqI].size(); ++rI )
			{
				if( (l < ranges[seqI][rI].first + 1 && ranges[seqI][rI].first + 1 <= r) ||		// left overlap and complete contains
					(l <= ranges[seqI][rI].second - 1 && ranges[seqI][rI].first + 1 <= r) )		// right overlap and inside
				{
					intergenic_segs[seqI].set(rI);
					break;
				}
			}
		}

		for( size_t bbI = 0; bbI < alterna_ptrs.size(); ++bbI )
		{
			gnLocation bb_loc;
			if( (*alterna_ptrs[bbI])[seqI].first > 0 )
				bb_loc = gnLocation((*alterna_ptrs[bbI])[seqI].first, (*alterna_ptrs[bbI])[seqI].second);
			else
				bb_loc = gnLocation(-(*alterna_ptrs[bbI])[seqI].first, -(*alterna_ptrs[bbI])[seqI].second);
			if( (*alterna_ptrs[bbI])[0].first > 2302400 && (*alterna_ptrs[bbI])[0].second < 2303211 )
			{
				cerr << "debugme\n";
			}
			for( size_t featI = 0; featI < bb_features[bbI].size(); ++featI )
			{
				gnBaseFeature* feat = seq_table[seqI]->getFeature( bb_features[bbI][featI] );
				gnLocation feat_loc = feat->GetLocation(0);
				gnLocation intersect = feat_loc.GetIntersection( bb_loc, gnLocation::determinedRegions );
				string name = feat->GetName();
				if( intersect.GetFirst() == bb_loc.GetFirst() &&
					intersect.GetLast() == bb_loc.GetLast() && 
					name == "CDS" )
				{
					if( intersect.GetLast() - intersect.GetFirst() > ALTERNALOG_MIN_SIZE ||
						intersect.GetFirst() - ALTERNALOG_MIN_SIZE > feat_loc.GetFirst() ||
						intersect.GetLast() + ALTERNALOG_MIN_SIZE < feat_loc.GetLast() )
					{
						// alternalog completely contained by CDS, at least ALTERNALOG_MIN_SIZE inside the CDS
						genic[seqI].set( ptr_map[bbI] );
					}else{
						genic_fudge[seqI].set( ptr_map[bbI] );	// small and close to the edge
					}
				}
				else if( (intersect.GetFirst() == bb_loc.GetFirst() ||
					intersect.GetLast() == bb_loc.GetLast()) && 
					name == "CDS" )
				{
					bool up = false;
					// overlaps a CDS by at least ALTERNALOG_MIN_SIZE nucleotides,
					// but does not contain the CDS, nor is it contained by the CDS
					if( intersect.GetFirst() == bb_loc.GetFirst() )
					{
						if( feat->GetLocationType() != gnLocation::LT_Standard )
							up = true;
					}else if( feat->GetLocationType() == gnLocation::LT_Standard )
						up = true;

					if( !up && intersect.GetLast() - intersect.GetFirst() > ALTERNALOG_MIN_SIZE )
						overlaps_cds_downstream[seqI].set( ptr_map[bbI] );
					else if( !up )
						overlaps_cds_downstream_fudge[seqI].set( ptr_map[bbI] );
					else if( up && intersect.GetLast() - intersect.GetFirst() > ALTERNALOG_MIN_SIZE )
						overlaps_cds_upstream[seqI].set( ptr_map[bbI] );
					else
						overlaps_cds_upstream_fudge[seqI].set( ptr_map[bbI] );

				}else if( intersect.GetLast() - intersect.GetFirst() > ALTERNALOG_MIN_SIZE &&
					name == "CDS" )
				{
					// spans CDS
					spanner[seqI].set( ptr_map[bbI] );
				}
				if( intersect.GetLast() - intersect.GetFirst() > ALTERNALOG_MIN_SIZE &&
					name == "rRNA" )
				{
					// overlaps a rRNA by at least ALTERNALOG_MIN_SIZE nucleotides
					rrna[seqI].set( ptr_map[bbI] );
				}
				if( intersect.GetLast() - intersect.GetFirst() > ALTERNALOG_MIN_SIZE &&
					name == "tRNA" )
				{
					// overlaps a tRNA by at least ALTERNALOG_MIN_SIZE nucleotides
					trna[seqI].set( ptr_map[bbI] );
					variable_trna[seqI].set(bb_features[bbI][featI]);
				}
				if( intersect.GetLast() - intersect.GetFirst() > ALTERNALOG_MIN_SIZE &&
					name == "misc_RNA" )
				{
					// overlaps a misc_RNA by at least ALTERNALOG_MIN_SIZE nucleotides
					miscrna[seqI].set( ptr_map[bbI] );
					variable_miscrna[seqI].set(bb_features[bbI][featI]);
				}
				delete feat;
			}
		}
		intergenic[seqI] = genic[seqI] | overlaps_cds_upstream[seqI] | overlaps_cds_downstream[seqI] | rrna[seqI] | trna[seqI];
		intergenic[seqI].flip();
	}

	// identify pseudogenizing segments as intergenic segments in one genome that
	// are genic in other genomes
	size_t seqI = 0;
	for( seqI = 0; seqI < seq_count; ++seqI )
	{
		bitset_t pseudo = bbclass_tmp;
		for( size_t seqJ = 0; seqJ < seq_count; ++seqJ )
		{
			if( seqJ == seqI )
				continue;
			pseudo |= genic[seqJ] | overlaps_cds_upstream[seqJ] | overlaps_cds_downstream[seqJ];
		}
		bitset_t fudge = genic_fudge[seqI] | overlaps_cds_upstream_fudge[seqI] | overlaps_cds_downstream_fudge[seqI];
		fudge.flip();	// if it's questionably within a gene then don't let it be a pseudogene.  we want to be sure
						// about these
		pseudogenized[seqI] = intergenic[seqI] & pseudo & fudge;
	}
}

void analyzeVariableSegments( ostream& os, const vector< bb_entry_t >& bb_list, const vector< double >& avg_lens, uint anno_seqI, vector< gnSequence* >& seq_table, string site_class_name = "alternalog", bool analyze_alternalogs = true )
{
	gnSequence* annotated_seq = seq_table[anno_seqI];
	vector< pair< size_t, size_t > > alternalogs;
	vector< bb_seqentry_t > alternabb_list;
	findVariableSegmentsWithFlankingBB( bb_list, avg_lens, alternalogs, ALTERNALOG_MIN_SIZE, ALTERNALOG_MIN_SIZE, analyze_alternalogs );
	makeVariableSegmentsCoordinateList( bb_list, alternalogs, alternabb_list );

	os << "There are " << alternalogs.size() << " " << site_class_name << " sites\n";

	// count genic vs. intergenic alternalogs
	// classify alternalogs as genic, intergenic, etc.
	vector< bitset_t > alt_genic, alt_overlaps_cds_upstream, alt_overlaps_cds_downstream;
	vector< bitset_t > alt_intergenic, alt_spanner, alt_trna, alt_rrna, alt_pseudogenized;
	vector< bitset_t > alt_genic_fudge, alt_overlaps_cds_upstream_fudge, alt_overlaps_cds_downstream_fudge;
	vector< bitset_t > alt_miscrna, v_miscrna, v_trna;
	vector< bitset_t > intergenic_segs;

	classifyCoordinateRanges( 
		alternabb_list, annotated_seq, seq_table, alt_genic, alt_genic_fudge, alt_overlaps_cds_upstream,
		alt_overlaps_cds_upstream_fudge, alt_overlaps_cds_downstream, alt_overlaps_cds_downstream_fudge,
		alt_intergenic, alt_spanner, alt_trna, alt_rrna, alt_miscrna, alt_pseudogenized, v_miscrna, v_trna, 
		intergenic_segs 
		);

	// find alternalogs that are always inside annotated genes
	bitset_t bbclass_tmp( alternabb_list.size() );
	bitset_t alt_multi_allelic_genes( bbclass_tmp );
	alt_multi_allelic_genes.flip();
	// alternalogs that are always outside genes
	bitset_t alt_multi_allelic_intergenic( bbclass_tmp );
	bitset_t alt_multi_allelic_entirely_intergenic( bbclass_tmp );
	alt_multi_allelic_intergenic.flip();
	alt_multi_allelic_entirely_intergenic.flip();
	uint seq_count = bb_list.front().bb_seq.size();
	for( uint seqI = 0; seqI < seq_count; ++seqI )
	{
		alt_multi_allelic_genes &= alt_genic[seqI];
		alt_multi_allelic_intergenic &= alt_intergenic[seqI];
		bitset_t spanner_flip = alt_spanner[seqI];
		spanner_flip.flip();
		alt_multi_allelic_entirely_intergenic &= alt_intergenic[seqI] & spanner_flip;
	}


	os << " There are " << alt_multi_allelic_genes.count() << " apparently multi-allelic genes (" << site_class_name << ")\n";
	os << " There are " << alt_multi_allelic_intergenic.count() << " apparently multi-allelic regions with intergenic endpoints (" << site_class_name << ")\n";
	os << " Of those, " << alt_multi_allelic_entirely_intergenic.count() << " contain no annotated CDS (" << site_class_name << ")\n";
	os << " The remaining segments span gene boundaries, but are not entirely contained in annotated genes\n";

	bitset_t trna_neighbor;
	bitset_t miscrna_neighbor;
	bitset_t converging_cds;
	bitset_t diverging_cds;
	bitset_t inline_cds;
	bitset_t vv_miscrna;
	bitset_t vv_trna;
	classifyIntergenic( os, alternabb_list, alt_multi_allelic_intergenic, anno_seqI, 
		annotated_seq, trna_neighbor, miscrna_neighbor, converging_cds, diverging_cds, inline_cds, vv_miscrna, vv_trna );


	os << "There are " << trna_neighbor.count() << " intergenic segments with a tRNA nearest neighbor\n";
	os << "There are " << miscrna_neighbor.count() << " intergenic segments with a miscRNA nearest neighbor\n";
	os << "There are " << converging_cds.count() << " intergenic segments surrounded by converging CDS\n";
	os << "There are " << diverging_cds.count() << " intergenic segments surrounded by diverging CDS\n";
	os << "There are " << inline_cds.count() << " intergenic segments surrounded by inline CDS\n";
	bitset_t miscrna_inter = v_miscrna[anno_seqI] | vv_miscrna;
	os << "There are " << miscrna_inter.count() << " annotated misc_RNA associated with variable segments\n";
	os << "There are " << intergenic_segs[anno_seqI].size() << " intergenic sites in the ref genome, of which " << intergenic_segs[anno_seqI].count() << " exhibit variability\n";
	bitset_t trna_inter = v_trna[anno_seqI] | vv_trna;
	os << "There are " << trna_inter.count() << " annotated tRNA associated with variable segments\n";

	if( miscrna_neighbor.count() > 0 )
	{
		os << "coordinates of variable segs with misc_RNA neighboring:\n";
		printFilteredBbSeqList( os, alternabb_list, miscrna_neighbor );
	}
	if( diverging_cds.count() > 0 )
	{
		os << "coordinates of variable segs with diverging_cds neighboring:\n";
		printFilteredBbSeqList( os, alternabb_list, diverging_cds );
	}
	bitset_t total_miscrna = alt_miscrna[anno_seqI] | miscrna_neighbor;
	os << "Total variable intergenic segs that neighbor or contain miscRNA: " << total_miscrna.count() << endl;

	os << "coordinates of multi-allelic genes:\n";
	printFilteredBbSeqList( os, alternabb_list, alt_multi_allelic_genes );

	os << "coordinates of multi-allelic intergenic regions without CDS:\n";
	printFilteredBbSeqList( os, alternabb_list, alt_multi_allelic_entirely_intergenic );

	for( size_t seqI = 0; seqI < seq_count; ++seqI )
	{
		os << "genome " << seqI << " has " << alt_genic[seqI].count() << " " << site_class_name << " within CDS\n";
		os << "genome " << seqI << " has " << alt_spanner[seqI].count() << " " << site_class_name << " that span CDS boundaries\n";
		os << "genome " << seqI << " has " << alt_intergenic[seqI].count() << " " << site_class_name << " that lie entirely in intergenic regions\n";
		os << "genome " << seqI << " has " << alt_rrna[seqI].count() << " " << site_class_name << " that contain rRNA\n";
		os << "genome " << seqI << " has " << alt_trna[seqI].count() << " " << site_class_name << " that contain tRNA\n";
		os << "genome " << seqI << " has " << alt_miscrna[seqI].count() << " " << site_class_name << " that contain misc_RNA\n";
		os << "genome " << seqI << " has " << alt_pseudogenized[seqI].count() << " apparent recent pseudogenes in " << site_class_name << "\n";
		os.flush();

/*
		os << "coordinates of genic alternalogs:\n";
		printFilteredBbSeqList( os, alternabb_list, genic[seqI] );
*/

		if( alt_trna[seqI].count() > 0 )
		{
			os << "coordinates of tRNA " << site_class_name << ":\n";
			printFilteredBbSeqList( os, alternabb_list, alt_trna[seqI] );
		}

		if( alt_rrna[seqI].count() > 0 )
		{
			os << "coordinates of rRNA " << site_class_name << ":\n";
			printFilteredBbSeqList( os, alternabb_list, alt_rrna[seqI] );
		}

		if( alt_miscrna[seqI].count() > 0 )
		{
			os << "coordinates of misc_RNA " << site_class_name << ":\n";
			printFilteredBbSeqList( os, alternabb_list, alt_miscrna[seqI] );
		}

		os << "coordinates of possible pseudogenes:\n";
		printFilteredBbSeqList( os, alternabb_list, alt_pseudogenized[seqI] );
		os.flush();
	}
}

const uint INTERNAL_NODE = (std::numeric_limits<uint>::max)();
const uint INTERVAL_UNKNOWN = (std::numeric_limits<uint>::max)();

int main( int argc, char* argv[] )
{
#if	WIN32
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif

	if( argc < 7 )
	{
		cerr << "bbAnalyze <xmfa file> <guide tree> <backbone seqpos file> <backbone col file> <annotated seq index> <output file>\n";
		cerr << "annotated seq index starts at 0.\n";
		return -1;
	}
	string aln_fname( argv[1] );
	string guide_tree_fname( argv[2] );
	string bbseq_fname( argv[3] );
	string bbcol_fname( argv[4] );
	int gff_seq_index = atoi( argv[5] );
	string output_fname( argv[6] );

	ifstream aln_input( aln_fname.c_str() );
	if( !aln_input.is_open() ){
		cerr << "Error opening \"" << aln_fname << "\"" << endl;
		return -2;
	}
	ifstream tree_input( guide_tree_fname.c_str() );
	if( !tree_input.is_open() ){
		cerr << "Error opening \"" << guide_tree_fname << "\"" << endl;
		return -3;
	}
	ifstream bbseq_input( bbseq_fname.c_str() );
	if( !bbseq_input.is_open() ){
		cerr << "Error opening \"" << bbseq_fname << "\"" << endl;
		return -4;
	}
	ifstream bbcol_input( bbcol_fname.c_str() );
	if( !bbcol_input.is_open() ){
		cerr << "Error opening \"" << bbcol_fname << "\"" << endl;
		return -4;
	}
	ofstream anal_output( output_fname.c_str() );
	if( !anal_output.is_open() ){
		cerr << "Error opening \"" << output_fname << "\" for writing" << endl;
		return -6;
	}
	
	// read the guide tree
	PhyloTree< TreeNode > tree;
	tree.readTree( tree_input );

	// read the backbone column file	
	vector< bb_seqentry_t > bb_seq_list;
	vector< pair< size_t, ULA > > bb_col_list;
	readBackboneSeqFile( bbseq_input, bb_seq_list );
	readBackboneColsFile( bbcol_input, bb_col_list );

	// read the alignment
	IntervalList iv_list;
	iv_list.ReadStandardAlignment( aln_input );

	LoadSequences(iv_list, &cout);



	const size_t seq_count = iv_list.seq_table.size();

	vector< bb_entry_t > bb_list( bb_seq_list.size() );
	for( size_t i = 0; i < bb_seq_list.size(); ++i )
	{
		bb_list[i].bb_seq = bb_seq_list[i];
		bb_list[i].bb_cols = bb_col_list[i].second;
		bb_list[i].iv = bb_col_list[i].first;
		// awful hack: homogenize the parity inside intervals.  this is a bug in progressiveMauve
		for( size_t seqI = 0; seqI < seq_count; ++seqI )
		{
			AbstractMatch::orientation o = iv_list[bb_list[i].iv].Orientation(seqI);
			if( o == AbstractMatch::undefined )
				continue;
			if( bb_list[i].bb_cols.LeftEnd(seqI) != NO_MATCH )
				bb_list[i].bb_cols.SetOrientation( seqI, o );
			if( (bb_list[i].bb_seq[seqI].first < 0 && o == AbstractMatch::forward) ||
				(bb_list[i].bb_seq[seqI].first > 0 && o == AbstractMatch::reverse) )
				bb_list[i].bb_seq[seqI].first *= -1;
			if( (bb_list[i].bb_seq[seqI].second < 0 && o == AbstractMatch::forward) ||
				(bb_list[i].bb_seq[seqI].second > 0 && o == AbstractMatch::reverse) )
				bb_list[i].bb_seq[seqI].second *= -1;
			if( genome::absolut( bb_list[i].bb_seq[seqI].first ) > genome::absolut( bb_list[i].bb_seq[seqI].second ) )
				swap( bb_list[i].bb_seq[seqI].first, bb_list[i].bb_seq[seqI].second );
		}
	}


	// make faux single-genome bb segments for anything not contained in
	// real backbone
	for( uint seqI = 0; seqI < seq_count; seqI++ )
	{
		vector< AbstractMatch* > seq_beeb;
		ULA single_ula(1);
		for( size_t i = 0; i < bb_seq_list.size(); ++i )
		{
			if( bb_seq_list[i][seqI].first == 0 )
				continue;
			single_ula.SetStart( 0, genome::absolut(bb_seq_list[i][seqI].first) );
			single_ula.SetLength( genome::absolut(bb_seq_list[i][seqI].second - bb_seq_list[i][seqI].first) + 1 );
			seq_beeb.push_back( single_ula.Copy() );
		}
		SingleStartComparator<AbstractMatch> ssc(0);
		sort( seq_beeb.begin(), seq_beeb.end(), ssc );
		// HACK!!
		// trim single base pair overlaps in seq_beeb that arise due to an off-by-one bug in the backbone output...
		EliminateOverlaps_v2( seq_beeb );
		sort( seq_beeb.begin(), seq_beeb.end(), ssc );
		list< AbstractMatch* > seq_beeb_list( seq_beeb.begin(), seq_beeb.end() );
		AddGapMatches( seq_beeb_list, seq_beeb_list.begin(), seq_beeb_list.end(), 
			   0, 1, iv_list.seq_table[seqI]->length()+1, AbstractMatch::forward, 1 );
		sort( seq_beeb.begin(), seq_beeb.end() );
		vector< AbstractMatch* > tmp_list( seq_beeb_list.begin(), seq_beeb_list.end() );
		sort( tmp_list.begin(), tmp_list.end() );
		vector< AbstractMatch* > new_beeb( seq_beeb_list.size() - seq_beeb.size() );
		std::set_difference(  tmp_list.begin(), tmp_list.end(),
				seq_beeb.begin(), seq_beeb.end(), new_beeb.begin() );

		// add each new_beeb to the backbone list
		size_t bbI = bb_list.size();
		bb_list.resize( bbI + new_beeb.size() );
		for( size_t i = 0; i < new_beeb.size(); ++i )
		{
			bb_list[bbI].bb_seq.resize( seq_count );
			bb_list[bbI].bb_seq[seqI] = make_pair( new_beeb[i]->LeftEnd(0), new_beeb[i]->RightEnd(0) );
			bb_list[bbI].iv = INTERVAL_UNKNOWN;
			ULA cols(seq_count);
			cols.SetLeftEnd(seqI, 1);
			cols.SetLength(new_beeb[i]->Length(0));
			bb_list[bbI].bb_cols = cols;
			bbI++;
		}
		for( size_t i = 0; i < tmp_list.size(); ++i )
			tmp_list[i]->Free();
	}

	// create a map between tree nodes and sequences
	vector< uint > node_sequence_map( tree.size(), -1 );
	for( uint seqI = 0; seqI < seq_count; seqI++ )
	{
		stringstream seq_name;
		seq_name << "seq" << seqI + 1;
		node_id_t nodeI = 0;
		for( ; nodeI < tree.size(); nodeI++ )
		{
			if( seq_name.str() == tree[nodeI].name )
			{
				node_sequence_map[nodeI] = seqI;
				break;
			}
		}
		if( nodeI == tree.size() )
			throw "Phylogenetic tree names unrecognized.  Should follow seqN naming format\n";
	}

	// mark small backbone segments
	bitset_t too_small( bb_list.size(), false );
	for( size_t bbI = 0; bbI < bb_list.size(); ++bbI )
		if( bb_list[bbI].bb_cols.Length() < DISCARD_SEGMENT )
			too_small.set(bbI, true);
	bitset_t not_small = too_small;
	not_small.flip();

	vector< double > avg_lens( bb_list.size(), 0 );
	for( size_t bbI = 0; bbI < bb_list.size(); ++bbI )
	{
		double ct = 0;
		for( size_t seqI = 0; seqI < bb_list[bbI].bb_seq.size(); ++seqI )
		{
			if( bb_list[bbI].bb_seq[seqI].first != 0 )
			{
				ct++;
				avg_lens[bbI] += genome::absolut( bb_list[bbI].bb_seq[seqI].second - bb_list[bbI].bb_seq[seqI].first ) + 1;
			}
		}
		avg_lens[bbI] /= ct;
	}


	// got the backbone.  now do something with it.
	// at each node of the tree, count the total amount backbone contained in nodes
	// below that tree, both inside genes and outside genes


	vector< node_id_t > all_leaves;
	getLeaves( tree, tree.root, all_leaves );
	sort( all_leaves.begin(), all_leaves.end() );

	bitset_t true_temper( bb_list.size() );
	true_temper.reset();
	true_temper.flip();
	bitset_t false_temper( bb_list.size() );
	false_temper.reset();

	vector< bitset_t > unique( tree.size(), true_temper );
	// partial contains bb segs that have representation among two or more genomes below a given node
	vector< bitset_t > partial( tree.size(), true_temper );
	// conserved have representation in all genomes below a node, and possibly others
	vector< bitset_t > conserved( tree.size(), true_temper );
	// child partial have representation in one or more genomes below a node
	vector< bitset_t > c1_partial( tree.size(), false_temper );
	vector< bitset_t > c2_partial( tree.size(), false_temper );
	vector< bitset_t > c1_complete( tree.size(), false_temper );
	vector< bitset_t > c2_complete( tree.size(), false_temper );

	// calculate which segments have heterogenous occurrence at each node
	vector< bitset_t > hop_one( tree.size(), false_temper );
	vector< bitset_t > hop_two( tree.size(), false_temper );
	vector< double > pan_genome_size( tree.size(), 0 );
	// hop_two if(c1_partial && c2_partial) && !c1_complete && !c2_complete
	// hop_one if !hop_two && (!c1_complete || !c2_complete) && (c1_partial && c2_partial) && !(hop_one at incomplete child)

	stack< node_id_t > node_stack;
	node_stack.push( tree.root );
	bitset_t visited( tree.size(), false );
	while( node_stack.size() > 0 )
	{
		node_id_t nI = node_stack.top();
		if( !visited[nI] && tree[nI].children.size() > 0 )
		{
			node_stack.push( tree[nI].children[0] );
			node_stack.push( tree[nI].children[1] );
			visited.set(nI,true);
			continue;	// visit post-order
		}
		node_stack.pop();

		vector< node_id_t > leaves;
		getLeaves( tree, nI, leaves );
		sort( leaves.begin(), leaves.end() );

		vector< node_id_t > not_leaves( all_leaves.size() - leaves.size() );
		std::set_difference( all_leaves.begin(), all_leaves.end(), 
			leaves.begin(), leaves.end(), 
			not_leaves.begin() );


		vector< node_id_t > c1_leaves;
		vector< node_id_t > c2_leaves;
		if( tree[nI].children.size() > 0 )
		{
			getLeaves( tree, tree[nI].children[0], c1_leaves );
			getLeaves( tree, tree[nI].children[1], c2_leaves );
		}

		for( size_t bbI = 0; bbI < bb_list.size(); ++bbI )
		{
			// do all the leaves have this segment?
			size_t lI = 0;
			size_t ct = 0;
			for( lI = 0; lI < leaves.size(); ++lI )
			{
				if( bb_list[bbI].bb_seq[ node_sequence_map[ leaves[lI] ] ].first != 0 )
					ct++;
			}
			unique[nI].set(bbI, ct == leaves.size());

			// was this conserved in more than one?
			partial[nI].set(bbI, ct > 1);
			conserved[nI].set(bbI, ct == leaves.size());

			// if this one was represented at all then it's part of the pan-genome
			if( ct > 0 )
				pan_genome_size[nI] += avg_lens[bbI];

			// do only the leaves below this node have this segment?
			for( lI = 0; lI < not_leaves.size(); ++lI )
			{
				if( bb_list[bbI].bb_seq[ node_sequence_map[ not_leaves[lI] ] ].first != 0 )
					unique[nI].set(bbI, false);
			}

			// is the segment present in both children?
			bool c1 = false;
			bool c2 = false;
			uint c1_ct = 0;
			uint c2_ct = 0;
			for( lI = 0; lI < c1_leaves.size(); ++lI )
			{
				if( bb_list[bbI].bb_seq[ node_sequence_map[ c1_leaves[lI] ] ].first != 0 )
					c1_ct++;
			}
			for( lI = 0; lI < c2_leaves.size(); ++lI )
			{
				if( bb_list[bbI].bb_seq[ node_sequence_map[ c2_leaves[lI] ] ].first != 0 )
					c2_ct++;
			}
			c1_partial[nI].set(bbI, c1_ct > 0);
			c2_partial[nI].set(bbI, c2_ct > 0);
			c1_complete[nI].set(bbI, c1_ct == c1_leaves.size());
			c2_complete[nI].set(bbI, c2_ct == c2_leaves.size());
		}
	}

	node_stack.push( tree.root );
	visited = bitset_t( tree.size(), false );
	vector< bitset_t > all_unique( tree.size(), false_temper );
	while( node_stack.size() > 0 )
	{
		node_id_t nI = node_stack.top();
		if( !visited[nI] && tree[nI].children.size() > 0 )
		{
			node_stack.push( tree[nI].children[0] );
			node_stack.push( tree[nI].children[1] );
			visited.set(nI,true);
			continue;	// visit post-order
		}
		node_stack.pop();

		all_unique[nI] = unique[nI];

		if( tree[nI].children.size() == 0 )
			continue;	// hop concept doesn't apply to leaf nodes
		bitset_t not_c1_comp = c1_complete[nI];
		not_c1_comp.flip();
		bitset_t not_c2_comp = c2_complete[nI];
		not_c2_comp.flip();
		hop_two[nI] = c1_partial[nI] & c2_partial[nI] & not_c1_comp & not_c2_comp;
		bitset_t not_hop_two_nI = hop_two[nI];
		not_hop_two_nI.flip();
		bitset_t not_child_hop = hop_one[ tree[nI].children[0] ] | hop_one[ tree[nI].children[1] ];
		not_child_hop.flip();
		hop_one[nI] = not_hop_two_nI & (not_c1_comp | not_c2_comp) & c1_partial[nI] & c2_partial[nI] & not_child_hop;

		// don't count small segments in anything
		hop_two[nI] &= not_small;
		hop_one[nI] &= not_small;
		unique[nI] &= not_small;
		conserved[nI] &= not_small;
		partial[nI] &= not_small;
		all_unique[nI] = unique[nI] | all_unique[ tree[nI].children[0] ] | all_unique[ tree[nI].children[1] ];
	}

	// compute length statistics for various types of backbone
	vector< double > conserved_len( tree.size(), 0 );
	vector< double > unique_len( tree.size(), 0 );
	vector< double > hop_one_len( tree.size(), 0 );
	vector< double > hop_two_len( tree.size(), 0 );

	for( size_t nI = 0; nI < tree.size(); nI++ )
	{
		// count up avg lengths
		for( size_t bbI = 0; bbI < bb_list.size(); ++bbI )
		{
			if( conserved[nI].test(bbI) )
				conserved_len[nI] += avg_lens[bbI];
			if( unique[nI].test(bbI) )
				unique_len[nI] += avg_lens[bbI];
			if( hop_one[nI].test(bbI) )
				hop_one_len[nI] += avg_lens[bbI];
			if( hop_two[nI].test(bbI) )
				hop_two_len[nI] += avg_lens[bbI];
		}
	}

	// print a general summary of how clustered variable segments are...
	bitset_t uni_root = unique[0] & not_small;
	anal_output << "There are " << uni_root.count() << " segments conserved among all genomes\n";
	anal_output << "and " << not_small.count()-uni_root.count() << " variable segments fall in between these\n";


	// prepare to analyze distribution of gene functions in backbone
	vector< bb_seqentry_t > m_bbseq_list( bb_list.size() );
	for( size_t bbI = 0; bbI < bb_list.size(); ++bbI )
		m_bbseq_list[bbI] = bb_list[bbI].bb_seq;
	multifun_map_t all_mf_count;
	multifun_names_t all_mf_names;
	size_t cds_count = getCDScount( iv_list.seq_table[gff_seq_index] );
	bitset_t all_features( iv_list.seq_table[gff_seq_index]->getFeatureListLength() );
	all_features.flip();
	makeMultiFunCount( iv_list.seq_table[gff_seq_index], all_mf_count, all_mf_names, all_features );

	// print summaries for each node
	anal_output << "#\n";
	anal_output << "# Alignment tree summary\n";
	anal_output << "#\n";
	for( size_t nI = 0; nI < tree.size(); nI++ )
	{
		anal_output << "Node " << nI << endl;
		vector< node_id_t > leaves;
		getLeaves( tree, nI, leaves );
		anal_output << "Genomes at or below this node:\n";
		for( size_t lI = 0; lI < leaves.size(); ++lI )
			anal_output << '\t' << iv_list.seq_filename[ node_sequence_map[ leaves[ lI ] ] ] << endl;

		anal_output << "\tNumber of unique segments at this node: " << unique[nI].count() << endl; 
		anal_output << "\tNumber of hop one (single deletion) segments at this node: " << hop_one[nI].count() << endl; 
		anal_output << "\tNumber of hop two (multiple deletion or lgt) segments at this node: " << hop_two[nI].count() << endl; 

		anal_output << "total avg. \"core-genome\" size at this node: " << conserved_len[nI] << endl;
		anal_output << "total avg. unique length at this node: " << unique_len[nI] << endl;
		anal_output << "total avg. hop one length at this node: " << hop_one_len[nI] << endl;
		anal_output << "total avg. hop two length at this node: " << hop_two_len[nI] << endl;
		anal_output << "total \"pan-genome\" size at this node: " << pan_genome_size[nI] << endl;

		// if this node has the annotated genome below it then analyze the distribution of
		// backbone content
		vector< uint > leaf_seqids( leaves.size() );
		for( size_t i = 0; i < leaves.size(); ++i )
			leaf_seqids[i] = node_sequence_map[leaves[i]];
		vector< uint >::iterator id_iter =std::find( leaf_seqids.begin(), leaf_seqids.end(), gff_seq_index );
		if( id_iter != leaf_seqids.end() )
		{
			vector< vector< size_t > > intersecting;
			featureIntersect( m_bbseq_list, gff_seq_index, intersecting, iv_list.seq_table[gff_seq_index] );
			bitset_t features_hit;
			getFeatureHits( intersecting, conserved[nI], features_hit );
			multifun_map_t bb_mf_count;
			multifun_names_t bb_mf_names;
			double expect_freq = (double)features_hit.count() / (double)cds_count;
			makeMultiFunCount( iv_list.seq_table[gff_seq_index], bb_mf_count, bb_mf_names, features_hit );
			anal_output << "#\n#Conserved gene content distribution\n#\n";
			anal_output << "Avg percent conserved " << setprecision(3) << expect_freq * 100 << endl;
			mfAnalyze( anal_output, all_mf_count, bb_mf_count, all_mf_names, expect_freq );	

			// analyze hop_one distributions
			intersecting.clear();
			featureIntersect( m_bbseq_list, gff_seq_index, intersecting, iv_list.seq_table[gff_seq_index] );
			features_hit.clear();
			getFeatureHits( intersecting, hop_one[nI], features_hit );
			bb_mf_count.clear();
			bb_mf_names.clear();
			expect_freq = (double)features_hit.count() / (double)cds_count;
			makeMultiFunCount( iv_list.seq_table[gff_seq_index], bb_mf_count, bb_mf_names, features_hit );
			anal_output << "#\n#Hop one gene content distribution\n#\n";
			anal_output << "Avg percent in hop_one " << setprecision(3) << expect_freq * 100 << endl;
			mfAnalyze( anal_output, all_mf_count, bb_mf_count, all_mf_names, expect_freq );


			// analyze hop_two distributions
			intersecting.clear();
			featureIntersect( m_bbseq_list, gff_seq_index, intersecting, iv_list.seq_table[gff_seq_index] );
			features_hit.clear();
			getFeatureHits( intersecting, hop_two[nI], features_hit );
			bb_mf_count.clear();
			bb_mf_names.clear();
			expect_freq = (double)features_hit.count() / (double)cds_count;
			makeMultiFunCount( iv_list.seq_table[gff_seq_index], bb_mf_count, bb_mf_names, features_hit );
			anal_output << "#\n#Hop two gene content distribution\n#\n";
			anal_output << "Avg percent in hop_two " << setprecision(3) << expect_freq * 100 << endl;
			mfAnalyze( anal_output, all_mf_count, bb_mf_count, all_mf_names, expect_freq );


			// analyze distributions of segments unique to this clade
			intersecting.clear();
			featureIntersect( m_bbseq_list, gff_seq_index, intersecting, iv_list.seq_table[gff_seq_index] );
			features_hit.clear();
			getFeatureHits( intersecting, all_unique[nI], features_hit );
			bb_mf_count.clear();
			bb_mf_names.clear();
			expect_freq = (double)features_hit.count() / (double)cds_count;
			makeMultiFunCount( iv_list.seq_table[gff_seq_index], bb_mf_count, bb_mf_names, features_hit );
			anal_output << "#\n#Unique to this clade gene content distribution\n#\n";
			anal_output << "Avg percent in unique_to_clade " << setprecision(3) << expect_freq * 100 << endl;
			mfAnalyze( anal_output, all_mf_count, bb_mf_count, all_mf_names, expect_freq );

		}
	}

	// first analyze all variable segments
	analyzeVariableSegments( anal_output, bb_list, avg_lens, gff_seq_index, iv_list.seq_table, "variable segments", false );

	// then analyze "alternalogs": variable segments with at least two non-null alleles
	analyzeVariableSegments( anal_output, bb_list, avg_lens, gff_seq_index, iv_list.seq_table, "alternalogs", true );
	anal_output.flush();
}

