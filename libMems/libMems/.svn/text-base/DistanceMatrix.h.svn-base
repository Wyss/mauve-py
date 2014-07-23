/*******************************************************************************
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __DistanceMatrix_h__
#define __DistanceMatrix_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/SubstitutionMatrix.h"
#include "libMems/IntervalList.h"
#include "libMems/MatchList.h"
#include "libMems/GappedAlignment.h"
#include "libMems/CompactGappedAlignment.h"


namespace mems {


void TransformDistanceIdentity( NumericMatrix<double>& identity );

void DistanceMatrix( const MatchList& mlist, NumericMatrix<double>& identity );


template< class AbstractMatchVectorType >
void IdentityMatrix( const AbstractMatchVectorType& matches, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity );
template<class AbstractMatchType>
void MatchIdentityMatrix( const AbstractMatchType& amt, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity);

void DistanceMatrix( uint seq_count, const std::vector< std::pair< uint64, uint64 > >& detail_list, NumericMatrix<double>& distance );

void IdentityMatrix( const IntervalList& iv_list, NumericMatrix<double>& identity );
inline
void IdentityMatrix( const IntervalList& iv_list, NumericMatrix<double>& identity )
{
	std::vector< const AbstractMatch* > am_list;
	for( size_t ivI = 0; ivI < iv_list.size(); ivI++ )
		am_list.push_back( &iv_list[ivI] );
	IdentityMatrix( am_list, iv_list.seq_table, identity );
}

template< class AbstractMatchVectorType >
void IdentityMatrix( const AbstractMatchVectorType& matches, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity ){
	if( matches.size() == 0 )
		return;

	uint seq_count = seq_table.size();
	identity = NumericMatrix<double>( seq_count, seq_count );
	identity.init( 0 );
	NumericMatrix<double> possible( seq_count, seq_count );
	possible.init( 0 );
	
	for( uint ivI = 0; ivI < matches.size(); ivI++ ){
		AddToMatchIdentityMatrix( *matches[ ivI ], seq_table, identity );
	}
	for( uint seqI = 0; seqI < seq_count; seqI++ ){
		for( uint seqJ = 0; seqJ < seq_count; seqJ++ ){
			gnSeqI shorter_len = seq_table[seqI]->length() < seq_table[seqJ]->length() ? seq_table[seqI]->length() : seq_table[seqJ]->length();
			possible( seqI, seqJ ) += shorter_len;
		}
	}
	identity /= possible;
}


template< class AbstractMatchVectorType >
void BackboneIdentityMatrix( const AbstractMatchVectorType& matches, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity ){
	if( matches.size() == 0 )
		return;

	size_t seq_count = seq_table.size();
	identity = NumericMatrix<double>( seq_count, seq_count );
	identity.init( 0 );
	
	for( uint ivI = 0; ivI < matches.size(); ivI++ ){
		AddToMatchIdentityMatrix( *matches[ ivI ], seq_table, identity );
	}

	NumericMatrix<double> possible( seq_count, seq_count );
	possible.init( 0 );

	for( size_t mI = 0; mI < matches.size(); ++mI ){
		std::vector< std::string > alignment;
		GetAlignment( *(matches[mI]), seq_table, alignment );
		for( gnSeqI charI = 0; charI < matches[mI]->AlignmentLength(); charI++ ){
			for( size_t seqI = 0; seqI < seq_count; seqI++ ){
				for( size_t seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
					if( alignment[ seqI ][ charI ] != '-' &&
						alignment[ seqJ ][ charI ] != '-' ){
							possible( seqI, seqJ ) += 1;
					}
				}
			}
		}
	}

	identity /= possible;
}


template<class AbstractMatchType>
void MatchIdentityMatrix( const AbstractMatchType& amt, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity)
{
	if( amt.SeqCount() == 0 )
		return;
	uint seq_count = amt.SeqCount();
	identity = NumericMatrix<double>( seq_count, seq_count );
	identity.init( 0 );
	uint seqI;
	uint seqJ;

	std::vector< std::string > alignment;
	GetAlignment( amt, seq_table, alignment );
	for( gnSeqI charI = 0; charI < amt.AlignmentLength(); charI++ ){
		for( seqI = 0; seqI < seq_count; seqI++ ){
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
				if( ( toupper( alignment[ seqI ][ charI ] ) == 
					toupper( alignment[ seqJ ][ charI ] ) ) &&
					alignment[ seqI ][ charI ] != '-' ){
					
						identity( seqI, seqJ ) += 1;
				}
			}
		}
	}

	for( seqI = 0; seqI < seq_count; seqI++ ){
		for( seqJ = seq_count; seqJ > 0; seqJ-- ){
			if( seqI == seqJ - 1 )
				// set the diagonal to identical
				identity( seqI, seqJ - 1 ) = 1;
			else if( seqI < seqJ - 1 ){
				// determine the length of the shorter sequence
				gnSeqI shorter_len = amt.Length( seqI ) < amt.Length( seqJ - 1 ) ? amt.Length( seqI ) : amt.Length( seqJ - 1 );
				// divide through
				identity( seqI, seqJ - 1 ) /= (double)shorter_len;
				// maxes out at 1
				if( identity( seqI, seqJ - 1 ) > 1 )
					identity( seqI, seqJ - 1 ) = 1;
			}else	// copy the other one
				identity( seqI, seqJ - 1 ) = identity( seqJ - 1, seqI );
		}
	}
}



template<class AbstractMatchType>
void AddToMatchIdentityMatrix( const AbstractMatchType& amt, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity)
{
	if( amt.SeqCount() == 0 )
		return;
	uint seq_count = amt.SeqCount();
	uint seqI;
	uint seqJ;

	std::vector< std::string > alignment;
	GetAlignment( amt, seq_table, alignment );
	for( gnSeqI charI = 0; charI < amt.AlignmentLength(); charI++ ){
		for( seqI = 0; seqI < seq_count; seqI++ ){
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
				if( ( toupper( alignment[ seqI ][ charI ] ) == 
					toupper( alignment[ seqJ ][ charI ] ) ) &&
					alignment[ seqI ][ charI ] != '-' ){
					
						identity( seqI, seqJ ) += 1;
				}
			}
		}
	}
}

/*
// template specialization for (exact) matches
inline
void AddToMatchIdentityMatrix( const Match& m, const std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& identity)
{
	if( m.SeqCount() == 0 )
		return;
	for( uint seqI = 0; seqI < m.SeqCount(); seqI++ )
		if( m.LeftEnd(seqI) != NO_MATCH )
			for( uint seqJ = seqI + 1; seqJ < m.SeqCount(); seqJ++ )
				if( m.LeftEnd(seqJ) != NO_MATCH )
					identity(seqI,seqJ) += m.Length();
}
*/

template< typename MatchVector >
void SingleCopyDistanceMatrix( MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, NumericMatrix<double>& distance )
{
	uint seq_count = seq_table.size();
	distance = NumericMatrix<double>( seq_count, seq_count );
	distance.init( 0 );
	uint seqI;
	uint seqJ;
	std::vector< std::pair< bitset_t, bitset_t > > tmp_comp( seq_count );
	std::vector< std::vector< std::pair< bitset_t, bitset_t > > > pair_comp( seq_count, tmp_comp );
	for( uint seqI = 0; seqI < seq_count; ++seqI )
	{
		for( uint seqJ = seqI+1; seqJ < seq_count; ++seqJ )
		{
			pair_comp[seqI][seqJ].first.resize( seq_table[seqI]->length(), false );
			pair_comp[seqI][seqJ].second.resize( seq_table[seqJ]->length(), false );
		}
	}
#pragma omp parallel for
	for( int ivI = 0; ivI < iv_list.size(); ++ivI )
	{
		std::vector< bitset_t > aln_table;
#pragma omp critical
{
		iv_list[ivI]->GetAlignment(aln_table);
}
		for( uint seqI = 0; seqI < seq_count; ++seqI )
		{
			for( uint seqJ = seqI+1; seqJ < seq_count; ++seqJ )
			{
				gnSeqI seqI_pos = iv_list[ivI]->LeftEnd(seqI);
				gnSeqI seqJ_pos = iv_list[ivI]->LeftEnd(seqJ);
				AbstractMatch::orientation o_i = iv_list[ivI]->Orientation(seqI);
				AbstractMatch::orientation o_j = iv_list[ivI]->Orientation(seqJ);
				if( o_i == AbstractMatch::reverse )
					seqI_pos = iv_list[ivI]->RightEnd(seqI);
				if( o_j == AbstractMatch::reverse )
					seqJ_pos = iv_list[ivI]->RightEnd(seqJ);
				if( seqI_pos == NO_MATCH || seqJ_pos == NO_MATCH )
					continue;
				for( size_t colI = 0; colI < aln_table[seqI].size(); ++colI )
				{
					if( aln_table[seqI].test(colI) && aln_table[seqJ].test(colI) )
					{
						pair_comp[seqI][seqJ].first.set(seqI_pos-1,true);
						pair_comp[seqI][seqJ].second.set(seqJ_pos-1,true);
					}
					if( aln_table[seqI].test(colI) )
						if( o_i == AbstractMatch::forward )
							seqI_pos++;
						else
							seqI_pos--;
					if( aln_table[seqJ].test(colI) )
						if( o_j == AbstractMatch::forward )
							seqJ_pos++;
						else
							seqJ_pos--;
				}
			}
		}
	}
	for( uint seqI = 0; seqI < seq_count; ++seqI )
	{
		distance(seqI,seqI) = 1;
		for( uint seqJ = seqI+1; seqJ < seq_count; ++seqJ )
		{
			double pI = ((double)pair_comp[seqI][seqJ].first.count())/((double)pair_comp[seqI][seqJ].first.size());
			double pJ = ((double)pair_comp[seqI][seqJ].second.count())/((double)pair_comp[seqI][seqJ].second.size());
			distance(seqI,seqJ) = (pI + pJ) / 2.0;
			distance(seqJ,seqI) = (pI + pJ) / 2.0;
		}
	}
	TransformDistanceIdentity(distance);
}

inline
void DistanceMatrix( const MatchList& mlist, NumericMatrix<double>& distance ){
	IdentityMatrix(mlist, mlist.seq_table, distance );
	TransformDistanceIdentity( distance );
}

inline
void TransformDistanceIdentity( NumericMatrix<double>& identity ){
	for( int i = 0; i < identity.cols(); i++ ){
		for( int j = 0; j < identity.rows(); j++ ){
			identity( i, j ) = 1 - identity( i, j );
		}
	}
}

inline
void DistanceMatrix( uint seq_count, const std::vector< std::pair< uint64, uint64 > >& detail_list, NumericMatrix<double>& distance ){
	distance = NumericMatrix<double>( seq_count, seq_count );
	distance.init( 0 );
	uint seqI;
	uint seqJ;
	for( seqI = 0; seqI < seq_count; seqI++ ){
		uint64 seqI_mask = 1;
		seqI_mask <<= seq_count - seqI - 1;
		for( seqJ = 0; seqJ < seq_count; seqJ++ ){
			uint64 seqJ_mask = 1;
			seqJ_mask <<= seq_count - seqJ - 1;
			for( uint pairI = 0; pairI < detail_list.size(); pairI++ ){
				if( (detail_list[ pairI ].first & seqI_mask) != 0 &&
					(detail_list[ pairI ].first & seqJ_mask) != 0 ){
					distance( seqI, seqJ ) += detail_list[ pairI ].second;
				}
			}
		}
	}
	
	for( seqI = 0; seqI < seq_count; seqI++ ){
		for( seqJ = 0; seqJ < seq_count; seqJ++ ){
			if( seqI == seqJ )
				continue;
			double avg_length = ( distance( seqI, seqI ) + distance( seqJ, seqJ ) ) / 2;
			distance( seqI, seqJ ) = 1.0 - ( distance( seqI, seqJ ) / avg_length );
			if( !(distance( seqI, seqJ ) == distance( seqI, seqJ )) ){
				distance( seqI, seqJ ) = 1.0;
			}
		}
	}

	// set the diagonal identical to itself
	for( seqI = 0; seqI < seq_count; seqI++ )
		distance( seqI, seqI ) = 0;
}


}	// namespace mems


#endif	// __DistanceMatrix_h__

