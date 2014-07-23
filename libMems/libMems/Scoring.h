/*******************************************************************************
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __Scoring_h__
#define __Scoring_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/SubstitutionMatrix.h"
#include <string>
#include <vector>

namespace mems {

static const score_t INVALID_SCORE = (std::numeric_limits<score_t>::max)();

//tjtaed: function to compute the SP column score, and cumulative SP score from an alignment
void computeSPScore( const std::vector<std::string>& alignment, const PairwiseScoringScheme& pss, std::vector<score_t>& scores, score_t& score );
//tjt: function to compute the consensus column score, consensus sequence, and cumulative consensus score from an alignment
void computeConsensusScore( const std::vector<std::string>& alignment, const PairwiseScoringScheme& pss, std::vector<score_t>& scores, std::string& consensus, score_t& score );
void computeMatchScores( const std::string& seq1, const std::string& seq2, const PairwiseScoringScheme& scoring, std::vector<score_t>& scores );
void computeGapScores( const std::string& seq1, const std::string& seq2, const PairwiseScoringScheme& scoring, std::vector<score_t>& scores );


//tjt: function to compute the consensus column score, consensus sequence, and cumulative consensus score from an alignment 
inline
void computeConsensusScore( const std::vector<std::string>& alignment, const PairwiseScoringScheme& pss, 
						   std::vector<score_t>& scores, std::string& consensus, score_t& score )
{

	consensus.clear();
	std::vector< std::vector< score_t > > allscores;

	scores.resize( alignment.at(0).size() );
	std::fill(scores.begin(), scores.end(), INVALID_SCORE);

	score =	INVALID_SCORE;

	std::vector< string > nucleotides;
	nucleotides.push_back(std::string(alignment.at(0).size(),'A'));
	nucleotides.push_back(std::string(alignment.at(0).size(),'G'));
	nucleotides.push_back(std::string(alignment.at(0).size(),'C'));
	nucleotides.push_back(std::string(alignment.at(0).size(),'T'));
	
	for( size_t i = 0; i < nucleotides.size(); i++)
	{
		//tjt: score alignment!
		//for each row in the alignment, compare to string of A,G,C,T and build consensus
		std::vector< score_t > consensus_scores(alignment.at(0).size(), 0);
		
		for( gnSeqI j = 0; j < alignment.size(); j++)
		{
			std::vector< score_t > tscores(alignment.at(0).size(), 0);
		
			computeMatchScores( alignment.at(j), nucleotides.at(i), pss, tscores );
			
			for( gnSeqI k = 0; k < alignment.at(j).size(); k++)
				if( tscores.at(k) != INVALID_SCORE )
					consensus_scores.at(k) += tscores.at(k);

			computeGapScores( alignment.at(j), nucleotides.at(i), pss, tscores );

			for( gnSeqI k = 0; k < alignment.at(j).size(); k++)
				if( tscores.at(k) != INVALID_SCORE )
					consensus_scores.at(k) += tscores.at(k);
			
		}
		allscores.push_back(consensus_scores);
	}
	
	//tjt: find maxvalue for each column
	// 0 = A, 1 = G, 2 = C, 3 = T
	
	std::vector< int > columnbp( alignment.at(0).size(), (std::numeric_limits<int>::min)());
	
	//for A,G,C,T
	for( size_t i = 0; i < nucleotides.size(); i++)
	{
		//for each column
		for( size_t j = 0; j < alignment.at(0).size(); j++)
		{
			if( allscores.at(i).at(j) == INVALID_SCORE )
				continue;
			if( i == 0  )
			{				
				scores.at(j) = allscores.at(i).at(j);
				columnbp.at(j) = 0;
			}
			else if (allscores.at(i).at(j) > scores.at(j))
			{
				scores.at(j) = allscores.at(i).at(j);
				columnbp.at(j) = i;
			}
		}
	}
	//update score with maxvalue from each column
	for( size_t j = 0; j < alignment.at(0).size(); j++)
	{
		if( scores.at(j) != INVALID_SCORE )
			score += scores.at(j);
		if (columnbp.at(j) == 0)
			consensus.append("A");
		else if (columnbp.at(j) == 1)
			consensus.append("G");
		else if (columnbp.at(j) == 2)
			consensus.append("C");
		else if (columnbp.at(j) == 3)
			consensus.append("T");
	
	}
}

inline
void computeMatchScores( const std::string& seq1, const std::string& seq2, 
						const PairwiseScoringScheme& scoring, std::vector<score_t>& scores )
{
	scores.resize( seq1.size() );
	std::fill(scores.begin(), scores.end(), INVALID_SCORE);
	const uint8* table = SortedMerList::BasicDNATable();

	for (unsigned uColIndex = 0; uColIndex < seq1.size(); ++uColIndex)
	{
		char c1 = seq1[uColIndex];
		char c2 = seq2[uColIndex];
		if( c1 == '-' || c2 == '-' )
			continue;
		unsigned uLetter1 = table[c1];
		unsigned uLetter2 = table[c2];

		score_t scoreMatch = scoring.matrix[uLetter1][uLetter2];
		scores[uColIndex] = scoreMatch;
	}
}

inline
void computeGapScores( const std::string& seq1, const std::string& seq2, const PairwiseScoringScheme& scoring, 
					  std::vector<score_t>& scores )
{
	scores.resize(seq1.size());

	bool bGapping1 = false;
	bool bGapping2 = false;
	score_t gap_open_score = scoring.gap_open;
	score_t gap_extend_score = scoring.gap_extend;
	score_t term_gap_score = gap_open_score;

	unsigned uColCount = seq1.size();
	unsigned uColStart = 0;
	bool bLeftTermGap = false;
	for (unsigned uColIndex = 0; uColIndex < seq1.size(); ++uColIndex)
	{
		bool bGap1 = seq1[uColIndex] == '-';
		bool bGap2 = seq2[uColIndex] == '-';
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bLeftTermGap = true;
			uColStart = uColIndex;
			break;
			}
		}

	unsigned uColEnd = uColCount - 1;
	bool bRightTermGap = false;
	for (int iColIndex = (int) uColCount - 1; iColIndex >= 0; --iColIndex)
		{
		bool bGap1 = seq1[iColIndex] == '-';
		bool bGap2 = seq2[iColIndex] == '-';
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bRightTermGap = true;
			uColEnd = (unsigned) iColIndex;
			break;
			}
		}

	unsigned gap_left_col = 0;
	score_t cur_gap_score = 0;
	for (unsigned uColIndex = uColStart; uColIndex <= uColEnd; ++uColIndex)
		{
		bool bGap1 = seq1[uColIndex] == '-';
		bool bGap2 = seq2[uColIndex] == '-';

		if (bGap1 && bGap2)
			continue;

		if (bGap1)
			{
			if (!bGapping1)
				{
				gap_left_col = uColIndex;
				if (uColIndex == uColStart)
					{
					cur_gap_score += term_gap_score;
				}else{
					cur_gap_score += gap_open_score;
					}
				bGapping1 = true;
				}
			else
				{
				cur_gap_score += gap_extend_score;
				}
			continue;
			}

		else if (bGap2)
			{
			if (!bGapping2)
				{
				gap_left_col = uColIndex;
				if (uColIndex == uColStart)
					{
					cur_gap_score += term_gap_score;
				}else{
					cur_gap_score += gap_open_score;
					}
				bGapping2 = true;
				}
			else
				{
				cur_gap_score += gap_extend_score;
				}
			continue;
			}

		if( (bGapping1 || bGapping2) )
		{
			score_t valid_cols = 0;
			for( unsigned uGapIndex = gap_left_col; uGapIndex < uColIndex; ++uGapIndex )
				if( seq1[uGapIndex] != '-' || seq2[uGapIndex] != '-' )
					valid_cols++;
			// spread the total gap penalty evenly across all columns
			score_t per_site_penalty = cur_gap_score / valid_cols;
			score_t extra = cur_gap_score - (per_site_penalty * valid_cols);
			for( unsigned uGapIndex = gap_left_col; uGapIndex < uColIndex; ++uGapIndex )
			{
				if( seq1[uGapIndex] == '-' && seq2[uGapIndex] == '-' )
					continue;
				if( scores[uGapIndex] != INVALID_SCORE )
				{
					genome::breakHere();
					cerr << "asdgohasdoghasodgh\n";
				}
				scores[uGapIndex] = per_site_penalty;
			}
			if( scores[gap_left_col] == INVALID_SCORE )
			{
				cerr << "crap!\n";
				genome::breakHere();
			}
			scores[gap_left_col] += extra;
			gap_left_col = (std::numeric_limits<unsigned>::max)();
			cur_gap_score = 0;
		}
		bGapping1 = false;
		bGapping2 = false;
		}

	if (bGapping1 || bGapping2)
		{
		cur_gap_score -= gap_open_score;
		cur_gap_score += term_gap_score;

		score_t valid_cols = 0;
		for( unsigned uGapIndex = gap_left_col; uGapIndex < uColCount; ++uGapIndex )
			if( seq1[uGapIndex] != '-' || seq2[uGapIndex] != '-' )
				valid_cols++;
		// spread the total gap penalty evenly across all columns
		score_t per_site_penalty = cur_gap_score / valid_cols;
		score_t extra = cur_gap_score - (per_site_penalty * valid_cols);
		for( unsigned uGapIndex = gap_left_col; uGapIndex < uColCount; ++uGapIndex )
		{
			if( seq1[uGapIndex] == '-' && seq2[uGapIndex] == '-' )
				continue;
			scores[uGapIndex] = per_site_penalty;
		}
		if( valid_cols > 0 )
		{
			if( scores[gap_left_col] == INVALID_SCORE )
			{
				cerr << "crap!\n";
				genome::breakHere();
			}
			scores[gap_left_col] += extra;
		}
	}
}

inline
void computeSPScore( const std::vector<string>& alignment, const PairwiseScoringScheme& pss, 
					std::vector<score_t>& scores, score_t& score )
{
	std::vector< score_t > cur_m_scores( alignment[0].size(), INVALID_SCORE );
	std::vector< score_t > cur_g_scores( alignment[0].size(), INVALID_SCORE );
	scores.resize(alignment[0].size());
	std::fill(scores.begin(), scores.end(), 0);
	score = 0;
	double w = 1;	// weight, to be determined later...
	for( size_t i = 0; i < alignment.size(); ++i )
	{
		for( size_t j = i+1; j < alignment.size(); ++j )
		{
			std::fill( cur_m_scores.begin(), cur_m_scores.end(), INVALID_SCORE );
			std::fill( cur_g_scores.begin(), cur_g_scores.end(), INVALID_SCORE );
			computeMatchScores( alignment.at(i), alignment.at(j), pss, cur_m_scores );
			computeGapScores( alignment.at(i), alignment.at(j), pss, cur_g_scores );
			for( size_t k = 0; k < cur_m_scores.size(); ++k )
			{
				score_t s = 0;
				if( cur_m_scores[k] != INVALID_SCORE )
					s += cur_m_scores[k];
				if( cur_g_scores[k] != INVALID_SCORE )
					s += cur_g_scores[k];
				scores[k] += (score_t)(w * (double)s);
			}
		}
	}
	for( size_t k = 0; k < scores.size(); ++k )
		score += scores[k];
}


}	// namespace mems


#endif	// __Scoring_h__

