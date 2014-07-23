#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "libMems/SuperInterval.h"

using namespace std;
using namespace genome;

namespace mems {
// working in mems

bool debug_aligner = false;

SuperInterval::SuperInterval() :
length(0),
left_end(0),
c1_siv((std::numeric_limits<size_t>::max)()),
c2_siv((std::numeric_limits<size_t>::max)()),
parent_siv((std::numeric_limits<size_t>::max)())
{}

SuperInterval::SuperInterval( const Interval& reference_iv ) :
reference_iv(reference_iv),
length(0),
left_end(0),
c1_siv((std::numeric_limits<size_t>::max)()),
c2_siv((std::numeric_limits<size_t>::max)()),
parent_siv((std::numeric_limits<size_t>::max)())
{
}

SuperInterval::SuperInterval(const SuperInterval& siv) :
left_end(siv.left_end),
length( siv.length ),
reference_iv( siv.reference_iv ),
c1_siv(siv.c1_siv),
c2_siv(siv.c2_siv),
parent_siv(siv.parent_siv)
{
}
SuperInterval& SuperInterval::operator=(const SuperInterval& siv)
{
	left_end = siv.left_end;
	length = siv.length;
	reference_iv = siv.reference_iv;
	c1_siv = siv.c1_siv;
	c2_siv = siv.c2_siv;
	parent_siv = siv.parent_siv;
	return *this;
}



/** Sets the length of this match to @param len */
void SuperInterval::SetLength( gnSeqI len )
{
	length = len;
}

void SuperInterval::CropLeft( gnSeqI amount )
{
	reference_iv.CropStart(amount);

	left_end += amount;
	length -= amount;

	if(debug_aligner)
		ValidateSelf();
}

void SuperInterval::CropRight( gnSeqI amount )
{
	reference_iv.CropEnd(amount);
	length -= amount;

	if(debug_aligner)
		ValidateSelf();
}

void SuperInterval::ValidateSelf() const
{
	vector< bitset_t > aln_mat;
	reference_iv.GetAlignment(aln_mat);
	if( aln_mat[0].size() != reference_iv.AlignmentLength() )
	{
		breakHere();
		cerr << "trouble! aln_mat[0].size() is: " << aln_mat[0].size() << " while reference_iv.AlignmentLength() is: " << reference_iv.AlignmentLength() << endl;
		cerr << "mult: " << reference_iv.Multiplicity() << endl;
		cerr << "matches.size(): " << reference_iv.GetMatches().size() << endl;
	}
	for( size_t i = 0; i < aln_mat.size(); i++ )
	{
		gnSeqI lenny = 0;
		for( size_t j = 0; j < aln_mat[i].size(); j++ )
			if( aln_mat[i][j] )
				lenny++;
		if( lenny != reference_iv.Length(i) )
		{
			cerr << "krudunkle, ref_iv.Length(" << i << "): " << reference_iv.Length(i) << "\n";
			cerr << "should be: " << lenny << endl;
			breakHere();
		}
	}
	if( reference_iv.LeftEnd(0) != NO_MATCH && reference_iv.Length(0) == 0 )
	{
		cerr << "brozooka\n";
		breakHere();
	}
	if( reference_iv.LeftEnd(1) != NO_MATCH && reference_iv.Length(1) == 0 )
	{
		cerr << "brokazooka\n";
		breakHere();
	}

	if( Length() != reference_iv.AlignmentLength() )
	{
		breakHere();
		cerr << "crapola\n";
	}
}

} // namespace mems
