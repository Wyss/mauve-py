#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __SuperInterval_h__
#define __SuperInterval_h__

#include "libMems/Interval.h"

namespace mems {

/**
 * A class that stores an alignment and coordinate mapping between collinear segments of an ancestral genome and two
 * descendant genomes.
 */
class SuperInterval
{
public:

	SuperInterval();
	/**
	 * Creates a new SuperInterval.
	 */
	SuperInterval( const mems::Interval& reference_iv );
	SuperInterval(const SuperInterval& siv);
	SuperInterval& operator=(const SuperInterval& siv);
	~SuperInterval(){}
	
	/** Returns the length */
	virtual gnSeqI Length() const { return length; }

	/** Sets the length to @param len */
	virtual void SetLength( gnSeqI len );

	virtual int64 LeftEnd() const { return left_end; }

	virtual void SetLeftEnd( const int64& left_end ) { this->left_end = left_end; }

	mems::Interval reference_iv;

	/** the index of the SuperInterval this is aligned to in c1 */
	size_t c1_siv;
	/** the index of the SuperInterval this is aligned to in c2 */
	size_t c2_siv;
	/** the index of the SuperInterval this is aligned to in the parent */
	size_t parent_siv;

	void CropLeft( gnSeqI amount );
	void CropRight( gnSeqI amount );

	bool operator<( const SuperInterval& si ) const{ return left_end < si.left_end; }

	void ValidateSelf() const;

	void swap( SuperInterval& other )
	{
		reference_iv.swap(other.reference_iv);
		std::swap(c1_siv, other.c1_siv);
		std::swap(c2_siv, other.c2_siv);
		std::swap(parent_siv, other.parent_siv);
		std::swap(left_end, other.left_end);
		std::swap(length, other.length);
	}

protected:
	int64 left_end;
	int64 length;
};


} // namespace mems

namespace std {
template<> inline
void swap( mems::SuperInterval& a, mems::SuperInterval& b )
{
	a.swap(b);
}
}

#endif //__SuperInterval_h__
