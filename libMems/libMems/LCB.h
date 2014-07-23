#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __LCB_h__
#define __LCB_h__

#include <vector>
#include <libGenome/gnDefs.h>

namespace mems {

/** 
 * This class is used to track relationships between LCBs during the LCB determination process.
 */
class LCB{
public:
	LCB() : lcb_id(0), weight(0), to_be_deleted(false) {};
	std::vector< int64 > left_end;	/**< The left end position of the LCB in each sequence */
	std::vector< int64 > right_end;  /**< The right end position of the LCB in each sequence */
	std::vector< uint > left_adjacency;	/**< 'Pointers' (actually IDs) to the LCBs on the left in each sequence */
	std::vector< uint > right_adjacency;	/**< 'Pointers' (actually IDs) to the LCBs on the right in each sequence */
	int lcb_id;			/**< A numerical ID that can be assigned to this LCB */
	double weight;		/**< The weight (or coverage) of this LCB */
	bool to_be_deleted;	/**< set to true if this LCB is about to be deleted, but the deletion hasn't yet been processed */
};

/**
 * Compares LCBs.
 * Used by LCB construction algorithm 
 */
class LCBLeftComparator {
public:
	LCBLeftComparator( uint seq ) : m_seq(seq){};
	bool operator()(const LCB& a, const LCB& b) const{
		
		int64 a_start = a.left_end[ m_seq ], b_start = b.left_end[ m_seq ];
		if( a_start == NO_MATCH || b_start == NO_MATCH ){
			if( b_start != NO_MATCH )
				return true;
			return false;
		}
		if(a_start < 0)
			a_start = -a_start;
		if(b_start < 0)
			b_start = -b_start;

		int64 diff = a_start - b_start;
		return diff < 0;
	}
protected:
	uint m_seq;
private:
	LCBLeftComparator();
};

class LCBIDComparator {
public:
	bool operator()(const LCB& a, const LCB& b) const
	{
		return a.lcb_id < b.lcb_id;
	}
};


} // namespace mems


#endif  // __LCB_h__

