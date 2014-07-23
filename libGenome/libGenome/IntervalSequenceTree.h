#ifndef __IntervalSequenceTree_h__
#define __IntervalSequenceTree_h__

#include <vector>

/*
class interval {
	gnSeqI length();
	cropEnd();
	cropStart();
}
*/

// select idea shamelessly ripped from Metrowerks
template <bool Condition, class If, class Then>
class type_select
{
public:
	typedef If type;
};

template <class If, class Then>
class type_select<false, If, Then>
{
public:
	typedef Then type;
};

static const uint64 IST_END = ((0xFFFFFFFF << 31) << 1 ) + 0xFFFFFFFF;

/**
 * class implementing an Interval Sequence Tree
 * this is a tree for storing a changing sequence of intervals 
 * that was invented one rainy afternoon by Aaron Darling and
 * Michael Rusch.  
 * Important features are that stabbing queries,
 * stabbing insertions, and stabbing deletions are O( log n )
 * assuming a uniform distribution of stab sites.
 */
template< class Key, class Allocator = std::allocator<Key> >
class IntervalSequenceTree
{
public:
	typedef Key value_type;
//	typedef unsigned long long size_type;

	typedef Allocator                             allocator_type;
	typedef typename Allocator::reference         reference;
	typedef typename Allocator::const_reference   const_reference;
	typedef typename Allocator::size_type         size_type;
	typedef typename Allocator::difference_type   difference_type;
	typedef typename Allocator::pointer           pointer;
	typedef typename Allocator::const_pointer     const_pointer;

// node types and iterator definitions
protected:
	/**
	 * This class represents nodes of an Interval Sequence Tree.  Internal
	 * nodes define any of left, center, and right to be non-null and key to
	 * be null.  Leaf nodes define left, center, and right as null and key points
	 * to an interval.  The length field in an internal node is always the sum
	 * of lengths of the leaf nodes in its subtree.  The subtree_size field is
	 * defined as the number of nodes (leaf and internal) below the node.
	 */
	class IstNode {
	public:			
		IstNode* parent;
		IstNode* left;
		IstNode* right;
		size_type subtree_size;
		size_type length;
		Key* key;

		IstNode() :
			parent( NULL ),
			left( NULL ),
			right( NULL ),
			subtree_size( 0 ),
			length( 0 ),
			key( NULL ) {}
	};
	typedef typename Allocator::template rebind<IstNode>::other node_allocator_type;
	typedef typename node_allocator_type::pointer node_pointer;
	typedef typename node_allocator_type::const_pointer const_node_pointer;

//	typedef typename IstNode* node_pointer;
//	typedef typename const IstNode* const_node_pointer;

public:	

	// generic bidirectional iterator interface ripped from MSL, thanks guys
	template <bool is_const>
	class __generic_iterator
	{
	public:
		typedef typename IntervalSequenceTree::value_type       value_type;
//		typedef typename IntervalSequenceTree::difference_type  difference_type;
		typedef typename type_select<is_const, typename IntervalSequenceTree::const_pointer,
		                                  typename IntervalSequenceTree::pointer>::type pointer;
		typedef typename type_select<is_const, typename IntervalSequenceTree::const_reference,
		                                  typename IntervalSequenceTree::reference>::type reference;
		typedef std::bidirectional_iterator_tag        iterator_category;
		
		__generic_iterator() {}
		__generic_iterator(const __generic_iterator<false>& i) : ptr_(i.ptr_) {}
		reference operator * () const {return ptr_->data_;}
		pointer operator -> () const  {return &ptr_->data_;}
		__generic_iterator& operator ++ () {increment((const IstNode*&)ptr_); return *this;}
		__generic_iterator operator ++ (int) {__generic_iterator tmp(*this); ++(*this); return tmp;}
		__generic_iterator& operator -- () {decrement((const IstNode*&)ptr_); return *this;}
		__generic_iterator operator -- (int) {__generic_iterator tmp(*this); --(*this); return tmp;}
		friend bool operator ==(const __generic_iterator& x, const __generic_iterator& y) {return x.ptr_ == y.ptr_;}
		friend bool operator !=(const __generic_iterator& x, const __generic_iterator& y) {return x.ptr_ != y.ptr_;}
	private:
		typedef typename type_select<is_const, typename IntervalSequenceTree::node_pointer,
		                                  typename IntervalSequenceTree::const_node_pointer>::type node_pointer;

		node_pointer ptr_;

		explicit __generic_iterator(node_pointer n) : ptr_(n) {}

		friend class __generic_iterator<true>;
		friend class IntervalSequenceTree;
	};

	friend class __generic_iterator<false>;
	friend class __generic_iterator<true>;
	typedef __generic_iterator<false> iterator;
	typedef __generic_iterator<true>  const_iterator;
	typedef std::reverse_iterator< iterator > reverse_iterator;
	typedef std::reverse_iterator< const_iterator > const_reverse_iterator;


// constructor related methods
	IntervalSequenceTree();
	template< class InputIterator >
	IntervalSequenceTree( InputIterator first, InputIterator last );
	IntervalSequenceTree( const IntervalSequenceTree& ist );
	IntervalSequenceTree& operator=( const IntervalSequenceTree& ist );
	~IntervalSequenceTree();

// standard container methods
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
	reverse_iterator rbegin();
	const_reverse_iterator rbegin() const;
	reverse_iterator rend();
	const_reverse_iterator rend() const;
	size_type max_size() const;
	bool empty() const;
	
// insertion and erasure
	iterator insert( const value_type& val, size_type point = IST_END );
	template <class InputIterator>
	void insert(InputIterator first, InputIterator last, size_type point = IST_END );

	size_type erase( size_type point, size_type length );
	void erase( iterator first, iterator last );

// search
	iterator find( size_type point );
	const_iterator find( size_type point ) const;

// interval sequence specific:
	/**
	 * Returns the total length of intervals contained in this interval sequence
	 */
	size_type length() const;
	size_type nodeCount() const;
	size_type countNodes( IstNode* x = NULL ) const;

protected:
	IstNode *root;		/**< Root of the tree */
	IstNode *leftmost;	/**< Left most tree node, for begin() method */
	IstNode *rightmost;	/**< Right most tree node, for end() method */
	
	static void propogateChanges( IstNode* cur_node, int64 length_diff, int64 subtree_diff );
	static IstNode* recursiveFind( size_type& point, IstNode* node );
	static void increment( IstNode*& x);
	void decrement( IstNode*& x) const;
	static void deleteSubtree( IstNode*& istn );
	static void checkTree( node_pointer cur_node );
};



//template< class Key, class Allocator >
//IntervalSequenceTree< Key, Allocator >::IST_END = -1;

template< class Key, class Allocator >
inline
IntervalSequenceTree< Key, Allocator >::IntervalSequenceTree(){
	root = NULL;
	leftmost = NULL;
	rightmost = NULL;
//	IST_END = -1;	// wraps to UINT64_MAX because IST_END is unsigned
}

template< class Key, class Allocator >
template< class InputIterator >
IntervalSequenceTree< Key, Allocator >::IntervalSequenceTree( InputIterator first, InputIterator last ){
	insert( first, last );
}

template< class Key, class Allocator >
IntervalSequenceTree< Key, Allocator >::IntervalSequenceTree( const IntervalSequenceTree& ist ){
	insert( ist.begin(), ist.end() );
}

template< class Key, class Allocator >
IntervalSequenceTree< Key, Allocator >& IntervalSequenceTree< Key, Allocator >::operator=( const IntervalSequenceTree& ist ){
	insert( ist.begin(), ist.end() );
}

template< class Key, class Allocator >
IntervalSequenceTree< Key, Allocator >::~IntervalSequenceTree(){
	deleteSubtree( root );
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::iterator 
IntervalSequenceTree< Key, Allocator >::begin(){
	return iterator( leftmost );
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::const_iterator 
IntervalSequenceTree< Key, Allocator >::begin() const{
	return const_iterator( leftmost );
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::iterator 
IntervalSequenceTree< Key, Allocator >::end(){
	return iterator( NULL );
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::const_iterator 
IntervalSequenceTree< Key, Allocator >::end() const{
	return const_iterator( NULL );
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::reverse_iterator 
IntervalSequenceTree< Key, Allocator >::rbegin(){
	return reverse_iterator( end() );
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::const_reverse_iterator 
IntervalSequenceTree< Key, Allocator >::rbegin() const{
	return const_reverse_iterator( end() );
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::reverse_iterator 
IntervalSequenceTree< Key, Allocator >::rend(){
	return reverse_iterator( begin() );
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::const_reverse_iterator 
IntervalSequenceTree< Key, Allocator >::rend() const{
	return const_reverse_iterator( begin() );
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::size_type 
IntervalSequenceTree< Key, Allocator >::max_size() const{
	return IST_END - 1;
}

template< class Key, class Allocator >
bool IntervalSequenceTree< Key, Allocator >::empty() const{
	return root == NULL ? true : false;
}

template< class Key, class Allocator >
void IntervalSequenceTree< Key, Allocator >::checkTree( 
//	IntervalSequenceTree< Key, Allocator >::IstNode*
	node_pointer cur_node ){
	if( cur_node ){
		if( cur_node->parent == cur_node )
			std::cerr << "freakout\n";
		checkTree( cur_node->left );
		checkTree( cur_node->right );
	}
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::iterator 
IntervalSequenceTree< Key, Allocator >::insert( 
	const Key& val, 
	typename IntervalSequenceTree< Key, Allocator >::size_type point )
{
	size_type iv_offset = point;
	IstNode* ins_node = recursiveFind( iv_offset, root );
	IstNode* new_node = new IstNode();
	new_node->key = new Key( val );
	new_node->length = val.GetLength();
	new_node->subtree_size = 0;
	if( ins_node == NULL ){
		// end insert
		rightmost = new_node;
		if( root == NULL ){
			root = new_node;
			leftmost = new_node;
			return iterator( new_node );
		}
		// find the shallowest right insertion point
		ins_node = NULL;
		decrement( ins_node );
		// make a new parent node
		IstNode* new_parent = new IstNode();
		new_parent->left = ins_node;
		new_parent->right = new_node;
		new_parent->parent = ins_node->parent;
		if( new_parent->parent == NULL )
			root = new_parent;
		else
			new_parent->parent->right = new_parent;
		ins_node->parent = new_parent;
		new_parent->length = ins_node->length;

		// update lengths and subtree_sizes along the path to the root
//		checkTree( root );
//		propogateChanges( new_node, 0, 0 );
//		propogateChanges( ins_node, 0, 0 );
		propogateChanges( new_parent, new_node->length, 2 );
		return iterator( new_node );
	}

	// iv_offset is the distance into the node that the leaf should be split
	// 0 is a special case (left insert)
	if( iv_offset == 0 ){
		IstNode* new_parent = new IstNode();
		new_parent->left = new_node;
		new_parent->right = ins_node;
		new_parent->parent = ins_node->parent;
		if( new_parent->parent->right == ins_node )
			new_parent->parent->right = new_parent;
		else
			new_parent->parent->left = new_parent;
		new_parent->length = ins_node->length;

		ins_node->parent = new_parent;
		new_node->parent = new_parent;

		if( point == 0 )
			leftmost = new_node;
		// update lengths and subtree_sizes along the path to the root
//		checkTree( root );
//		propogateChanges( new_node, 0, 0 );
//		propogateChanges( ins_node, 0, 0 );
		propogateChanges( new_parent, new_node->length, 2 );
	}else{
		// need to split a leaf node
		IstNode* new_gp = new IstNode();
		IstNode* new_parent = new IstNode();
		new_gp->parent = ins_node->parent;
		new_gp->right = new_parent;
		new_gp->left = new IstNode();
		new_gp->left->key = new Key( *ins_node->key );
		new_gp->left->key->CropEnd( ins_node->length - iv_offset );
		new_gp->left->length = new_gp->left->key->GetLength();
		new_gp->left->parent = new_gp;
		
		ins_node->key->CropStart( iv_offset );
		ins_node->length = ins_node->key->GetLength();
		ins_node->parent = new_parent;
		new_node->parent = new_parent;
		new_parent->left = new_node;
		new_parent->right = ins_node;
		new_parent->parent = new_gp;
		new_parent->length = new_node->length + ins_node->length;
		new_parent->subtree_size = 2;

		new_gp->length = ins_node->length + new_gp->left->length;
		new_gp->subtree_size = 1;
		if( new_gp->parent == NULL ){
			root = new_gp;
			leftmost = new_gp->left;
			rightmost = ins_node;
		}else if( new_gp->parent->right == ins_node )
			new_gp->parent->right = new_gp;
		else
			new_gp->parent->left = new_gp;

		// update lengths and subtree_sizes along the path to the root
//		checkTree( root );
//		propogateChanges( new_node, 0, 0 );
//		propogateChanges( ins_node, 0, 0 );
//		propogateChanges( new_gp->left, 0, 0 );
//		propogateChanges( new_parent, 0, 0 );
		new_gp->subtree_size = -1;
		propogateChanges( new_gp, new_node->length, 4 );
	}
	return iterator( new_node );
}

template< class Key, class Allocator >
template <class InputIterator>
void IntervalSequenceTree< Key, Allocator >::insert( 
	InputIterator first, 
	InputIterator last, 
	typename IntervalSequenceTree< Key, Allocator >::size_type point )
{

}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::size_type 
IntervalSequenceTree< Key, Allocator >::erase( 
	typename IntervalSequenceTree< Key, Allocator >::size_type point, 
	typename IntervalSequenceTree< Key, Allocator >::size_type length )
{
	size_type iv_offset = point;
	IstNode* ins_node = recursiveFind( iv_offset, root );

	// iv_offset is the distance into the node that the leaf should be split
	// 0 is a special case (left delete)
	size_type deleted_nodes = 0;
	while( length > 0 ){
		if( ins_node == NULL ){
			// end delete?  that's illegal
			return deleted_nodes;
		}
		if( iv_offset == 0 ){
			if( length >= ins_node->length ){
				// delete the whole thing
				length -= ins_node->length;
				if( ins_node->parent == NULL ){
					// deleting the root
					delete ins_node;
					root = NULL;
					leftmost = NULL;
					rightmost = NULL;
					return deleted_nodes + 1;
				}

				IstNode* other_child, *del_node;
				if( ins_node->parent->left == ins_node ){
					other_child = ins_node->parent->right;
				}else if( ins_node->parent->right == ins_node ){
					other_child = ins_node->parent->left;
				}
				del_node = ins_node;
				increment( ins_node );

				// update tree structure
				IstNode* tmp_parent = other_child->parent;
				IstNode* tmp_gp = tmp_parent->parent;
				*tmp_parent = *other_child;
				tmp_parent->parent = tmp_gp;
				if( tmp_parent->left )
					tmp_parent->left->parent = tmp_parent;
				if( tmp_parent->right )
					tmp_parent->right->parent = tmp_parent;
				if( ins_node == other_child )
					ins_node = tmp_parent;
				delete other_child;
				
				// propogate deletion length thru root
				tmp_parent = tmp_parent->parent;
//				checkTree( root );
				propogateChanges( tmp_parent, -del_node->length, -2 );

				// finally delete ins_node
				delete del_node;
				++deleted_nodes;
			}else{
				// crop from start
				ins_node->key->CropStart( length );
//				checkTree( root );
				propogateChanges( ins_node, -length, 0 );
				return deleted_nodes;
			}
		}else if( length >= ins_node->length - iv_offset ){
			// crop from end
			ins_node->key->CropEnd( ins_node->length - iv_offset );
			length -= ins_node->length - iv_offset;
//			checkTree( root );
			propogateChanges( ins_node, -(ins_node->length - iv_offset), 0 );
			increment( ins_node );
			iv_offset = 0;
		}else{
			// delete from middle (nastee part)
			IstNode* new_parent = new IstNode();
			new_parent->left = ins_node;
			new_parent->length = ins_node->length;
			new_parent->right = new IstNode();
			new_parent->right->key = new Key( *ins_node->key );
			new_parent->right->length = ins_node->length - length - iv_offset;
			new_parent->right->key->CropStart( length + iv_offset );
			new_parent->left->key->CropEnd( ins_node->length - iv_offset );
			new_parent->left->length = iv_offset;
			new_parent->parent = ins_node->parent;
			if( new_parent->parent == NULL ){
				root = new_parent;
				rightmost = new_parent->right;
			}else if( new_parent->parent->left == ins_node )
				new_parent->parent->left = new_parent;
			else if( new_parent->parent->right == ins_node )
				new_parent->parent->right = new_parent;
			
			ins_node->parent = new_parent;
			new_parent->right->parent = new_parent;
//			checkTree( root );
//			propogateChanges( ins_node, 0, 0 );
//			propogateChanges( new_parent->right, 0, 0 );
			propogateChanges( new_parent, -length, 2 );
			return deleted_nodes;
		}
	}
	return deleted_nodes;
}

template< class Key, class Allocator >
void IntervalSequenceTree< Key, Allocator >::propogateChanges( 
	IstNode* cur_node,
	int64 length_diff, 
	int64 subtree_diff )
{
	std::vector< IstNode* > node_stack;
	while( cur_node != NULL ){
		if( cur_node->parent == cur_node )
			std::cerr << "when I say oh, you say shit!";
		cur_node->length += length_diff;
		cur_node->subtree_size += subtree_diff;
		node_stack.push_back( cur_node );
		cur_node = cur_node->parent;
	}
}

template< class Key, class Allocator >
void IntervalSequenceTree< Key, Allocator >::erase( 
	typename IntervalSequenceTree< Key, Allocator >::iterator first, 
	typename IntervalSequenceTree< Key, Allocator >::iterator last )
{

}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::iterator 
IntervalSequenceTree< Key, Allocator >::find( 
	typename IntervalSequenceTree< Key, Allocator >::size_type point ) 
{
	return iterator( IntervalSequenceTree< Key, Allocator >::recursiveFind( point, root ) );
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::const_iterator 
IntervalSequenceTree< Key, Allocator >::find( 
	typename IntervalSequenceTree< Key, Allocator >::size_type point ) const
{
	return const_iterator( IntervalSequenceTree< Key, Allocator >::recursiveFind( point, root ) );
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::IstNode* 
IntervalSequenceTree< Key, Allocator >::recursiveFind( 
	typename IntervalSequenceTree< Key, Allocator >::size_type& point, 
	IstNode* node ) {

	if( node == NULL )
		return NULL;

	// return this node if it's a leaf
	if( node->key != NULL )
		return node;
	// look for the next node to recurse to
	if( point < node->length ){
		if( node->left ){
			if( point < node->left->length )
				return recursiveFind( point, node->left );
			point -= node->left->length;
		}
		return recursiveFind( point, node->right );
	}
	point -= node->length;
	// out of range
	return NULL;
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::size_type 
IntervalSequenceTree< Key, Allocator >::length() const{
	return root == NULL ? 0 : root->length;
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::size_type 
IntervalSequenceTree< Key, Allocator >::nodeCount() const{
	return root == NULL ? 0 : root->subtree_size + 1;
}

template< class Key, class Allocator >
typename IntervalSequenceTree< Key, Allocator >::size_type 
IntervalSequenceTree< Key, Allocator >::countNodes( IstNode* x ) const{
	if( x == NULL )
		x = root;
	if( x->key == NULL )
		return countNodes( x->left ) + countNodes( x->right ) + 1;
	return 1;
}


template< class Key, class Allocator >
void IntervalSequenceTree< Key, Allocator >::increment( IstNode*& x) {
	// find the least-ancestor with another child
	// and set x to that child
	while( x->parent != NULL ){
		if( x == x->parent->left &&
			x->parent->right != NULL ){
			x = x->parent->right;
			break;
		}else
			x = x->parent;
	}
	// if there were no other children to the right then we're at the end
	if( x->parent == NULL ){
		x = NULL;
		return;
	}

	// find the left-most leaf node below x
	while( x->key == NULL ){
		if( x->left != NULL )
			x = x->left;
		else if( x->right != NULL )
			x = x->right;
	}
}

template< class Key, class Allocator >
void IntervalSequenceTree< Key, Allocator >::decrement( IstNode*& x) const{
	if( x != NULL ){
		// find the least-ancestor with another child to the left
		// and set x to that child
		while( x->parent != NULL ){
			if( x == x->parent->right &&
				x->parent->left != NULL){
				x = x->parent->left;
				break;
			}else
				x = x->parent;
		}
		// if there was no other children to the left then we're at the end
		// raise hell! (cause an access violation)
		if( x->parent == NULL )
			x = NULL;
	}else{
		x = root;
	}
	
	// find the right-most leaf node below x
	while( x->key == NULL ){
		if( x->right != NULL )
			x = x->right;
		else if( x->left != NULL )
			x = x->left;
	}
}

template< class Key, class Allocator >
void IntervalSequenceTree< Key, Allocator >::deleteSubtree( IstNode*& istn ) {
	if( istn->left != NULL )
		deleteSubtree( istn->left );
	if( istn->right != NULL )
		deleteSubtree( istn->right );
	if( istn->key != NULL )
		delete istn->key;
	delete istn;
}


#endif // __IntervalSequenceTree_h__
