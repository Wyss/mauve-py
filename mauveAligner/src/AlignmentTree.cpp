#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "PhyloTree.h"
#include <sstream>
#include <stack>
using namespace std;

typedef unsigned uint;

PhyloTree::PhyloTree() : vector< TreeNode >() {
	weight = 0;
	root = 0;
}

PhyloTree::PhyloTree( const PhyloTree& pt ) :
vector< TreeNode >( pt ),
weight( pt.weight ),
root( pt.root )
{}

PhyloTree& PhyloTree::operator=( const PhyloTree& pt )
{
	vector< TreeNode >::operator=( pt );
	weight = pt.weight;
	root = pt.root;
	return *this;
}

PhyloTree::~PhyloTree()
{}

void PhyloTree::clear()
{
	vector< TreeNode >::clear();
	weight = 0;
	root = 0;
}


/**
 *  readTree version 2.0: read in a phylogenetic tree in the Newick file format.
 *
 */
void PhyloTree::readTree( istream& tree_file ){
	string line;
	clear();
	if( !getline( tree_file, line ) )
		return;

	stringstream line_str( line );

	// look for a weight
	string::size_type open_bracket_pos = line.find( "[" );
	string::size_type bracket_pos = line.find( "]" );
	if( open_bracket_pos != string::npos && bracket_pos != string::npos && 
		open_bracket_pos < bracket_pos && bracket_pos < line.find( "(" ) ){
		// read in a weight
		getline( line_str, line, '[' );
		getline( line_str, line, ']' );
		stringstream weight_str( line );
		weight_str >> weight;
	}
	
	// ready to begin parsing the tree data.
	string tree_line;
	getline( line_str, tree_line, ';' );
	uint read_state = 0;	/**< read_state of 0 indicates nothing has been parsed yet */
	uint section_start = 0;
	stack< node_id_t > node_stack;
	stringstream blen_str;
	TreeNode new_node;
	new_node.distance = 0;	// default the distance to 0
	for( uint charI = 0; charI < tree_line.size(); charI++ ){
		switch( tree_line[ charI ] ){
			// if this is an open parens then simply create a new
			// parent node and push it on the parent stack
			case '(':
				if( node_stack.size() > 0 ){
					new_node.parents.clear();
					new_node.parents.push_back( node_stack.top() );
					(*this)[ node_stack.top() ].children.push_back( (*this).size() );
				}
				node_stack.push( (*this).size() );
				push_back( new_node );
				read_state = 1;
				section_start = charI + 1;
				break;
			case ')':
				// read off a branch length
				blen_str.clear();
				blen_str.str( tree_line.substr( section_start, charI - section_start ) );
				blen_str >> (*this)[ node_stack.top() ].distance;
				if( read_state == 2 )
					node_stack.pop();
				section_start = charI + 1;
				// pop off the top of the node stack after its branch length is read:
				read_state = 2;
				break;
			case ',':
				// read off a branch length
				blen_str.clear();
				blen_str.str( tree_line.substr( section_start, charI - section_start ) );
				blen_str >> (*this)[ node_stack.top() ].distance;
				if( read_state == 2 )
					node_stack.pop();
				section_start = charI + 1;
				read_state = 1;	// indicates that we'll be creating a new node when we hit :
				break;
			case ':':
				// read off a name, if possible
				if( read_state == 1 ){
					new_node.parents.clear();
					new_node.parents.push_back( node_stack.top() );
					(*this)[ node_stack.top() ].children.push_back( (*this).size() );
					node_stack.push( (*this).size() );
					push_back( new_node );
					read_state = 2;	// pop this node after reading its branch length
				}
				(*this)[ node_stack.top() ].name = tree_line.substr( section_start, charI - section_start );
				section_start = charI + 1;
				break;
			default:
				break;
		}
	}

}


void PhyloTree::writeTree( ostream& os ) const{
	stack< node_id_t > node_stack;
	stack< uint > child_stack;
	node_stack.push( root );
	child_stack.push( 0 );

	if( (*this).weight != 0 )
		os << "[" << weight << "]";
	os << "(";

	while( node_stack.size() > 0 ) {
		if( (*this)[ node_stack.top() ].children.size() != 0 ){
			// this is a parent node
			// if we have scanned all its children then pop it
			if( child_stack.top() == (*this)[ node_stack.top() ].children.size() ){
				os << ")";
				if( node_stack.size() > 1 )
					os << ":" << (*this)[ node_stack.top() ].distance;
				node_stack.pop();
				child_stack.pop();
				continue;
			}
			// try to recurse to its children
			// if the child is a parent as well spit out a paren
			node_id_t child = (*this)[ node_stack.top() ].children[ child_stack.top() ];
			node_stack.push( child );
			child_stack.top()++;
			// print a comma to separate multiple children
			if( child_stack.top() > 1 )
				os << ",";
			if( (*this)[ child ].children.size() > 0 ){
				child_stack.push( 0 );
				os << "(";
			}
			continue;
		}
		
		// this is a leaf node
		os << (*this)[ node_stack.top() ].name << ":" << (*this)[ node_stack.top() ].distance;
		
		// pop the child
		node_stack.pop();
	}
	os << ";" << endl;
}


double PhyloTree::getHeight() const
{
	return getHeight( root );
}
double PhyloTree::getHeight( node_id_t nodeI ) const
{
	if( (*this)[ nodeI ].children.size() == 0 )
		return (*this)[ nodeI ].distance;
	return (*this)[ nodeI ].distance + getHeight( (*this)[ nodeI ].children[ 0 ] );
}
