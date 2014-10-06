#include "libMems/PhyloTree.h"
#include "libMems/TreeUtilities.h"
#include <vector>
#include <sstream>
#include <algorithm>
#include <utility>
#include <fstream>

using namespace std;

typedef unsigned int uint;

/**
 * Depth first search to check whether a subtree contains a given node
 */
bool containsNode( PhyloTree< TreeNode >& t, node_id_t subtree_nodeI, node_id_t query_nodeI )
{
	stack< node_id_t > node_stack;
	node_stack.push( subtree_nodeI );
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		if( cur_node == query_nodeI )
			return true;
		if( t[cur_node].children.size() > 0 )
		{
			for( size_t childI = 0; childI < t[cur_node].children.size(); childI++ )
				node_stack.push( t[cur_node].children[childI] );
		}
	}
	return false;
}

/** place a root on the branch with endpoints root_left and root_right
 */
void rerootTree( PhyloTree< TreeNode >& t, node_id_t new_root )
{
	// new root must be an internal node
	if( t[new_root].children.size() == 0 )
		throw "Can't root on a leaf node";
	if( new_root == t.root )
		return;	// idiot caller didn't realize it's already rooted here

	// change the old root node to an internal node
	uint childI = 0;
	for( ; childI < t[t.root].children.size(); childI++ ){
		if( containsNode( t, t[t.root].children[childI], new_root ) )
		{
			t[t.root].parents.push_back( t[t.root].children[childI] );
			std::vector<node_id_t>::iterator last = std::remove( t[t.root].children.begin(), t[t.root].children.end(), t[t.root].children[childI] );
			t[t.root].children.erase(last,t[t.root].children.end());
			break;
		}
	}
	// shake the tree out on the new root node
	t.root = new_root;
	t[t.root].children.insert( t[t.root].children.end(), t[t.root].parents.begin(), t[t.root].parents.end() );

	stack<node_id_t> node_stack;
	node_stack.push(t.root);
	while( node_stack.size() > 0 )
	{
		// delete the current node from all of its child nodes lists 
		// and insert it as a parent
		// make all other nodes reference by the child grandchildren
		// recurse on each child
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		for( uint childI = 0; childI < t[cur_node].children.size(); childI++ )
		{
			TreeNode& child_n = t[t[cur_node].children[childI]]; 
			std::vector<node_id_t>::iterator last = std::remove( child_n.children.begin(), child_n.children.end(), cur_node );
			child_n.children.erase(last,child_n.children.end());
			last = std::remove( child_n.parents.begin(), child_n.parents.end(), cur_node );
			child_n.parents.erase(last,child_n.parents.end());
			child_n.children.insert( child_n.children.end(), child_n.parents.begin(), child_n.parents.end() );
			child_n.parents.clear();
			child_n.parents.push_back(cur_node);
			node_stack.push(t[cur_node].children[childI]);
		}
	}
}


int main( int argc, char* argv[] )
{
	if( argc < 3 )
	{
		cerr << "Usage: rootTrees <nexus input file> <nexus output file>\n";
	}
	string input_filename = argv[1];
	string output_filename = argv[2];
	ifstream input_file( input_filename.c_str() );
	if( !input_file.is_open() )
	{
		cerr << "Error opening \"" << input_filename << "\"\n";
		return -1;
	}
	ofstream output_file( output_filename.c_str() );
	if( !output_file.is_open() )
	{
		cerr << "Error opening \"" << output_filename << "\"\n";
		return -1;
	}
	
	uint tree_count = 0;
	vector< string > tree_list;
	while( true )
	{
		PhyloTree< TreeNode > t;
		t.readTree( input_file );
		if( t.size() == 0 )
			break;
		vector< PhyloTree< TreeNode > > rooted_trees;
//		rootAtEachNode( t, rooted_trees );
		for( size_t treeI = 0; treeI < rooted_trees.size(); treeI++ )
		{
			rooted_trees[treeI].writeTree( output_file );
		}
		tree_count++;
		if( tree_count % 100 == 0 )
			cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
			cout << "Processed " << tree_count << " trees";
	}
	cerr << "Wrote rooted trees to \"" << output_filename << "\"\n";
	return 0;
}