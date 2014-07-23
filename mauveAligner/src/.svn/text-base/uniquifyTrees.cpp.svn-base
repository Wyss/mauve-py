#include "libMems/PhyloTree.h"
#include <vector>
#include <sstream>
#include <algorithm>
#include <utility>
#include <fstream>

using namespace std;

typedef unsigned int uint;

bool taxonNameLessThan( string name1, string name2 )
{
	stringstream n1_str( name1 );
	stringstream n2_str( name2 );
	int n1, n2;
	n1_str >> n1;
	n2_str >> n2;
	return n1 < n2;
}

template<class T, class S>
void findAndErase( T& container, S& item )
{
	T new_container;
	for( typename T::iterator t_iter = container.begin(); t_iter != container.end(); t_iter++ )
		if( *t_iter != item )
			new_container.push_back( *t_iter );
	container = new_container;
};

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
			findAndErase( t[t.root].children, t[t.root].children[childI] );
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
			findAndErase( child_n.children, cur_node );
			findAndErase( child_n.parents, cur_node );
			child_n.children.insert( child_n.children.end(), child_n.parents.begin(), child_n.parents.end() );
			child_n.parents.clear();
			child_n.parents.push_back(cur_node);
			node_stack.push(t[cur_node].children[childI]);
		}
	}
}

/**
 * Find the leaf node lexicographically least taxon name in the 
 * subtree below nodeI
 */
node_id_t getRepresentativeTaxon( PhyloTree< TreeNode >& t, node_id_t nodeI )
{
	stack< node_id_t > node_stack;
	node_stack.push( nodeI );
	string least_name = "";
	node_id_t least_node = nodeI;
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		if( t[cur_node].children.size() > 0 )
		{
			for( size_t childI = 0; childI < t[cur_node].children.size(); childI++ )
				node_stack.push( t[cur_node].children[childI] );
		}
		else
		{
			if( least_name == "" )
			{
				least_name = t[cur_node].name;
				least_node = cur_node;
			}
			if( taxonNameLessThan( t[cur_node].name, least_name ) )
			{
				least_name = t[cur_node].name;
				least_node = cur_node;
			}
		}
	}
	return least_node;
}

class TaxonNamePairComparator
{
public:
	bool operator()( const pair<string, size_t>& p1, const pair<string, node_id_t>& p2 )
	{
		return taxonNameLessThan( p1.first, p2.first );
	}
};

void sortTaxa( PhyloTree< TreeNode >& t )
{
	for( node_id_t nodeI = 0; nodeI < t.size(); nodeI++ )
	{
		if( t[nodeI].children.size() == 0 )
			continue;
		// get the "representative" of each subtree
		vector< pair<string, node_id_t> > representatives = vector< pair<string, node_id_t> >( t[nodeI].children.size() );
		for( size_t repI = 0; repI < representatives.size(); repI++ )
		{
			node_id_t rep_node = getRepresentativeTaxon( t, t[nodeI].children[ repI ] );
			representatives[ repI ] = make_pair( t[rep_node].name, repI );
		}
		// sort children on their representative taxon names
		TaxonNamePairComparator tnc;
		sort( representatives.begin(), representatives.end(), tnc );
		// repopulate the children array with the sorted order
		vector< node_id_t > sorted_children;
		for( size_t repI = 0; repI < representatives.size(); repI++ )
			sorted_children.push_back( t[nodeI].children[representatives[repI].second] );
		t[nodeI].children = sorted_children;
	}
}

/**
 * Assumes that taxa have numeric labels starting at 1 and simply
 * subtracts 1 from each node label
 */
void relabelTaxaToStartWithZero( PhyloTree< TreeNode >& t )
{
	for( node_id_t nodeI = 0; nodeI < t.size(); nodeI++ )
	{
		if( t[nodeI].name == "" )
			continue;
		stringstream name_str( t[nodeI].name );
		uint number;
		name_str >> number;
		number--;
		stringstream new_name_str;
		new_name_str << number;
		t[nodeI].name = new_name_str.str();
	}
}

int main( int argc, char* argv[] )
{
	if( argc < 3 )
	{
		cerr << "Usage: uniquifyTrees <nexus input file> <nexus output file>\n";
		cerr << "All trees in the input file must have the same number of taxa and the same taxon labels\n";
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
	
	size_t tree_sizes = 0;
	uint tree_count = 0;
	vector< string > tree_list;
	while( true )
	{
		PhyloTree< TreeNode > t;
		t.readTree( input_file );
		if( t.size() == 0 )
			break;
		if( tree_sizes == 0 )
			tree_sizes = t.size();
		if( t.size() != tree_sizes )
		{
			cerr << "Error: tree " << tree_count + 1 << " has a different number of taxa\n";
			return -2;
		}
 		sortTaxa( t );
		relabelTaxaToStartWithZero( t );
		stringstream ss;
		t.writeTree(ss);
		tree_list.push_back(ss.str());
		cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
		cout << "Read " << tree_list.size() << " trees";
	}
	cout << endl;
	cout << "Writing unique trees to \"" << output_filename << "\"\n";
	sort(tree_list.begin(), tree_list.end() );
	size_t unique_count = 0;
	for( size_t treeI = 0; treeI < tree_list.size(); treeI++ )
	{
		if( treeI > 0 && tree_list[treeI] == tree_list[treeI - 1] )
			continue;
		output_file << tree_list[treeI] << endl;
		unique_count++;
	}
	cerr << "There are " << unique_count << " unique trees\n";
	return 0;
}
