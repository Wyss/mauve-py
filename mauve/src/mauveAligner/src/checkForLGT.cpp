#include "libMems/PhyloTree.h"
#include <vector>
#include <sstream>
#include <algorithm>
#include <utility>
#include <fstream>
#include <set>

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

void setTaxonNames( PhyloTree< TreeNode >& t, char** taxon_names )
{
	for( node_id_t nI = 0; nI < t.size(); nI++ )
	{
		if( t[nI].name.size() == 0 )
			continue;
		stringstream ss( t[nI].name );
		uint num;
		ss >> num;
		t[nI].name = taxon_names[num];
	}
}

int main( int argc, char* argv[] )
{
	if( argc != 3 )
	{
		cerr << "Usage: checkForLGT <newick input file> <newick output file>\n";
		return -1;
	}
	string input_filename = argv[1];
	string tree_outfname = argv[2];
	vector< string > group_1;
	vector< string > group_2;
	for( uint taxonI = 0; taxonI < 16; taxonI++ )
	{
		stringstream ss;
		ss << taxonI;
		group_1.push_back( ss.str() );
	}
	for( uint taxonI = 16; taxonI < 21; taxonI++ )
	{
		stringstream ss;
		ss << taxonI;
		group_2.push_back( ss.str() );
	}
	char* taxon_names[] = {
		"E. coli 53638",
		"E. coli b171",
		"E. coli b7a",
		"E. coli e110019",
		"E. coli e22",
		"E. coli e24377a",
		"E. coli f11",
		"E. coli HS",
		"S. boydii BS512",
		"S. sonnei Ss046",
		"S. flexneri 2457T",
		"S. flexneri 301",
		"E. coli CFT073",
		"E. coli O157_H7 RIMD",
		"E. coli O157_H7 EDL933",
		"E. coli K-12 MG1655",
		"S. enterica B67",
		"S. enterica CT18",
		"S. enterica LT2",
		"S. enterica PA9150",
		"S. enterica Ty2",
	};

	ifstream input_file( input_filename.c_str() );
	if( !input_file.is_open() )
	{
		cerr << "Error opening \"" << input_filename << "\"\n";
		return -1;
	}
	
	uint tree_count = 0;
	vector< PhyloTree< TreeNode > > tree_list;
	while( true )
	{
		PhyloTree< TreeNode > new_t;
		tree_list.push_back( new_t );
		PhyloTree< TreeNode >& t = tree_list[tree_list.size() - 1];
		t.readTree( input_file );
		if( t.size() == 0 )
			break;
		tree_count++;
	}
	tree_list.erase( tree_list.end() - 1 );

	for( size_t treeI = 0; treeI < tree_list.size(); treeI++ )
	{
		PhyloTree< TreeNode >& t = tree_list[treeI];

		if( t[t.root].children.size() != 2 )
		{
			cout << treeI << "\t1\n";
			continue;
		}

		vector< node_id_t > group1_id;
		vector< node_id_t > group2_id;
		node_id_t nI = 0;
		size_t gI = 0;
		for( gI = 0; gI < group_1.size(); gI++ )
		{
			nI = 0;
			for( ; nI < t.size(); nI++ )
			{
				if( t[nI].name == group_1[gI] )
				{
					group1_id.push_back( nI );
					break;
				}
			}
			if( nI == t.size() )
			{
				cerr << "Couldn't find node " << group_1[gI] << " in tree " << treeI << endl;
				return -1;
			}
		}
		for( gI = 0; gI < group_2.size(); gI++ )
		{
			nI = 0;
			for( ; nI < t.size(); nI++ )
			{
				if( t[nI].name == group_2[gI] )
				{
					group2_id.push_back( nI );
					break;
				}
			}
			if( nI == t.size() )
			{
				cerr << "Couldn't find node " << group_2[gI] << " in tree " << treeI << endl;
				return -1;
			}
		}


		node_id_t g1_subtree;
		if( containsNode( t, t[t.root].children[0], group1_id[0] ) )
			g1_subtree = t[t.root].children[0];
		else
			g1_subtree = t[t.root].children[1];

		bool g1_monophyletic = true;
		bool g2_monophyletic = true;

		node_id_t cur_parent = group1_id[0];
		set<node_id_t> g1_remaining;
		g1_remaining.insert( group1_id.begin(), group1_id.end() );
		// find the least common ancestor of all g1 nodes
		while(g1_remaining.size() > 0)
		{
			// go to parent
			cur_parent = t[cur_parent].parents[0];
			set<node_id_t>::iterator iter = g1_remaining.begin();
			while( iter != g1_remaining.end() )
			{
				if( containsNode( t, cur_parent, *iter ) )
				{
					set<node_id_t>::iterator erase_iter = iter;
					iter++;
					g1_remaining.erase( erase_iter );
				}else
					iter++;
			}
		}
		// check none of group 2 is below the group 1 LCA
		for( gI = 0; gI < group2_id.size(); gI++ )
			if( containsNode( t, cur_parent, group2_id[gI] ) )
				break;
		if( gI < group2_id.size() )
			g1_monophyletic = false;


		cur_parent = group2_id[0];
		set<node_id_t> g2_remaining;
		g2_remaining.insert( group2_id.begin(), group2_id.end() );
		// find the least common ancestor of all g1 nodes
		while(g2_remaining.size() > 0)
		{
			// go to parent
			cur_parent = t[cur_parent].parents[0];
			set<node_id_t>::iterator iter = g2_remaining.begin();
			while( iter != g2_remaining.end() )
			{
				if( containsNode( t, cur_parent, *iter ) )
				{
					set<node_id_t>::iterator erase_iter = iter;
					iter++;
					g2_remaining.erase( erase_iter );
				}else
					iter++;
			}
		}

		// check none of group 1 is below the group 2 LCA
		for( gI = 0; gI < group1_id.size(); gI++ )
			if( containsNode( t, cur_parent, group1_id[gI] ) )
				break;
		if( gI < group1_id.size() )
			g2_monophyletic = false;

		if( !g1_monophyletic && !g2_monophyletic )
			cout << treeI << "\t2\n"; // found something interesting?
		else if( !g1_monophyletic )
			cout << treeI << "\t3\n";
		else if( !g2_monophyletic )
			cout << treeI << "\t4\n";
		else
			cout << treeI << "\t0\n"; // nothing to see here
	}

	ofstream tree_out( tree_outfname.c_str() );
	if( !tree_out.is_open() )
	{
		cerr << "Error opening \"" << tree_outfname << "\"\n";
		return -1;
	}
	for( size_t treeI = 0; treeI < tree_list.size(); treeI++ )
	{
		setTaxonNames(tree_list[treeI], taxon_names);
		tree_list[treeI].writeTree(tree_out);
	}
	return 0;
}