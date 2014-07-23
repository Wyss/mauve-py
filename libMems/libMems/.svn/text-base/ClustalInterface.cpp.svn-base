/*******************************************************************************
 * $Id: ClustalInterface.cpp,v 1.27 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/ClustalInterface.h"
#include <sstream>
#include "libGenome/gnFilter.h"

#include <fstream>

extern "C" {
#include "libClustalW/clustalw.h"

extern sint max_names;
extern Boolean usemenu, dnaflag, explicit_dnaflag;
extern Boolean interactive;
extern char *seqname;
extern sint nseqs;
extern sint *seqlen_array;
extern char **names,**titles;
extern char **seq_array;
//extern Boolean profile1_empty, profile2_empty;
extern sint max_aln_length;
//extern char *gap_penalty_mask, *sec_struct_mask;
//extern sint struct_penalties;
extern float    gap_open,      gap_extend;
extern float  	dna_gap_open,  dna_gap_extend;
//extern char *gap_penalty_mask1,*gap_penalty_mask2;
//extern char *sec_struct_mask1,*sec_struct_mask2;
//extern sint struct_penalties1,struct_penalties2;
//extern char *ss_name1,*ss_name2;
extern float    pw_go_penalty,      pw_ge_penalty;
extern float  	dna_pw_go_penalty,  dna_pw_ge_penalty;
//extern sint    wind_gap,ktup,window,signif;
//extern sint    dna_wind_gap, dna_ktup, dna_window, dna_signif;

extern Boolean 	output_clustal, output_nbrf, output_phylip, output_gcg, output_gde, output_nexus;
extern FILE 	*clustal_outfile, *gcg_outfile, *nbrf_outfile, *phylip_outfile, *nexus_outfile;
//extern char 	clustal_outname[FILENAMELEN+1], gcg_outname[FILENAMELEN+1];
extern char* amino_acid_codes;
extern sint max_aa;

//extern short   blosum45mt[];
//extern short   def_aa_xref[];
extern sint gap_pos1;
extern double** tmat;

extern Boolean		use_endgaps;
extern Boolean		endgappenalties;

extern sint output_order;
extern Boolean 	no_weights;

}

using namespace std;
using namespace genome;
namespace mems {

/** 
 * When performing progressive alignment, clustalW misaligns the first sequence for
 * some reason.  define MISALIGNMENT_WORKAROUND to enable a workaround for this bug.
 * The workaround adds an additional copy of the first sequence to each alignment
 * then removes the misaligned copy of the first sequence.
 */
#define MISALIGNMENT_WORKAROUND

lint get_aln_score(void);


ClustalInterface& ClustalInterface::getClustalInterface()
{
        static ClustalInterface m_ci;

        return m_ci;
}

ClustalInterface::ClustalInterface(){
	// some defaults can't hurt
	max_alignment_length = 10000;
	min_flank_size = 3;
	clustal_score_cutoff = 0;

	// shut off end gaps...
	use_endgaps = FALSE;
	// enable end gap penalties
	endgappenalties = TRUE;
	// force same input/output order
	output_order = INPUT;
	no_weights = FALSE; // TRUE;


	init_amenu();
	init_interface();
	init_matrix();

	fill_chartab();
	allocated_aln = false;
	
	// shut off end gaps...
	use_endgaps = FALSE;
	// enable end gap penalties
	endgappenalties = TRUE;
}

ClustalInterface& ClustalInterface::operator=( const ClustalInterface& ci )
{
	GappedAligner::operator=( ci );
	min_flank_size = ci.min_flank_size;
	clustal_score_cutoff = ci.clustal_score_cutoff;
	distance_matrix = ci.distance_matrix;
	allocated_aln = ci.allocated_aln;
	return *this;
}

void ClustalInterface::SetDistanceMatrix( NumericMatrix< double >& distance_matrix, string& tree_filename ){
	SetDistanceMatrix( distance_matrix, tree_filename, false );
}

void ClustalInterface::SetDistanceMatrix( NumericMatrix< double >& distance_matrix, string& tree_filename, boolean reread_tree ){
	char* phylip_name;
	uint seqI, seqJ;
#ifdef MISALIGNMENT_WORKAROUND
	if( reread_tree == false ){
		NumericMatrix< double > dist_plus_matrix( distance_matrix.cols() + 1, distance_matrix.cols() + 1 );
		for( seqI = 0; seqI < dist_plus_matrix.cols(); seqI++ ){
			for( seqJ = 0; seqJ < dist_plus_matrix.cols(); seqJ++ ){
				double new_val = 0;
				if( seqI == 0 ){
					if( seqJ == 0 )
						new_val = 0;
					else
						new_val = distance_matrix( seqI, seqJ - 1 );
				}else{
					if( seqJ == 0 )
						new_val = distance_matrix( seqI - 1, seqJ );
					else
						new_val = distance_matrix( seqI - 1, seqJ - 1 );
				}
				dist_plus_matrix( seqI, seqJ ) = new_val;
			}
		}
		SetDistanceMatrix( dist_plus_matrix, tree_filename, true );
	}
#else
	reread_tree = true;
#endif
	if( reread_tree )
		this->distance_matrix = distance_matrix;
	free_aln( nseqs );
	nseqs = distance_matrix.cols();
	alloc_aln( nseqs );
	allocated_aln = true;

	for( seqI = 1; seqI <= distance_matrix.cols(); seqI++ ){
		ostringstream ss;
		ss << "seq" << seqI;
		int namelen = MAXNAMES < ss.str().size() ? MAXNAMES : ss.str().size();
		strncpy( names[ seqI ], ss.str().c_str(), namelen);		/*    "   "  name   */
		strncpy( titles[ seqI ], ss.str().c_str(), namelen);		/*    "   "  title  */

		alloc_seq( seqI, 1 );
		// set max_names and max_aln_length
		if( strlen( names[ seqI ] ) > max_names )
			max_names = strlen( names[ seqI ] );
	}
	// copy phylo tree name
	phylip_name = (char * ) ckalloc( tree_filename.length() + 1);
	strcpy( phylip_name, tree_filename.c_str() );
	
	// copy tmat entries
	for( seqI = 0; seqI < nseqs; seqI++ )
		for( uint seqJ = 0; seqJ < nseqs; seqJ++ )
			tmat[ seqI + 1][ seqJ + 1 ] = distance_matrix( seqI, seqJ );
	
	FILE* tree;
	if((tree = open_explicit_file( phylip_name ))==NULL) return;
	if (nseqs >= 2) {
		guide_tree(tree,1,nseqs);
	}


// read the tree back in
	if( reread_tree )
		int status = read_tree(phylip_name, (sint)0, nseqs);
	phylip_name = (char*)ckfree( phylip_name );
	allocated_aln = false;

}

// tries to read in a guide tree from a particular file,
// throws an exception if it doesn't work out
void ClustalInterface::setGuideTree( string& tree_filename, NumericMatrix< double >& dist_mat, uint seq_count ){
#ifdef MISALIGNMENT_WORKAROUND
	seq_count++;
#endif
	distance_matrix = dist_mat;
	// check whether the file exists
	ifstream guide_file( tree_filename.c_str() );
	if( guide_file.is_open() )
		guide_file.close();	// success
	else
		throw( "Unable to open guide tree file" );

	char* phylip_name;
	uint seqI;

	free_aln( nseqs );
	nseqs = seq_count;
	alloc_aln( nseqs );
	allocated_aln = true;

	for( seqI = 1; seqI <= seq_count; seqI++ ){
		ostringstream ss;
		ss << "seq" << seqI;
		int namelen = MAXNAMES < ss.str().size() ? MAXNAMES : ss.str().size();
		strncpy( names[ seqI ], ss.str().c_str(), namelen);		/*    "   "  name   */
		strncpy( titles[ seqI ], ss.str().c_str(), namelen);		/*    "   "  title  */

		alloc_seq( seqI, 1 );
		// set max_names and max_aln_length
		if( strlen( names[ seqI ] ) > max_names )
			max_names = strlen( names[ seqI ] );
	}

	// copy tmat entries
	for( seqI = 0; seqI < nseqs; seqI++ )
		for( uint seqJ = 0; seqJ < nseqs; seqJ++ )
			tmat[ seqI + 1][ seqJ + 1 ] = 1 - distance_matrix( seqI, seqJ );

	// copy phylo tree name
	phylip_name = (char * ) ckalloc( tree_filename.length() + 1);
	strcpy( phylip_name, tree_filename.c_str() );
	int success = read_tree(phylip_name, (sint)0, nseqs);
	phylip_name = (char*)ckfree( phylip_name );
	allocated_aln = false;
	if( !success )
		throw "Error loading guide tree\n";
}

boolean ClustalInterface::Align( GappedAlignment& cr, Match* r_begin, Match* r_end, vector< gnSequence* >& seq_table ){
	boolean flank = false;
	gnSeqI gap_size = 0;
	boolean create_ok = true;
	uint seq_count = seq_table.size();
	uint seqI;
	uint align_seqs = 0;
//	
//  get the size of the largest intervening gap
//  also do some sanity checking while we're at it.
//
try{
	for( seqI = 0; seqI < seq_count; seqI++ ){
		int64 gap_start = 0;
		int64 gap_end = 0;
		create_ok = getInterveningCoordinates( seq_table, r_begin, r_end, seqI, gap_start, gap_end );
		// skip this sequence if it's undefined
		if( gap_start == NO_MATCH || gap_end == NO_MATCH )
			continue;
		if( !create_ok )
			break;

		int64 diff = gap_end - gap_start;
		if( diff <= 0 ){
			continue;	// can't align nothing
		}
		if( diff > max_alignment_length ){
			cout << "gap from " << gap_start << " to " << gap_end << " is too big for ClustalW\n";
			continue;	// can't align if it's too big
		}
		gap_size = diff < gap_size ? gap_size : diff;
		align_seqs++;
	}

	if( align_seqs <= 1 )
		create_ok = false;
// 
//	Get the sequence in the intervening gaps between these two matches
//  Include a flank of matching sequence on either side
//
	vector< string > seq_data;
	vector< int64 > starts;
	gnSeqI left_flank, right_flank;
	const gnFilter* rc_filter = gnFilter::DNAComplementFilter();

	if( create_ok ){
//		left_flank = min( r_begin->Length(), max( gap_size, min_flank_size ) );
//		right_flank = min( r_end->Length(), max( gap_size, min_flank_size ) );
		left_flank = 0;
		right_flank = 0;
		
		for( seqI = 0; seqI < seq_count; seqI++ ){
			// cheap hack to avoid mysterious clustalW misalignment
#ifdef MISALIGNMENT_WORKAROUND
			if( seqI == 1 )
				seq_data.push_back( seq_data[ 0 ] );
#endif
			// skip this sequence if it's undefined
			if( (r_end != NULL && r_end->Start( seqI ) == NO_MATCH ) ||
				(r_begin != NULL && r_begin->Start( seqI ) == NO_MATCH) ){
				starts.push_back( NO_MATCH );
				seq_data.push_back( "" );
				continue;
			}

			// determine the size of the gap
			int64 gap_start = 0;
			int64 gap_end = 0;
			getInterveningCoordinates( seq_table, r_begin, r_end, seqI, gap_start, gap_end );

			int64 diff = gap_end - gap_start;
			if( diff <= 0 || diff > max_alignment_length ){
				starts.push_back( NO_MATCH );
				seq_data.push_back( "" );
				continue;
			}
			// calculate flank size and extract sequence data
			if( r_end == NULL || r_end->Start( seqI ) > 0 ){
				starts.push_back( gap_start );
				seq_data.push_back( seq_table[ seqI ]->ToString( left_flank + diff + right_flank, gap_start - left_flank ) );
			}else{
				// reverse complement the sequence data.
				starts.push_back( -gap_start );
				string cur_seq_data = seq_table[ seqI ]->ToString( left_flank + diff + right_flank, gap_start - right_flank );
				rc_filter->ReverseFilter( cur_seq_data );
				seq_data.push_back( cur_seq_data );
			}
		}
	}

	if( create_ok ){
		if( !CallClustal( seq_data ) ){
			cout << "Clustal was unable to align:\n";
			cout << "Left match: " << *r_begin << endl;
			cout << "Right match: " << *r_end << endl;
			return false;
		}

		// ensure that the flanks were successfully aligned
		boolean good_alignment = true;
		gnSeqI flankI = 0;
		gnSeqI align_length=0;
		for( seqI = 1; seqI <= seq_count; seqI++ )
			if( align_length < ( seqlen_array[seqI] < 0 ? 0 : (gnSeqI)seqlen_array[seqI] ))
				align_length = seqlen_array[seqI];

		if( !good_alignment ){
			// just align without the flanking regions for now??
			return false;
		}else{
			// now extract the alignment from clustal's global variables
			cr = GappedAlignment( seq_count, align_length );
			vector< string > align_array;
			int64 last_residue = -1;	// tracks the right-most residue in the alignment
			int64 first_residue = align_length + 2;	// tracks the left-most residue in the alignment
#ifdef MISALIGNMENT_WORKAROUND
			for( seqI = 2; seqI <= seq_count + 1; seqI++ ){
#else
			for( seqI = 1; seqI <= seq_count; seqI++ ){
#endif
				string new_seq = string( seqlen_array[ seqI ] - left_flank - right_flank, '-' );
				uint new_seq_charI = 0;
				uint cur_seq_len = 0;
				for( uint charJ = left_flank + 1; charJ <= seqlen_array[ seqI ] - right_flank; charJ++ ){
					char val = seq_array[ seqI ][ charJ ];
					if( val >= 0 && val <= max_aa ){
						if( charJ > last_residue )
							last_residue = charJ;
						if( charJ < first_residue )
							first_residue = charJ;
						new_seq[ new_seq_charI ]= amino_acid_codes[ val ];
						cur_seq_len++;
					}
					new_seq_charI++;
				}
				align_array.push_back( new_seq );
//				cerr << "new_seq.size() is: " << new_seq.size() << endl;
#ifdef MISALIGNMENT_WORKAROUND
				cr.SetStart( seqI - 2, starts[ seqI - 2 ] );
				cr.SetLength( cur_seq_len, seqI - 2 );
#else
				cr.SetStart( seqI - 1, starts[ seqI - 1 ] );
				cr.SetLength( cur_seq_len, seqI - 1 );
#endif
			}
			int64 end_gap_count = align_array[ 0 ].size() - (last_residue - left_flank);
			if( last_residue != -1 && end_gap_count > 0 ){
				for( seqI = 0; seqI < align_array.size(); seqI++ ){
					align_array[ seqI ] = align_array[ seqI ].substr( 0, align_array[ seqI ].size() - end_gap_count );
				}
			}
			int64 start_gap_count = left_flank + 1 - first_residue;
			if( first_residue != align_length && start_gap_count > 0 ){
				for( seqI = 0; seqI < align_array.size(); seqI++ ){
					align_array[ seqI ] = align_array[ seqI ].substr( start_gap_count, align_array[ seqI ].size() - start_gap_count );
				}
			}
			cr.SetAlignment( align_array );
		}
		return true;
	}
}catch(exception& e){
	cerr << "At: " << __FILE__ << ":" << __LINE__ << endl;
	cerr << e.what();
}
	return false;
}


boolean ClustalInterface::CallClustal( vector< string >& seq_table ){
	char* phylip_name;
	
//	if( allocated_aln )
		free_aln( nseqs );
	alloc_aln( seq_table.size() );
	allocated_aln = true;

	if( distance_matrix.cols() == seq_table.size() ){
		// copy tmat entries
		for( uint seqI = 0; seqI < nseqs; seqI++ )
			for( uint seqJ = 0; seqJ < nseqs; seqJ++ )
				tmat[ seqI + 1][ seqJ + 1 ] = 1 - distance_matrix( seqI, seqJ );
	}else{
		// prepare to infer a phylo tree
		phylip_name = (char * ) ckalloc( strlen( "tmp_tree.txt" ) + 1);
		strcpy( phylip_name, "tmp_tree.txt" );
	}

	uint seqI;
	max_aln_length = 0;
	max_names = 0;
	for( seqI = 1; seqI <= seq_table.size(); seqI++ ){
		seqlen_array[ seqI ] = seq_table[ seqI - 1 ].length();		/* store the length */
		ostringstream ss;
		ss << "seq" << seqI;
		int namelen = ss.str().size();
		names[ seqI ] = (char * ) ckalloc( namelen + 1 );
		titles[ seqI ] = (char * ) ckalloc( namelen + 1 );
		strcpy( names[ seqI ], ss.str().c_str());		/*    "   "  name   */
		strcpy( titles[ seqI ], ss.str().c_str());		/*    "   "  title  */

		// set max_names and max_aln_length
		if( (int)strlen( names[ seqI ] ) > max_names )
			max_names = strlen( names[ seqI ] );
		if( seqlen_array[ seqI ] > max_aln_length )
			max_aln_length = seqlen_array[ seqI ];
	}

	for( seqI = 1; seqI <= seq_table.size(); seqI++ ){

		alloc_seq( seqI, max_aln_length );
		char* seq_char_array = new char[ seq_table[ seqI - 1 ].length() + 2];
		uint copyI = 0;
		string& dna_seq = seq_table[ seqI - 1 ];
		for( ; copyI < dna_seq.length(); copyI++ )
			seq_char_array[ copyI + 1 ] = toupper( dna_seq[ copyI ] );
		seq_char_array[ 0 ] = '-';	// silly clustal ignores the first character.
		seq_char_array[ copyI + 1 ] = 0;
		n_encode( seq_char_array, seq_array[ seqI ], dna_seq.length() );
		delete[] seq_char_array;
	}
	max_aln_length *= 2;

/*	struct_penalties1 = struct_penalties2 = NONE;
	if (sec_struct_mask1 != NULL) sec_struct_mask1=( char* )ckfree(sec_struct_mask1);
	if (sec_struct_mask2 != NULL) sec_struct_mask2=( char* )ckfree(sec_struct_mask2);
	if (gap_penalty_mask1 != NULL) gap_penalty_mask1=( char* )ckfree(gap_penalty_mask1);
	if (gap_penalty_mask2 != NULL) gap_penalty_mask2=( char* )ckfree(gap_penalty_mask2);
	if (ss_name1 != NULL) ss_name1=( char* )ckfree(ss_name1);
	if (ss_name2 != NULL) ss_name2=( char* )ckfree(ss_name2);
*/	
	nseqs = seq_table.size();
    gap_open   = dna_gap_open;
    gap_extend = dna_gap_extend;
    pw_go_penalty  = dna_pw_go_penalty;
    pw_ge_penalty  = dna_pw_ge_penalty;
/*    ktup       = dna_ktup;
    window     = dna_window;
    signif     = dna_signif;
    wind_gap   = dna_wind_gap;
*/	dnaflag = TRUE;
	output_clustal = FALSE;
	
	int retval = 0;
	if( distance_matrix.cols() == seq_table.size() ){
//		char* dump_file = "clustalout.txt";
//		output_clustal = TRUE;
//		if((clustal_outfile = open_explicit_file( dump_file ))==NULL) return false;
		
		retval = malign_nofiles( 0, false );
//		create_alignment_output( 1, nseqs );
//		fclose( clustal_outfile );	// this is done by the clustal output function

	}else{
		pairalign((sint)0,nseqs,(sint)0,nseqs);
		
		FILE* tree;
		if((tree = open_explicit_file( phylip_name ))==NULL) return false;
		if (nseqs >= 2) {
			guide_tree(tree,1,nseqs);
		}

//		char* dump_file = "clustalout.txt";
//		if((clustal_outfile = open_explicit_file( dump_file ))==NULL) return false;
		
		retval = malign( 0, phylip_name );

		phylip_name = (char*)ckfree( phylip_name );
	//	fclose( clustal_outfile );	// this is done by the clustal output function
	}

	if( retval <= 0 )
		return false;
	return true;

}

/*
lint get_aln_score(void)
{
  static short  *mat_xref, *matptr;
  static sint maxres;
  static sint  s1,s2,c1,c2;
  static sint    ngaps;
  static sint    i,l1,l2;
  static lint    score;
  static sint   matrix[NUMRES][NUMRES];


  matptr = blosum45mt;
  mat_xref = def_aa_xref;
  maxres = get_matrix(matptr, mat_xref, matrix, TRUE, 100);
  if (maxres == 0)
    {
       fprintf(stdout,"Error: matrix blosum30 not found\n");
       return -1;
    }

  score=0;
  for (s1=1;s1<=nseqs;s1++)
   {
    for (s2=1;s2<s1;s2++)
      {

        l1 = seqlen_array[s1];
        l2 = seqlen_array[s2];
        for (i=1;i<l1 && i<l2;i++)
          {
            c1 = seq_array[s1][i];
            c2 = seq_array[s2][i];
            if ((c1>=0) && (c1<=max_aa) && (c2>=0) && (c2<=max_aa))
                score += matrix[c1][c2];
          }

        ngaps = count_gaps(s1, s2, l1);

        score -= (int)(100 * gap_open * ngaps);

      }
   }

  score /= 100;

  return score;
}
*/

}
