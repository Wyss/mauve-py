/*******************************************************************************
 * $Id: GenericIntervalList.h,v 1.6 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _IntervalList_h_
#define _IntervalList_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <list>
#include <sstream>

#include "libMems/SortedMerList.h"
#include "libGenome/gnSequence.h"
#include "libMems/Interval.h"
#include "libMems/MemHash.h"
#include "libMems/CompactGappedAlignment.h"
#include "libGenome/gnSourceFactory.h"
#include "libGenome/gnFASSource.h"
#include "libGenome/gnSEQSource.h"
#include "libGenome/gnGBKSource.h"
#include "libGenome/gnRAWSource.h"

namespace mems {

/**
 * This class represents a set Intervals, each of which is a collinear aligned region
 * There are functions to read and write an GenericIntervalList.
 * @see Interval
 */
template< class MatchType = Interval >
class GenericIntervalList : public std::vector< MatchType > {
public:
	GenericIntervalList(){};
	GenericIntervalList( const GenericIntervalList& ml );
	GenericIntervalList& operator=( const GenericIntervalList& ml );
	
	/**
	 * Deletes the objects associated
	 * with this GenericIntervalList.
	 */
	void Clear();

	/**
	 * Reads a GenericIntervalList from an input stream
	 * Sequence and SML file names are read into the seq_filename
	 * and sml_filename vectors, but the actual files are not
	 * opened.  The calling function should load them after
	 * using this method.
	 * @param match_stream The input stream to read from
	 */
	void ReadList( std::istream& match_stream );

	/**
	 *  Writes a GenericIntervalList to the designated output stream
	 * @param match_stream The outptu stream to write to
	 */
	void WriteList( std::ostream& match_stream ) const;
	
	/**
	 *  Writes a gapped alignment of sequences to the output stream
	 */
	void WriteAlignedSequences(std::ostream& match_file) const;
	
	/**
	 *	Writes a gapped alignment of sequences in a standard format
	 */
	void WriteStandardAlignment( std::ostream& out_file ) const;

    /**
	 *	Writes a gapped alignment of sequences in xml format
	 */
    void WriteXMLAlignment( std::ostream& out_file ) const;
	
	/**
	 * Reads in a set of intervals that are in xmfa (eXtended MultiFastA) format
	 */
	void ReadStandardAlignment( std::istream& in_stream );

	/**
	 * Reads in a set of intervals that are in xmfa (eXtended MultiFastA) format
	 * and stores them in CompactGappedAlignments<>
	 */
	void ReadStandardAlignmentCompact( std::istream& in_stream );
	
	std::vector<std::string> seq_filename;	/**< The names of files associated with the sequences used by this alignment */
	std::vector<genome::gnSequence*> seq_table;	/**< The actual sequences used in this alignment */

	std::string backbone_filename;	/**< The name of an associated backbone file, or empty if none exists */
protected:

};


typedef GenericIntervalList<> IntervalList;

template< class MatchType >
GenericIntervalList<MatchType>::GenericIntervalList( const GenericIntervalList<MatchType>& ml )
{
	*this = ml;
}

template< class MatchType >
GenericIntervalList<MatchType>& GenericIntervalList<MatchType>::operator=( const GenericIntervalList<MatchType>& ml )
{
	std::vector< MatchType >::operator=( ml );
	seq_filename = ml.seq_filename;
	seq_table = ml.seq_table;
	return *this;
}

template< class MatchType >
void GenericIntervalList<MatchType>::Clear() 
{
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ ){
		if( seq_table[ seqI ] != NULL )
			delete seq_table[ seqI ];
	}
	seq_filename.clear();
	this->clear();
}

template< class MatchType >
void GenericIntervalList<MatchType>::ReadList(std::istream& match_file)
{
	std::string tag;
	gnSeqI len;
	int64 start;
	unsigned int seq_count;
	uint seqI;
	
	match_file >> tag;	//format version tag
	if( tag != "FormatVersion" ){
		Throw_gnEx(InvalidFileFormat());
	}
	match_file >> tag;	//format version
	if( tag != "4" ){
		Throw_gnEx(InvalidFileFormat());
	}
	match_file >> tag;	//sequence count tag
	if( tag != "SequenceCount" ){
		Throw_gnEx(InvalidFileFormat());
	}
	match_file >> seq_count;	//sequence count
	if(seq_count < 2){
		Throw_gnEx(InvalidFileFormat());
	}
	
	std::vector< std::string > alignment;
	// read the sequence file names and lengths
	for( seqI = 0; seqI < seq_count; seqI++ ){
		match_file >> tag;	// name tag
		getline( match_file, tag );
		// skip the tab character
		tag = tag.substr( 1 );
		seq_filename.push_back(tag);
//		try{
//			gnSequence *new_seq = new gnSequence();
//			new_seq->LoadSource(tag);
//			seq_table.push_back( new_seq );
//		}catch( gnException& gne );
		match_file >> tag;	// length tag
		match_file >> tag;	// length

		alignment.push_back( "" );	// initialize alignment vector
	}
	uint interval_count;
	match_file >> tag;	// interval count tag
	match_file >> interval_count;	// interval count
	
	
	// read the matches
	std::string cur_line;
	Interval* cur_iv = NULL;
	boolean clustal_match;
	std::vector< AbstractMatch* > iv_matches;
	bool parsing = false;
	
	while( std::getline( match_file, cur_line ) ){
		if( cur_line.find( "Interval" ) != std::string::npos ){
			// end the old interval
			if( iv_matches.size() > 0 )
			{
				this->push_back( Interval(iv_matches.begin(), iv_matches.end()) );
//				for( size_t mI = 0; mI < iv_matches.size(); mI++ )
//					delete iv_matches[mI];
				iv_matches.clear();
			}
			parsing = true;
			continue;
		}
		if( !parsing )
			continue;
		if( cur_line.length() == 0 )
			continue;
		
		clustal_match = false;
		if( cur_line == "GappedAlignment" ){
			clustal_match = true;
			getline( match_file, cur_line );

			std::stringstream line_stream( cur_line );
			line_stream >> len;
			GappedAlignment* cr = new GappedAlignment( seq_count, len );
			
			for( seqI = 0; seqI < seq_count; seqI++ ){
				line_stream >> start;
				cr->SetStart( seqI, start );
				std::getline( match_file, alignment[ seqI ] );
				int64 seq_len = 0;
				for( uint charI = 0; charI < alignment[ seqI ].length(); charI++ )
					if( alignment[ seqI ][ charI ] != '-' )
						seq_len++;
				cr->SetLength( seq_len, seqI );
			}
			cr->SetAlignment( alignment );
			iv_matches.push_back( cr );
		}
		else{

			Match mmhe( seq_count );
			Match* mhe = mmhe.Copy();
			std::stringstream line_stream( cur_line );
			
			line_stream >> len;
			mhe->SetLength(len);

			for( seqI = 0; seqI < seq_count; seqI++){
				line_stream >> start;
				mhe->SetStart(seqI, start);
			}
		
			iv_matches.push_back( mhe );
		}
	}
	if( iv_matches.size() > 0 )
		this->push_back( Interval(iv_matches.begin(), iv_matches.end()) );
	if( interval_count != this->size() ){
		Throw_gnEx(InvalidFileFormat());
	}
	
}

template< class MatchType >
void GenericIntervalList<MatchType>::WriteList(std::ostream& match_file) const
{

	unsigned int seq_count = seq_table.size();
	uint seqI;
	
	match_file << "FormatVersion" << '\t' << 4 << "\n";
	match_file << "SequenceCount" << '\t' << seq_count << "\n";
	for(seqI = 0; seqI < seq_count; seqI++){
		match_file << "Sequence" << seqI << "File" << '\t';
		if( seq_filename.size() > seqI )
			match_file << seq_filename[seqI];
		else
			match_file << "null";
		match_file << "\n";
		match_file << "Sequence" << seqI << "Length" << '\t';
		if( seq_table.size() > seqI )
			match_file << seq_table[seqI]->length();
		else
			match_file << "0";
		match_file << "\n";
	}

	match_file << "IntervalCount" << '\t' << this->size() << std::endl;
	
	for( uint ivI = 0; ivI < this->size(); ivI++ ){
		match_file << "Interval " << ivI << std::endl;
		const std::vector<AbstractMatch*>& matches = (*this)[ ivI ].GetMatches();
		for( uint matchI = 0; matchI < matches.size(); matchI++ ){
			const AbstractMatch* m = matches[ matchI ];
			const GappedAlignment* cr = dynamic_cast< const GappedAlignment* >( m );
			const Match* match = dynamic_cast< const Match* >( m );
			if( match != NULL ){
				match_file << *match << std::endl;
			}
			else if( cr != NULL ){
				match_file << "GappedAlignment\n";
				match_file << cr->Length();
				for( seqI = 0; seqI < seq_count; seqI++ )
					match_file << '\t' << cr->Start( seqI );
				match_file << std::endl;

				const std::vector< std::string >& align_matrix = GetAlignment( *cr, seq_table );
				for( seqI = 0; seqI < seq_count; seqI++ )
					match_file << align_matrix[ seqI ] << std::endl;
			}
		}
		match_file << std::endl;
	}
}

//stub for now, later use a XML library to write/read alignments in xml format..
template< class MatchType >
void GenericIntervalList<MatchType>::WriteXMLAlignment( std::ostream& out_file ) const 
{
	if( this->size() == 0 )
		return;
    // write source sequence filenames and formats
	// to make Paul happy
	boolean single_input = true;
    uint seqI = 0;
	for( seqI = 1; seqI < seq_filename.size(); seqI++ ){
		if( seq_filename[ 0 ] != seq_filename[ seqI ] ){
			single_input = false;
			break;
		}
	}
//	unsigned int seq_count = seq_table.size();
	
    out_file << "<procrastAlignment sequence=\"" << seq_filename[ 0 ] << "\">" << std::endl;
	for( uint ivI = 0; ivI < this->size(); ivI++ ){
		if( (*this)[ ivI ].AlignmentLength() == 0 ){
			continue;
		}
        out_file << "\t<localAlignment id = \"" << ivI+1 << "\" length = \"" << (*this)[ ivI ].AlignmentLength() << "\" multiplicity = \"" << (*this)[ ivI ].Multiplicity() << "\" spscore=\"" << (*this)[ ivI ].spscore << "\">" << std::endl;
    
		std::vector<std::string> alignment;
		if( seq_table.size() == 1 && seq_table.size() != (*this)[ ivI ].SeqCount() )
		{
			GetAlignment( (*this)[ ivI ], std::vector<genome::gnSequence*>((*this)[ ivI ].SeqCount(), seq_table[0]), alignment );
		}else
			GetAlignment( (*this)[ ivI ], seq_table, alignment );
		for( seqI = 0; seqI < (*this)[ ivI ].SeqCount(); seqI++ ){
			int64 startI = (*this)[ ivI ].Start( seqI );
			gnSeqI length = (*this)[ ivI ].Length( seqI );
			// if this genome doesn't have any sequence in this
			// interval then skip it...
			if( startI == 0 &&ivI > 0)	// kludge: write all seqs into the first interval so java parser can read it
				continue;
   
		    out_file << "\t\t<component id=\"" << seqI+1 << "\" seqid=\"1\" leftend=\"" << (*this)[ ivI ].LeftEnd( seqI ) << "\" length=\"" << (*this)[ ivI ].Length( seqI ) << "\" orientation=\"" <<  (*this)[ ivI ].Orientation( seqI) << "\">" << alignment[ seqI ].data();
            out_file << "\t\t</component> " << std::endl;


		}
		out_file << "\t</localAlignment>" << std::endl;
	}
	out_file << "</procrastAlignment>" << std::endl;
}

template< class MatchType >
void GenericIntervalList<MatchType>::WriteStandardAlignment( std::ostream& out_file ) const 
{
	if( this->size() == 0 )
		return;

//	unsigned int seq_count = seq_table.size();
	uint seqI = 0;
	
	// write out the format version
	out_file << "#FormatVersion Mauve1" << std::endl;
	
	// write source sequence filenames and formats
	// to make Paul happy
	boolean single_input = true;
	for( seqI = 1; seqI < seq_filename.size(); seqI++ ){
		if( seq_filename[ 0 ] != seq_filename[ seqI ] ){
			single_input = false;
			break;
		}
	}
	for( seqI = 0; seqI < seq_filename.size(); seqI++ ){
		out_file << "#Sequence" << seqI + 1 << "File\t" << seq_filename[ seqI ] << std::endl;
		if( single_input )
			out_file << "#Sequence" << seqI + 1 << "Entry\t" << seqI + 1 << std::endl;
		
		genome::gnSourceFactory* sf = genome::gnSourceFactory::GetSourceFactory();
		genome::gnBaseSource* gnbs = sf->MatchSourceClass( seq_filename[ seqI ] );
		genome::gnFASSource* gnfs = dynamic_cast< genome::gnFASSource* >(gnbs);
		genome::gnRAWSource* gnrs = dynamic_cast< genome::gnRAWSource* >(gnbs);
		genome::gnSEQSource* gnss = dynamic_cast< genome::gnSEQSource* >(gnbs);
		genome::gnGBKSource* gngs = dynamic_cast< genome::gnGBKSource* >(gnbs);
		if( gnfs != NULL )
			out_file << "#Sequence" << seqI + 1 << "Format\tFastA" << std::endl;
		else if( gnrs != NULL )
			out_file << "#Sequence" << seqI + 1 << "Format\traw" << std::endl;
		else if( gnss != NULL ){
			out_file << "#Sequence" << seqI + 1 << "Format\tDNAstar" << std::endl;
			out_file << "#Annotation" << seqI + 1 << "File\t" << seq_filename[ seqI ] << std::endl;
			out_file << "#Annotation" << seqI + 1 << "Format\tDNAstar" << std::endl;
		}else if( gngs != NULL ){
			out_file << "#Sequence" << seqI + 1 << "Format\tGenBank" << std::endl;
			out_file << "#Annotation" << seqI + 1 << "File\t" << seq_filename[ seqI ] << std::endl;
			out_file << "#Annotation" << seqI + 1 << "Format\tGenBank" << std::endl;
		}
	}

	if( this->backbone_filename != "" )
		out_file << "#BackboneFile\t" << this->backbone_filename << std::endl;
	
	for( uint ivI = 0; ivI < this->size(); ivI++ ){
		if( (*this)[ ivI ].AlignmentLength() == 0 ){
			continue;
		}
		std::vector<std::string> alignment;
		if( seq_table.size() == 1 && seq_table.size() != (*this)[ ivI ].SeqCount() )
		{
			GetAlignment( (*this)[ ivI ], std::vector<genome::gnSequence*>((*this)[ ivI ].SeqCount(), seq_table[0]), alignment );
		}else
			GetAlignment( (*this)[ ivI ], seq_table, alignment );
		for( seqI = 0; seqI < (*this)[ ivI ].SeqCount(); seqI++ ){
			int64 startI = (*this)[ ivI ].Start( seqI );
			gnSeqI length = (*this)[ ivI ].Length( seqI );
			// if this genome doesn't have any sequence in this
			// interval then skip it...
			if( startI == 0 &&ivI > 0)	// kludge: write all seqs into the first interval so java parser can read it
				continue;
			out_file << "> " << seqI + 1 << ":";
			if( startI > 0 ){
				out_file << genome::absolut( startI ) << "-" << genome::absolut( startI ) + length - 1 << " + ";
			}else if(startI == 0){
				out_file << 0 << "-" << 0 << " + ";
			}else{
				out_file << genome::absolut( startI ) << "-" << genome::absolut( startI ) + length - 1 << " - ";
			}
			if( single_input )
				out_file << seq_filename[ 0 ];	// write the sequence filename as the seq name
			else				
				out_file << seq_filename[ seqI ];	// write the sequence filename as the seq name
			out_file << std::endl;
			gnSeqI cur_pos = 0;
			for( ; cur_pos < alignment[ seqI ].length(); cur_pos += 80 ){
				gnSeqI cur_len = cur_pos + 80 < alignment[ seqI ].length() ? 80 : alignment[ seqI ].length() - cur_pos;
				out_file.write( alignment[ seqI ].data() + cur_pos, cur_len );
				out_file << std::endl;
			}
		}
		out_file << "=" << std::endl;
	}
	
}

template< class MatchType >
void GenericIntervalList<MatchType>::ReadStandardAlignment( std::istream& in_stream ) 
{
	uint seq_count = 0;
	gnSeqI max_len = 0;
	std::string cur_line;
	if( !std::getline( in_stream, cur_line ) )
	{
		Clear();	// if we can't read from the file then just return an empty interval list
		return;
	}
	uint seqI = 0;
	std::vector< gnSeqI > lengths;
	std::vector< GappedAlignment* > ga_list;
	GappedAlignment cr;
	std::string empty_line;
	std::vector< std::string > aln_mat;
	uint line_count = 1;
	while( true ){
		
		while( cur_line[0] == '#' ){
			// hit a comment or metadata.  try to parse it if it's a filename
			std::getline( in_stream, cur_line );
			line_count++;
			std::stringstream ss( cur_line );
			std::string token;
			std::getline( ss, token, '\t' );
			if( token.substr(1, 8) != "Sequence" || token.find( "File" ) == std::string::npos )
				continue;
			std::getline( ss, token );
			seq_filename.push_back( token );
		}
		
		// read and parse the def. line
		std::stringstream line_str( cur_line );
		std::getline( line_str, cur_line, '>' );
		std::getline( line_str, cur_line, ':' );
		// take off leading whitespace
		std::stringstream parse_str( cur_line );

		parse_str >> seqI;	// the sequence number
				
		int64 start, stop;
		std::getline( line_str, cur_line, '-' );
		parse_str.clear();
		parse_str.str( cur_line );
		parse_str >> start;
		line_str >> stop;
		std::string strand;
		line_str >> strand;

		std::string name;	// anything left is the name
		std::getline( line_str, name );

		// read and parse the sequence
		while( aln_mat.size() < seqI )
			aln_mat.push_back( empty_line );

		gnSeqI chars = 0;
		while( std::getline( in_stream, cur_line ) ){
			line_count++;
			if( (cur_line[ 0 ] == '>' ) || (cur_line[ 0 ] == '=' ))
				break;
			for( uint charI = 0; charI < cur_line.length(); charI++ )
				if( cur_line[ charI ] != '-' )
					chars++;
			aln_mat[ seqI - 1 ] += cur_line;
		}
		while( lengths.size() < seqI )
			lengths.push_back(0);

		lengths[ seqI - 1 ] = chars;

// temporary workaround for file format inconsistency
		if( strand == "+" )
			cr.SetStart( seqI - 1, start );
		else if( start < stop ){
			if( chars == 0 )
				cr.SetStart( seqI - 1, 0 );
			else
				cr.SetStart( seqI - 1, -start );
			if( chars != stop - start + 1 && !(chars == 0 && stop - start == 1) ){
				std::cerr << "Error in XMFA file format\n";
				std::cerr << "Before line " << line_count << std::endl;
				std::cerr << "Expecting " << stop - start + 1 << " characters based on defline\n";
				std::cerr << "Actually read " << chars << " characters of sequence\n";
				Throw_gnEx(InvalidFileFormat());
			}
		}else{
			if( chars == 0 )
				cr.SetStart( seqI - 1, 0 );
			else
				cr.SetStart( seqI - 1, -stop );
			if( chars != start - stop + 1 && !(chars == 0 && stop - start == 1) ){
				std::cerr << "Error in XMFA file format\n";
				std::cerr << "Before line " << line_count << std::endl;
				std::cerr << "Expecting " << start - stop + 1 << " characters based on defline\n";
				std::cerr << "Actually read " << chars << " characters of sequence\n";
				Throw_gnEx(InvalidFileFormat());
			}
		}

		if( chars > max_len )
			max_len = aln_mat[ seqI - 1 ].length();
			
		if( cur_line.size() == 0 )
			break;
		// did we finish an aligned region?
		if( cur_line[ 0 ] != '>' ){
			GappedAlignment *new_cr = new GappedAlignment( aln_mat.size(), max_len );
			for( uint seqJ = 0; seqJ < seqI; seqJ++ ){
				new_cr->SetStart( seqJ, cr.Start( seqJ ) );
				new_cr->SetLength( lengths[ seqJ ], seqJ );
				cr.SetStart( seqJ, NO_MATCH );
			}
			for( uint seqJ = 0; seqJ < seqI; seqJ++ )
				aln_mat[seqJ].resize( max_len, '-' );

			new_cr->SetAlignment(aln_mat);
			lengths.clear();
			if( seq_count < seqI )
				seq_count = seqI;

			ga_list.push_back( new_cr );

			max_len = 0;	// reset length for the next interval
			aln_mat.clear();	// reset cr for next interval

			// bail out on EOF or corruption
			if( cur_line[ 0 ] != '=' )
				break;
			// otherwise read up to the next def. line
			while( std::getline( in_stream, cur_line ) ){
				line_count++;
				if( cur_line[ 0 ] == '>' )
					break;
			}
			if( cur_line[ 0 ] != '>' )
				break;
		}
	}

	// now process all GappedAlignments into Intervals
	for( uint ivI = 0; ivI < ga_list.size(); ivI++ ){
		GappedAlignment* cr = ga_list[ ivI ];
		GappedAlignment* new_cr = new GappedAlignment( seq_count, cr->AlignmentLength() );

		const std::vector< std::string >& align_matrix = GetAlignment( *cr, seq_table );
		std::vector< std::string > new_aln_mat(seq_count);
		for( seqI = 0; seqI < align_matrix.size(); seqI++ ){
			new_cr->SetLength( cr->Length( seqI ), seqI );
			new_cr->SetStart( seqI, cr->Start(seqI) );
			new_aln_mat[ seqI ] = align_matrix[ seqI ];
			if( new_aln_mat[ seqI ].length() == 0 )
				new_aln_mat[ seqI ] = std::string( new_cr->AlignmentLength(), '-' );
		}
		for( ; seqI < seq_count; seqI++ ){
			new_cr->SetLength( 0, seqI );
			new_cr->SetStart( seqI, 0 );
			new_aln_mat[ seqI ] = std::string( new_cr->AlignmentLength(), '-' );
		}
		new_cr->SetAlignment(new_aln_mat);
		delete cr;
		cr = new_cr;
		ga_list[ ivI ] = new_cr;

		std::vector<AbstractMatch*> asdf(1, cr);
		Interval iv( asdf.begin(), asdf.end() );
		this->push_back( iv );
	}
}

template< class MatchType >
void GenericIntervalList<MatchType>::ReadStandardAlignmentCompact( std::istream& in_stream ) 
{
	uint seq_count = 0;
	gnSeqI max_len = 0;
	std::string cur_line;
	std::getline( in_stream, cur_line );
	uint seqI = 0;
	std::vector< gnSeqI > lengths;
	std::vector< GappedAlignment* > ga_list;
	GappedAlignment cr;
	std::string empty_line;
	std::vector< std::string > aln_mat;
	uint line_count = 1;
	while( true ){
		
		while( cur_line[0] == '#' ){
			// hit a comment or metadata.  try to parse it if it's a filename
			std::getline( in_stream, cur_line );
			line_count++;
			std::stringstream ss( cur_line );
			std::string token;
			std::getline( ss, token, '\t' );
			if( token.substr(1, 8) != "Sequence" || token.find( "File" ) == std::string::npos )
				continue;
			std::getline( ss, token );
			seq_filename.push_back( token );
		}
		
		// read and parse the def. line
		std::stringstream line_str( cur_line );
		std::getline( line_str, cur_line, '>' );
		std::getline( line_str, cur_line, ':' );
		// take off leading whitespace
		std::stringstream parse_str( cur_line );

		parse_str >> seqI;	// the sequence number
				
		int64 start, stop;
		std::getline( line_str, cur_line, '-' );
		parse_str.clear();
		parse_str.str( cur_line );
		parse_str >> start;
		line_str >> stop;
		std::string strand;
		line_str >> strand;

		std::string name;	// anything left is the name
		std::getline( line_str, name );

		// read and parse the sequence
		while( aln_mat.size() < seqI )
			aln_mat.push_back( empty_line );

		gnSeqI chars = 0;
		while( std::getline( in_stream, cur_line ) ){
			line_count++;
			if( (cur_line[ 0 ] == '>' ) || (cur_line[ 0 ] == '=' ))
				break;
			for( uint charI = 0; charI < cur_line.length(); charI++ )
				if( cur_line[ charI ] != '-' )
					chars++;
			aln_mat[ seqI - 1 ] += cur_line;
		}
		while( lengths.size() < seqI )
			lengths.push_back(0);

		lengths[ seqI - 1 ] = chars;

// temporary workaround for file format inconsistency
		if( strand == "+" )
			cr.SetStart( seqI - 1, start );
		else if( start < stop ){
			if( chars == 0 )
				cr.SetStart( seqI - 1, 0 );
			else
				cr.SetStart( seqI - 1, -start );
			if( chars != stop - start + 1 && !(chars == 0 && stop - start == 1) ){
				std::cerr << "Error in XMFA file format\n";
				std::cerr << "Before line " << line_count << std::endl;
				std::cerr << "Expecting " << stop - start + 1 << " characters based on defline\n";
				std::cerr << "Actually read " << chars << " characters of sequence\n";
				Throw_gnEx(InvalidFileFormat());
			}
		}else{
			if( chars == 0 )
				cr.SetStart( seqI - 1, 0 );
			else
				cr.SetStart( seqI - 1, -stop );
			if( chars != start - stop + 1 && !(chars == 0 && stop - start == 1) ){
				std::cerr << "Error in XMFA file format\n";
				std::cerr << "Before line " << line_count << std::endl;
				std::cerr << "Expecting " << start - stop + 1 << " characters based on defline\n";
				std::cerr << "Actually read " << chars << " characters of sequence\n";
				Throw_gnEx(InvalidFileFormat());
			}
		}

		if( chars > max_len )
			max_len = aln_mat[ seqI - 1 ].length();
			
		if( cur_line.size() == 0 )
			break;
		// did we finish an aligned region?
		if( cur_line[ 0 ] != '>' ){
			GappedAlignment *new_cr = new GappedAlignment( aln_mat.size(), max_len );
			for( uint seqJ = 0; seqJ < seqI; seqJ++ ){
				new_cr->SetStart( seqJ, cr.Start( seqJ ) );
				new_cr->SetLength( lengths[ seqJ ], seqJ );
				cr.SetStart( seqJ, NO_MATCH );
			}
			for( uint seqJ = 0; seqJ < seqI; seqJ++ )
				aln_mat[seqJ].resize( max_len, '-' );

			new_cr->SetAlignment(aln_mat);
			lengths.clear();
			if( seq_count < seqI )
				seq_count = seqI;

			ga_list.push_back( new_cr );

			max_len = 0;	// reset length for the next interval
			aln_mat.clear();	// reset cr for next interval

			// bail out on EOF or corruption
			if( cur_line[ 0 ] != '=' )
				break;
			// otherwise read up to the next def. line
			while( std::getline( in_stream, cur_line ) ){
				line_count++;
				if( cur_line[ 0 ] == '>' )
					break;
			}
			if( cur_line[ 0 ] != '>' )
				break;
		}
	}

	// now process all GappedAlignments into Intervals
	//cerr << "Stuffing all GappedAlignments into Intervals" << endl;
	for( uint ivI = 0; ivI < ga_list.size(); ivI++ )
	{	
		GappedAlignment* cr = ga_list[ ivI ];
		uint compact_seq_count =  cr->SeqCount();
		CompactGappedAlignment<>* new_cr = new CompactGappedAlignment<>(compact_seq_count, cr->AlignmentLength() );
		const std::vector< std::string > align_matrix = GetAlignment( *cr, seq_table );
		//cout << cr->SeqCount() << " " << seq_count << " "  << align_matrix.size() << endl;
		
		std::vector< std::string > new_aln_mat(compact_seq_count);
		for( seqI = 0; seqI < compact_seq_count; seqI++ ){
			new_cr->SetLength( cr->Length( seqI ), seqI );
			new_cr->SetStart( seqI, cr->Start(seqI) );
			new_aln_mat[ seqI ] = align_matrix[ seqI ];
			if( new_aln_mat[ seqI ].length() == 0 )
				new_aln_mat[ seqI ] = std::string( new_cr->AlignmentLength(), '-' );
		}
		
		for( ; seqI < compact_seq_count; seqI++ ){
			new_cr->SetLength( 0, seqI );
			new_cr->SetStart( seqI, 0 );
			new_aln_mat[ seqI ] = std::string( new_cr->AlignmentLength(), '-' );
		}
		
		new_cr->SetAlignment( new_aln_mat );
		delete cr;

		//CompactGappedAlignment<>* cga =  new_cr;
		//ga_list[ ivI ] = dynamic_cast<GappedAlignment*>(cga);
		Interval iv;
		this->push_back( iv );
		std::vector< AbstractMatch* > matches(1, new_cr);
		this->back().SetMatches( matches );
	}
}


template< class MatchType >
void GenericIntervalList<MatchType>::WriteAlignedSequences(std::ostream& match_file) const
{

	unsigned int seq_count = seq_table.size();
	uint seqI;
	
	match_file << "mauveAligner data\n";
	match_file << "FormatVersion" << '\t' << 5 << "\n";
	match_file << "SequenceCount" << '\t' << seq_count << "\n";
	for(seqI = 0; seqI < seq_count; seqI++){
		match_file << "Sequence" << seqI << "File" << '\t';
		if( seq_filename.size() > seqI )
			match_file << seq_filename[seqI];
		else
			match_file << "null";
		match_file << "\n";
		match_file << "Sequence" << seqI << "Length" << '\t';
		if( seq_table.size() > seqI )
			match_file << seq_table[seqI]->length();
		else
			match_file << "0";
		match_file << "\n";
	}

	match_file << "AlignmentCount" << '\t' << this->size() << std::endl;

	if( this->size() == 0 )
		return;
	
	for( uint ivI = 0; ivI < this->size(); ivI++ ){
		
		match_file << (*this)[ ivI ].AlignmentLength();
		for( seqI = 0; seqI < seq_count; seqI++ )
			match_file << '\t' << (*this)[ ivI ].Start( seqI );
		match_file << std::endl;

		std::vector<std::string> alignment;
		GetAlignment( (*this)[ ivI ], this->seq_table, alignment );
		for( seqI = 0; seqI < seq_count; seqI++ )
			match_file << alignment[ seqI ] << std::endl;
		match_file << std::endl;
	}
	
}


}

#endif	//_IntervalList_h_
