/*******************************************************************************
 * $Id: gnAlignedSequences.cpp,v 1.11 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/gnAlignedSequences.h"
#include <sstream>

using namespace std;
using namespace genome;
namespace mems {

gnAlignedSequences::gnAlignedSequences()
{
	alignedSequenceFileName = "";
	
	
}


gnAlignedSequences::gnAlignedSequences(const gnAlignedSequences &toCopy)
{
	alignedSequenceFileName = toCopy.alignedSequenceFileName;
	consensus = toCopy.consensus;
	
	names = toCopy.names;
	sequences = toCopy.sequences;
	positions = toCopy.positions;
}


void gnAlignedSequences::constructFromClustalW(string alignedFileName)
{
	alignedSequenceFileName = alignedFileName;
	
	readClustalWAlignment();
	buildConsensus();
	
	indexPositions.resize(consensus.size());
	for (int i=0; i<consensus.size(); i++)
		indexPositions[i] = i+1;
}


void gnAlignedSequences::constructFromPhylip(string alignedFileName)
{
	alignedSequenceFileName = alignedFileName;
	
	readPhylipAlignment();
	buildConsensus();
	
	indexPositions.resize(consensus.size());
	for (int i=0; i<consensus.size(); i++)
		indexPositions[i] = i+1;
}


void gnAlignedSequences::constructFromMSF(string alignedFileName)
{
	alignedSequenceFileName = alignedFileName;
	
	readMSFAlignment();
	buildConsensus();
	
	indexPositions.resize(consensus.size());
	for (int i=0; i<consensus.size(); i++)
		indexPositions[i] = i+1;
}


void gnAlignedSequences::constructFromRelaxedNexus( istream& align_stream ){
	readRelaxedNexusAlignment( align_stream );
//	buildConsensus();
	
//	indexPositions.resize(consensus.size());
//	for (int i=0; i<consensus.size(); i++)
//		indexPositions[i] = i+1;
}

void gnAlignedSequences::constructFromNexus(string alignedFileName)
{
	alignedSequenceFileName = alignedFileName;
	
	readNexusAlignment();
	buildConsensus();
	
	indexPositions.resize(consensus.size());
	for (int i=0; i<consensus.size(); i++)
		indexPositions[i] = i+1;
}


void gnAlignedSequences::constructFromMega(string alignedFileName)
{
	alignedSequenceFileName = alignedFileName;
	
	readMegaAlignment();
	buildConsensus();
	
	indexPositions.resize(consensus.size());
	for (int i=0; i<consensus.size(); i++)
		indexPositions[i] = i+1;
}

const vector< string >& gnAlignedSequences::getSupportedFormats()
{
	static vector< string > formats;
	if( formats.size() == 0 ){
		formats.push_back( "phylip" );
		formats.push_back( "clustal" );
		formats.push_back( "msf" );
		formats.push_back( "nexus" );
		formats.push_back( "mega" );
		formats.push_back( "codon" );
	}
	return formats;
}

boolean gnAlignedSequences::isSupportedFormat( const string& format_name )
{
	const vector< string >& formats = getSupportedFormats();
	for( int formatI = 0; formatI < formats.size(); formatI++ ){
		if( formats[ formatI ] == format_name )
			return true;
	}
	return false;
}
void gnAlignedSequences::output( const string& format_name, ostream& os ) const
{
	bool rval = false;

	if( format_name == "phylip" )
		rval = outputPhylip( os );

	if( format_name == "clustal" )
		rval = outputClustalW( os );

	if( format_name == "msf" )
		rval = outputMSF( os );

	if( format_name == "nexus" )
		rval = outputNexus( os );

	if( format_name == "mega" )
		rval = outputMega( os );
	
	if( format_name == "codon" )
		rval = outputCodon( os );
	
	if( !rval )
		throw "Error writing alignment\n";

}

bool gnAlignedSequences::outputPhylip(ostream& os) const
{
	
	os << "Sequences in Alignment: " << sequences.size()
		<< "  Bases in Each Aligned Sequence: " << sequences[0].length() << endl;
	
	int offset = 10;
	uint seqI;
	for( seqI = 0; seqI < sequences.size(); seqI++ )
	{
		int position = 0;
		const string& seq = sequences[ seqI ];
		string seqName = names[ seqI ].substr( 0, offset );
		seqName.append( offset - seqName.length() + 1, ' ' ); 
		
		os << seqName;

		for ( position=0; position + offset < seq.size(); position += offset){
			if ( position % 50 == 0)
				os << endl;
			os.write( seq.data() + position, offset );
			os << ' ';
		}

		if ( position % 50 == 0)
			os << endl;

		os.write( seq.data() + position, seq.size() - position );
		os << endl;
	}
	
	return true;
}

uint64 countGaps( string& seq );
uint64 countGaps( string& seq ){
	uint gap_count = 0;
	for( uint charI = 0; charI < seq.length(); charI++ )
		if( seq[ charI ] == '-' )
			gap_count++;
	return gap_count;
}

bool gnAlignedSequences::outputClustalW(ostream& os) const
{
	boolean output_positions = true;
	
	os << "Clustal W multiple sequence alignment" << endl;
	
	vector< int64 > seq_pos( sequences.size(), 0 );
	if( positions.size() == sequences.size() )
		seq_pos = positions;
	vector< string > seq_names;
	int pos;
	uint seqI = 0;
	int longestNameSize = 0;
	for( ; seqI < sequences.size(); seqI++ )
	{
		seq_names.push_back( names[ seqI ].substr( 0, 30 ) );
		if ( seq_names[ seq_names.size() - 1 ].length() > longestNameSize)
			longestNameSize=seq_names[ seq_names.size() - 1 ].length();
	}
	// add space padding to the names
	for( seqI = 0; seqI < seq_names.size(); seqI++ )
		seq_names[ seqI ] += string( (longestNameSize - seq_names[ seqI ].length()) + 6, ' ' ); 
	for (pos=0; pos+60 < alignedSeqsSize(); pos+=60)
	{
		os << endl
		   << endl;
		for( seqI = 0; seqI < sequences.size(); seqI++ )
		{
			os << seq_names[ seqI ];
			const string& seq = sequences[ seqI ];
			string cur_seq = seq.substr( pos, 60 );
			os << cur_seq;
			if( output_positions ){
				seq_pos[ seqI ] += 60 - countGaps( cur_seq );
				os << " " << seq_pos[ seqI ];
			}
			os << endl;
		}
	}
	
	if (pos<alignedSeqsSize())
	{
		os << endl
		   << endl;
	
		for( seqI = 0; seqI < sequences.size(); seqI++ )
		{
			os << seq_names[ seqI ];
			const string& seq = sequences[ seqI ];
			string cur_seq = seq.substr( pos, 60 );
			os << cur_seq;
			if( output_positions ){
				seq_pos[ seqI ] += 60 - countGaps( cur_seq );
				os << " " << seq_pos[ seqI ];
			}
			os << endl;
		}
	}		   
	return true;
}


bool gnAlignedSequences::outputMSF(ostream& os) const
{
	os << "//" << endl;
	
	list <pair <string*, string*> >::const_iterator sequenceItr = alignedSequences.begin();
	int longestSeqNameLength = 0;
	for ( ; sequenceItr!=alignedSequences.end(); sequenceItr++)
	{
		if ((*(*sequenceItr).first).length() > longestSeqNameLength)
			longestSeqNameLength = (*(*sequenceItr).first).length();
	}
	
	int pos = 0;
	for ( ; pos+60<(*(*alignedSequences.begin()).second).size(); pos+=60)
	{
		// output spaces until sequence ordinates
		for (int i=0; i<longestSeqNameLength+2; i++)
			os << " ";
		
		os << pos+1;
		for (int i=0; i<54; i++) // output appropriate number of spaces on ordinate line
			os << " ";
		os << pos+60 << endl;
		
		for (sequenceItr=alignedSequences.begin(); sequenceItr!=alignedSequences.end(); sequenceItr++)
		{
			int spaces = longestSeqNameLength-(*(*sequenceItr).first).length();
			for (int i=0; i<spaces; i++)
				os << " ";
				
			os << (*(*sequenceItr).first) << "  ";
			
			string seq = (*(*sequenceItr).second).substr(pos, 60);
			for (int i=0; i<60; i++)
			{
				if (seq[i]=='-')
					os << ".";
				else
					os << seq[i];
			}
			os << endl;
		}
		
		os << endl;
	}
	
	if (pos<(*(*alignedSequences.begin()).second).size())
	{
		// output spaces until sequence ordinates
		for (int i=0; i<longestSeqNameLength+2; i++)
			os << " ";
		
		os << pos+1;
		for (int i=0; i<(*(*alignedSequences.begin()).second).size()-pos; i++) // output appropriate number of spaces on ordinate line
			os << " ";
		os << (*(*alignedSequences.begin()).second).size() << endl;
		
		for (sequenceItr=alignedSequences.begin(); sequenceItr!=alignedSequences.end(); sequenceItr++)
		{
			int spaces = longestSeqNameLength-(*(*sequenceItr).first).length();
			for (int i=0; i<spaces; i++)
				os << " ";
				
			os << (*(*sequenceItr).first) << "  ";
			
			string seq = (*(*sequenceItr).second).substr(pos, (*(*alignedSequences.begin()).second).size()-pos );
			for (int i=0; i<seq.length(); i++)
			{
				if (seq[i]=='-')
					os << ".";
				else
					os << seq[i];
			}
			os << endl;
		}
		
		os << endl;
	}
	
	return false;
}



bool gnAlignedSequences::outputNexus(ostream& os) const
{
	os << "begin data;" << endl
	   << "  dimensions ntax=" << sequences.size();
	if( sequences.size() == 0 )
		return true;
	os << " nchar=" 
	   << sequences[0].length() << ";" << endl
	   << "  ;" << endl
	   << "  matrix" << endl;
	   
	list <pair <string*, string*> >::const_iterator sequenceItr = alignedSequences.begin();
	int i;
	int seqI;
	int longestSeqNameLength = 0;
	for( seqI = 0; seqI < sequences.size(); seqI++ ){
		if( names[ seqI ].length() > longestSeqNameLength )
			longestSeqNameLength = names[ seqI ].length();
	}
	
	int pos = 1;
	for ( ; pos+59 < sequences[0].size(); pos+=60)
	{
		os << "[";
		// output spaces until sequence ordinates
		for (i = 0; i < longestSeqNameLength+2; i++)
			os << " ";
		
		os << pos;
		for (i = 0; i < 54; i++) // output appropriate number of spaces on ordinate line
			os << " ";
		os << pos+59 << "]" << endl;
		
		for( seqI = 0; seqI < sequences.size(); seqI++ )
		{
			os << names[ seqI ];
			
			int spaces = longestSeqNameLength - names[ seqI ].length();
			for (i = 0; i < spaces + 2; i++)
				os << " ";
				
			string seq = sequences[ seqI ].substr( pos, 60 );
			os << seq << endl;
		}
		
		os << endl;
	}
	
	// write out the last little bit
	if (pos - 1 < sequences[0].size())
	{
		// output spaces until sequence ordinates
		os << "[";
		// output spaces until sequence ordinates
		for (i = 0; i < longestSeqNameLength + 2; i++)
			os << " ";
		
		os << pos;
		for (i=0; i < sequences[0].size() - pos + 1; i++) // output appropriate number of spaces on ordinate line
			os << " ";
		os << pos+59 << "]" << endl;
		
		for (sequenceItr = alignedSequences.begin(); sequenceItr != alignedSequences.end(); sequenceItr++)
		for( seqI = 0; seqI < sequences.size(); seqI++ )
		{
			os << names[ seqI ];
			
			int spaces = longestSeqNameLength - names[ seqI ].length();
			for (i=0; i<spaces+2; i++)
				os << " ";
				
			string seq = sequences[seqI].substr( pos, sequences[seqI].size()-pos+1 );
			os << seq << endl;
		}
		
		os << endl;
	}
	
	return true;
}
/*
bool gnAlignedSequences::outputNexus(ostream& os) const
{
	os << "begin data;" << endl
	   << "  dimensions ntax=" << alignedSequences.size() << " nchar=" 
	   << alignedSequences.begin()->second->size() << ";" << endl
	   << "  ;" << endl
	   << "  matrix" << endl;
	   
	list <pair <string*, string*> >::const_iterator sequenceItr = alignedSequences.begin();
	int i;
	int longestSeqNameLength = 0;
	for ( ; sequenceItr != alignedSequences.end(); sequenceItr++)
	{
		if ( sequenceItr->first->length() > longestSeqNameLength )
			longestSeqNameLength = sequenceItr->first->length();
	}
	
	int pos = 1;
	for ( ; pos+59 < alignedSequences.begin()->second->size(); pos+=60)
	{
		os << "[";
		// output spaces until sequence ordinates
		for (i = 0; i < longestSeqNameLength+2; i++)
			os << " ";
		
		os << pos;
		for (i = 0; i < 54; i++) // output appropriate number of spaces on ordinate line
			os << " ";
		os << pos+59 << "]" << endl;
		
		for (sequenceItr=alignedSequences.begin(); sequenceItr != alignedSequences.end(); sequenceItr++)
		{
			os << (*(*sequenceItr).first);
			
			int spaces = longestSeqNameLength - sequenceItr->first->length();
			for (i = 0; i < spaces + 2; i++)
				os << " ";
				
			string seq = sequenceItr->second->substr( pos, 60 );
			for (i = 0; i < 60; i++)
				os << seq[i];
			os << endl;
		}
		
		os << endl;
	}
	
	if (pos - 1 < alignedSequences.begin()->second->size())
	{
		// output spaces until sequence ordinates
		os << "[";
		// output spaces until sequence ordinates
		for (i = 0; i < longestSeqNameLength + 2; i++)
			os << " ";
		
		os << pos;
		for (i=0; i < alignedSequences.begin()->second->size() - pos + 1; i++) // output appropriate number of spaces on ordinate line
			os << " ";
		os << pos+59 << "]" << endl;
		
		for (sequenceItr = alignedSequences.begin(); sequenceItr != alignedSequences.end(); sequenceItr++)
		{
			os << *(sequenceItr->first);
			
			int spaces = longestSeqNameLength-(*(*sequenceItr).first).length();
			for (i=0; i<spaces+2; i++)
				os << " ";
				
			string seq = (*(*sequenceItr).second).substr( pos, (*(*alignedSequences.begin()).second).size()-pos+1 );
			for (i=0; i<seq.length(); i++)
				os << seq[i];
			os << endl;
		}
		
		os << endl;
	}
	
	return false;
}
*/
bool gnAlignedSequences::outputMega(ostream& os) const
{
	os << "#MEGA" << endl
	   << "TITLE:" << endl;
	   
	list <pair <string*, string*> >::const_iterator sequenceItr = alignedSequences.begin();
	int longestSeqNameLength = 0;

	for ( ; sequenceItr!=alignedSequences.end(); sequenceItr++){
		if (sequenceItr->first->length() > longestSeqNameLength)
			longestSeqNameLength = sequenceItr->first->length();
	}
	
	gnSeqI pos = 1;
	gnSeqI remaining_len = alignedSequences.begin()->second->size();	//determine the amount to be written
	// loop while there is more to write
	while(remaining_len > 0){
		os << endl;
		gnSeqI write_chars = MEGA_ALIGN_COLUMNS < remaining_len ? MEGA_ALIGN_COLUMNS : remaining_len;

		//write each sequence's line
		for (sequenceItr = alignedSequences.begin(); sequenceItr != alignedSequences.end(); sequenceItr++){
			os << "#" << *(sequenceItr->first);
			
			int spaces = longestSeqNameLength - sequenceItr->first->length();
			for (int i = 0; i < spaces + 5; i++)
				os << " ";

			string seq = sequenceItr->second->substr( pos, write_chars );
			for (int i = 0; i < write_chars; i++)
				os << seq[i];
			os << endl;
		}
		os << endl;

		pos += write_chars;
		remaining_len -= write_chars;
	}
	return true;
}


bool gnAlignedSequences::outputCodon(ostream& os) const
{
	list <pair <string*, string*> >::const_iterator sequenceItr = alignedSequences.begin();
	
	os << '\t' << alignedSequences.size() << '\t' << (*(*sequenceItr).second).size() << endl;
	
	int offset = 10;
	for ( ; sequenceItr != alignedSequences.end(); sequenceItr++)
	{
		int position = 0;
		string seq = (*(*sequenceItr).second);
		string seqName = (*(*sequenceItr).first);
		if (seqName.size() <= offset) 
		{
			for (int i=seqName.size(); i<offset; i++)
				seqName += " ";
		}
		
		else 
		{
			string temp = seqName;
			seqName = "";
			for (int i=0; i<offset; i++)
				seqName += temp[i];
		}
		
		os << seqName;
		int count = 0;
		for ( ; position+3<seq.size(); position+=3)
		{
			if (count == 20)
			{
				count = 0;
				os << endl;
			}
			for (int i=position; i<position+3; i++)
			   os << seq[i];
			   
			os << ' ';
			count++;
		}
		
		for ( ; position < seq.size(); position++)
			os << seq[position];
		
		os << endl;
	}
	
	return false;
}


bool gnAlignedSequences::outputWithConsensus(ostream& os)
{
	list <pair <string*, string*> >::iterator sequenceItr = alignedSequences.begin();
	
	os << '\t' << alignedSequences.size() << '\t' << (*(*sequenceItr).second).size() << endl;
	
	int offset = 10;
	for ( ; sequenceItr != alignedSequences.end(); sequenceItr++)
	{
		int position = 0;
		string seq = (*(*sequenceItr).second);
		string seqName = (*(*sequenceItr).first);
		if (seqName.size() <= offset) 
		{
			for (int i=seqName.size(); i<offset; i++)
				seqName += " ";
		}
		
		else 
		{
			string temp = seqName;
			seqName = "";
			for (int i=0; i<offset; i++)
				seqName += temp[i];
		}
		
		os << seqName;
		int count = 0;
		for ( ; position+10<seq.size(); position+=10)
		{
			if (count == 5)
			{
				count = 0;
				os << endl;
			}
			for (int i=position; i<position+10; i++)
			   os << seq[i];
			   
			os << ' ';
			count++;
		}
		
		for ( ; position < seq.size(); position++)
			os << seq[position];
		
		os << endl;
	}
	
	int position = 0;
	int count = 0;
	os << "Consensus:";
	for ( ; position+10<consensus.size(); position+=10)
	{
		if (count == 5)
		{
			count = 0;
			os << endl;
		}
		for (int i=position; i<position+10; i++)
		   os << consensus[i];
		   
		os << ' ';
		count++;
	}
	
	for ( ; position < consensus.size(); position++)
		os << consensus[position];
	
	os << endl;
	
	return false;
}


gnAlignedSequences gnAlignedSequences::getAlignedSegment(unsigned start, unsigned stop)
{
	gnAlignedSequences newAlignment;
	
	addAllSegments(newAlignment, start, stop);
	newAlignment.buildConsensus();
	
	return newAlignment;
}


gnAlignedSequences gnAlignedSequences::getCodons(int readingFrame, int startCodon, int codonMultiple)
{
	gnAlignedSequences toReturn;
	int startBase = ((startCodon*3)-2)+(readingFrame-1);
	
	for (int index=startBase; (index+2)<(*(*alignedSequences.begin()).second).size(); index+=(codonMultiple*3))
		addAllSegmentsReplaceGaps(toReturn, index, index+2);
		
	toReturn.buildConsensus();
	
	return toReturn;
}


gnSeqI gnAlignedSequences::alignedSeqsSize() const
{
	if( sequences.size() > 0 )
		return sequences[ 0 ].size();
	return 0;
}


bool gnAlignedSequences::removeAlignedSeq(string seqName)
{
	list <pair <string*, string*> >::iterator sequenceItr = alignedSequences.begin();
	
	for ( ; sequenceItr != alignedSequences.end(); sequenceItr++)
	{
		if ((*(*sequenceItr).first) == seqName)
		{
			alignedSequences.erase(sequenceItr);
			return true;
		}
	}	
	
	return false;
}


bool gnAlignedSequences::removeAlignedSeq(unsigned index)
{
	list <pair <string*, string*> >::iterator sequenceItr = alignedSequences.begin();
	int i = 0;
	
	for ( ; sequenceItr != alignedSequences.end(); sequenceItr++)
	{
		if (i == index)
		{
			alignedSequences.erase(sequenceItr);
			return true;
		}
		
		i++;
	}	
	
	return false;
}


void gnAlignedSequences::concatenateAlignedSequences(gnAlignedSequences toConcat)
{
	list <pair <string*, string*> >::iterator toConcatItr = toConcat.alignedSequences.begin();
	list <pair <string*, string*> >::iterator originalItr;
	
	unsigned largestSeqSize = 0;
	
	for ( ; toConcatItr != toConcat.alignedSequences.end(); toConcatItr++)
	{
		for (originalItr = alignedSequences.begin(); originalItr != alignedSequences.end(); originalItr++)
		{
			if ((*(*originalItr).second).size() > largestSeqSize)
				largestSeqSize = (*(*originalItr).second).size();
			
			if ((*(*toConcatItr).first) == (*(*originalItr).first)) {
				string seq = (*(*originalItr).second);
				seq += (*(*toConcatItr).second);
				(*(*originalItr).second) = seq;
				break;
			}
		}
	}
	
	for (originalItr = alignedSequences.begin(); originalItr != alignedSequences.end(); originalItr++)
	{
		while ((*(*originalItr).second).size() < largestSeqSize)
			(*(*originalItr).second).append("-");
	}
	
	buildConsensus();
}


void gnAlignedSequences::extractVariableSites(gnAlignedSequences &variableSites, bool countGapsAsMismatches)
{
	list <pair <string*, string*> >::iterator originalItr = alignedSequences.begin();
	
	int alignedSeqSize = (*((*originalItr).second)).size();
	
	char positionBase;
	int matchStart = alignedSeqSize,
		matchStop = alignedSeqSize;
		
	bool mismatch = false;
	
	indexPositions.resize(0);
	
	for (int position=alignedSeqSize; position > 0; position--)
	{
		originalItr = alignedSequences.begin();
		positionBase = (*((*originalItr).second))[position-1];
		while (!countGapsAsMismatches && (*((*originalItr).second))[position-1] == '-')
		{
			originalItr++;
			positionBase = (*((*originalItr).second))[position-1];
			if (originalItr == alignedSequences.end()) break;
		}
		
		if (originalItr == alignedSequences.end()) break;
			
		for ( ; originalItr != alignedSequences.end(); originalItr++)
		{
			// extend matched segment before adding match to variableSites
			// much less expensive to add blocks of sites rather than a single site at a time
			if (positionBase != (*((*originalItr).second))[position-1])// && matchStop==position)
			{
				if (!(!countGapsAsMismatches && (*((*originalItr).second))[position-1] == '-'))
				{
					mismatch = true;
					break;
				}
			}
		}
		
		if (!mismatch)
			matchStart--;
		
		else
		{
			matchStart--;
			matchStop = matchStart;
			
			//variableSites.indexPositions.resize(variableSites.indexPositions.size()+1);
			variableSites.indexPositions.push_back(position);//[indexPositions.size()-1]=position;
		}
		
		mismatch = false;
	}
	
	for (int i=variableSites.indexPositions.size()-1; i>=0; i--)
		addAllSegments(variableSites, variableSites.indexPositions[i], variableSites.indexPositions[i]);
		
	variableSites.buildConsensus();
}


bool gnAlignedSequences::collapseIdenticalSequences()
{
	list <pair <string*, string*> >::iterator itr1 = alignedSequences.begin();
	list <pair <string*, string*> >::iterator itr2;
	bool toReturn = false;
	
	for ( ; itr1!=alignedSequences.end(); itr1++)
	{
		itr2=alignedSequences.begin();
		for (itr2++; itr2!=alignedSequences.end(); itr2++)
		{
			if (((*(*itr1).second)==(*(*itr2).second)) && itr1!=itr2)
			{
				list <pair <string*, string*> >::iterator itrTemp = itr2;
				itr2--;
				alignedSequences.erase(itrTemp);
				toReturn = true;
			}
		}
	}

	return toReturn;
}


vector <char> gnAlignedSequences::operator[]( const int offset ) //const
{
	vector <char> toReturn;
	list <pair <string*, string*> >::iterator itr;
	
	for (itr=alignedSequences.begin(); itr!=alignedSequences.end(); itr++)
		toReturn.push_back((*(*itr).second)[offset]);
	
	return toReturn;
}


bool gnAlignedSequences::readClustalWAlignment()
{
	ifstream alignmentFile;
	
	alignmentFile.open(alignedSequenceFileName.c_str(), ios::in | ios::binary);
	
	if (!(alignmentFile.is_open()))
	{
		cout << "Unable to open " << alignedSequenceFileName << ".\n"
			 << "Exiting.\n";
		
		exit(-1);
	}
	
	string line;
	
	// REMOVE 1st 3 LINES FROM .ALN FILE - SEQUENCE BEGINS ON LINE 4
	getline(alignmentFile, line);
	getline(alignmentFile, line);
	getline(alignmentFile, line);
	
	bool constructSuccess = constructClustalWAlignedSequenceList(alignmentFile);
	
	alignmentFile.close();
	
	if (constructSuccess) return true;
	
	return false;
}


bool gnAlignedSequences::readPhylipAlignment()
{
	ifstream alignmentFile;
	
	alignmentFile.open(alignedSequenceFileName.c_str(), ios::in | ios::binary);
	
	if (!(alignmentFile.is_open()))
	{
		cout << "Unable to open " << alignedSequenceFileName << ".\n"
			 << "Exiting.\n";
		
		exit(-1);
	}
	
	string line;
	
	// REMOVE 1st LINE FROM PHYLIP FILE - SEQUENCE NUMBER AND LENGTH OF SEQUENCES
	getline(alignmentFile, line);
	
	bool constructSuccess = constructPhylipAlignedSequenceList(alignmentFile);
	
	alignmentFile.close();
	
	if (constructSuccess) return true;
	
	return false;
}


bool gnAlignedSequences::readMSFAlignment()
{
	ifstream alignmentFile;
	
	alignmentFile.open(alignedSequenceFileName.c_str(), ios::in | ios::binary);
	
	if (!(alignmentFile.is_open()))
	{
		cout << "Unable to open " << alignedSequenceFileName << ".\n"
			 << "Exiting.\n";
		
		exit(-1);
	}
	
	string line;
	getline(alignmentFile, line);
	
	// remove format's initial annotation
	while (line.find("//")<0 || line.find("//")>line.size())
		getline(alignmentFile, line);
		
	bool constructSuccess = constructMSFAlignedSequenceList(alignmentFile);
	
	alignmentFile.close();
	
	if (constructSuccess) return true;
	
	return false;
}


/**
 * This function assumes that the #NEXUS at the beginning of the file has
 * been read off already.  It will read a single aligned sequences entry.
 */
bool gnAlignedSequences::readRelaxedNexusAlignment( istream& align_stream ){
	
	string line;
	string comments;
	getline( align_stream, line );
	if( line == "#NEXUS" ){
		getline( align_stream, line );
	}
	if( line[0] == '[' ){
		getline( align_stream, line );
		while( line[0] != ']' ){
			comments += line + "\n";
			getline( align_stream, line );
		}
		getline( align_stream, line );	// possibly empty line
		if( line.size() == 0 )
			getline( align_stream, line );
	}
	while( line.length() == 0 )
		getline( align_stream, line );
	// this is the alignment info line
	stringstream align_info( line );
	uint seq_count;
	gnSeqI align_len;
	align_info >> seq_count;
	align_info >> align_len;
	sequences = vector< string >( seq_count );
	// now read in each alignment line
	for( uint seqI = 0; seqI < seq_count; seqI++ ){
		align_stream >> line;
		names.push_back( line );
//		getline( align_stream, line );
		align_stream >> sequences[ seqI ];
//		Array< char > seq_data( align_len );
//		align_stream.read( seq_data.data, align_len );
//		sequences.push_back( seq_data.data );
	}
	
	// read off the trailing newline
	getline( align_stream, line );
	return true;
}


bool gnAlignedSequences::readNexusAlignment()
{
	ifstream alignmentFile;
	
	alignmentFile.open(alignedSequenceFileName.c_str(), ios::in | ios::binary);
	
	if (!(alignmentFile.is_open()))
	{
		cout << "Unable to open " << alignedSequenceFileName << ".\n"
			 << "Exiting.\n";
		
		exit(-1);
	}
	
	string line;
	getline(alignmentFile, line);
	
	// remove format's initial annotation
	while (line.find("begin data;")<0 || line.find("begin data;")>line.length()) // searching for "begin data;"
		getline(alignmentFile, line);
		
	bool constructSuccess = constructNexusAlignedSequenceList(alignmentFile);
	
	alignmentFile.close();
	
	if (constructSuccess) return true;
	
	return false;
}


bool gnAlignedSequences::readMegaAlignment()
{
	ifstream alignmentFile;
	
	alignmentFile.open(alignedSequenceFileName.c_str(), ios::in | ios::binary);
	
	if (!(alignmentFile.is_open()))
	{
		cout << "Unable to open " << alignedSequenceFileName << ".\n"
			 << "Exiting.\n";
		
		exit(-1);
	}
	
	string line;
	// remove first three lines from mega file - prior to begining of sequence data
	getline(alignmentFile, line);
	getline(alignmentFile, line);
	getline(alignmentFile, line);
		
	bool constructSuccess = constructMegaAlignedSequenceList(alignmentFile);
	
	alignmentFile.close();
	
	if (constructSuccess) return true;
	
	return false;
}


bool gnAlignedSequences::constructClustalWAlignedSequenceList(ifstream& alignmentFile)
{
	string line;
	
	// GET THE 1st LINE OF SEQUENCE
	getline(alignmentFile, line);
	
	while (alignmentFile.good())
	{
		while (line[0] != ' ' && line[0] != '\0')
		{
			string sequenceName;
			int i;
			for (i=0; line[i] != ' '; i++)
				sequenceName += line[i];
				
			const gnFilter* newFilter = gnFilter::fullDNASeqFilter();
			string sequenceBases;
			for(int i=sequenceName.size(); i < line.length(); i++){
				if ((*newFilter).IsValid(line[i]))
			    	sequenceBases += line[i];
			}
			
			list <pair <string*, string*> >::iterator sequenceItr; 
			if (!(sequenceNameInList(sequenceName, sequenceItr)))
			{
				pair <string*, string*> sequence;
				sequence.first = new string(sequenceName);
				
				sequence.second = new string( sequenceBases );
				
				alignedSequences.push_back(sequence);
			}
			
			else
				(*(*sequenceItr).second).append(sequenceBases);
	
			getline(alignmentFile, line);
		}
		
		getline(alignmentFile, line);
		getline(alignmentFile, line);
	}
	
	
	return false;
}


bool gnAlignedSequences::constructPhylipAlignedSequenceList(ifstream& alignmentFile)
{
	string line;
	
	// GET THE 1st LINE OF SEQUENCE
	getline(alignmentFile, line);
	
	while (alignmentFile.good())
	{
		if (line[10]!=' ')
		{
			string sequenceName = line.substr(0,10);
			cout << sequenceName << endl;
			const gnFilter* newFilter = gnFilter::fullDNASeqFilter();
			string sequenceBases;
			for(int i=10; i < line.length(); i++)
			{
				if ((*newFilter).IsValid(line[i]))
			    	sequenceBases += line[i];
			}
		
			pair <string*, string*> sequence;
			sequence.first = new string(sequenceName);
			
			sequence.second = new string( sequenceBases );
			
			alignedSequences.push_back(sequence);
			
			getline(alignmentFile, line);
		}

		// NOT THE 1st LINE IN SEQUENCE (CONTAINS SEQ NAME)
		else
		{
			string sequenceBases;
			while (line[10]==' ' && line[0]!='\0' && line.length()>0)
			{
				const gnFilter* newFilter = gnFilter::fullDNASeqFilter();
				for(int i=0; i < line.length(); i++)
				{
					if ((*newFilter).IsValid(line[i]))
				    	sequenceBases += line[i];
				}
				
				getline(alignmentFile, line);
			}
			
			list <pair <string*, string*> >::iterator sequenceItr = alignedSequences.end();
			sequenceItr--;
			(*(*sequenceItr).second) += sequenceBases;
		}
	}
	
	return false;
}


bool gnAlignedSequences::constructMSFAlignedSequenceList(ifstream& alignmentFile)
{
	string line;
	
	// clear coordinate line
	getline(alignmentFile, line);

	while (alignmentFile.good())//line[0] != '\0')
	{
		getline(alignmentFile, line); // 1st line of sequence
		while (!coordinates(line))
		{
			string sequenceName;
			int i;

			for (i=0; line[i] == ' '; i++) {}

			for (; line[i] != ' '; i++)
				sequenceName += line[i];
				
			const gnFilter* newFilter = gnFilter::fullDNASeqFilter();
			string sequenceBases;
			for( ; i < line.length(); i++){
				if ((*newFilter).IsValid(line[i]))
			    	sequenceBases += line[i];
			    else if (line[i] == '.' || line[i]=='~')
			    	sequenceBases += '-';
			}
			
			list <pair <string*, string*> >::iterator sequenceItr; 
			if (!(sequenceNameInList(sequenceName, sequenceItr)))
			{
				pair <string*, string*> sequence;
				sequence.first = new string(sequenceName);
				
				sequence.second = new string( sequenceBases );
				
				alignedSequences.push_back(sequence);
			}
			
			else
				(*(*sequenceItr).second).append(sequenceBases);
	
			getline(alignmentFile, line);
		}
	}	
	
	return false;
}


bool gnAlignedSequences::constructNexusAlignedSequenceList(ifstream& alignmentFile)
{
	string line;
	
	// GET THE 1st LINE OF SEQUENCE
	getline(alignmentFile, line);
	
	// searching for "endblock;"
	while (alignmentFile.good() && (line.find("endblock;")<0 || line.find("endblock;")>line.length())) 
	{
		while (line[0]!='[' && line[0]!=' ' && line[0]!='\n' && line[0]!='\r' && alignmentFile.good())
		{
			string sequenceName;
			int i;
			for (i=0; line[i] != ' '; i++)
				sequenceName += line[i];
				
			const gnFilter* newFilter = gnFilter::fullDNASeqFilter();
			string sequenceBases;
			for(int i=sequenceName.size(); i < line.length(); i++){
				if ( line[i] != '\r' && line[i] != '\n' && line[i] != ' ' )
			    	sequenceBases += line[i];
			}
			
			list <pair <string*, string*> >::iterator sequenceItr; 
			if (!(sequenceNameInList(sequenceName, sequenceItr)))
			{
				pair <string*, string*> sequence;
				sequence.first = new string(sequenceName);
				
				sequence.second = new string( sequenceBases );
				
				alignedSequences.push_back(sequence);
			}
			
			else
				(*(*sequenceItr).second).append(sequenceBases);
	
			getline(alignmentFile, line);
		}
		
		getline(alignmentFile, line);
	}
	
	
	return false;
}


bool gnAlignedSequences::constructMegaAlignedSequenceList(ifstream& alignmentFile)
{
	string line;
	string consensusSequenceBases;
	list <pair <string*, string*> >::iterator alignedSequencesItr;
	
	// GET THE 1st LINE OF SEQUENCE
	getline(alignmentFile, line);
	
	int previousLineLength = 0;
	
	// searching for "endblock;"
	while (alignmentFile.good()) 
	{
		while (line.length()>0 && line[0]=='#')
		{
			string sequenceName;
			for (int i=1; line[i] != ' '; i++)
				sequenceName += line[i];
				
			const gnFilter* newFilter = gnFilter::fullDNASeqFilter();
			string sequenceBases;
			bool isInSeqName = true;
			if (alignedSequences.size()>0)// && consensusSequenceBases.size()>0)
				consensusSequenceBases = (*(*alignedSequencesItr).second);
				
			for(int i=sequenceName.size(); i < line.length(); i++)
			{
				// allow only valid characters to be placed - if '.' replace leter
				// with consensus data
				if ((*newFilter).IsValid(line[i]) && !isInSeqName)
			    	sequenceBases += line[i];
			    	
			    else if (line[i] == ' ') isInSeqName=false;
			    	
			    else if (line[i]=='.' && alignedSequences.size()>0 && !isInSeqName) // a reference to the consensus
			    	sequenceBases += consensusSequenceBases[sequenceBases.size()+previousLineLength];
			}
			
			list <pair <string*, string*> >::iterator sequenceItr; 
			if (!(sequenceNameInList(sequenceName, sequenceItr)))
			{
				pair <string*, string*> sequence;
				sequence.first = new string(sequenceName);
				
				sequence.second = new string( sequenceBases );
				
				alignedSequences.push_back(sequence);
				alignedSequencesItr = alignedSequences.begin();
			}
			
			else
				(*(*sequenceItr).second).append(sequenceBases);
	
			getline(alignmentFile, line);
		}
		
		if (alignedSequences.size() > 0)
			previousLineLength = (*(*alignedSequences.begin()).second).size();
		
		getline(alignmentFile, line);
	}
	
	
	return false;
}


int gnAlignedSequences::sequenceNameInList( string& sequenceName ){
	for( uint nameI = 0; nameI < names.size(); nameI++ ){
		if( sequenceName == names[ nameI ] )
			return nameI;
	}
	return -1;
}

bool gnAlignedSequences::sequenceNameInList(string sequenceName, list <pair <string*, string*> >::iterator &sequenceItr)
{
	for (sequenceItr = alignedSequences.begin(); sequenceItr != alignedSequences.end(); sequenceItr++)
	{
		if (sequenceName == (*(*sequenceItr).first))
			return true;
	}
	
	return false;
}


bool gnAlignedSequences::buildConsensus()
{
	char consensusBase = '-';

	consensus = "";
	
	vector <char> crossAlignmentBases;
	for (int index=0; index<(*(*alignedSequences.begin()).second).size(); index++)
	{
		vector <int> baseCounts(26, 0);
		crossAlignmentBases = (*this)[index];
		/*list <pair <string*, string*> >::iterator itr = alignedSequences.begin();
		itr++;*/
		for (int i=0; i<crossAlignmentBases.size(); i++)
		{
			// to hold knowledge of consensus if MEGA '.' format employed
			// ('.'==same as base in 1st sequence)
			if (i == 0)
				consensusBase=crossAlignmentBases[i];
			
			// consensus already established if in MEGA '.' format - the 1st seq	
			if (i>0 && crossAlignmentBases[i]=='.')
				break;
		
			else if (crossAlignmentBases[i] != '-')
			{
				int baseIndex = determineBaseIndex(crossAlignmentBases[i]);
				baseCounts[baseIndex]++;
			}
		}
		
		int toAppendToConsensus = 0;
		for (int i=1; i<baseCounts.size(); i++)
		{
			// strictly alphabetic - count ties are broken lexigraphically
			if (baseCounts[i] > baseCounts[toAppendToConsensus])
				toAppendToConsensus = i;

			/* nearly functional code for replacing '.'s w/ consensus data
			if (crossAlignmentBases[i]=='.')
			{
				(*(*itr).second).erase(index, 1);
				string toInsert;
				toInsert += crossAlignmentBases[0];
				(*(*itr).second).insert(index, toInsert);
			}

			itr++;*/
		}
		
		consensus += (toAppendToConsensus+65);
	}

	return false;
}


void gnAlignedSequences::addSequence(string& seqToAdd, string& seqName)
{
	sequences.push_back( seqToAdd );
	names.push_back( seqName );
}


void gnAlignedSequences::addSequence(gnSequence& seqToAdd, string& seqName)
{

	ErrorMsg( "Fix gnAlignedSequences::addSequence()" );
	sequences.push_back( seqToAdd.ToString() );
	names.push_back( seqName );

/*	list <pair <string*, string*> >::iterator itr;
	if (!sequenceNameInList(seqName, itr))
	{
		pair <string*, string*> toAdd;
		toAdd.first = new string(seqName);
		toAdd.second = new string( seqToAdd.ToString() );
		
		alignedSequences.push_back(toAdd);
	}

	else
	{
		(*((*itr).second)) += seqToAdd.ToString();
	}
*/
}


void gnAlignedSequences::addSequence(gnSequence seqToAdd, string seqName, int consensusStart, string originalConsensus)
{
	list <pair <string*, string*> >::iterator itr;
	if (!sequenceNameInList(seqName, itr))
	{
		pair <string*, string*> toAdd;
		toAdd.first = new string(seqName);
		string seq = seqToAdd.ToString();
		toAdd.second = new string( seq );
		(*toAdd.second).erase();
		
		for (int i=0; i<(*toAdd.second).size(); i++)
		{
			if (seq[i] == '-')
				seq[i] = originalConsensus[consensusStart+i-1];
		}
		
		(*toAdd.second) = seq;
		
		alignedSequences.push_back(toAdd);
	}

	else
	{
		string seq = (*((*itr).second));
		seq += seqToAdd.ToString();
		for (int i=0; i<seq.size(); i++)
		{
			if (seq[i+(*((*itr).second)).size()]=='-' && originalConsensus.size()>0)
				seq[i+(*((*itr).second)).size()] = originalConsensus[consensusStart+i-1];
		}
		
		(*((*itr).second)) = seq;
	}
}


void gnAlignedSequences::addAllSegments(gnAlignedSequences &alignment, unsigned start, unsigned stop)
{
	for ( uint seqI = 0; seqI < alignment.sequences.size(); seqI++ ){
		if (stop == 0 || stop == alignment.sequences[ seqI ].size()-1)
			stop = alignment.sequences[ seqI ].size();
		string seq = alignment.sequences[ seqI ].substr(start, stop-start+1);
		alignment.addSequence( seq, alignment.names[ seqI ] );

	}
}


void gnAlignedSequences::addAllSegmentsReplaceGaps(gnAlignedSequences &alignment, unsigned start, unsigned stop)
{
	list <pair <string*, string*> >::iterator alignedItr = alignedSequences.begin();
	for ( ; alignedItr != alignedSequences.end(); alignedItr++)
	{
		if (stop == 0 || stop == (*(*alignedItr).second).size()-1)
			stop = (*(*alignedItr).second).size();
			
		alignment.addSequence(((*(*alignedItr).second).substr(start, stop-start+1)), 
							  ((*(*alignedItr).first)), start, consensus);
	}
}


void gnAlignedSequences::removeAllSegments(unsigned start, unsigned stop)
{
	list <pair <string*, string*> >::iterator alignedItr = alignedSequences.begin();
	for ( ; alignedItr != alignedSequences.end(); alignedItr++)
	{
		if (stop == 0)
			stop = (*(*alignedItr).second).size();
			
		(alignedItr->second)->erase(start, stop-start+1);
	}

	cout << start << " " << stop << ": " << stop-start+1 << endl;
}


int gnAlignedSequences::determineBaseIndex(char base)
{
	if (base < 91) // Upper Case
		return (base-65);
		
	// Lower Case
	return (base-97);
}


bool gnAlignedSequences::coordinates(string line)
{
	bool toReturn = true;
	
	for (int i=0; i<line.length(); i++)
	{
		if (line[i]!=' ' && line[i]!='\r' && line[i]!='\n' && (line[i]<48 || line[i]>57))
		{
			toReturn = false;
			break;
		}
	}
	
	return toReturn;
}

}
