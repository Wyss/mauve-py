/////////////////////////////////////////////////////////////////////////////
// File:            gnSequence.h
// Purpose:         Sequence class
// Description:     Provides a high level sequence interface to all types of
//					sequence data.
// Changes:        
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Modified by:     
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details 
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnSequence_h_
#define _gnSequence_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <cstdlib>
#include <iostream>
#include <list>
#include "libGenome/gnBaseSource.h"
#include "libGenome/gnGenomeSpec.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnStringSpec.h"
#include "libGenome/gnBaseFilter.h"
#include "libGenome/gnCompare.h"

namespace genome {


/**
 * gnSequence is the most commonly used class in libGenome.  It provides a simple
 * and general way to manipulate genetic sequences, regardless of what kind of
 * database, web site, or file they are stored in.  Sequence data can be
 * manipulated like a c++ std::string by using subseq() and erase().  gnSequence
 * updates annotated sequences with each change, breaking features if necessary.
 */
class GNDLLEXPORT gnSequence : public gnClone
{
public:
	/**
	 * Empty Constructor, creates an empty gnSequence.
	 */
	gnSequence();
	/**
	 * Creates a gnSequence with a single fragment containing the bases in "seq".
	 * @param seq The null terminated character array of base pairs to use.
	 */
	gnSequence( const gnSeqC* seq );
	/**
	 * Creates a gnSequence with a single fragment containing the bases in "str".
	 * @param str The base pairs to use.
	 */
	gnSequence( const std::string& str );
	/**
	 * Creates a gnSequence with the genome stored in the gnGenomeSpec "gngs".
	 * @param gngs the gnGenomeSpec to get contigs from.
	 */
	gnSequence( const gnGenomeSpec& gngs );
	/**
	 * Creates a gnSequence with the sequence fragment stored in the gnFragmentSpec "gnfs".
	 * @param gnfs the gnFragmentSpec to get contigs from.
	 */
	gnSequence( const gnFragmentSpec& gnfs );
	/**
	 * Creates a gnSequence with the contig stored in the gnContigSpec "gncs".
	 * Usually gncs will be a gnStringSpec or a gnSourceSpec.
	 * @param gncs the gnContigSpec to get contigs from.
	 */
	gnSequence( const gnContigSpec& gncs );
	/**
	 * Creates a gnSequence with a single fragment containing the bases in "bases".
	 * @param bases The base pairs to use
	 * @param length The length of the base pair array.
	 */
	gnSequence( gnSeqC *bases, const gnSeqI length);
	/**
	 * Copies the gnSequence "seq".
	 * @param seq The gnSequence to copy.
	 */
	gnSequence( const gnSequence& seq);
	/**
	 * Destructor, frees the memory used by this sequence.
	 */
	~gnSequence();

	/**
	 * Copies the gnSequence.
	 */
	gnSequence& operator=( const gnSequence& seq);

	gnSequence* Clone() const;

	/**
	 * Returns the number of sequence fragments in this sequence.
	 * @return the number of sequence fragments in this sequence.
	 */
	gnSeqI contigListSize() const;
	/**
	 * Returns the number of sequence fragments in this sequence.
	 * @return the number of sequence fragments in this sequence.
	 */
	gnSeqI contigListLength() const;
	/**
	 * Returns the index of the contig which contains the specified base.
	 * @param baseI A base pair in the contig.
	 * @throws SeqIndexOutOfBounds will be propagated if baseI is invalid
	 * @return The contig index which contains the base.
	 */
	uint32 contigIndexByBase( const gnSeqI baseI) const;
	/**
	 * Returns a gnSequence containing the specified fragment.
	 * @param contigI The index of the fragment to get.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 * @return A gnSequence containing only the fragment.
	 */
	gnSequence contig( const uint32 contigI) const;
	/**
	 * Returns a gnSequence containing only the fragment which contains the specified base.
	 * @param baseI A base pair index in the contig
	 * @throws SeqIndexOutOfBounds will be propagated if baseI is invalid
	 * @return A gnSequence containing only the contig.
	 */
	gnSequence contigByBase( const gnSeqI baseI) const;
	/**
	 * Returns the index of the base pair where the specified contig starts in this sequence.
	 * @param contigI The index of the contig.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 * @return The starting base pair index.
	 */
	virtual gnSeqI contigStart( const uint32 contigI) const;
	/**
	 * Returns the length in base pairs of the specified contig in this sequence.
	 * @param contigI The index of the contig.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 * @return The length of the contig.
	 */
	virtual gnSeqI contigLength( const uint32 contigI) const;
	/**
	 * Returns the index of the contig with the given name.  
	 * If two contigs have the same name, contigIndexByName() will return the first.
	 * @param contigName The name of the contig.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigName is not found
	 * @return The index of the contig.
	 */
	virtual uint32 contigIndexByName( std::string& contigName) const;
	/**
	 * Returns the name of the specified contig.  
	 * @param contigI The index of the contig.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 * @return The name of the contig.
	 */
	virtual std::string contigName( const uint32 contigI) const;
	/**
	 * Returns a gnSequence containing only the named contig.
	 * If two contigs have the same name, contigByName() will return the first.
	 * @param contigName The name of the contig.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigName is not found
	 * @return A gnSequence containing the contig.
	 */
	virtual gnSequence contigByName( std::string& contigName) const;
	/**
	 * Merges the bases starting at base index "startI" and ending at "endI" into one contig, splitting existing contigs.
	 * @param startI The starting base pair of the new contig.
	 * @param endI The ending base pair of the new contig.
	 * @throws SeqIndexOutOfBounds exception is propagated if startI or endI are invalid
	 */
	virtual void merge(const gnSeqI startI, const gnSeqI endI);
	/**
	 * Merges the contigs starting at the contig index "startC" and ending at "endC" into one contig.
	 * @param startC The starting contig of the new contig.
	 * @param endC The ending contig of the new contig.
	 * @throws FragmentIndexOutOfBounds exception is propagated if startC or endC are invalid
	 */
	virtual void mergeContigs(const uint32 startC, const uint32 endC);
	/**
	 * Splits the specified contig after splitI
	 * @param splitI The base pair to split after.
	 * @param contigI The index of the contig to split or ALL_CONTIGS by default.
	 */
	virtual void splitContig(const gnSeqI splitI, const uint32 contigI=ALL_CONTIGS);

	virtual void setContigName( const uint32 contigI, const std::string& contig_name);

	/**
	 * Returns the size of the feature list for the specified contig.
	 * @return The feature list size.
	 */
	virtual uint32 getFeatureListLength() const;
	/**
	 * Returns the feature specified by featureI.
	 * @param featureI The index of the feature to return.
	 * @throws FeatureIndexOutOfBounds exception is thrown if featureI does not reference a valid feature
	 * @return A copy of the feature allocated on the heap.  You must delete the feature when finished with it.
	 */
	virtual gnBaseFeature* getFeature(const uint32 featureI) const;
	/**
	 * Creates a list of all features which are contained by coordinates specified.
	 * @param lt The coordinates containing the features to return.
	 * @param feature_vector The std::vector to store features in.  Each feature pointer points to a copy of the feature on the heap which must be deleted.
	 * @param index_vector A std::vector of indices which correspond to the features.
	 */
	virtual void getContainedFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const;
	/**
	 * Creates a list of all features which intersect the coordinates specified.
	 * @param lt The coordinates intersecting the features to return.
	 * @param feature_vector The std::vector to store features in.  Each feature pointer points to a copy of the feature on the heap which must be deleted.
	 * @param index_vector A std::vector of indices which correspond to the features.
	 */
	virtual void getIntersectingFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const;
	/**
	 * Adds the specified feature to the feature list.
	 * @param feature The feature to add.
	 * @return The feature index.
	 */
	virtual uint32 addFeature(gnBaseFeature* feature);
	/**
	 * Removes the specified feature.
	 * @param featureI The index in the feature list of the feature to remove.
	 * @throws FeatureIndexOutOfBounds exception is thrown if featureI does not reference a valid feature
	 */
	virtual void removeFeature(const uint32 featureI);
	/**
	 * Creates a list of features which may have been broken by an edit.
	 * @param lt The coordinates containing the features to return.
	 * @param feature_vector The std::vector to store features in.  Each feature pointer points to a copy of the feature on the heap which must be deleted.
	 */
	virtual void getBrokenFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector) const;
	/**
	 * Returns the size of the header list for the specified contig.
	 * @param contigI The index of the contig to use, or ALL_CONTIGS to add a general header.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 * @return The header list size.
	 */
	virtual uint32 getHeaderListLength(const uint32 contigI) const;
	/**
	 * Returns the feature specified by featureI.
	 * @param contigI The index of the contig to use, or ALL_CONTIGS to get a general header.
	 * @param headerI The index of the header to return.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 * @throws HeaderIndexOutOfBounds will be propagated if headerI is invalid
	 * @return The header, or NULL if headerI is out of range.
	 */
	virtual gnBaseHeader* getHeader(const uint32 contigI, const uint32 headerI) const;
	/**
	 * Adds header information to a specified contig.
	 * @param contigI The index of the contig to use, or ALL_CONTIGS to add a general header.
	 * @param header The header to add.
	 * @param headerI The index in the header list where this header should be inserted.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 */
	virtual void addHeader(const uint32 contigI, gnBaseHeader* header, const uint32 headerI);
	/**
	 * Removes header information from a specified contig.
	 * @param contigI The index of the contig to use, or ALL_CONTIGS to remove a general header.
	 * @param headerI The index in the header list of the header to remove.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 */
	virtual void removeHeader(const uint32 contigI, const uint32 headerI);
	/**
	 * Reverse complements a specified contig, or the entire sequence if ALL_CONTIGS is specified.
	 * @param revComp True if the area should be reverse complement, false otherwise.
	 * @param contigI The index of the contig to use, or ALL_CONTIGS.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 */
	virtual void setReverseComplement( const boolean revComp, const uint32 contigI=ALL_CONTIGS);
	/**
	 * Returns true if a specified contig, or the entire sequence is reverse complement.
	 * @param contigI The index of the contig to use, or ALL_CONTIGS for the whole sequence.
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 * @return True if the sequence is reverse complement, false otherwise.
	 */
	virtual boolean isReverseComplement( const uint32 contigI=ALL_CONTIGS );
	/**
	 * Returns true if this sequence is circular.
	 * @return True if this sequence is circular.
	 */
	virtual boolean isCircular() const;
	/**
	 * Sets whether this sequence should be read circular.
	 * If circular is set, reads beyond the end of the sequence will pick up
	 * at the beginning.
	 * @param value True for circular, false otherwise.
	 */
	virtual void setCircular( const boolean value );
	
	/**
	 * Converts the global sequence coordinate baseI to a contig local coordinate.
	 * @param contigI This is set to the index of the contig containing baseI
	 * @param baseI This is the global coordinate to be converted to contig local.
	 * @throws SeqIndexOutOfBounds will be thrown if baseI is out of range
	 */
	virtual void globalToLocal(uint32& contigI, gnSeqI& baseI) const;

	/**
	 * Converts the local contig coordinate baseI to a global sequence coordinate.
	 * @param contigI The index of the contig containing baseI.
	 * @param baseI The local contig coordinate to be converted to global.
	 * @throws SeqIndexOutOfBounds will be thrown if baseI is out of range
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 */
	virtual void localToGlobal(const uint32 contigI, gnSeqI& baseI) const;
	/**
	 * Converts the global sequence coordinate baseI to a local coordinate in the original
	 * data source.  globalToSource() will overwrite any values passed to it!
	 * @param contigI This will be set to contig index in the original source.
	 * @param baseI This will be set to the contig local base index in the original source.
	 * @throws SeqIndexOutOfBounds will be thrown if baseI is out of range
	 */
	virtual void globalToSource(uint32& contigI, gnSeqI& baseI) const;
	/**
	 * Converts the contig local sequence coordinate baseI in contig contigI to a
	 * local coordinate in the original data source.  localToSource() will overwrite any values
	 * passed to it!
	 * @param contigI This will be set to contig index in the original source.
	 * @param baseI This will be set to the contig local base index in the original source.
	 * @throws SeqIndexOutOfBounds will be thrown if baseI is out of range
	 * @throws FragmentIndexOutOfBounds will be propagated if contigI is invalid
	 */
	virtual void localToSource(uint32& contigI, gnSeqI& baseI) const;
	/**
	 * Loads the sequence located at the URL in "sourcename".
	 * Possible URLs currently include only "file:///" URLs.
	 * If no URL prefix is found then LoadSource assumes that "sourcename"
	 * contains the name of a local file.
	 * @param sourcename The location of the data to load.
	 * @return True if successful.  False otherwise.
	 */
	virtual bool LoadSource(const std::string sourcename);

	/** Assigns a filter which all sequence data must pass through when
	 *  read from the object.
	 *  @param filt The filter to use.
	 */
	virtual void setFilter(const gnBaseFilter* filt);

	/** Assigns a list of filters which all sequence data passes through
	 *  in order when read from the object.  There may not be any NULL
	 *  pointers in the list.
	 *  @param filt_list The ordered list of filters to use.
	 */
	virtual void setFilterList(std::list<const gnBaseFilter*>& filt_list);

	/** Returns the list of filters currently being used
	 *  @return The list of filters in use.
	 */
	virtual std::list<const gnBaseFilter*> getFilterList() const;

	/**
	 * Assigns the sequence "seq" to this sequence.
	 * @param seq The seqence to assign to this sequence.
	 */
	virtual void assign(gnSequence& seq);

// Comparison Operators
	/**
	 * Assigns the sequence "seq" to this sequence.
	 */
	void operator=(gnSequence& seq);
	boolean operator==(const gnSequence& seq) const;
	boolean operator!=(const gnSequence& seq) const;
	boolean operator<(const gnSequence& seq) const;
	boolean operator<=(const gnSequence& seq) const;
	boolean operator>(const gnSequence& seq) const;
	boolean operator>=(const gnSequence& seq) const;
	/**
	 * Appends the bases in "seq" to this sequence.
	 */
	gnSequence& operator+=(const gnSequence& seq);
    /**
	 * Compares the bases in "seq" to this sequence.
	 * @param seq The sequence to compare this sequence to.
	 * @return Negative if this sequence is lesser, 0 if the two sequences are
	 * equal, and positive if this sequence is greater.
	 */
	virtual int compare(const gnSequence& seq) const;
	virtual int compare(const std::string& str) const;
// Append and Insert Operators
	/**
	 * Appends the bases in "seq" to this sequence.
	 * @param seq The sequence to append to this sequence.
	 */
	virtual void append( const gnSequence& seq);
	/**
	 * Inserts the first "len" bases in "bases" into this sequence.at "offset".
	 * insert() will update the locations of all affected features.
	 * @param offset The base pair to insert before.
	 * @param bases The character array of bases to insert.
	 * @param length The length of the character array.
	 */
	virtual void insert( const gnSeqI offset, const gnSeqC *bases, const gnSeqI length);
	/**
	 * Inserts the annotated sequence in "seq" into this sequence.at "offset".
	 * insert() will update the locations of all affected features.
	 * @param offset The base pair to insert before.
	 * @param seq The sequence to insert.
	 */
	virtual void insert( const gnSeqI offset, const gnSequence& seq);
	/**
	 * Inserts the annotated sequence in "gnbs" into this sequence.at "offset".
	 * insert() will update the locations of all affected features.
	 * @param offset The base pair to insert before.
	 * @param gnbs The spec to insert.
	 */
	virtual void insert( const gnSeqI offset, const gnGenomeSpec& gnbs);
// Concatenation Operators
	/**
	 * Concatenates this sequence with the annotated sequence in "seq".
	 */
	gnSequence const operator+(const gnSequence& seq) const;
// Substrings and base deletion
	/**
	 * Creates a sequence containing the "length" bases starting at "offset".
	 * @param offset The base pair index to start the subsequence.
	 * @param length The length of the subsequence.
	 * @return The subsequence.
	 */
	gnSequence subseq(const gnSeqI offset, const gnSeqI length) const;
	/**
	 * Deletes the "len" bases starting at "offset".
	 * @param offset The base pair index to start erasing.
	 * @param length The length of the sequence to erase
	 */
	virtual void erase( const gnSeqI offset=0, const gnSeqI length=GNSEQI_END );
// IO operators
	/**
	 * Reads bases from the specified input stream (e.g. cin).
	 */
	friend std::istream& operator>>(std::istream& is, gnSequence& gns);	//read from source.
	/**
	 * Writes the bases in this sequence to the specified output stream (e.g. cout).
	 */
	friend std::ostream& operator<<(std::ostream& os, const gnSequence& gns); //write to source.
// Size functions
	/**
	 * Returns the length of this sequence.
	 * @return the length of this sequence.
	 */
	virtual gnSeqI length() const;
	/**
	 * Returns the length of this sequence.
	 * @return the length of this sequence.
	 */
	virtual gnSeqI size() const;
// Raw Sequence Access
	/**
	 * Returns the "length" bases starting at "offset" as a std::string.
	 * @param length The length of the sequence to convert
	 * @param offset The base pair index of the sequence to convert
	 * @return The std::string of base pairs.
	 */
	virtual std::string ToString( const gnSeqI length=GNSEQI_END, const gnSeqI offset=1 ) const;
	/**
	 * Converts the "length" bases starting at "offset" into the std::string "str".
	 * @param str The std::string to store bases in.
	 * @param length The length, in base pairs, to convert.
	 * @param offset The base pair index to start converting.
	 * @return True if successful, false otherwise.
	 */
	virtual boolean ToString( std::string& str, const gnSeqI length=GNSEQI_END, const gnSeqI offset=1 ) const;
	/**
	 * Converts the "length" bases starting at "offset" into the character array "pSeqC"..
	 * After converting, "length" will be set to the actual length of the sequence.
	 * Be sure to null terminate the character array if you are going to print it!
	 * @param pSeqC The character array of bases to store bases in.
	 * @param length The length, in base pairs, to convert.
	 * @param offset The base pair index to start converting.
	 * @return True if successful, false otherwise.
	 */
	virtual boolean ToArray( gnSeqC* pSeqC, gnSeqI length, const gnSeqI offset=1 ) const;
	/**
	 * Returns the base at "offset".
	 * @param offset The index of the base to get.
	 * @return The base.
	 */
	virtual gnSeqC GetSeqC( const gnSeqI offset ) const; // return gnSeqC illegal char.
	/**
	 * Returns the base at the specified index.
	 * @return The base at the specified index.
	 */
	gnSeqC operator[]( const gnSeqI offset ) const;

	/**
	 * Get the spec (annotated sequence) which this sequence represents.
	 * @return The spec represented by this gnSequence.
	 */
	virtual gnGenomeSpec* GetSpec() const;
	/**
	 * Set the spec (annotated sequence) which this sequence represents.
	 * @param s The new spec.
	 */
	virtual void SetSpec(gnGenomeSpec* s);
	
	/**
	 * Find looks for the search sequence within this gnSequence and returns
	 * the position of the first match if any exists.
	 * @param search The gnSequence to be found
	 * @param offset The position where searching will begin (default is 0)
	 * @return The first position where the search std::string is found or
	 * GNSEQI_ERROR if the search sequence is not found.
	 */
	virtual gnSeqI find(const gnSequence& search, const gnSeqI offset=0) const;
	
private:
	gnGenomeSpec *spec;
	std::list<const gnBaseFilter*> filter_list;
	const gnCompare* comparator;
}; // class gnSequence


GNDLLEXPORT
std::istream& operator>>(std::istream& is, gnSequence& gns);	//read from source.
GNDLLEXPORT
std::ostream& operator<<(std::ostream& os, const gnSequence& gns); //write to source.

// Assignment Operators
inline
void gnSequence::operator=(gnSequence& seq){
	spec = seq.spec->Clone();
}
inline
void gnSequence::assign(gnSequence& seq){
	spec = seq.spec->Clone();
}
inline
boolean gnSequence::operator==(const gnSequence& seq) const{
	return (compare(seq)==0);
}
inline
boolean gnSequence::operator!=(const gnSequence& seq) const{
	return (compare(seq)!=0);
}
inline
boolean gnSequence::operator<(const gnSequence& seq) const{
	return (compare(seq)<0);
}
inline
boolean gnSequence::operator<=(const gnSequence& seq) const{
	return (compare(seq)<=0);
}
inline
boolean gnSequence::operator>(const gnSequence& seq) const{
	return (compare(seq)>0);
}
inline
boolean gnSequence::operator>=(const gnSequence& seq) const{
	return (compare(seq)>=0);
}
// Append and Insert Operators
inline
gnSequence& gnSequence::operator+=(const gnSequence& seq){
	insert(GNSEQI_END, *seq.spec);
	return *this;
}
inline
void gnSequence::append( const gnSequence& seq){
	insert(GNSEQI_END, *seq.spec);
}
//Feature functions
inline
uint32 gnSequence::getFeatureListLength() const{
	return spec->GetFeatureListLength();
}
inline
gnBaseFeature* gnSequence::getFeature(const uint32 featureI) const{
	return spec->GetFeature(featureI);
}
inline
uint32 gnSequence::addFeature(gnBaseFeature* feature){
	return spec->AddFeature(feature);
}
inline
void gnSequence::removeFeature(const uint32 featureI){
	spec->RemoveFeature(featureI);
}
inline
void gnSequence::getContainedFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const{
	spec->GetContainedFeatures(lt, feature_vector, index_vector);
}
inline
void gnSequence::getIntersectingFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const{
	spec->GetIntersectingFeatures(lt, feature_vector, index_vector);
}
inline
void gnSequence::getBrokenFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector) const{
	spec->GetBrokenFeatures(lt, feature_vector);
}

inline
boolean gnSequence::isCircular() const{
	return spec->IsCircular();
}

inline
void gnSequence::setCircular( const boolean value ){
	spec->SetCircular(value);
}
inline
void gnSequence::insert( const gnSeqI offset, const gnSequence& seq){
	insert(offset, *seq.spec);
}

// Size functions
inline
gnSeqI gnSequence::length() const{
	return spec->GetLength();
}
inline
gnSeqI gnSequence::size() const{
	return spec->GetLength();
}
inline
gnGenomeSpec* gnSequence::GetSpec() const{
	return spec;
}
inline
void gnSequence::SetSpec(gnGenomeSpec* s){
	if(spec != NULL)
		delete spec;
	spec = s;
}

inline
void gnSequence::setFilter(const gnBaseFilter* filt){
	filter_list.clear();
	if(filt != NULL)
		filter_list.push_back(filt);
}
inline
void gnSequence::setFilterList(std::list<const gnBaseFilter*>& filt_list){
	filter_list = filt_list;
}
inline
std::list<const gnBaseFilter*> gnSequence::getFilterList() const{
	return filter_list;
}



}	// end namespace genome

#endif
	// _gnSequence_h_
