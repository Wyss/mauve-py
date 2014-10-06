#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libGenome/gnGBKSource.h"
#include "libGenome/gnFASSource.h"
#include <algorithm>

static boolean LocationLessthan(const gnLocation& a, const gnLocation& b);

static boolean LocationLessthan(const gnLocation& a, const gnLocation& b){
	return a.GetStart() < b.GetStart();
}

static boolean LocationEndLessthan(const gnLocation& a, const gnLocation& b);

static boolean LocationEndLessthan(const gnLocation& a, const gnLocation& b){
	return a.GetEnd() < b.GetStart();
}

static boolean LocationSizeLessthan(const gnLocation& a, const gnLocation& b);

static boolean LocationSizeLessthan(const gnLocation& a, const gnLocation& b){
	return a.GetEnd() - a.GetStart() < b.GetEnd() - b.GetStart();
}

// create a location list which contains every intervening region
// in l_list
void ComplementLocation(vector<gnLocation>& l_list){
	if(l_list.size() < 2)
		return;	//need 2 or more locations to get the in between regions
	
	vector<gnLocation> comp_list;
	
	sort(l_list.begin(), l_list.end(), &LocationLessthan);
	
	for(uint32 locationI = 1; locationI < l_list.size(); locationI++){

		//if they don't intersect there is an intervening region
		if( !l_list[locationI].Intersects( l_list[locationI-1] ) ){

			//create the region between locationI and locationI - 1
			gnLocation between_loc;
			between_loc.SetStart(l_list[locationI-1].GetEnd());
			between_loc.SetEnd(l_list[locationI].GetStart());
			//add it to the complement list
			comp_list.push_back(between_loc);
		}
		
	}

	l_list = comp_list;
}

void DeleteRepeats(list<gnLocation>& loc_list){

	loc_list.sort(&LocationLessthan);

	list<gnLocation>::iterator loc_iter = loc_list.begin();
	list<gnLocation>::iterator cur_loc;
	for(; loc_iter != loc_list.end(); loc_iter++){
		cur_loc = loc_iter;
		cur_loc++;
		while(cur_loc != loc_list.end() && !LocationEndLessthan(*loc_iter, *cur_loc)){
			if(loc_iter->Intersects(*cur_loc)){
				*loc_iter = loc_iter->GetUnion(*cur_loc);
				list<gnLocation>::iterator to_del = cur_loc;
				cur_loc--;
				loc_list.erase(to_del);
			}
			cur_loc++;
		}
	}
}

void WriteData(gnSequence& file_sequence, list<gnLocation>& loc_list, string& filename){
	gnSequence loc_seq;
	cout << "Gathering sequence data to write...\n";
	list<gnLocation>::iterator loc_iter = loc_list.begin();
	uint32 loc_count = loc_list.size();
	uint32 notice_interval = loc_count / 10;
	uint32 cur_locI = 0;
	for(;loc_iter != loc_list.end(); loc_iter++){
		loc_seq += file_sequence.ToString( loc_iter->GetEnd() - loc_iter->GetStart(), loc_iter->GetStart());
		if(cur_locI % notice_interval == 0){
			cout << (cur_locI / loc_count) * 10 << "%..";
		}
		cur_locI++;
	}
	cout << "\nWriting sequence data...\n";
	gnFASSource::Write(loc_seq, filename);
}

void WriteList(list<gnLocation>& loc_list, string& filename){
	ofstream list_out(filename.c_str());
	list<gnLocation>::iterator loc_iter = loc_list.begin();
	for(;loc_iter != loc_list.end(); loc_iter++){
		list_out << loc_iter->GetStart() << '\t' << loc_iter->GetEnd() << '\n';
	}
	list_out.close();
}

void print_usage(char* pname){
	cout << "Usage: " << pname << " <genbank file> <exon_list> <intron_list> <exon_seq> <intron_seq>\n";
}

int main( int argc, char* argv[] )
{
	boolean run_interactive = false;
	string seq_filename;
	string exon_list_filename;
	string intron_list_filename;
	string exon_seq_filename;
	string intron_seq_filename;
	// define a gnSequence to store the sequence
	gnSequence file_sequence;

	if(!run_interactive){
		// check for correct calling semantics
		if(argc != 6){
			print_usage(argv[0]);
			return -1;
		}

		seq_filename = argv[1];
		exon_list_filename = argv[2];
		intron_list_filename = argv[3];
		exon_seq_filename = argv[4];
		intron_seq_filename = argv[5];
	}else{

		// Get the name of the sequence to load
		cout << "Enter the name of the sequence file to load:\n" ;
		cin >> seq_filename;
	}

	// Load the sequence and tell the user if it loaded successfully
	if(file_sequence.LoadSource(seq_filename))
	{
		cout << seq_filename << " loaded successfully.  " << file_sequence.length() << " base pairs.\n";
	}else{
		cout << "Error loading file.\n";
		return -1;
	}
	
	uint32 feature_count = file_sequence.getFeatureListLength();
	uint32 cds_count = 0;

	//define a sequence for each type of data
	gnSequence exon_seq;
	gnSequence intron_seq;
	
	list<gnLocation> exon_list;
	list<gnLocation> intron_list;

	cout << "There are " << feature_count << " known features.\n";
	for(uint32 featureI = 0; featureI < feature_count; featureI++){
		gnBaseFeature* cur_feat = file_sequence.getFeature(featureI);
		if(cur_feat == NULL){
			cout << "cur_feat is NULL, trace me!\n";
			cur_feat = file_sequence.getFeature(featureI);
			continue;
		}
		if(cur_feat->GetName() == "CDS"){
			cds_count++;
			uint32 loc_count = cur_feat->GetLocationListLength();
			vector<gnLocation> cur_exon_list, cur_intron_list;

			for(uint32 locationI = 0; locationI < loc_count; locationI++){
				cur_exon_list.push_back(cur_feat->GetLocation(locationI));
				exon_list.push_back(cur_exon_list[locationI]);
			}
			if(loc_count >= 2){
				cur_intron_list = cur_exon_list;
				ComplementLocation(cur_intron_list);
				for(uint32 locationI = 0; locationI < cur_intron_list.size(); locationI++)
					intron_list.push_back(cur_intron_list[locationI]);
			}
/*  Code for exon size tracking  */
/*			sort(cur_exon_list.begin(), cur_exon_list.end(), &LocationSizeLessthan);
			gnSeqI cur_size = cur_exon_list[loc_count - 1].GetEnd() - cur_exon_list[loc_count - 1].GetStart();
			size_out << cur_size << '\n';
			if(cur_size > 5000){
				cout << "psycho bob is at it again\n";
				cout << "Found " << cur_size << " length exon in feature: " << cur_feat->GetName() << "\n";
				cout << "Feature Index: " << featureI << "\n";
				for(uint32 i=0; i < cur_feat->GetQualifierListLength(); i++){
					cout << cur_feat->GetQualifierName(i) << "\t" << cur_feat->GetQualifierValue(i) << "\n";
				}
				cout << "\n";
			}		*/
			delete cur_feat;
		}
	}

	cout << "There are " << cds_count << " cds features.\n";
	
	//now sort the lists and delete repeats and fix overlaps
	uint32 exon_count = exon_list.size();
	cout << "Deleting exon repeats...  ";
	DeleteRepeats(exon_list);
	cout << "Eliminated " << exon_count - exon_list.size() << " repeats for a total of ";
	cout << exon_list.size() << " repeats\n";
	uint32 intron_count = intron_list.size();
	cout << "Deleting intron repeats...  ";
	DeleteRepeats(intron_list);
	cout << "Eliminated " << intron_count - intron_list.size() << " repeats for a total of ";
	cout << intron_list.size() << " repeats\n";
	WriteList(exon_list, exon_list_filename);
	WriteList(intron_list, intron_list_filename);
	WriteData(file_sequence, exon_list, exon_seq_filename);
	WriteData(file_sequence, intron_list, intron_seq_filename);
	if(run_interactive)
		cin >> seq_filename;
	return 0;
};
