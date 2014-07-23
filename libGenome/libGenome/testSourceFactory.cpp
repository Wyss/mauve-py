#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSourceFactory.h"
#include "libGenome/gnFASSource.h"
#include "libGenome/gnDNXSource.h"
#include <iostream>

int32 main( int32 argc, char* argv[])
{
	argc; argv;
	
	gnSourceFactory* sf = gnSourceFactory::GetSourceFactory();

	// ########################################################
	// ## Path
	// ########################################################
	cout << " ***** TEST Path methods ***** " << endl;
	  // add
	string path = "c:/tmp";
	cout << "AddPath( " << path  << " )" << endl;
	if( !sf->AddPath( path ) )
		cout << "   bad Unable to add path: " << path << endl;
	else
		cout << "   good Path List Size: " << sf->GetPathListSize() << endl;
	  // insert
	path = "c:\\temp";
	cout << "InsPath( " << path  << ", 0 )"<< endl;
	if( !sf->InsPath( path, 0 ) )
		cout << "   bad Unable to Insert path at 0: " << path << endl;
	else
		cout << "   good Path List Size: " << sf->GetPathListSize() << endl;
	  // has
	path = "c:\\tmp";
	cout << "HasPath( " << path << " )" << endl;
	if( sf->HasPath( path ) )
		cout << "   good Does have " << path << endl;
	else
		cout << "   bad Does NOT have " << path << endl;
	path = "c:\\temp";
	cout << "HasPath( " << path << " )" << endl;
	if( sf->HasPath( path ) )
		cout << "   bad Does have " << path << endl;
	else
		cout << "   good Does NOT have " << path << endl;
	  // get
	cout << "GetPath(0) = " << endl;
	if( !(sf->GetPath(0).empty()) )
		cout << "   good" << endl;
	else
		cout << "   bad" << endl;
	cout << "GetPath(300)" <<  endl;
	if( sf->GetPath(300).empty() )
		cout << "   good" << endl;
	else
		cout << "   bad" << endl;
	  // del
	cout << "DelPath( 0 )" << endl;
	if( !sf->DelPath( 0 ) )
		cout << "   bad Unable to delete 0" << endl;
	else
		cout << "   good Path List Size: " << sf->GetPathListSize() << endl;
	cout << "DelPath( 0 )" << endl;
	if( !sf->DelPath( 0 ) )
		cout << "   bad Unable to delete 0" << endl;
	else
		cout << "   good Path List Size: " << sf->GetPathListSize() << endl;
	cout << "DelPath( 0 )" << endl;
	if( !sf->DelPath( 0 ) )
		cout << "   good Unable to delete 0" << endl;
	else
		cout << "   bad Path List Size: " << sf->GetPathListSize() << endl;
	

	// ########################################################
	// ## Source Class
	// ########################################################
	cout << " ***** TEST Source Class methods ***** " << endl;
	string extStr = ".fas";
	  // Set
	cout << "SetSourceClass( " << extStr << " , gnFASSource )" << endl;
	gnFASSource fasSource;
	if( !sf->SetSourceClass( extStr, fasSource ) )
		cout << "   bad Unable to add " << extStr << " , gnFASSource" << endl;
	else
		cout << "   good SourceClass List Size: " << sf->GetSourceClassListSize() << endl;
	  // Has
	cout << "HasSourceClass( " << extStr << " )" << endl;
	if( sf->HasSourceClass( extStr ) )
		cout << "   good Does have " << extStr << endl;
	else
		cout << "   bad Does NOT have " << extStr << endl;
	extStr = ".seq";
	cout << "HasSourceClass( " << extStr << " )" << endl;
	if( sf->HasSourceClass( extStr ) )
		cout << "   bad Does have " << extStr << endl;
	else
		cout << "   good Does NOT have " << extStr << endl;
	  // default
	cout << "SetDefaultSourceClass( new gnFASSource() )" << endl;
	if( !sf->SetDefaultSourceClass( fasSource  ) )
		cout << "   bad Cannot set Default Source Class" << endl;
	else
		cout << "   good DEFAULT source class set" << endl;
	  // Get
	extStr = ".fas";
	cout << "GetSourceClass( " << extStr << " )" << endl;
	if( sf->GetSourceClass( extStr ) != 0 )
		cout << "   good Does have " << extStr << endl;
	else
		cout << "   bad Does NOT have " << extStr << endl;
	extStr = ".seq";
	cout << "GetSourceClass( " << extStr << " )" << endl;
	if( sf->GetSourceClass( extStr ) != 0 )
		cout << "   good Does have " << extStr << endl;
	else
		cout << "   bad Does NOT have " << extStr << endl;
	  // Match
	string sourceStr = "c:/tests.fas" ;
	cout << "MatchSourceClass( " << sourceStr << " )" << endl;
	if( sf->MatchSourceClass( sourceStr ) != 0 )
		cout << "   good Does match" << endl;
	else
		cout << "   bad Does NOT match" << endl;
	sourceStr = "c:/tests.dnx";
	cout << "MatchSourceClass( " << sourceStr << " )" << endl;
	if( sf->MatchSourceClass( sourceStr ) != sf->GetDefaultSourceClass() )
		cout << "   bad Does match" << endl;
	else
		cout << "   good Does NOT match" << endl;
	  // Del
	extStr = ".fas";
	cout << "DelSourceClass( " << extStr << " )" << endl;
	if( !sf->DelSourceClass( extStr ) )
		cout << "   bad Cannot Delete " << extStr << endl;
	else
		cout << "   good SourceClass List Size: " << sf->GetSourceClassListSize() << endl;
	extStr = ".fas";
	cout << "DelSourceClass( " << extStr << " )" << endl;
	if( !sf->DelSourceClass( extStr ) )
		cout << "   good Cannot Delete " << extStr << endl;
	else
		cout << "   bad SourceClass List Size: " << sf->GetSourceClassListSize() << endl;

	// ########################################################
	// ## Source
	// ########################################################
	extStr = ".fas";
	  // Set
	cout << "SetSourceClass( " << extStr << " , gnFASSource )" << endl;
	if( !sf->SetSourceClass( extStr, fasSource ) )
		cout << "   bad Unable to add " << extStr << " , gnFASSource" << endl;
	else
		cout << "   good SourceClass List Size: " << sf->GetSourceClassListSize() << endl;

	extStr = ".dnx";
	  // Set
	gnDNXSource dnxSource;
	cout << "SetSourceClass( " << extStr << " , gnDNXSource )" << endl;
	if( !sf->SetSourceClass( extStr, dnxSource ) )
		cout << "   bad Unable to add " << extStr << " , gnDNXSource" << endl;
	else
		cout << "   good SourceClass List Size: " << sf->GetSourceClassListSize() << endl;


	cout << " ***** TEST Source Class methods ***** " << endl;
	  // add
	path = "W:\\Developm\\_Data\\_Testing_data\\fas\\";
	cout << "AddPath( " << path  << " )" << endl;
	if( !sf->AddPath( path ) )
		cout << "   bad Unable to add path: " << path << endl;
	else
		cout << "   good Path List Size: " << sf->GetPathListSize() << endl;

	  // Add
	sourceStr = "W:/Developm/_Data/_Testing_data/fas/test.fas";
	cout << "AddSource( " << sourceStr << ", true )" << endl;
	if( !sf->AddSource( sourceStr ) )
		cout << "   Source " << sourceStr << " NOT added" << endl;
	else
		cout << "   Source List Size: " << sf->GetSourceListSize() << endl;

	sourceStr = "test.abi";
	cout << "AddSource( " << sourceStr << ", true )" << endl;
	if( !sf->AddSource( sourceStr ) )
		cout << "   Source " << sourceStr << " NOT added" << endl;
	else
		cout << "   Source List Size: " << sf->GetSourceListSize() << endl;
	

	char bubba[50];
	
	cout << "Opening W:/Developm/_Data/_Testing_data/fas/V13anooz.fas...\n"; 
	gpSourceSeq *asplundh = new gpSourceSeq("W:/Developm/_Data/_Testing_data/fas/test.fas");
	cout << "Reading some bases.\n";
	string mystring = asplundh->ToString(130, 2480);
	cout << mystring << "\n";
	cout << "If that worked then congratulations!  You're one step closer than you realized.\n";
	cin >> bubba;
	delete asplundh;

// test out DNX
	cout << "Opening W:/Developm/_Data/_Testing_data/fas/deenex.dnx...\n"; 
	asplundh = new gpSourceSeq("W:/Developm/_Data/_Testing_data/fas/deenex.dnx");
	cout << "Reading some bases.\n";
	mystring = asplundh->ToString();
	cout << mystring << "\n";
	cout << "If that worked then congratulations!  You're one step closer than you realized.\n";
	cin >> bubba;
}