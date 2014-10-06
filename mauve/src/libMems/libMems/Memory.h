#ifndef __libMems_Memory_h__
#define __libMems_Memory_h__


void printMemUsage();
static bool debugging_memory = false;
#include <iostream>

#ifdef WIN32
#include <windows.h>
#include <PSAPI.h>
inline
void printMemUsage()
{
//	if(!debugging_memory)
//		return;

	DWORD proclist[500];
	DWORD cbNeeded;
	BOOL rval;
	rval = EnumProcesses( proclist, sizeof(proclist), &cbNeeded );
	int p_count = cbNeeded / sizeof(DWORD);
	HANDLE phand;
	HMODULE hMod;
	char szFileName[MAX_PATH];
	for( int p = 0; p < p_count; p++ )
	{
		phand = OpenProcess( PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, 0, proclist[p] );
		DWORD dwSize2;
		if (EnumProcessModules(phand, &hMod, sizeof(hMod), &dwSize2)) 
		{

			// Get the module name
			if ( !GetModuleBaseName(phand, hMod, szFileName, sizeof(szFileName)) )
				szFileName[0] = 0;
			if( strncmp( szFileName, "progressiveMauve", 16 ) == 0 )
				break;	// found the right module
		}
		CloseHandle(phand);
	}

	PROCESS_MEMORY_COUNTERS mem_info;

	if( GetProcessMemoryInfo( phand, &mem_info, sizeof(mem_info) ) )
	{
			std::cout << "Working set size: " << mem_info.WorkingSetSize / (1024 * 1024) << " Mb\n";
//		cout << "Paged pool usage: " << mem_info.QuotaPagedPoolUsage << endl;
//		cout << "Non-Paged pool usage: " << mem_info.QuotaNonPagedPoolUsage << endl;
			std::cout << "Pagefile usage: " << mem_info.PagefileUsage / (1024 * 1024) << " Mb\n";
			std::cout.flush();
	}
}
#else
inline
void printMemUsage()
{};
#endif

#endif	//__libMems_Memory_h__

