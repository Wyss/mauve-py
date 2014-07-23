;NSIS Script For libGenome

;Backgound Colors
BGGradient 8080FF 0000A8 FFFFFF
BrandingText " "

;Title Of Your Application
Name "libGenome 0.5.2"

;Do A CRC Check
CRCCheck On

;Output File Name
OutFile "libGenome_installer.exe"

;License Page Introduction
LicenseText "You must agree to this license before installing."

;License Data
LicenseData "C:\Development\libGenome\COPYING"

ComponentText "Select libGenome components to install" "All are recommended"

CompletedText "Completed successfully.  You may need to log out or restart the computer."

;The Default Installation Directory
InstallDir "$DESKTOP\libGenome"

;The text to prompt the user to enter a directory
DirText "Please select a libGenome development folder below"

AllowRootDirInstall true

SubSection "libGenome Development"

Section "libGenomeBase"
  ;Install Files
  SetOutPath $INSTDIR
  SetCompress Auto
  SetOverwrite IfNewer

  ; Top level directory files
  File "C:\Development\libGenome\AUTHORS"
  File "C:\Development\libGenome\ChangeLog"
  File "C:\Development\libGenome\configure.in"
  File "C:\Development\libGenome\COPYING"
  File "C:\Development\libGenome\INSTALL"
  File "C:\Development\libGenome\Makefile.am"
  File "C:\Development\libGenome\NEWS"
  File "C:\Development\libGenome\README"
SectionEnd

Section "Header Files"
  ;Install Files
  SetOutPath $INSTDIR
  SetCompress Auto
  SetOverwrite IfNewer
  ; /include files
SetOutPath "$INSTDIR\include"
File "C:\Development\libGenome\libGenome\gn_cw.h"
File "C:\Development\libGenome\libGenome\gn_cw.pch"
File "C:\Development\libGenome\libGenome\gn_cw.pch++"
File "C:\Development\libGenome\libGenome\gn_cw_d.h"
File "C:\Development\libGenome\libGenome\gn_cw_d.pch"
File "C:\Development\libGenome\libGenome\gn_cw_d.pch++"
File "C:\Development\libGenome\libGenome\gn_cw_dll.h"
File "C:\Development\libGenome\libGenome\gn_cw_dll.pch"
File "C:\Development\libGenome\libGenome\gn_cw_dll.pch++"
File "C:\Development\libGenome\libGenome\gn_cw_dlld.h"
File "C:\Development\libGenome\libGenome\gn_cw_dlld.pch"
File "C:\Development\libGenome\libGenome\gn_cw_dlld.pch++"
File "C:\Development\libGenome\libGenome\Makefile.am"

  ; /include/gn files
SetOutPath "$INSTDIR\include"
File "C:\Development\libGenome\libGenome\gn\gnABISource.h"
File "C:\Development\libGenome\libGenome\gn\gnBaseFeature.h"
File "C:\Development\libGenome\libGenome\gn\gnBaseFilter.h"
File "C:\Development\libGenome\libGenome\gn\gnBaseHeader.h"
File "C:\Development\libGenome\libGenome\gn\gnBaseQualifier.h"
File "C:\Development\libGenome\libGenome\gn\gnBaseSource.h"
File "C:\Development\libGenome\libGenome\gn\gnBaseSpec.h"
File "C:\Development\libGenome\libGenome\gn\gnClone.h"
File "C:\Development\libGenome\libGenome\gn\gnCompare.h"
File "C:\Development\libGenome\libGenome\gn\gnContigSpec.h"
File "C:\Development\libGenome\libGenome\gn\gnDataBaseSource.h"
File "C:\Development\libGenome\libGenome\gn\gnDebug.h"
File "C:\Development\libGenome\libGenome\gn\gnDefs.h"
File "C:\Development\libGenome\libGenome\gn\gnDNASequence.h"
File "C:\Development\libGenome\libGenome\gn\gnDNXSource.h"
File "C:\Development\libGenome\libGenome\gn\gnException.h"
File "C:\Development\libGenome\libGenome\gn\gnExceptionCode.h"
File "C:\Development\libGenome\libGenome\gn\gnFASSource.h"
File "C:\Development\libGenome\libGenome\gn\gnFastTranslator.h"
File "C:\Development\libGenome\libGenome\gn\gnFeature.h"
File "C:\Development\libGenome\libGenome\gn\gnFileContig.h"
File "C:\Development\libGenome\libGenome\gn\gnFileSource.h"
File "C:\Development\libGenome\libGenome\gn\gnFilter.h"
File "C:\Development\libGenome\libGenome\gn\gnFragmentSpec.h"
File "C:\Development\libGenome\libGenome\gn\gnGBKSource.h"
File "C:\Development\libGenome\libGenome\gn\gnGenomeSpec.h"
File "C:\Development\libGenome\libGenome\gn\gnLocation.h"
File "C:\Development\libGenome\libGenome\gn\gnMultiSpec.h"
File "C:\Development\libGenome\libGenome\gn\gnPosSpecificTranslator.h"
File "C:\Development\libGenome\libGenome\gn\GNPREC.H"
File "C:\Development\libGenome\libGenome\gn\gnProteinSequence.h"
File "C:\Development\libGenome\libGenome\gn\gnRAWSource.h"
File "C:\Development\libGenome\libGenome\gn\gnRNASequence.h"
File "C:\Development\libGenome\libGenome\gn\gnSEQSource.h"
File "C:\Development\libGenome\libGenome\gn\gnSequence.h"
File "C:\Development\libGenome\libGenome\gn\gnSetup.h"
File "C:\Development\libGenome\libGenome\gn\gnSourceFactory.h"
File "C:\Development\libGenome\libGenome\gn\gnSourceHeader.h"
File "C:\Development\libGenome\libGenome\gn\gnSourceQualifier.h"
File "C:\Development\libGenome\libGenome\gn\gnSourceSpec.h"
File "C:\Development\libGenome\libGenome\gn\gnStringHeader.h"
File "C:\Development\libGenome\libGenome\gn\gnStringQualifier.h"
File "C:\Development\libGenome\libGenome\gn\gnStringSpec.h"
File "C:\Development\libGenome\libGenome\gn\gnStringTools.h"
File "C:\Development\libGenome\libGenome\gn\gnTranslator.h"
File "C:\Development\libGenome\libGenome\gn\gnVersion.h"
File "C:\Development\libGenome\libGenome\gn\Makefile.am"
SectionEnd

Section "Precompiled Libraries"
; lib/ files
SetOutPath "$INSTDIR\lib"
File "C:\Development\libGenome\lib\genome.dll"
File "C:\Development\libGenome\lib\genome.lib"
File "C:\Development\libGenome\lib\genomed.dll"
File "C:\Development\libGenome\lib\genomed.lib"
File "C:\Development\libGenome\lib\genomeddll.lib"
File "C:\Development\libGenome\lib\genomedll.lib"
SectionEnd

Section "libGenome Project Files"
; projects/ files
SetOutPath "$INSTDIR\projects"
File "C:\Development\libGenome\projects\libGenome.mcp"
File "C:\Development\libGenome\projects\libGenome.nsi"
SectionEnd

Section "Source Files"
; src/ files
SetOutPath "$INSTDIR\src"
File "C:\Development\libGenome\libGenome\gnABISource.cpp"
File "C:\Development\libGenome\libGenome\gnBaseFeature.cpp"
File "C:\Development\libGenome\libGenome\gnCompare.cpp"
File "C:\Development\libGenome\libGenome\gnContigSpec.cpp"
File "C:\Development\libGenome\libGenome\gnDNXSource.cpp"
File "C:\Development\libGenome\libGenome\gnException.cpp"
File "C:\Development\libGenome\libGenome\gnExceptionCode.cpp"
File "C:\Development\libGenome\libGenome\gnFASSource.cpp"
File "C:\Development\libGenome\libGenome\gnFastTranslator.cpp"
File "C:\Development\libGenome\libGenome\gnFeature.cpp"
File "C:\Development\libGenome\libGenome\gnFileContig.cpp"
File "C:\Development\libGenome\libGenome\gnFileSource.cpp"
File "C:\Development\libGenome\libGenome\gnFilter.cpp"
File "C:\Development\libGenome\libGenome\gnFragmentSpec.cpp"
File "C:\Development\libGenome\libGenome\gnGBKSource.cpp"
File "C:\Development\libGenome\libGenome\gnGenomeSpec.cpp"
File "C:\Development\libGenome\libGenome\gnLocation.cpp"
File "C:\Development\libGenome\libGenome\gnMultiSpec.cpp"
File "C:\Development\libGenome\libGenome\gnPosSpecificTranslator.cpp"
File "C:\Development\libGenome\libGenome\gnRAWSource.cpp"
File "C:\Development\libGenome\libGenome\gnSEQSource.cpp"
File "C:\Development\libGenome\libGenome\gnSeqStringTest.cpp"
File "C:\Development\libGenome\libGenome\gnSequence.cpp"
File "C:\Development\libGenome\libGenome\gnSourceFactory.cpp"
File "C:\Development\libGenome\libGenome\gnSourceHeader.cpp"
File "C:\Development\libGenome\libGenome\gnSourceQualifier.cpp"
File "C:\Development\libGenome\libGenome\gnSourceSpec.cpp"
File "C:\Development\libGenome\libGenome\gnStringHeader.cpp"
File "C:\Development\libGenome\libGenome\gnStringQualifier.cpp"
File "C:\Development\libGenome\libGenome\gnStringSpec.cpp"
File "C:\Development\libGenome\libGenome\gnStringTools.cpp"
File "C:\Development\libGenome\libGenome\gnTranslator.cpp"
File "C:\Development\libGenome\libGenome\Makefile.am"
SectionEnd

Section "libGenome Stationary"
; stationary/ files
SetOutPath "$INSTDIR\stationary\libGenome_App"
File "C:\Development\libGenome\stationary\libGenome_App\libGenome_App.mcp"
SetOutPath "$INSTDIR\stationary\libGenome_App\src"
File "C:\Development\libGenome\stationary\libGenome_App\src\genomeApp.cpp"

push $0
push $1
push $2
push $3

push 0
pop $0
EnumRegKey $1 HKLM "SOFTWARE\Metrowerks\CodeWarrior\Product Versions" $0
ReadRegStr $2 HKLM "SOFTWARE\Metrowerks\CodeWarrior\Product Versions\$1" VERSION
ReadRegStr $3 HKLM "SOFTWARE\Metrowerks\CodeWarrior\Product Versions\$1" PATH

DetailPrint "Installing libGenome stationary for $1 version $2 in directory $3"
SetOutPath "$3\Stationery\libGenome_App"
File "C:\Development\libGenome\stationary\libGenome_App\libGenome_App.mcp"
SetOutPath "$3\Stationery\libGenome_App\src"
File "C:\Development\libGenome\stationary\libGenome_App\src\genomeApp.cpp"
CreateDirectory "$3\Stationery\libGenome_App\include"
CreateDirectory "$3\Stationery\libGenome_App\doc"
CreateDirectory "$3\Stationery\libGenome_App\bin"

pop $3
pop $2
pop $1
pop $0

; Add a registry key for the GenomeDev source tree
push $0
push $1
StrLen $0 $INSTDIR
IntOp $0 $0 - 10
StrCpy $1 $INSTDIR $0
WriteRegStr HKLM "SOFTWARE\libGenome\" "GenomeDev" "$1"
pop $1
pop $0

SectionEnd

SubSectionEnd

SubSection "Documentation"
Section "HTML documentation"
; copy documentation directories
SetOutPath "$INSTDIR\doc"
File "C:\Development\libGenome\doc\Doxyfile.txt"
File "C:\Development\libGenome\doc\Makefile.am"
File /r "C:\Development\libGenome\doc\html"
SectionEnd

Section "UNIX man page documentation"
SetOutPath "$INSTDIR\doc"
File "C:\Development\libGenome\doc\Doxyfile.txt"
File "C:\Development\libGenome\doc\Makefile.am"
File /r "C:\Development\libGenome\doc\man"
SectionEnd

Section "Documentation Start Menu Shortcut"
  ;Add Shortcuts
  CreateShortCut "$SMPROGRAMS\libGenome Documentation.lnk" "$INSTDIR\doc\html\index.html"
SectionEnd

SubSectionEnd

Section "Uninstaller"
  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\libGenome" "DisplayName" "libGenome (remove only)"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\libGenome" "UninstallString" "$INSTDIR\Uninst.exe"
WriteUninstaller "Uninst.exe"
SectionEnd


UninstallText "This will uninstall libGenome from your system"

Section Uninstall

  ;Delete Files
; Top level directory files
Delete /REBOOTOK "$INSTDIR\AUTHORS"
Delete /REBOOTOK "$INSTDIR\ChangeLog"
Delete /REBOOTOK "$INSTDIR\configure.in"
Delete /REBOOTOK "$INSTDIR\COPYING"
Delete /REBOOTOK "$INSTDIR\INSTALL"
Delete /REBOOTOK "$INSTDIR\Makefile.am"
Delete /REBOOTOK "$INSTDIR\NEWS"
Delete /REBOOTOK "$INSTDIR\README"

; /include/gn files
Delete /REBOOTOK "$INSTDIR\libGenome\gnABISource.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnBaseFeature.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnBaseFilter.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnBaseHeader.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnBaseQualifier.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnBaseSource.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnBaseSpec.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnClone.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnCompare.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnContigSpec.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnDataBaseSource.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnDebug.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnDefs.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnDNASequence.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnDNXSource.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnException.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnExceptionCode.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFASSource.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFastTranslator.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFeature.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFileContig.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFileSource.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFilter.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFragmentSpec.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnGBKSource.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnGenomeSpec.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnLocation.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnMultiSpec.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnPosSpecificTranslator.h"
Delete /REBOOTOK "$INSTDIR\libGenome\GNPREC.H"
Delete /REBOOTOK "$INSTDIR\libGenome\gnProteinSequence.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnRAWSource.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnRNASequence.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSEQSource.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSequence.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSetup.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSourceFactory.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSourceHeader.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSourceQualifier.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSourceSpec.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnStringHeader.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnStringQualifier.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnStringSpec.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnStringTools.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnTranslator.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gnVersion.h"
Delete /REBOOTOK "$INSTDIR\libGenome\Makefile.am"
RMDir "$INSTDIR\include"

  ; /include files
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw.pch"
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw.pch++"
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw_d.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw_d.pch"
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw_d.pch++"
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw_dll.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw_dll.pch"
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw_dll.pch++"
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw_dlld.h"
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw_dlld.pch"
Delete /REBOOTOK "$INSTDIR\libGenome\gn_cw_dlld.pch++"
Delete /REBOOTOK "$INSTDIR\libGenome\Makefile.am"
RMDir "$INSTDIR\include"


; lib/ files
Delete /REBOOTOK "$INSTDIR\lib\genome.dll"
Delete /REBOOTOK "$INSTDIR\lib\genome.lib"
Delete /REBOOTOK "$INSTDIR\lib\genomed.dll"
Delete /REBOOTOK "$INSTDIR\lib\genomed.lib"
Delete /REBOOTOK "$INSTDIR\lib\genomeddll.lib"
Delete /REBOOTOK "$INSTDIR\lib\genomedll.lib"
RMDir "$INSTDIR\lib"

; projects/ files
Delete /REBOOTOK "$INSTDIR\projects\libGenome.mcp"
Delete /REBOOTOK "$INSTDIR\projects\libGenome.nsi"
RMDir "$INSTDIR\projects"

; src/ files
Delete /REBOOTOK "$INSTDIR\libGenome\gnABISource.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnBaseFeature.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnCompare.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnContigSpec.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnDNXSource.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnException.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnExceptionCode.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFASSource.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFastTranslator.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFeature.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFileContig.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFileSource.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFilter.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnFragmentSpec.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnGBKSource.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnGenomeSpec.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnLocation.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnMultiSpec.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnPosSpecificTranslator.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnRAWSource.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSEQSource.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSeqStringTest.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSequence.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSourceFactory.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSourceHeader.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSourceQualifier.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnSourceSpec.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnStringHeader.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnStringQualifier.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnStringSpec.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnStringTools.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\gnTranslator.cpp"
Delete /REBOOTOK "$INSTDIR\libGenome\Makefile.am"
RMDir "$INSTDIR\src"

; stationary/ files
Delete /REBOOTOK "$INSTDIR\stationary\libGenome_App\src\genomeApp.cpp"
RMDir "$INSTDIR\stationary\libGenome_App\src"
Delete /REBOOTOK "$INSTDIR\stationary\libGenome_App\libGenome_App.mcp"
RMDir "$INSTDIR\stationary\libGenome_App"
RMDir "$INSTDIR\stationary\include"
RMDir "$INSTDIR\stationary\doc"
RMDir "$INSTDIR\stationary\bin"
RMDir "$INSTDIR\stationary\"

; copy documentation directories
RMDir /r "$INSTDIR\doc\html"
RMDir /r "$INSTDIR\doc\man"
RMDir /r "$INSTDIR\doc"

; Delete documentation shortcut
Delete "$SMPROGRAMS\libGenome Documentation.lnk"

; Delete the stationery
push $0
push $1
push $2
push $3

push 0
pop $0
EnumRegKey $1 HKLM "SOFTWARE\Metrowerks\CodeWarrior\Product Versions" $0
ReadRegStr $2 HKLM "SOFTWARE\Metrowerks\CodeWarrior\Product Versions\$1" VERSION
ReadRegStr $3 HKLM "SOFTWARE\Metrowerks\CodeWarrior\Product Versions\$1" PATH

DetailPrint "Removing libGenome stationary for $1 version $2 in directory $3"
Delete /REBOOTOK "$3\Stationery\libGenome_App\libGenome_App.mcp"
Delete /REBOOTOK "$3\Stationery\libGenome_App\src\genomeApp.cpp"
RMDir "$3\Stationery\libGenome_App\src"
RMDir "$3\Stationery\libGenome_App\include"
RMDir "$3\Stationery\libGenome_App\doc"
RMDir "$3\Stationery\libGenome_App\bin"
RMDir "$3\Stationery\libGenome_App"

pop $3
pop $2
pop $1
pop $0


  ;Delete Uninstaller And Unistall Registry Entries
  Delete "$INSTDIR\Uninst.exe"
  DeleteRegKey HKEY_LOCAL_MACHINE "SOFTWARE\libGenome"
  DeleteRegKey HKEY_LOCAL_MACHINE "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\libGenome"
  RMDir "$INSTDIR"
SectionEnd


