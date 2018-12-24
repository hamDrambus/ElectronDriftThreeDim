> Release/ThreeDimSimulation 42 7.0 100 "Output/v01.0/eData_7.0Td.root"
			<random seed> <Td> <N of electrons> <output file>
====================================================================================
WINDOWS:
====================================================================================
Setups:
1) Add 'root\include' to Project Preferences->Configuration Properties->Include Directories,
	'root\lib' and 'root\bin' to Configuration Properties->Library Directories,
	source directory for user code (if is in not default vs folder) to Configuration Properties->Source Directories

2) Add 'root\include' to c/c++->General->Additional Include Directories
	
3) Add _CRT_SECURE_NO_WARNINGS and __WIN32__ to c/c++->Preprocessor->Preprocessor Definitions

4) Add 'root\include\w32pragma.h' to c/c++->Advanced->Forced Include Files

5) Add 'root\lib' and 'root\bin' to Linker->General->Additional Library Directories

6) Add 'root_v5.34.34\lib\*.lib' to Linker->Input->Additional Dependencies

7) Debuging->Working Directory - set to root of project folder (which contatins 'include' and 'src')

Notes:
#include "Windows4Root.h" instead of #include "windows.h"

====================================================================================
LINUX (Eclipse):
====================================================================================