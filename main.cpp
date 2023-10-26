#include <iostream>
#include <fstream>
#include <unistd.h>
#include "HalfEdgeMesh.h"
#include "AreaExpansionRate.h"

void FrontMatter(int argc, char **argv){
	if(argc != 2){
		std::cout << "Error: the number of values in console command\n";
		std::cout << "Please put the directory name of workspace.\n"; 
		std::cout << std::endl;
		exit(0);
	}
	char* dirName = argv[1];
	std::cout << "Workspace name : " << dirName << std::endl;
	
	if(chdir(dirName) != 0){ // ディレクトリの移動
		std::cout << "Error: moving the directory of workspace" << std::endl;
		exit(1);
	}
	char dirPath[256];
	if(getcwd(dirPath,256) == NULL){ // (移動後の)ディレクトリを取得
		std::cout << "Error: getting workspace" << std::endl;
		exit(1);
	}
	std::cout << "Workspace path : " << dirPath << "\n" << std::endl;
}

int main(int argc, char** argv){
	FrontMatter(argc, argv);

	AreaExpansionRate AER;
	if(AER.PreSurfaceMode == 0){
		AER.HarmonicMapPost();
		AER.CalcAreaExpRatePerFacePostFromPlane();
	}else{
		AER.HarmonicMaps();
		AER.SetCorrespPostOnPre();
		AER.CalcAreaExpRatePerFace();
	}
	AER.WriteAreaExpRate();

	return 0;
}
