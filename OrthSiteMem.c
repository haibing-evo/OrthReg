#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "OrthSiteInfo.h"

int main(int argc, char *argv[])
{
	OrthSiteInfo orth;
	long vectorpos1,vectorpos2;
	
	if(argc!=5)
	{
		printf("Usage:\n");
		printf("%s chainfile tspeices qspecies min_chainsize\n",argv[0]);
		return 0;
	}
	
	orth.SetChainFile(argv[1]);
	orth.SetTargetSpecies(argv[2]);
	orth.SetQuerySpecies(argv[3]);
	orth.SetMinChainSize(atol(argv[4]));
	
	orth.LoadChainData();
	
	while(orth.ReadChainHeader())
	{
		if(orth.skip)
		{
			cout << "Skipped a chain ... " << endl;
			orth.SeekToNextChain();
			continue;
		}
			
		//cout << "----------new chain------------" << endl;
		orth.SetRelatedFileNames();

		vectorpos1 = orth.FindPositionInOrthSiteVector(orth.running_chr_query_target_name);
		if(vectorpos1>=0)
		{
			//cout << "Load Query from OrthSite Vector ... " << endl;
			orth.LoadQueryFromVector(vectorpos1);
		}
		else
		{
			//cout << "Generate a new Query OrthSite ... " << endl;
			orth.GenerateQueryOrthSite();
			orth.LoadQuerySeq();
		}
		
		vectorpos2 = orth.FindPositionInOrthSiteVector(orth.running_chr_target_query_name);
		if(vectorpos2>=0)
		{
			//cout << "Load Target from OrthSite Vector ... " << endl;
			orth.LoadTargetFromVector(vectorpos2);
		}
		else
		{
			//cout << "Generate a new Target OrthSite ... " << endl;
			orth.GenerateTargetOrthSite();
			orth.LoadTargetSeq();
		}
		
		if(orth.GetConflicts())
		{
			cout << "Skipped a chain by conflicts ... " << endl;
			orth.SeekToNextChain();
			continue;
		}
		
		//cout << "Read Chain Details ... " << endl;
		orth.LoadChainDetails();
		
		//cout << "Save to OrthSite Vector ... " << endl;
		orth.SaveOrthSiteToVector();
	}
	//cout << "Save OrthSite Vector to Disk Files ... " << endl;
	orth.SaveOrthSiteVectorToFile();
}