#ifndef __OrthSiteInfo_H
#define __OrthSiteInfo_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include "MemoryCache.h"

using namespace std;

#ifdef  __cplusplus
extern "C" {
#endif

typedef struct
{
	char chr[6];
	long orthsite;
	char nt;
	char strand;
	float encode;
} orthinfo;

class OrthSiteInfo
{
public:
string workingdir;
string queryspecies;
string targetspecies;
string chainfilename;
string qfilename;
string tfilename;
string qchrfile;
string tchrfile;

MemoryCache mc;
unsigned long long chainpos;
char *linestr;

long qchrsize;
long tchrsize;
long qstart;
long qend;
long tstart;
long tend;
string qchr;
string tchr;
string qstrand;
string tstrand;
long score;
char savestrand;
bool skip;

long alignsize;
long conflictsite;
long tspansize;
long qspansize;
long minchainsize;

vector <string> chr_query_target_name;
vector <orthinfo *> chr_query_target_data;
vector <long> chr_genome_size;

string running_chr_query_target_name;
orthinfo * running_chr_query_target_data;
char * running_chr_query_seq;

string running_chr_target_query_name;
orthinfo * running_chr_target_query_data;
char * running_chr_target_seq;

orthinfo * tmp_chr_data;

OrthSiteInfo();
virtual ~OrthSiteInfo();

void SetChainFile(char *chainfile);
int ReadChainHeader();
int LoadChainDetails();
int SeekToNextChain();
int LoadChainData();
int GetConflicts();
int SetMinChainSize(long pminsize);

void SetQuerySpecies(char *qs);
void SetTargetSpecies(char *ts);
void SetRelatedFileNames();

long LoadOrthSiteFromFile(string filename);
int SaveOrthSiteVectorToFile();

int SaveOrthSiteToVector();
long FindPositionInOrthSiteVector(string search_query_target);

int GenerateQueryOrthSite();
int GenerateTargetOrthSite();

void LoadQueryFromVector(long index);
void LoadTargetFromVector(long index);

int GetFastaFileLineLength(string filename);
void LoadQuerySeq();
void LoadTargetSeq();

void FreeRunningMemory();
void FreeMemory();

};

#ifdef  __cplusplus
}
#endif

#endif