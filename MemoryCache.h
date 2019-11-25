#ifndef __MemoryCache_H
#define __MemoryCache_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;

#ifdef  __cplusplus
extern "C" {
#endif

class MemoryCache
{
public:
MemoryCache();
~MemoryCache();

vector <string> name;
vector <char *> data;
vector <unsigned long long> datasize;

void FreeResource();

long Find(string query);
long LoadDataFromFile(string pname,string filename);
int SaveVectorToFiles(string dir);
long long ReadLineFromData(string pname,char **linedata,long long curpos);
long AddData(string pname,char *pdata,unsigned long long pdatasize);
void DeleteData(string pname);

char * GetData(long item);
unsigned long long GetDatasize(long item);
char * GetData(string pname);
unsigned long long GetDatasize(string pname);
};

#ifdef  __cplusplus
}
#endif

#endif