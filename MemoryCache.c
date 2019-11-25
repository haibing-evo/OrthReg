#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <memory.h>
#include <math.h>

#include "MemoryCache.h"

using namespace std;

MemoryCache::MemoryCache()
{
}

MemoryCache::~MemoryCache()
{
	FreeResource();
}

void MemoryCache::FreeResource()
{
	long elements;
	long elementid;
	
	elements = data.size();
	for(elementid=0;elementid<elements;elementid++)
	{
		if(data[elementid]!=NULL) 
			free(data[elementid]);
		data[elementid] = NULL;
	}
}

long MemoryCache::Find(string query)
{
	long elements;
	long elementid;
	
	elements = name.size();
	for(elementid=0;elementid<elements;elementid++)
	{
		if(name[elementid] == query)
			return elementid;
	}
	return -1;
}

long MemoryCache::AddData(string pname,char *pdata,unsigned long long pdatasize)
{
	char *tmpdata;
	long pos;
	
	pos = Find(pname);
	
	if(pos>=0)
	{
		if(data[pos]!=NULL) free(data[pos]);
		data[pos] = pdata;
		datasize[pos]=pdatasize;
	}
	else
	{
		name.push_back(pname);
		data.push_back(pdata);
		datasize.push_back(pdatasize);
	}
}

long MemoryCache::LoadDataFromFile(string pname,string filename)
{
	FILE *pfile;
	long pos;
	unsigned long long localdatasize;
	char *tmpdata;
	
	pfile = fopen(filename.c_str(),"rb");
	if(pfile==NULL) return 0;
	
	fseek(pfile,0,SEEK_END);
	localdatasize = ftell(pfile);
	fseek(pfile,0,SEEK_SET);
	
	pos = Find(pname);
	if(pos>=0)
	{
		data[pos] = (char*)realloc(data[pos],localdatasize);
		fread(data[pos],sizeof(char),localdatasize,pfile);
		datasize[pos] = localdatasize;
	}
	else
	{
		tmpdata = (char*)malloc(localdatasize);
		fread(tmpdata,sizeof(char),localdatasize,pfile);
		name.push_back(pname);
		data.push_back(tmpdata);
		datasize.push_back(localdatasize);
	}
	fclose(pfile);
	return 1;
}

int MemoryCache::SaveVectorToFiles(string dir)
{
	FILE *pfile;
	long elementid,elementcount;
	string filename;
	elementcount = name.size();
	for(elementid=0; elementid<elementcount; elementid++)
	{
		filename = dir + name[elementid];
		pfile = fopen(filename.c_str(),"wb");
		fwrite(data[elementid],sizeof(char),datasize[elementid],pfile);
		fclose(pfile);		
	}
	return 1;
}


long long MemoryCache::ReadLineFromData(string pname,char **linestr,long long curpos)
{
	long long pos1,pos2;
	long long pos;
	char *pMem;
	long long pdatasize;
	
	pMem = GetData(pname);
	if(pMem==NULL)	return -1;

	pdatasize = GetDatasize(pname);

	if(curpos>=pdatasize)
		return -1;
	
	pos1 = curpos;
	*linestr = (char*)realloc(*linestr,1);
	(*linestr)[0]=0;
		
	if(pos1>=0)
		if(pMem[pos1]=='\n')
			return pos1;
		 
	pos2 = curpos+1;

	while(pMem[pos2]!='\n')
	{
		pos2++;
		if(pos2==pdatasize)
			break;
	}
	*linestr = (char*)realloc(*linestr,pos2-pos1+2);
	memcpy(*linestr,pMem+pos1,pos2-pos1+1);
	(*linestr)[pos2-pos1+1]=0;
	
	for(pos=pos2-pos1; pos>=0; pos--)
	{
		if((*linestr)[pos]=='\n')
		{
			(*linestr)[pos] = 0;
			if(pos>0)
			{
				if((*linestr)[pos-1]=='\r')
					(*linestr)[pos-1] = 0;
			}
			break;
		}
	}
	
	return pos2;
	
}

void MemoryCache::DeleteData(string pname)
{
	vector<string>::iterator it1,it2;
	vector<char *>::iterator it3;
	vector<unsigned long long>::iterator it5;
	
	it2 = name.end();
	it1 = name.begin();
	it3 = data.begin();
	it5 = datasize.begin();
	for(;it1!=it2;)
	{
		if(*it1 == pname)
		{
			it1 = name.erase(it1);
			free(*it3);
			it3 = data.erase(it3);
			it5 = datasize.erase(it5);
		}
		else
		{
			++it1;
			++it3;
			++it5;
		}
	}
}

char* MemoryCache::GetData(string pname)
{
	long pos;
	pos = Find(pname);
	if(pos>=0)
	{
		return(data[pos]);
	}
	else
	{
		return NULL;
	}
}

unsigned long long MemoryCache::GetDatasize(string pname)
{
	long pos;
	pos = Find(pname);
	if(pos>=0)
	{
		return(datasize[pos]);
	}
	else
	{
		return 0;
	}
}

char * MemoryCache::GetData(long item)
{
	if(data.size()>=item)
		return data[item];
	return NULL;
}

unsigned long long MemoryCache::GetDatasize(long item)
{
	if(datasize.size()>=item)
		return datasize[item];
	return 0;
}

