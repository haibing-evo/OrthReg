#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory.h>
#include <unistd.h>
#include <math.h>

#include "OrthSiteInfo.h"
#include "xiesequence.h"
#include "MemoryCache.h"

using namespace std;
#define LINE_SIZE 10000

OrthSiteInfo::OrthSiteInfo()
{
	char dir[1024];
	getcwd(dir,1024);
	workingdir = dir;
	chainpos = 0;
	linestr = NULL;
}

OrthSiteInfo::~OrthSiteInfo()
{
	if(linestr!=NULL)	free(linestr);
	FreeMemory();
}

void OrthSiteInfo::SetChainFile(char *chainfile)
{
	chainfilename = chainfile;
}

void OrthSiteInfo::SetQuerySpecies(char *qs)
{
	queryspecies = qs;
}

void OrthSiteInfo::SetTargetSpecies(char *ts)
{
	targetspecies = ts;
}

int OrthSiteInfo::LoadChainData()
{
	mc.LoadDataFromFile("chain",chainfilename);
	chainpos=0;
}

int OrthSiteInfo::SetMinChainSize(long pminsize)
{
	minchainsize = pminsize;
}

int OrthSiteInfo::SeekToNextChain()
{
	while(chainpos>0)
	{
		chainpos = mc.ReadLineFromData("chain",&linestr,chainpos);
		chainpos++;
		if(strlen(linestr)==0)
			break;
	}
}

int OrthSiteInfo::ReadChainHeader()
{
	chainpos = mc.ReadLineFromData("chain",&linestr,chainpos);
	chainpos++;
	if(linestr[0]==0)
		return 0;
	
	string s;
	s = linestr;
	stringstream ss(s);
	
	string header;
	long tmpvar;
	
	long chrstrpos;
	string chrheader="chr";
	
	ss >> header;
	ss >> score;
	ss >> tchr;
	ss >> tchrsize;
	ss >> tstrand;
	ss >> tstart;
	ss >> tend;
	ss >> qchr;
	ss >> qchrsize;
	ss >> qstrand;
	ss >> qstart;
	ss >> qend;
	
	if(qstrand==tstrand)
		savestrand='+';
	else
		savestrand='-';

	skip = false;
	
	alignsize = 0;
	conflictsite = 0;
	
	/*
	chrstrpos = tchr.find(chrheader);
	if(chrstrpos>=0)
		tchr = tchr.substr(chrstrpos+chrheader.length());
	else
		skip=true;

	chrstrpos = qchr.find(chrheader);
	if(chrstrpos>=0)
		qchr = qchr.substr(chrstrpos+chrheader.length());
	else
		skip=true;
	*/
	
	cout << header << "\t" << score << "\t" << tchr  << "\t" << tchrsize << "\t" << tstrand << "\t" << tstart << "\t" << tend << "\t" << qchr << "\t" << qchrsize << "\t" << qstrand << "\t" << qstart << "\t" << qend << endl;
	if(tchr.find("-")!=string::npos || tchr.find("Un")!=string::npos || tchr.find("Y")!=string::npos || tchr.find("M")!=string::npos || tchr.find("_")!=string::npos || tchr.find("random")!=string::npos || tchr.find("hap")!=string::npos || tchr.find(".")!=string::npos || qchr.find("-")!=string::npos || qchr.find("random")!=string::npos || qchr.find("hap")!=string::npos || qchr.find(".")!=string::npos || qchr.find("Un")!=string::npos || qchr.find("_")!=string::npos || qchr.find("Y")!=string::npos || qchr.find("M")!=string::npos)
	{
		skip = true;
		cout << "Skipped by chromosome name" << endl;
	}
	if(tstrand=="+")
		tstart=tstart+1;
	else 
	{
		tmpvar = tchrsize-tend;
		tend = tchrsize-tstart;
		tstart = tmpvar;
	}
	if(qstrand=="+")
		qstart=qstart+1;
	else 
	{
		tmpvar = qchrsize-qend;
		qend = qchrsize-qstart;
		qstart = tmpvar;
	}
	tspansize = tend-tstart+1;
	qspansize = qend-qstart+1;

	if(tspansize<minchainsize||qspansize<minchainsize)
	{
		cout << "Skipped by chain size" << endl;
		skip=true;
	}

	return 1;
}

int OrthSiteInfo::LoadChainDetails()
{
	long match;
	long qmatchpos,tmatchpos;
	long compareid;
	long t_insert=1;
	long q_insert=1;
	conflictsite = 0;
	alignsize = 0;
	
	while(t_insert>=0)
	{
		chainpos = mc.ReadLineFromData("chain",&linestr,chainpos);
		chainpos++;
		string s;
		s = linestr;
		stringstream ss(s);

		ss >> match;
		alignsize += match;
		t_insert = -1;
		q_insert = -1;
		if(alignsize<tspansize)
		{
			ss >> t_insert;
			ss >> q_insert;
			alignsize += t_insert;
		}
		for(compareid=0;compareid<match;compareid++)
		{
			if(tstrand=="+")
			{
				tmatchpos = tstart+compareid;
			}
			else
			{
				tmatchpos = tend-compareid;
			}
			if(qstrand=="+")
			{
				qmatchpos = qstart+compareid;
			}
			else
			{
				qmatchpos = qend-compareid;
			}
			if(running_chr_query_target_data[qmatchpos-1].chr[0]!=0)
			{
				conflictsite++;
				continue;
			}
			
			strcpy(running_chr_target_query_data[tmatchpos-1].chr,qchr.c_str());
			strcpy(running_chr_query_target_data[qmatchpos-1].chr,tchr.c_str());
			running_chr_target_query_data[tmatchpos-1].orthsite = qmatchpos;
			running_chr_target_query_data[tmatchpos-1].strand = savestrand;
			running_chr_query_target_data[qmatchpos-1].orthsite = tmatchpos;
			running_chr_query_target_data[qmatchpos-1].strand = savestrand;
		}
		if(t_insert==1&&q_insert==1)
		{
			if(tstrand=="+")
			{
				tmatchpos = tstart+match;
			}
			else
			{
				tmatchpos = tend-match;
			}
			if(qstrand=="+")
			{
				qmatchpos = qstart+match;
			}
			else
			{
				qmatchpos = qend-match;
			}
			if(running_chr_query_target_data[qmatchpos-1].chr[0]==0)
			{
				strcpy(running_chr_target_query_data[tmatchpos-1].chr,qchr.c_str());
				strcpy(running_chr_query_target_data[qmatchpos-1].chr,tchr.c_str());
				running_chr_target_query_data[tmatchpos-1].orthsite = qmatchpos;
				running_chr_target_query_data[tmatchpos-1].strand = savestrand;
				running_chr_query_target_data[qmatchpos-1].orthsite = tmatchpos;
				running_chr_query_target_data[qmatchpos-1].strand = savestrand;
			}
			else
				conflictsite++;
		}
		
		if(t_insert<LINE_SIZE)
		{
			for(compareid=0;compareid<t_insert;compareid++)
			{
				if(tstrand=="+")
				{
					tmatchpos = tstart+match+compareid;
				}
				else
				{
					tmatchpos = tend-match-compareid;
				}
				if(running_chr_target_query_data[tmatchpos-1].chr[0]==0)
					strcpy(running_chr_target_query_data[tmatchpos-1].chr,qchr.c_str());
			}
		}

		if(q_insert<LINE_SIZE)
		{
			for(compareid=0;compareid<q_insert;compareid++)
			{
				if(qstrand=="+")
				{
					qmatchpos = qstart+match+compareid;
				}
				else
				{
					qmatchpos = qend-match-compareid;
				}
				if(running_chr_query_target_data[qmatchpos-1].chr[0]==0)
					strcpy(running_chr_query_target_data[qmatchpos-1].chr,tchr.c_str());
			}
		}
		
		if(t_insert>=0)
		{
			if(tstrand=="+")
				tstart = tstart + match + t_insert;
			else
				tend = tend - match - t_insert;
			if(qstrand=="+")
				qstart = qstart + match + q_insert;
			else
				qend = qend - match - q_insert;
		}
	}
	chainpos = mc.ReadLineFromData("chain",&linestr,chainpos);
	chainpos++;
}

int OrthSiteInfo::GetConflicts()
{
	long qmatchpos;
	conflictsite = 0;
	
	for(qmatchpos=qstart;qmatchpos<qend;qmatchpos++)
	{
		if(running_chr_query_target_data[qmatchpos-1].chr[0]!=0)
			conflictsite++;
	}

	if((double)conflictsite/(double)tspansize<=0.2||(tspansize>=1000000&&qspansize>=1000000&&(tspansize-conflictsite)>=100000))
	{
		return 0;
	}
	else
	{
		skip = true;
		return 1;
	}
}

void OrthSiteInfo::SetRelatedFileNames()
{
	running_chr_query_target_name = "" + qchr + "." + queryspecies + "." + targetspecies;
	running_chr_target_query_name = "" + tchr + "." + targetspecies + "." + queryspecies;
	qfilename = workingdir + "/OrthSiteFile/" + running_chr_query_target_name;
	tfilename = workingdir + "/OrthSiteFile/" + running_chr_target_query_name;
	qchrfile = workingdir + "/"+queryspecies+"/" + qchr+".fa";
	tchrfile = workingdir + "/"+targetspecies+"/" + tchr+".fa";
}

void OrthSiteInfo::FreeMemory()
{
	long elements;
	long elementid;
	
	elements = chr_query_target_data.size();
	for(elementid=0;elementid<elements;elementid++)
	{
		free(chr_query_target_data[elementid]);
	}
}

long OrthSiteInfo::LoadOrthSiteFromFile(string filename)
{
	FILE *pfile;
	long ntsize;
	pfile = fopen(filename.c_str(),"rb");
	
	tmp_chr_data = NULL;
	
	if(pfile!=NULL)
	{
		fseek(pfile,0,SEEK_END);
		ntsize = ftell(pfile)/sizeof(orthinfo);
		tmp_chr_data = (orthinfo *)malloc(sizeof(orthinfo)*ntsize);
		memset(tmp_chr_data,0,ntsize*sizeof(orthinfo));
		fread(tmp_chr_data,sizeof(orthinfo),ntsize,pfile);
		fclose(pfile);
		return ntsize;
	}
	else
		return 0;
}

long OrthSiteInfo::FindPositionInOrthSiteVector(string search_query_target)
{
	long elements;
	long elementid;
	
	elements = chr_query_target_name.size();
	for(elementid=0;elementid<elements;elementid++)
	{
		if(chr_query_target_name[elementid] == search_query_target)
			return elementid;
	}
	return -1;
}

int OrthSiteInfo::SaveOrthSiteToVector()
{
	long vectorpos;
	
	/*
	if(conflictsite/tspansize>0.01||tspansize<minchainsize)
		return 0;
	*/
	
	vectorpos = FindPositionInOrthSiteVector(running_chr_query_target_name);
	if(vectorpos==-1)
	{
		chr_query_target_name.push_back(running_chr_query_target_name);
		chr_query_target_data.push_back(running_chr_query_target_data);
		chr_genome_size.push_back(qchrsize);
	}
	/*
	else
	{
		//memcpy(chr_query_target_data[vectorpos],running_chr_query_target_data,qchrsize*sizeof(orthinfo));
		if(chr_query_target_data[vectorpos]!=NULL) free(chr_query_target_data[vectorpos]);
		chr_query_target_data[vectorpos] = running_chr_query_target_data;
		running_chr_query_target_data = NULL;
	}
	*/
	
	vectorpos = FindPositionInOrthSiteVector(running_chr_target_query_name);
	if(vectorpos==-1)
	{
		chr_query_target_name.push_back(running_chr_target_query_name);
		chr_query_target_data.push_back(running_chr_target_query_data);
		chr_genome_size.push_back(tchrsize);
	}
	/*
	else
	{
		//memcpy(chr_query_target_data[vectorpos],running_chr_target_query_data,tchrsize*sizeof(orthinfo));
		if(chr_query_target_data[vectorpos]!=NULL) free(chr_query_target_data[vectorpos]);
		chr_query_target_data[vectorpos] = running_chr_target_query_data;
		running_chr_target_query_data = NULL;
	}
	*/
}

int OrthSiteInfo::SaveOrthSiteVectorToFile()
{
	FILE *pfile;
	long elementid,elementcount;
	string filename;
	elementcount = chr_query_target_data.size();
	for(elementid=0; elementid<elementcount; elementid++)
	{
		filename = workingdir + "/OrthSiteFile/" + chr_query_target_name[elementid];
		pfile = fopen(filename.c_str(),"wb");
		fwrite(chr_query_target_data[elementid],sizeof(orthinfo),chr_genome_size[elementid],pfile);
		fclose(pfile);
	}
	return 1;
}

int OrthSiteInfo::GenerateQueryOrthSite()
{
	running_chr_query_target_data = (orthinfo*)malloc(sizeof(orthinfo)*qchrsize);
	running_chr_query_seq = (char*)malloc(sizeof(char)*(qchrsize+1));
	memset(running_chr_query_target_data,0,qchrsize*sizeof(orthinfo));
	memset(running_chr_query_seq,0,(qchrsize+1)*sizeof(char));
}

int OrthSiteInfo::GenerateTargetOrthSite()
{
	running_chr_target_query_data = (orthinfo*)malloc(sizeof(orthinfo)*tchrsize);
	running_chr_target_seq = (char*)malloc(sizeof(char)*(tchrsize+1));
	memset(running_chr_target_query_data,0,tchrsize*sizeof(orthinfo));
	memset(running_chr_target_seq,0,(tchrsize+1)*sizeof(char));
}

void OrthSiteInfo::FreeRunningMemory()
{
	if(running_chr_query_target_data!=NULL)
	{
		free(running_chr_query_target_data);
		running_chr_query_target_data=NULL;
	}
	if(running_chr_target_query_data!=NULL)
	{
		free(running_chr_target_query_data);
		running_chr_target_query_data=NULL;
	}
	if(running_chr_query_seq!=NULL)
	{
		free(running_chr_query_seq);
		running_chr_query_seq=NULL;
	}
	if(running_chr_target_seq!=NULL)
	{
		free(running_chr_target_seq);
		running_chr_target_seq=NULL;
	}
}

void OrthSiteInfo::LoadQueryFromVector(long index)
{
	//memcpy(running_chr_query_target_data,chr_query_target_data[index],qchrsize*sizeof(orthinfo));
	running_chr_query_target_data = (orthinfo*)chr_query_target_data[index];
}

void OrthSiteInfo::LoadTargetFromVector(long index)
{
	//memcpy(running_chr_target_query_data,chr_query_target_data[index],tchrsize*sizeof(orthinfo));
	running_chr_target_query_data = (orthinfo*)chr_query_target_data[index];
}

void OrthSiteInfo::LoadQuerySeq()
{
   int errorflag;
   long readstart,readend;
   long linelen;
   char tmpfilename[1024];
   strcpy(tmpfilename,qchrfile.c_str());
   linelen = GetFastaFileLineLength(qchrfile);
   ReadFasta(tmpfilename,"+",1,qchrsize,linelen,1,running_chr_query_seq,&readstart,&readend,&errorflag);
   SeqUpper(running_chr_query_seq);
   for(readstart=0;readstart<qchrsize;readstart++)
   {
		running_chr_query_target_data[readstart].nt = running_chr_query_seq[readstart];
   }
   free(running_chr_query_seq);
   running_chr_query_seq=NULL;
}

void OrthSiteInfo::LoadTargetSeq()
{
   int errorflag;
   long readstart,readend;
   long linelen;
   char tmpfilename[1024];
   strcpy(tmpfilename,tchrfile.c_str());
   linelen = GetFastaFileLineLength(tchrfile);
   ReadFasta(tmpfilename,"+",1,tchrsize,linelen,1,running_chr_target_seq,&readstart,&readend,&errorflag);
   SeqUpper(running_chr_target_seq);
   for(readstart=0;readstart<tchrsize;readstart++)
   {
		running_chr_target_query_data[readstart].nt = running_chr_target_seq[readstart];
   }
   free(running_chr_target_seq);
   running_chr_target_seq=NULL;
}

int OrthSiteInfo::GetFastaFileLineLength(string filename)
{
	ifstream in;
	string data;
	in.open(filename.c_str(),ios::in);
	getline(in,data);
	getline(in,data);
	in.close();
	return data.length();
}


