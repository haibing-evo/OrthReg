#include <iostream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

#include <algorithm>
#include <cctype>
#include <clocale>
#include <cstring>
#include <iterator>
//#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#ifdef _MAC
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

#include "MemoryCache.h"

char tspecies1[24];
char tspecies2[24];
char qspecies[24];
char workingdir[1024];

struct orthinfo
{
	char chr[6];
	long orthsite;
	char nt;
	char strand;
	float encode;
};

vector<string> loadchrinfo(string filename)
{
	vector<string> retchr;
	string linedata;
	ifstream infile;
	infile.open(filename.c_str(),ios::in);
	while(getline(infile, linedata))
	{
		if(linedata.length()>0)
			retchr.push_back(linedata);
	}
	infile.close();
	return(retchr);
}

int main(int argc, char* argv[])
{
  if(argc!=5)
  {
     printf("Usage:\n");
     printf("%s qspecies qchrfile tspeices1 tspecies2\n",argv[0]);
     return 0;
  }
  
  FILE *pfile1,*pfile2,*pfile3,*pfile4;
  long i,j;
  
  vector<string> qchr;
  char t1chr[3]="";
  char t2chr[3]="";
  
  long qchrsize,tchrsize1,tchrsize2;
  struct orthinfo *tflag1=NULL,*tflag2=NULL,*qflag1=NULL,*qflag2=NULL,*t2_nt_flag;
  string chrfile;
  
  string mcfilename,mcname;
  string dir;
  long long pos;
  
  char tflagfile1[1024],qflagfile1[1024],tflagfile2[1024],qflagfile2[1024];
  
  strcpy(qspecies,argv[1]);
  chrfile = argv[2];
  strcpy(tspecies1,argv[3]);
  strcpy(tspecies2,argv[4]);
  getcwd(workingdir,1024);
  
  qchr = loadchrinfo(argv[2]);

  MemoryCache qt;
  MemoryCache mc_tmp;
  
  for(i=0;i<qchr.size();i++)
  {
	   MemoryCache mc;

       sprintf(qflagfile1,"%s/OrthSiteFile/%s.%s.%s",workingdir,qchr[i].c_str(),qspecies,tspecies1);
       sprintf(tflagfile1,"%s/OrthSiteFile/%s.%s.%s",workingdir,qchr[i].c_str(),qspecies,tspecies2);

       pfile1 = fopen(qflagfile1,"rb");
       if(pfile1!=NULL)
       {
           fseek(pfile1,0,SEEK_END);
           qchrsize = ftell(pfile1)/sizeof(struct orthinfo);
           
           printf("qchrsize=%ld\n",qchrsize);
           
           qflag1 = (struct orthinfo*)malloc(sizeof(struct orthinfo)*qchrsize);
           tflag1 = (struct orthinfo*)malloc(sizeof(struct orthinfo)*qchrsize);
           
           memset(qflag1,0,qchrsize*sizeof(struct orthinfo));
           memset(tflag1,0,qchrsize*sizeof(struct orthinfo));
           
           fseek(pfile1,0,SEEK_SET);
           fread(qflag1,sizeof(struct orthinfo),qchrsize,pfile1);
           fclose(pfile1);
           printf("Read qflag1\n");
       }
       else
           continue;
       
       t1chr[0]=0;
       t2chr[0]=0;
           
       for(j=0;j<qchrsize;j++)
       {
           if(qflag1[j].chr[0]!=0 && qflag1[j].orthsite!=0)
           {
               if(strcmp(t1chr,qflag1[j].chr)!=0)
               {
                  strcpy(t1chr,qflag1[j].chr);
				  //mcname = "chr";
				  mcname = t1chr;
				  mcname += ".";
				  mcname += tspecies1;
				  mcname += ".";
				  mcname += tspecies2;
				  
				  mcfilename = workingdir;
				  mcfilename += "/OrthSiteFile/";
				  mcfilename += mcname;
				  
				  pos = mc.Find(mcname);
				  
				  if(pos<0)
				  {
					mc.LoadDataFromFile(mcname,mcfilename);
				  }
  				  pos = mc.Find(mcname);
				  qflag2 = (struct orthinfo *)mc.GetData(pos);
				  tchrsize1 = mc.GetDatasize(pos)/sizeof(struct orthinfo);
               }
               
               if(strcmp(t2chr,qflag2[abs(qflag1[j].orthsite)-1].chr)!=0 && qflag2[abs(qflag1[j].orthsite)-1].chr[0]!=0 && qflag2[abs(qflag1[j].orthsite)-1].orthsite!=0)
               {
                  strcpy(t2chr,qflag2[abs(qflag1[j].orthsite)-1].chr);
				  //mcname = "chr";
				  mcname = t2chr;
				  mcname += ".";
				  mcname += tspecies2;
				  mcname += ".";
				  mcname += qspecies;
				  pos = qt.Find(mcname);
				  if(pos<0)
				  {
					  sprintf(tflagfile2,"%s/OrthSiteFile/%s.%s.%s",workingdir,t2chr,tspecies2,qspecies);
					  mcfilename = workingdir;
					  mcfilename += "/OrthSiteFile/";
					  mcfilename += mcname;
					  
					  if(qt.LoadDataFromFile(mcname,mcfilename)==0)
					  {
						  sprintf(tflagfile2,"%s/OrthSiteFile/%s.%s.%s",workingdir,t2chr,tspecies2,tspecies1);
						  pfile4 = fopen(tflagfile2,"rb");
						  if(pfile4!=NULL)
						  {
							   fseek(pfile4,0,SEEK_END);
							   tchrsize2 = ftell(pfile4)/sizeof(struct orthinfo);
							   fclose(pfile4);
						  }
						  tflag2 = (struct orthinfo*)malloc(sizeof(struct orthinfo)*tchrsize2);
						  memset(tflag2,0,tchrsize2*sizeof(struct orthinfo));
						  
						  qt.AddData(mcname,(char*)tflag2,tchrsize2*sizeof(struct orthinfo));
					  }
					  else
					  {
						pos = qt.Find(mcname);
						tflag2 = (struct orthinfo*)qt.GetData(pos);
						tchrsize2 = qt.GetDatasize(pos)/sizeof(struct orthinfo);
					  }
				  }
				  else
				  {
					tflag2 = (struct orthinfo*)qt.GetData(pos);
					tchrsize2 = qt.GetDatasize(pos)/sizeof(struct orthinfo);
				  }
				  
				  //mcname = "chr";
				  mcname = t2chr;
				  mcname += ".";
				  mcname += tspecies2;
				  mcname += ".";
				  mcname += tspecies1;
				  pos = mc_tmp.Find(mcname);
				  if(pos<0)
				  {
					mcfilename = workingdir;
					mcfilename += "/OrthSiteFile/";
					mcfilename += mcname;
					mc_tmp.LoadDataFromFile(mcname,mcfilename);
					pos = mc_tmp.Find(mcname);
				  }
				  t2_nt_flag = (struct orthinfo*)mc_tmp.GetData(pos);
				  
               }
               
               if(qflag2[abs(qflag1[j].orthsite)-1].orthsite!=0)
               {
                   fflush(stdout);
				   strcpy(tflag1[j].chr,qflag2[abs(qflag1[j].orthsite)-1].chr);
				   tflag1[j].orthsite = qflag2[abs(qflag1[j].orthsite)-1].orthsite;
				   tflag1[j].nt = qflag1[j].nt;
				   
				   strcpy(tflag2[abs(tflag1[j].orthsite)-1].chr,qchr[i].c_str());
				   tflag2[abs(tflag1[j].orthsite)-1].orthsite=j+1;
				   tflag2[abs(tflag1[j].orthsite)-1].nt = t2_nt_flag[abs(tflag1[j].orthsite)-1].nt;
				   
				   if(qflag1[j].strand==qflag2[abs(qflag1[j].orthsite)-1].strand)
				   {
						tflag1[j].strand='+';
						tflag2[abs(tflag1[j].orthsite)-1].strand='+';
				   }
				   else
				   {
						tflag1[j].strand='-';
						tflag2[abs(tflag1[j].orthsite)-1].strand='-';
				   }
			   }
           }
       }
       
       pfile3 = fopen(tflagfile1,"wb");
       fwrite(tflag1,sizeof(struct orthinfo),qchrsize,pfile3);
       fclose(pfile3);
       
       free(qflag1);
       free(tflag1);
  }
  dir = workingdir;
  dir += "/OrthSiteFile/";
  qt.SaveVectorToFiles(dir);

  return 1;
}
