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
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

#include "xiestring.h"
#include "xiesequence.h"
#include "xiefile.h"

#ifdef _MAC
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

using namespace std;

int query_chromsome_count;
int target_chromsome_count;

vector<string> pre_vec_tchr;
vector<string> pre_vec_qchr;

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

int returntchrid(string chrname)
{
  for(int i=0; i<pre_vec_tchr.size();i++)
  {
	if(chrname == pre_vec_tchr[i])
		return i+1;
  }
  return 0;
}

int returnqchrid(string chrname)
{
  for(int i=0; i<pre_vec_qchr.size();i++)
  {
	if(chrname == pre_vec_qchr[i])
		return i+1;
  }
  return 0;
}

/*
int returntchrid(char *chrname)
{
  if(strcmp(chrname,"M")==0||strcmp(chrname,"Y")==0)
	return 0;
  if(strcmp(chrname,"X")==0)
     return target_chromsome_count;
  else
     return atoi(chrname);
}

int returnqchrid(char *chrname)
{
  if(strcmp(chrname,"M")==0||strcmp(chrname,"Y")==0)
	return 0;
  if(strcmp(chrname,"X")==0)
     return query_chromsome_count;
  else
     return atoi(chrname);
}
*/

char workingdir[1024];

int main(int argc, char* argv[])
{
  FILE *orth_segments,*bedfile,*bedfilelist;
  gzline gzbedfile;
  FILE *pqtflag,*ptqflag;

  char tmpstr[128],*linestr,bedfilename[1024],*bedfileline,experiment[1024],orth_segments_filename[1024];
  long startpos,endpos,checkpos;
  long qchrsize;
  vector<long> tchrsize;
  char qtfile[1024],tqfile[1024];
  char qchr[128]="",tchr[128]="";
  int pretchrid,preqchrid;
  int tchrid,qchrid;
  
  struct orthinfo **tqflag,*qtflag=NULL,**pretqflag,**preqtflag;
  
  char outputstr[102400];
  int flag;
  char Qchr[3];
  long chrsize,cycle,Tstart,Tend,Qstart,Qend,NegNum;
  int max_insertion;
  int insertionid;
  string qchrfilename;
  string tchrfilename;

  getcwd(workingdir,1024);
    
  if(argc!=7)
  {
     printf("Usage:\n");
     printf("%s bed_file_list qSpecies tSpecies qChrfile tChrfile Max_Insertion_Interval_To_Merge\n",argv[0]);
     return 0;
  }

  max_insertion = atoi(argv[6]);
  if(max_insertion<0)
	max_insertion = 5;

  qchrfilename = argv[4];
  tchrfilename = argv[5];
  pre_vec_qchr = loadchrinfo(qchrfilename);
  query_chromsome_count = pre_vec_qchr.size();
  pre_vec_tchr = loadchrinfo(tchrfilename);
  target_chromsome_count = pre_vec_tchr.size();
  
  bedfilelist = fopen(argv[1],"rb");
  if(bedfilelist==NULL)
      return 0;
      
  tqflag = (struct orthinfo **)malloc(sizeof(struct orthinfo*)*target_chromsome_count);
  pretqflag = (struct orthinfo **)malloc(sizeof(struct orthinfo*)*target_chromsome_count);
  preqtflag = (struct orthinfo **)malloc(sizeof(struct orthinfo*)*query_chromsome_count);

  //load tqflag into memory
  for(pretchrid=0;pretchrid<target_chromsome_count;pretchrid++)
  {
	cout << "pre_vec_tchr = " << pre_vec_tchr[pretchrid] << endl;
	cout << "argv[3] = " << argv[3] << endl;
	cout << "argv[2] = " << argv[2] << endl;
    sprintf(tqfile,"%s/OrthSiteFile/%s.%s.%s",workingdir,pre_vec_tchr[pretchrid].c_str(),argv[3],argv[2]);
    cout << tqfile << endl;
    tchrsize.push_back(0);
    
    ptqflag = fopen(tqfile,"rb");
    if(ptqflag!=NULL)
    {
       fseek(ptqflag,0,SEEK_END );
       tchrsize[pretchrid]=ftell(ptqflag)/sizeof(struct orthinfo); 
	   cout << "tchrsize[" << pretchrid << "] = " << tchrsize[pretchrid] << endl;
       fseek(ptqflag,0,SEEK_SET);
       pretqflag[pretchrid]=(struct orthinfo*)malloc(sizeof(struct orthinfo)*tchrsize[pretchrid]);
       tqflag[pretchrid]=(struct orthinfo*)malloc(sizeof(struct orthinfo)*tchrsize[pretchrid]);
       //memset(pretqflag[pretchrid],0,tchrsize[pretchrid]*sizeof(struct orthinfo));
       fread(pretqflag[pretchrid],sizeof(struct orthinfo),tchrsize[pretchrid],ptqflag);
       fclose(ptqflag);           
    }
    else
    {
       pretqflag[pretchrid]=NULL;
    }
    cout << "Done" << endl;
  }

  //load qtflag into memory
  for(preqchrid=0;preqchrid<query_chromsome_count;preqchrid++)
  {
    sprintf(qtfile,"%s/OrthSiteFile/%s.%s.%s",workingdir,pre_vec_qchr[preqchrid].c_str(),argv[2],argv[3]);
    cout << qtfile << endl;

    pqtflag = fopen(qtfile,"rb");
    if(pqtflag!=NULL)
    {
       fseek(pqtflag,0,SEEK_END );
       qchrsize=ftell(pqtflag)/sizeof(struct orthinfo);
       fseek(pqtflag,0,SEEK_SET);
       preqtflag[preqchrid]=(struct orthinfo*)malloc(sizeof(struct orthinfo)*qchrsize);
       memset(preqtflag[preqchrid],0,qchrsize*sizeof(struct orthinfo));
       fread(preqtflag[preqchrid],sizeof(struct orthinfo),qchrsize,pqtflag);
       fclose(pqtflag);           
    }
    else
    {
       preqtflag[preqchrid]=NULL;
    }
  }
  printf("Cached files has been loaded\n");

  //read bed file list
  bedfileline = (char*)malloc(1);
  while(!IsEndOfFile(bedfilelist))
  {
     ReadLine(bedfilelist,&bedfileline);
     if(bedfileline[0]==0) continue;

     ParseString(bedfileline,2,"\t",experiment);
     ParseString(bedfileline,1,"\t",bedfilename);

     //get next bed file from list
     //bedfile = fopen(bedfilename,"rb");
     //if(bedfile==NULL)
     //   continue;
        
     gz_open_linemode(bedfilename,&gzbedfile);

     fprintf(stdout,"start analysis on bed file: %s\n",bedfilename);
     fflush(stdout);        

    //init tqflag for each bed file
    for(pretchrid=0;pretchrid<target_chromsome_count;pretchrid++)
       memcpy(tqflag[pretchrid],pretqflag[pretchrid],tchrsize[pretchrid]*sizeof(struct orthinfo));
        
    linestr = (char*)malloc(1);
    while(gz_readline(&gzbedfile,&linestr)>0)  
    {
       if(linestr[0]==0) continue;
       
       ParseString(linestr,1,"\t",qchr);
       /*
       if(memcmp(qchr,"chr",3)==0)
          StringRight(qchr,strlen(qchr)-3,qchr);
       */

       //find qchr from bed file
       qchrid = returnqchrid(qchr);
       if(qchrid==0) continue;
       
       qtflag = preqtflag[qchrid-1];
       
       if(qtflag==NULL) continue;
       
       ParseString(linestr,2,"\t",tmpstr);
       startpos = atol(tmpstr);
       ParseString(linestr,3,"\t",tmpstr);
       endpos = atol(tmpstr);
       
       for(checkpos=startpos;checkpos<=endpos;checkpos++)
       {
          //mark encode flag in tqflag
          if(qtflag[checkpos-1].chr[0]!=0&&qtflag[checkpos-1].orthsite!=0)
          {
             tchrid = returntchrid(qtflag[checkpos-1].chr);
			 if(tchrid==0) continue;
				 tqflag[tchrid-1][abs(qtflag[checkpos-1].orthsite)-1].encode=1;
          }
       }
    }
    free(linestr);
    //fclose(bedfile); 
    gz_close_linemode(&gzbedfile);
    fprintf(stdout,"finish reading bedfile: %s\n",bedfilename);
    fflush(stdout);
    
    //output orth_segments
    sprintf(orth_segments_filename,"%s/OrthMapFile/%s.orth.segments",workingdir,experiment);
    orth_segments = fopen(orth_segments_filename,"wb");
    flag=0;
    NegNum=0;
    Tstart=0;
    Tend=0;
    Qstart=0;
    Qend=0;
    outputstr[0]=0;
    fprintf(stdout,"mapping and marking\n");
    fflush(stdout);

    for(pretchrid=0;pretchrid<target_chromsome_count;pretchrid++)
    {
       for(cycle=0;cycle<tchrsize[pretchrid];cycle++)
       {
         if(flag==0 && tqflag[pretchrid][cycle].encode==1)
         {
           flag=1;
           Tstart=cycle+1;
           Tend=cycle+1;
           strcpy(Qchr,tqflag[pretchrid][cycle].chr);
           Qstart=tqflag[pretchrid][cycle].orthsite;
           Qend=tqflag[pretchrid][cycle].orthsite;
           
           qchrid = returnqchrid(tqflag[pretchrid][cycle].chr);
           qtflag = preqtflag[qchrid-1];
		   if(tqflag[pretchrid][cycle].strand=='+')
		   {
			  if(tqflag[pretchrid][cycle].nt!=qtflag[Qend-1].nt)
				NegNum++;
		   }
		   else
		   {
			  if(!((tqflag[pretchrid][cycle].nt=='A'&&qtflag[Qend-1].nt=='T')||(tqflag[pretchrid][cycle].nt=='C'&&qtflag[Qend-1].nt=='G')||(tqflag[pretchrid][cycle].nt=='T'&&qtflag[Qend-1].nt=='A')||(tqflag[pretchrid][cycle].nt=='G'&&qtflag[Qend-1].nt=='C')))
				NegNum++;
		   }
         }
         else if(flag==1 && tqflag[pretchrid][cycle].encode==1)
         {
           Tend=cycle+1;
           Qend=tqflag[pretchrid][cycle].orthsite;

           qchrid = returnqchrid(tqflag[pretchrid][cycle].chr);
           qtflag = preqtflag[qchrid-1];
		   if(tqflag[pretchrid][cycle].strand=='+')
		   {
			  if(tqflag[pretchrid][cycle].nt!=qtflag[Qend-1].nt)
				NegNum++;
		   }
		   else
		   {
			  if(!((tqflag[pretchrid][cycle].nt=='A'&&qtflag[Qend-1].nt=='T')||(tqflag[pretchrid][cycle].nt=='C'&&qtflag[Qend-1].nt=='G')||(tqflag[pretchrid][cycle].nt=='T'&&qtflag[Qend-1].nt=='A')||(tqflag[pretchrid][cycle].nt=='G'&&qtflag[Qend-1].nt=='C')))
				NegNum++;
		   }
         }
         else if(flag==1 && tqflag[pretchrid][cycle].encode==0)
         {
		   //merge fragments with separating insertions <= 10bp
		   insertionid=1;
		   if(cycle+max_insertion<tchrsize[pretchrid])
		   {
				while(insertionid<=max_insertion)
				{
					if(tqflag[pretchrid][cycle+insertionid].encode==1)
					{
						cycle = cycle+insertionid-1;
						break;
					}
					insertionid++;
				}
				if(insertionid<=max_insertion)
					continue;
		   }
           if(Qend>Qstart)
             sprintf(outputstr,"%s%s\t+\t%ld\t%ld\t%s\t+\t%ld\t%ld\t%ld\n",outputstr,pre_vec_tchr[pretchrid].c_str(),Tstart,Tend,Qchr,Qstart,Qend,NegNum);
           else
             sprintf(outputstr,"%s%s\t+\t%ld\t%ld\t%s\t-\t%ld\t%ld\t%ld\n",outputstr,pre_vec_tchr[pretchrid].c_str(),Tstart,Tend,Qchr,Qend,Qstart,NegNum);
           
           if(strlen(outputstr)>=100000)
           {
             fprintf(orth_segments,"%s",outputstr);
             fflush(orth_segments);
             outputstr[0]=0;
           }
           
           flag=0;
           Tstart=0;
           Tend=0;
           Qstart=0;
           Qend=0;
           NegNum=0;
         }
       }
    }
    if(strlen(outputstr)>0)
    {
       fprintf(orth_segments,"%s",outputstr);
       fflush(orth_segments);
       outputstr[0]=0;
    }
    fclose(orth_segments);
    fprintf(stdout,"finished current bed file\n");
    fflush(stdout);

  }
  free(bedfileline);
  fclose(bedfilelist);

  for(pretchrid=0;pretchrid<target_chromsome_count;pretchrid++)
  {
     free(pretqflag[pretchrid]);
     free(tqflag[pretchrid]);
  }
  for(preqchrid=0;preqchrid<query_chromsome_count;preqchrid++)
  {
     free(preqtflag[preqchrid]);
  }
  free(pretqflag);
  free(tqflag);
  free(preqtflag);

  return 1;
}

