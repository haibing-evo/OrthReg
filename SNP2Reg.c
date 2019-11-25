#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <algorithm>
#include <cctype>
#include <clocale>
#include <cstring>
#include <iostream>
#include <iterator>
//#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "xiefile.h"
#include "xiestring.h"
#include "xiesequence.h"

using namespace std;

typedef struct geneinformation
{
  char gene[128];
  char mid[128];
  char strand;
  int exonid;
  long exonstart;
  long exonend;
  long cdsstart;
  long cdsend;
  int seqtype;
  int phase;
  struct geneinformation * pNext;
} geneinfo;

typedef struct idnamelist
{
  char id[128];
  char gene[128];
  char name[128];
  struct idnamelist *pNext;
} idname;

typedef struct
{
  int chrid;
  long snppos;
  char snpbase1[8];
  char snpbase2[8];
  int mt;
} mtsnp;

geneinfo **globalhead_mt[8],**curnode_mt[8],**tailnode_mt[8];

long maxseglen[8];
idname *idnamehead_mt[8],*pidname_mt[8],*pidtail_mt[8];
char retname[8][128]={"","","","","","","",""};
char mtgene1[8][128]={"","","","","","","",""};
char mtname1[8][128]={"","","","","","","",""};
char mtgene2[8][128]={"","","","","","","",""};
char mtname2[8][128]={"","","","","","","",""};

geneinfo *startnode[8][20];
vector<string> chr;

char **chrseq;

#define SNPSIZE 1000
#define SNPRETSIZE SNPSIZE*400

char threadsdata[8][SNPRETSIZE][512];
mtsnp snpdata[8][SNPSIZE];
int retstatus[8][SNPSIZE];
int pthreadid[8]={0,1,2,3,4,5,6,7};


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

int getchrindex(string chrname)
{
  for(int i=0; i<chr.size();i++)
  {
	if(chrname == chr[i])
		return i;
  }
  return -1;
}

int GetFastaFileLineLength(char *filename)
{
	ifstream in;
	string data;
	in.open(filename,ios::in);
    getline(infile, data);
    getline(infile, data);
	return data.length();
}


void load_genome_seq(char *dir)
{
  int i;
  char filename[1024];
  int linelen;
  long chrlength;
  long readstart,readend;
  int errorflag;
  
  chrseq = (char**)malloc(sizeof(char*)*chr.size());
  for(i=0;i<chr.size();i++)
  {
     sprintf(filename,"%s/%s.fa",dir,chr[i].c_str());
     fprintf(stderr,"Loading genome sequence for %s\n",filename);
     linelen = GetFastaFileLineLength(filename);
     chrlength = GetChrLength(filename,linelen,1);
     fprintf(stderr,"linelen=%ld  chrlength=%ld\n",linelen,chrlength);
     chrseq[i] = (char*)malloc(sizeof(char)*(chrlength+1));
     ReadFasta(filename,"+",1,chrlength,linelen,1,chrseq[i],&readstart,&readend,&errorflag);
     fprintf(stderr,"Loaded genome sequence for %s\n",filename);
  }
}

void unload_genome_seq()
{
  int i;
  for(i=0;i<chr.size();i++)
  {
     free(chrseq[i]);
  }
  free(chrseq);
}

void resetthreaddata(int ti)
{
  int i;
  for(i=0;i<SNPSIZE;i++)
     snpdata[ti][i].chrid=-1;
}

void resetsnpdata()
{
  int i,j;
  for(i=0;i<8;i++)
  {
   for(j=0;j<SNPSIZE;j++)
     retstatus[i][j]=0;
   for(j=0;j<SNPRETSIZE;j++)
     threadsdata[i][j][0]=0;
  }
}

void resetlistdata()
{
  int i;
  for(i=0;i<8;i++)
  {
    globalhead_mt[i]=NULL;
    curnode_mt[i]=NULL;
    tailnode_mt[i]=NULL;
    idnamehead_mt[i]=NULL;
    pidname_mt[i]=NULL;
    pidtail_mt[i]=NULL;
  }
}

void local_codon(char *codonseq,int mismatchpos,char base,char *aa)
{
  char localbase1,localbase2;
  char newcodon[4];
  strcpy(newcodon,codonseq);
  switch(base)
  {
    case 'A':
      newcodon[mismatchpos-1]=base;
      aa[0]=f_codon(newcodon);
      aa[1]=0;
      break;
    case 'C':
      newcodon[mismatchpos-1]=base;
      aa[0]=f_codon(newcodon);
      aa[1]=0;
      break;
    case 'T':
      newcodon[mismatchpos-1]=base;
      aa[0]=f_codon(newcodon);
      aa[1]=0;
      break;
    case 'G':
      newcodon[mismatchpos-1]=base;
      aa[0]=f_codon(newcodon);
      aa[1]=0;
      break;
    case 'Y':
      localbase1='C';
      localbase2='T';
      newcodon[mismatchpos-1]=localbase1;
      aa[0]=f_codon(newcodon);
      aa[1]=0;
      newcodon[mismatchpos-1]=localbase2;
      aa[1]=f_codon(newcodon);
      aa[2]=0;
      break;
    case 'M':
      localbase1='C';
      localbase2='A';
      newcodon[mismatchpos-1]=localbase1;
      aa[0]=f_codon(newcodon);
      aa[1]=0;
      newcodon[mismatchpos-1]=localbase2;
      aa[1]=f_codon(newcodon);
      aa[2]=0;
      break;
    case 'R':
      localbase1='A';
      localbase2='G';
      newcodon[mismatchpos-1]=localbase1;
      aa[0]=f_codon(newcodon);
      aa[1]=0;
      newcodon[mismatchpos-1]=localbase2;
      aa[1]=f_codon(newcodon);
      aa[2]=0;
      break;
    case 'S':
      localbase1='C';
      localbase2='G';
      newcodon[mismatchpos-1]=localbase1;
      aa[0]=f_codon(newcodon);
      aa[1]=0;
      newcodon[mismatchpos-1]=localbase2;
      aa[1]=f_codon(newcodon);
      aa[2]=0;
      break;
    case 'K':
      localbase1='G';
      localbase2='T';
      newcodon[mismatchpos-1]=localbase1;
      aa[0]=f_codon(newcodon);
      aa[1]=0;
      newcodon[mismatchpos-1]=localbase2;
      aa[1]=f_codon(newcodon);
      aa[2]=0;
      break;
    case 'W':
      localbase1='A';
      localbase2='T';
      newcodon[mismatchpos-1]=localbase1;
      aa[0]=f_codon(newcodon);
      aa[1]=0;
      newcodon[mismatchpos-1]=localbase2;
      aa[1]=f_codon(newcodon);
      aa[2]=0;
      break;
    default:
      break;
  }
}

void setnullnode(geneinfo *node)
{
  strcpy(node->gene,"");
  strcpy(node->mid,"");
  node->strand=0;
  node->exonid=-99;
  node->exonstart=0;
  node->exonend=0;
  node->cdsstart=0;
  node->cdsend=0;
  node->phase=-100;
  node->seqtype=-100;
  node->pNext=NULL;
}

void appendidname_mt(int mtid,char * id, char * gene, char * name)
{
   if(idnamehead_mt[mtid]==NULL)
   {
     idnamehead_mt[mtid] = (idname*)malloc(sizeof(idname));
     pidtail_mt[mtid] = idnamehead_mt[mtid];
   }
   else
   {
      //should not compare gene here, because multiple transcript variants of a gene would get error because of the being negelected annotation in the chain
      if(strcmp(pidtail_mt[mtid]->id,id)==0)
      {
         if(strlen(pidtail_mt[mtid]->name)>strlen(name))
            strcpy(pidtail_mt[mtid]->name,name);
         return;
      }
      pidtail_mt[mtid]->pNext=(idname*)malloc(sizeof(idname));
      pidtail_mt[mtid] = pidtail_mt[mtid]->pNext;
   }
   pidname_mt[mtid] = pidtail_mt[mtid];
   strcpy(pidname_mt[mtid]->id,id);
   strcpy(pidname_mt[mtid]->gene,gene);
   strcpy(pidname_mt[mtid]->name,name);
   pidname_mt[mtid]->pNext = NULL;
}

char * findidname_mt(int mtid,char *id,int savepos)
{
  idname * pidcur;
  char lid[3][128];
  int cols,colid;
  
  retname[mtid][0]=0;
  cols= StringColumns(id,",");
  if(cols>0)
  {
     for(colid=0;colid<=cols;colid++)
        ParseString(id,colid+1,",",lid[colid]);
  }

  if(cols==0)
  {
    if(strcmp(mtgene1[mtid],id)==0)
       return mtname1[mtid];
    strcpy(mtgene1[mtid],id);
    strcpy(mtname1[mtid],id);
    for(pidcur=idnamehead_mt[mtid];pidcur!=NULL;pidcur=pidcur->pNext)
    {
       if(strcmp(id,pidcur->id)==0 || strcmp(id,pidcur->gene)==0 )
       {
         strcpy(mtname1[mtid],pidcur->name);
         return mtname1[mtid];
       }
    }
  }
  else
  {
    if(savepos==0)
    {
      if(strcmp(mtgene1[mtid],id)==0)
         return mtname1[mtid];
      strcpy(mtgene1[mtid],id);
      strcpy(mtname1[mtid],id);
    }
    else
    {
      if(strcmp(mtgene2[mtid],id)==0)
         return mtname2[mtid];
      strcpy(mtgene2[mtid],id);
      strcpy(mtname2[mtid],id);
    }
     
    for(colid=0;colid<=cols;colid++)
    {
       for(pidcur=idnamehead_mt[mtid];pidcur!=NULL;pidcur=pidcur->pNext)
       {
          if(lid[colid][0]=='-') continue;
          if(strcmp(lid[colid],pidcur->id)==0 || strcmp(lid[colid],pidcur->gene)==0)
          {
             strcat(retname[mtid],pidcur->name);
             strcat(retname[mtid],",");
             break;
          }
       }
    }
    retname[mtid][strlen(retname[mtid])-1]=0;
    if(savepos==0)
       strcpy(mtname1[mtid],retname[mtid]);
    else
       strcpy(mtname2[mtid],retname[mtid]);
    
    return retname[mtid];
  }

  return id;
}

void deleteidname_mt (int mtid)
{
  idname * pidcur;
  for(pidcur=idnamehead_mt[mtid];pidcur!=NULL;pidcur=pidname_mt[mtid])
  {
     pidname_mt[mtid] = pidcur->pNext;
     free(pidcur);
  }
}

void appendnode_mt(int mtid,char *linestr)
{
  char gene[128],mid[128],name[128];
  char strand;
  long exonstart=0,exonend=0,cdsstart=0,cdsend=0,intronstart=0,intronend=0,intergenicstart=0,intergenicend=0,csstart=0,csend=0,regstart=0,regend=0;
  int phase=-100,exonid,intronid,identity;
  char tmpdata[128];
  string chrname;
  int iscds;
  int flag=0;
  int chrid;
  
  ParseString(linestr,7,"\t",tmpdata);
  strand=tmpdata[0];
  
  ParseString(linestr,3,"\t",tmpdata);
  if(strcmp(tmpdata,"CDS")==0)
    iscds=1;
  else
    if(strcmp(tmpdata,"exon")==0)
       iscds=0;
    else
        if(strcmp(tmpdata,"intron")==0)
          iscds=-1;
        else
          if(strcmp(tmpdata,"intergenic")==0)
            iscds=-2;
          else
             if(strcmp(tmpdata,"cs")==0)
                iscds=-3;
             else
                if(strcmp(tmpdata,"encode")==0)
                   iscds=-4;
                else
                   if(strcmp(tmpdata,"motif")==0)
                      iscds=-5;
                   else
                      return;
  

  ParseString(linestr,2,"\"",gene);
  ParseString(linestr,4,"\"",mid);

  if(iscds>=-1)
  {
     if(StringFind(linestr,"gene_name",0)>0)
     {
        MultiParseString(linestr,2,"gene_name",2,"\"",name);
        appendidname_mt(mtid,mid,gene,name);
     }
     else
        appendidname_mt(mtid,mid,gene,gene);
  }
  
  ParseString(linestr,4,"\t",tmpdata);
  if(iscds==1)
    cdsstart=atol(tmpdata);
  if(iscds==0)
    exonstart=atol(tmpdata);
  if(iscds==-1)
    intronstart=atol(tmpdata);
  if(iscds==-2)
    intergenicstart=atol(tmpdata);
  if(iscds==-3)
    csstart=atol(tmpdata);
  if(iscds==-4||iscds==-5)
    regstart=atol(tmpdata);
    
  ParseString(linestr,5,"\t",tmpdata);
  if(iscds==1)
    cdsend=atol(tmpdata);
  if(iscds==0)
    exonend=atol(tmpdata);
  if(iscds==-1)
    intronend=atol(tmpdata);
  if(iscds==-2)
    intergenicend=atol(tmpdata);
  if(iscds==-3)
    csend=atol(tmpdata);
  if(iscds==-4||iscds==-5)
    regend=atol(tmpdata);

  if(iscds>=0)
  {
     ParseString(linestr,6,"\"",tmpdata);
     exonid=atoi(tmpdata);
  }
  if(iscds==-1)
  {
     ParseString(linestr,6,"\"",tmpdata);
     intronid=atoi(tmpdata);
  }
  if(iscds==-3||iscds==-4||iscds==-5)
  {
     ParseString(linestr,6,"\"",tmpdata);
     identity=atoi(tmpdata);
  }

  ParseString(linestr,8,"\t",tmpdata);
  if(tmpdata[0]=='.')
    phase=-100;
  else
    phase=atoi(tmpdata);
  
  ParseString(linestr,1,"\t",tmpdata);
  chrname = tmpdata;
  chrid=getchrindex(chrname);
  if(chrid<0)
    return;
  
//  printf("%d %d append %ld\n",mtid, chrid, globalhead_mt[mtid][chrid]);
  
  if(globalhead_mt[mtid][chrid]==NULL)
  {
      globalhead_mt[mtid][chrid] = (geneinfo*)malloc(sizeof(geneinfo));
      tailnode_mt[mtid][chrid]=globalhead_mt[mtid][chrid];
      setnullnode(tailnode_mt[mtid][chrid]);
      strcpy(tailnode_mt[mtid][chrid]->gene,gene);
      strcpy(tailnode_mt[mtid][chrid]->mid,mid);
      tailnode_mt[mtid][chrid]->strand=strand;
      tailnode_mt[mtid][chrid]->seqtype=iscds;
      if(iscds>=0)
      {
         tailnode_mt[mtid][chrid]->exonid=exonid;
         tailnode_mt[mtid][chrid]->exonstart=exonstart;
         tailnode_mt[mtid][chrid]->exonend=exonend;
         tailnode_mt[mtid][chrid]->cdsstart=cdsstart;
         tailnode_mt[mtid][chrid]->cdsend=cdsend;
         tailnode_mt[mtid][chrid]->phase=phase;
      }
      if(iscds==-1)
      {
         tailnode_mt[mtid][chrid]->exonid=intronid;
         tailnode_mt[mtid][chrid]->exonstart=intronstart;
         tailnode_mt[mtid][chrid]->exonend=intronend;
      }
      if(iscds==-2)
      {
         tailnode_mt[mtid][chrid]->exonstart=intergenicstart;
         tailnode_mt[mtid][chrid]->exonend=intergenicend;
      }
      if(iscds==-3)
      {
         tailnode_mt[mtid][chrid]->exonid=identity;
         tailnode_mt[mtid][chrid]->exonstart=csstart;
         tailnode_mt[mtid][chrid]->exonend=csend;
      }
      if(iscds==-4||iscds==-5)
      {
         tailnode_mt[mtid][chrid]->exonid=identity;
         tailnode_mt[mtid][chrid]->exonstart=regstart;
         tailnode_mt[mtid][chrid]->exonend=regend;
      }
      
      tailnode_mt[mtid][chrid]->pNext=NULL;
  }
  else
  {
      if(iscds>=0)
      {
         /*
         for(curnode_mt[mtid][chrid]=globalhead_mt[mtid][chrid];curnode_mt[mtid][chrid]!=NULL;curnode_mt[mtid][chrid]=curnode_mt[mtid][chrid]->pNext)
         {
            if(curnode_mt[mtid][chrid]->exonid==exonid&&strcmp(curnode_mt[mtid][chrid]->gene,gene)==0&&strcmp(curnode_mt[mtid][chrid]->mid,mid)==0)
            {
                flag=1;
                if(iscds==1)
                {
                  curnode_mt[mtid][chrid]->cdsstart=cdsstart;
                  curnode_mt[mtid][chrid]->cdsend=cdsend;
                  curnode_mt[mtid][chrid]->phase=phase;
                }
                if(iscds==0)
                {
                  curnode_mt[mtid][chrid]->exonstart=exonstart;
                  curnode_mt[mtid][chrid]->exonend=exonend;
                }
                break;
            }
         }
         */
         if(tailnode_mt[mtid][chrid]->exonid==exonid&&strcmp(tailnode_mt[mtid][chrid]->gene,gene)==0&&strcmp(tailnode_mt[mtid][chrid]->mid,mid)==0)
         {
//            if(tailnode_mt[mtid][chrid]->cdsstart==0||tailnode_mt[mtid][chrid]->exonstart==0)
//            {
                flag=1;
                tailnode_mt[mtid][chrid]->seqtype=iscds;
                if(iscds==1)
                {
                  tailnode_mt[mtid][chrid]->cdsstart=cdsstart;
                  tailnode_mt[mtid][chrid]->cdsend=cdsend;
                  tailnode_mt[mtid][chrid]->phase=phase;
                }
                if(iscds==0)
                {
                  tailnode_mt[mtid][chrid]->exonstart=exonstart;
                  tailnode_mt[mtid][chrid]->exonend=exonend;
                }
//            }
         }
         
         /*
         if(iscds==-1)
         {
            if(curnode_mt[mtid][chrid]->intronid==intronid&&strcmp(curnode_mt[mtid][chrid]->gene,gene)==0&&strcmp(curnode_mt[mtid][chrid]->mid,mid)==0)
            {
               flag=1;
               curnode_mt[mtid][chrid]->intronstart=intronstart;
               curnode_mt[mtid][chrid]->intronend=intronend;
               break;
            }
         }
         */
      }
      if(flag==0)
      {
        tailnode_mt[mtid][chrid]->pNext=(geneinfo*)malloc(sizeof(geneinfo));
        tailnode_mt[mtid][chrid]=tailnode_mt[mtid][chrid]->pNext;
        setnullnode(tailnode_mt[mtid][chrid]);
        strcpy(tailnode_mt[mtid][chrid]->gene,gene);
        strcpy(tailnode_mt[mtid][chrid]->mid,mid);
        tailnode_mt[mtid][chrid]->strand=strand;
        tailnode_mt[mtid][chrid]->seqtype=iscds;
        if(iscds>=0)
        {
           tailnode_mt[mtid][chrid]->exonid=exonid;
           tailnode_mt[mtid][chrid]->exonstart=exonstart;
           tailnode_mt[mtid][chrid]->exonend=exonend;
           tailnode_mt[mtid][chrid]->cdsstart=cdsstart;
           tailnode_mt[mtid][chrid]->cdsend=cdsend;
           tailnode_mt[mtid][chrid]->phase=phase;
        }
        if(iscds==-1)
        {
           tailnode_mt[mtid][chrid]->exonid=intronid;
           tailnode_mt[mtid][chrid]->exonstart=intronstart;
           tailnode_mt[mtid][chrid]->exonend=intronend;
        }
        if(iscds==-2)
        {
           tailnode_mt[mtid][chrid]->exonstart=intergenicstart;
           tailnode_mt[mtid][chrid]->exonend=intergenicend;
        }
        if(iscds==-3)
        {
           tailnode_mt[mtid][chrid]->exonid=identity;
           tailnode_mt[mtid][chrid]->exonstart=csstart;
           tailnode_mt[mtid][chrid]->exonend=csend;
        }
        if(iscds==-4||iscds==-5)
        {
           tailnode_mt[mtid][chrid]->exonid=identity;
           tailnode_mt[mtid][chrid]->exonstart=regstart;
           tailnode_mt[mtid][chrid]->exonend=regend;
        }
        tailnode_mt[mtid][chrid]->pNext=NULL;
      }
  }
  if((tailnode_mt[mtid][chrid]->exonend-tailnode_mt[mtid][chrid]->exonstart+1)>maxseglen[mtid])
      maxseglen[mtid] = tailnode_mt[mtid][chrid]->exonend-tailnode_mt[mtid][chrid]->exonstart+1;
  //printnode(tailnode_mt[mtid][chrid]);
}


geneinfo * findnode_mt(int mtid,int chrid,char *gene,char *mid,int exonid)
{
  geneinfo *localnode=NULL;
  for(localnode=globalhead_mt[mtid][chrid];localnode!=NULL;localnode=localnode->pNext)
  {
    if(localnode->exonid==exonid&&localnode->cdsstart>0&&strcmp(localnode->gene,gene)==0&&strcmp(localnode->mid,mid)==0)
    {
      return localnode;
    }
  }
  return NULL;
}


void * findsnpincodon2(void * mt_id)
{
   int mtid;
   int retid=0;
   int chrid;
   int lastchrid=-1;
   long snppos;
   char snpbase1[32];
   char snpbase2[32];
   
   char leftid[128],rightid[128];
   
   char codon1[4],codon2[4];
   int snpincodonpos;
   long codonpos1,codonpos2,codonpos3;
   char filename[256];
   int findhit;
   geneinfo *nextnode;
   char strand[2],bases[2];
   char aa1[3],aa2[3];
   int snpcount;

   mtid = *(int *)mt_id;
   
   for(snpcount=0;snpcount<SNPSIZE;snpcount++)
   {
      chrid = snpdata[mtid][snpcount].chrid;
      if(chrid==-1)
      {
         retstatus[mtid][snpcount]=-1;
         break;
      }
      /*
      if(chrid!=lastchrid)
      {
          startnode=globalhead_mt[mtid][chrid];
          lastchrid=chrid;
      }
      */
      if(startnode[mtid][chrid]==NULL)
         startnode[mtid][chrid]=globalhead_mt[mtid][chrid];
         
      snppos = snpdata[mtid][snpcount].snppos;
//      memcpy(snpbase1,snpdata[mtid][snpcount].snpbase1,1);
//      memcpy(snpbase2,snpdata[mtid][snpcount].snpbase2,1);
      
      strand[1]=0;
      findhit=0;
      
      /*
      if(snppos==3857490)
      {
          fprintf(stderr, "3857490 startnode: %d  %d  %d \n",startnode->exonstart,startnode->exonend,startnode->seqtype);
      }
      */
      //sprintf(filename,"/lustre/user/xiehb/genome/Sscrofa10.2/chr%s.fa",chr[chrid]);
      for(curnode_mt[mtid][chrid]=startnode[mtid][chrid];curnode_mt[mtid][chrid]!=NULL;curnode_mt[mtid][chrid]=curnode_mt[mtid][chrid]->pNext)
      {
         strcpy(snpbase1,snpdata[mtid][snpcount].snpbase1);
         strcpy(snpbase2,snpdata[mtid][snpcount].snpbase2);
         //fprintf(stderr,"inside thread %d\n",mtid); 
         if(curnode_mt[mtid][chrid]->exonend<snppos)
         {
            startnode[mtid][chrid]=curnode_mt[mtid][chrid];
            continue;
         }
         if(curnode_mt[mtid][chrid]->exonend-snppos>maxseglen[mtid])
            break;
         
         if(curnode_mt[mtid][chrid]->exonend>=snppos&&curnode_mt[mtid][chrid]->exonstart<=snppos)
         {
             //printf("%s %s %s %ld %d %ld %ld %ld %ld\n",chr[chrid],curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,snppos,curnode_mt[mtid][chrid]->exonid,curnode_mt[mtid][chrid]->exonstart,curnode_mt[mtid][chrid]->exonend,curnode_mt[mtid][chrid]->cdsstart,curnode_mt[mtid][chrid]->cdsend);
             findhit=1;
             if(curnode_mt[mtid][chrid]->seqtype>=0)
             {
                  if(curnode_mt[mtid][chrid]->cdsend>=snppos&&curnode_mt[mtid][chrid]->cdsstart<=snppos)
                  {
                      snpincodonpos=0;
                      strcpy(codon1,"");
                      strcpy(codon2,"");
                  
                      if(curnode_mt[mtid][chrid]->strand=='+')
                      {
                         snpincodonpos = ((snppos-curnode_mt[mtid][chrid]->cdsstart+1)%3+curnode_mt[mtid][chrid]->phase)%3;
                         strand[0]='+';
                         if(snpincodonpos==0)
                            snpincodonpos=3;
                         switch(snpincodonpos)
                         {
                            case 1:
                                codonpos1=snppos;
                                codonpos2=codonpos1+1;
                                if(codonpos2>curnode_mt[mtid][chrid]->cdsend)
                                {
                                  //
                                  nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid+1);
                                  if(nextnode==NULL) pthread_exit((void *)0);
                                  codonpos2=nextnode->cdsstart;
                                  codonpos3=codonpos2+1;
                                }
                                else
                                {
                                  codonpos3=codonpos2+1;
                                  if(codonpos3>curnode_mt[mtid][chrid]->cdsend)
                                  {
                                    //
                                    nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid+1);
                                    if(nextnode==NULL) pthread_exit((void *)0);
                                    codonpos3=nextnode->cdsstart;
                                  }
                                }
                                break;
                            case 2:
                                codonpos2=snppos;
                                codonpos1=codonpos2-1;
                                if(codonpos1<curnode_mt[mtid][chrid]->cdsstart)
                                {
                                  //
                                  nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid-1);
                                  if(nextnode==NULL) pthread_exit((void *)0);
                                  codonpos1=nextnode->cdsend;
                                  codonpos3=codonpos2+1;
                                }
                                else
                                {
                                  codonpos3=codonpos2+1;
                                  if(codonpos3>curnode_mt[mtid][chrid]->cdsend)
                                  {
                                    //
                                    nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid+1);
                                    if(nextnode==NULL) pthread_exit((void *)0);
                                    codonpos3=nextnode->cdsstart;
                                  }
                                }
                                break;
                            case 3:
                                codonpos3=snppos;
                                codonpos2=codonpos3-1;
                                if(codonpos2<curnode_mt[mtid][chrid]->cdsstart)
                                {
                                  //
                                  nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid-1);
                                  if(nextnode==NULL) pthread_exit((void *)0);
                                  codonpos2=nextnode->cdsend;
                                  codonpos1=codonpos2-1;
                                }
                                else
                                {
                                  codonpos1=codonpos2-1;
                                  if(codonpos1<curnode_mt[mtid][chrid]->cdsstart)
                                  {
                                    //
                                    nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid-1);
                                    if(nextnode==NULL) pthread_exit((void *)0);
                                    codonpos1=nextnode->cdsend;
                                  }
                                }
                                break;
                            default:
                                break;
                         }
                      }
                      else
                      {
                         strand[0]='-';
                         snpincodonpos=((curnode_mt[mtid][chrid]->cdsend-snppos+1)%3+curnode_mt[mtid][chrid]->phase)%3; 
                         if(snpincodonpos==0)
                            snpincodonpos=3;
                  
                         switch(snpincodonpos)
                         {
                            case 1:
                                codonpos1=snppos;
                                codonpos2=codonpos1-1;
                                if(codonpos2<curnode_mt[mtid][chrid]->cdsstart)
                                {
                                  //
                                  nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid+1);
                                  if(nextnode==NULL) pthread_exit((void *)0);
                                  codonpos2=nextnode->cdsend;
                                  codonpos3=codonpos2-1;
                                }
                                else
                                {
                                  codonpos3=codonpos2-1;
                                  if(codonpos3<curnode_mt[mtid][chrid]->cdsstart)
                                  {
                                    //
                                    nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid+1);
                                    if(nextnode==NULL) pthread_exit((void *)0);
                                    codonpos3=nextnode->cdsend;
                                  }
                                }
                                break;
                            case 2:
                                codonpos2=snppos;
                                codonpos1=codonpos2+1;
                                if(codonpos1>curnode_mt[mtid][chrid]->cdsend)
                                {
                                   //
                                   nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid-1);
                                   if(nextnode==NULL) pthread_exit((void *)0);
                                   codonpos1=nextnode->cdsstart;
                                   codonpos3=codonpos2-1;
                                }
                                else
                                {
                                  codonpos3=codonpos2-1;
                                  if(codonpos3<curnode_mt[mtid][chrid]->cdsstart)
                                  {
                                    //
                                    nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid+1);
                                    if(nextnode==NULL) pthread_exit((void *)0);
                                    codonpos3=nextnode->cdsend;
                  
                                  }
                                }
                                break;
                            case 3:
                                codonpos3=snppos;
                                codonpos2=codonpos3+1;
                                if(codonpos2>curnode_mt[mtid][chrid]->cdsend)
                                {
                                  //
                                  nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid-1);
                                  if(nextnode==NULL) pthread_exit((void *)0);
                                  codonpos2=nextnode->cdsend;
                                  codonpos1=codonpos2+1;
                                }
                                else
                                {
                                  codonpos1=codonpos2+1;
                                  if(codonpos1>curnode_mt[mtid][chrid]->cdsend)
                                  {
                                    //
                                    nextnode=findnode_mt(mtid,chrid,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid-1);
                                    if(nextnode==NULL) pthread_exit((void *)0);
                                    codonpos1=nextnode->cdsstart;
                                 }
                                }
                                break;
                            default:
                                break;
                         }
                      }
                      
                      /*
                      ReadFasta(filename,strand,codonpos1,1,70,1,bases,&readstart,&readend,&errorflag);
                      codon1[0]=bases[0];
                      ReadFasta(filename,strand,codonpos2,1,70,1,bases,&readstart,&readend,&errorflag);
                      codon1[1]=bases[0];
                      ReadFasta(filename,strand,codonpos3,1,70,1,bases,&readstart,&readend,&errorflag);
                      codon1[2]=bases[0];
                      codon1[3]=0;
                      */
                      if(strand[0]=='+')
                      {
                         codon1[0]=chrseq[chrid][codonpos1-1];
                         codon1[1]=chrseq[chrid][codonpos2-1];
                         codon1[2]=chrseq[chrid][codonpos3-1];
                      }
                      else
                      {
                         codon1[0]=chrseq[chrid][codonpos3-1];
                         codon1[1]=chrseq[chrid][codonpos2-1];
                         codon1[2]=chrseq[chrid][codonpos1-1];
                      }
                      codon1[3]=0;
                      SeqUpper(codon1);
                      
                      if(strand[0]=='-')
                      {
                        ReverseSequence(codon1);
                        ComplementSequence(codon1);
                      }
                      
                      if(strand[0]=='-')
                      {
                        ReverseSequence(snpbase1);
                        ComplementSequence(snpbase1);
                      }
                      codon1[snpincodonpos-1]=snpbase1[0];
                      aa1[0]=f_codon(codon1);
                      aa1[1]=0;
                      
                      strcpy(codon2,codon1);
                      if(strand[0]=='-')
                      {
                        ReverseSequence(snpbase2);
                        ComplementSequence(snpbase2);
                      }
                      codon2[snpincodonpos-1]=snpbase2[0];
                      aa2[0]=f_codon(codon2);
                      aa2[1]=0;
                  
                  }
                  else
                  {
                      snpincodonpos=0;
                      strcpy(codon1,snpbase1);
                      strcpy(codon2,snpbase2); 
                      aa1[0]='-';
                      aa1[1]=0;
                      aa2[0]='-';
                      aa2[1]=0;
                  }
                  sprintf(threadsdata[mtid][retid],"%s\t%ld\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d\t%c\t%ld\t%ld\t%ld\t%ld\t%d\n",chr[chrid].c_str(),snppos,codon1,codon2,aa1,aa2,snpincodonpos,findidname_mt(mtid,curnode_mt[mtid][chrid]->mid,0),curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid,curnode_mt[mtid][chrid]->strand,curnode_mt[mtid][chrid]->exonstart,curnode_mt[mtid][chrid]->exonend,curnode_mt[mtid][chrid]->cdsstart,curnode_mt[mtid][chrid]->cdsend,curnode_mt[mtid][chrid]->phase);
                  retid++;
             }
             if(curnode_mt[mtid][chrid]->seqtype==-1)
             {
                  snpincodonpos=-1;
                  strcpy(codon1,snpbase1);
                  strcpy(codon2,snpbase2); 
                  aa1[0]='-';
                  aa1[1]=0;
                  aa2[0]='-';
                  aa2[1]=0;
                  sprintf(threadsdata[mtid][retid],"%s\t%ld\t%s\t%s\t-\t-\t%d\t%s\t%s\t%d\t%c\t%ld\t%ld\t%ld\t%ld\t-\n",chr[chrid].c_str(),snppos,codon1,codon2,snpincodonpos,findidname_mt(mtid,curnode_mt[mtid][chrid]->mid,0),curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid,curnode_mt[mtid][chrid]->strand,curnode_mt[mtid][chrid]->exonstart,curnode_mt[mtid][chrid]->exonend,snppos-curnode_mt[mtid][chrid]->exonstart+1,curnode_mt[mtid][chrid]->exonend-snppos+1);
                  retid++;
             }
             if(curnode_mt[mtid][chrid]->seqtype==-2)
             {
                  snpincodonpos=-10;
                  strcpy(codon1,snpbase1);
                  strcpy(codon2,snpbase2); 
                  aa1[0]='-';
                  aa1[1]=0;
                  aa2[0]='-';
                  aa2[1]=0;
                  strcpy(leftid,findidname_mt(mtid,curnode_mt[mtid][chrid]->gene,0));
                  strcpy(rightid,findidname_mt(mtid,curnode_mt[mtid][chrid]->mid,1));
                  sprintf(threadsdata[mtid][retid],"%s\t%ld\t%s\t%s\t-\t-\t%d\t%s\t%s\tX\t%c\t%ld\t%ld\t-\t-\t-\n",chr[chrid].c_str(),snppos,codon1,codon2,snpincodonpos,leftid,rightid,curnode_mt[mtid][chrid]->strand,curnode_mt[mtid][chrid]->exonstart,curnode_mt[mtid][chrid]->exonend);
                  retid++;
             }
             if(curnode_mt[mtid][chrid]->seqtype==-3)
             {
             	  snpincodonpos=-20;
                  strcpy(codon1,snpbase1);
                  strcpy(codon2,snpbase2); 
                  aa1[0]='-';
                  aa1[1]=0;
                  aa2[0]='-';
                  aa2[1]=0;
                  sprintf(threadsdata[mtid][retid],"%s\t%ld\t%s\t%s\t-\t-\t%d\t%s\t%s\t%d\t%c\t%ld\t%ld\t-\t-\t-\n",chr[chrid].c_str(),snppos,codon1,codon2,snpincodonpos,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid,curnode_mt[mtid][chrid]->strand,curnode_mt[mtid][chrid]->exonstart,curnode_mt[mtid][chrid]->exonend);
                  retid++;
             }
             if(curnode_mt[mtid][chrid]->seqtype==-4)
             {
                  snpincodonpos=-30;
                  strcpy(codon1,snpbase1);
                  strcpy(codon2,snpbase2); 
                  aa1[0]='-';
                  aa1[1]=0;
                  aa2[0]='-';
                  aa2[1]=0;
                  sprintf(threadsdata[mtid][retid],"%s\t%ld\t%s\t%s\t-\t-\t%d\t%s\t%s\t%d\t%c\t%ld\t%ld\t-\t-\t-\n",chr[chrid].c_str(),snppos,codon1,codon2,snpincodonpos,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid,curnode_mt[mtid][chrid]->strand,curnode_mt[mtid][chrid]->exonstart,curnode_mt[mtid][chrid]->exonend);
                  retid++;
             }
             if(curnode_mt[mtid][chrid]->seqtype==-5)
             {
                  snpincodonpos=-40;
                  strcpy(codon1,snpbase1);
                  strcpy(codon2,snpbase2); 
                  aa1[0]='-';
                  aa1[1]=0;
                  aa2[0]='-';
                  aa2[1]=0;
                  sprintf(threadsdata[mtid][retid],"%s\t%ld\t%s\t%s\t-\t-\t%d\t%s\t%s\t%d\t%c\t%ld\t%ld\t-\t-\t-\n",chr[chrid].c_str(),snppos,codon1,codon2,snpincodonpos,curnode_mt[mtid][chrid]->gene,curnode_mt[mtid][chrid]->mid,curnode_mt[mtid][chrid]->exonid,curnode_mt[mtid][chrid]->strand,curnode_mt[mtid][chrid]->exonstart,curnode_mt[mtid][chrid]->exonend);
                  retid++;
             }
         }
      }
     
      if(findhit==0)
      {
        strcpy(snpbase1,snpdata[mtid][snpcount].snpbase1);
        strcpy(snpbase2,snpdata[mtid][snpcount].snpbase2);
        sprintf(threadsdata[mtid][retid],"%s\t%ld\t%s\t%s\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n",chr[chrid].c_str(),snppos,snpbase1,snpbase2);
        retid++;
      }
      retstatus[mtid][snpcount]=1;
   }
   pthread_exit((void *)1);
}

void deletenodes_mt(int mtid)
{
  int chrid;
  for(chrid=0;chrid<chr.size();chrid++)
  {
    for(curnode_mt[mtid][chrid]=globalhead_mt[mtid][chrid];curnode_mt[mtid][chrid]!=NULL;)
    {
       tailnode_mt[mtid][chrid]=curnode_mt[mtid][chrid];
       curnode_mt[mtid][chrid]=curnode_mt[mtid][chrid]->pNext;
       free(tailnode_mt[mtid][chrid]);
    }
  }
  free(globalhead_mt[mtid]);
}


int main(int argc, char *argv[])
{
  char tmpstr[128];
  FILE *snpfile,*gtffile;
  int chrid;
  char *linedata;
  long snppos;
  char snpbase[10],snpbase1[10],snpbase2[10];
  int threadcount,threadseqid,retid;
  pthread_t *threadid;
  void *threadret;
  int snpcount=-1;
  geneinfo **sortlist,*pnode,*tmpnode;
  long maxpos,steppos;
  
  if(argc!=6)
  {
     printf("Usage: %s threads(<=8) gtf_file snp_file chrlistfile genome_directory\n",argv[0]);
     return 1;
  }
  else
  {
     snpfile = fopen(argv[3],"rb");
     if(snpfile==NULL)
     {
       return 0;
     }
     gtffile = fopen(argv[2],"rb");
     if(gtffile==NULL)
     {
       fclose(snpfile);
       return 0;
     }
     threadcount = atoi(argv[1]);
     if(threadcount>8||threadcount<1)
     {
       printf("Usage: %s threads(<=8) gtf_file snp_file chrlistfile genome_directory\n",argv[0]);
       return 1;
     }
     chr = loadchrinfo(argv[4]);
  }
  
  resetlistdata();
  
  for(threadseqid=0;threadseqid<threadcount;threadseqid++)
  {
     maxseglen[threadseqid]=0;
     globalhead_mt[threadseqid] = (geneinfo **)malloc(chr.size()*sizeof(geneinfo *));
     curnode_mt[threadseqid] = (geneinfo **)malloc(chr.size()*sizeof(geneinfo *));
     tailnode_mt[threadseqid] = (geneinfo **)malloc(chr.size()*sizeof(geneinfo *));
     for(chrid=0;chrid<chr.size();chrid++)
     {
       globalhead_mt[threadseqid][chrid]=NULL;
       curnode_mt[threadseqid][chrid]=NULL;
       tailnode_mt[threadseqid][chrid]=NULL;
       startnode[threadseqid][chrid]=NULL;
     }
  } 
  linedata = (char*)malloc(1);
  fprintf(stderr,"Read exon info ... ");
  fflush(stderr);
  /*
  while(ReadLine(gtffile,&linedata)>0)
  {
    for(threadseqid=0;threadseqid<threadcount;threadseqid++)
       appendnode_mt(threadseqid,linedata);
  }
  */ 
  threadseqid=0;
  while(ReadLine(gtffile,&linedata)>0)
  {
     appendnode_mt(threadseqid,linedata);
  }
  
  fclose(gtffile);
  fprintf(stderr,"Done.\n");
  fflush(stderr);
  
  
  //sort chains with the exonend variable
  fprintf(stderr,"Sort exon ... ");

  threadseqid=0;
  for(chrid=0;chrid<chr.size();chrid++)
  {
      maxpos = 0;
      for(pnode = globalhead_mt[threadseqid][chrid]; pnode!=NULL; pnode=pnode->pNext )
      {
         if(pnode->exonend > maxpos)
            maxpos = pnode->exonend;
      }
      
      sortlist = (geneinfo **)malloc(maxpos*sizeof(geneinfo *));
      memset(sortlist,0,maxpos*sizeof(geneinfo *));
      for(pnode = globalhead_mt[threadseqid][chrid]; pnode!=NULL; )
      {
         steppos = pnode->exonend;

         if(sortlist[steppos-1]==NULL)
         {
            sortlist[steppos-1] = pnode;
            pnode=pnode->pNext;
            sortlist[steppos-1]->pNext=NULL;
         }
         else
         {
            for(tmpnode=sortlist[pnode->exonend-1];tmpnode!=NULL;tmpnode=tmpnode->pNext)
            {
               if(tmpnode->pNext==NULL)
               {
                  tmpnode->pNext=pnode;
                  tmpnode=pnode;
                  pnode=pnode->pNext;
                  tmpnode->pNext=NULL;
                  break;
               }
            }
         }
      }
      
      globalhead_mt[threadseqid][chrid]=NULL;
      for(steppos=0;steppos<maxpos;steppos++)
      {
         if(sortlist[steppos]!=NULL)
         {
            if(globalhead_mt[threadseqid][chrid]==NULL)
            {
               globalhead_mt[threadseqid][chrid] = sortlist[steppos];
               curnode_mt[threadseqid][chrid]=globalhead_mt[threadseqid][chrid];
               curnode_mt[threadseqid][chrid]->pNext=NULL;
               for(pnode=sortlist[steppos]->pNext;pnode!=NULL;pnode=pnode->pNext)
               {
                  curnode_mt[threadseqid][chrid]->pNext=pnode;
                  curnode_mt[threadseqid][chrid]=curnode_mt[threadseqid][chrid]->pNext;
               }
            }
            else
            {
               curnode_mt[threadseqid][chrid]->pNext=sortlist[steppos];
               curnode_mt[threadseqid][chrid]=curnode_mt[threadseqid][chrid]->pNext;
               
               for(pnode=sortlist[steppos]->pNext;pnode!=NULL;pnode=pnode->pNext)
               {
                  curnode_mt[threadseqid][chrid]->pNext=pnode;
                  curnode_mt[threadseqid][chrid]=curnode_mt[threadseqid][chrid]->pNext;
               }
            }
         }
      }
      free(sortlist);
  }     
  
  fprintf(stderr,"Done.\n");
  
  fprintf(stderr,"Duplicate sorted_gtf chains ... ");
  for(threadseqid=1;threadseqid<threadcount;threadseqid++)
  {
     for(chrid=0;chrid<chr.size();chrid++)
     {
         for(pnode = globalhead_mt[0][chrid]; pnode!=NULL; pnode=pnode->pNext)
         {
            if(globalhead_mt[threadseqid][chrid]==NULL)
            {
               globalhead_mt[threadseqid][chrid] = (geneinfo*)malloc(sizeof(geneinfo));
               curnode_mt[threadseqid][chrid]=globalhead_mt[threadseqid][chrid];
            }
            else
            {
               curnode_mt[threadseqid][chrid]->pNext=(geneinfo*)malloc(sizeof(geneinfo));
               curnode_mt[threadseqid][chrid] = curnode_mt[threadseqid][chrid]->pNext;
            }
            strcpy(curnode_mt[threadseqid][chrid]->gene,pnode->gene);
            strcpy(curnode_mt[threadseqid][chrid]->mid,pnode->mid);
            curnode_mt[threadseqid][chrid]->strand=pnode->strand;
            curnode_mt[threadseqid][chrid]->exonid=pnode->exonid;
            curnode_mt[threadseqid][chrid]->exonstart=pnode->exonstart;
            curnode_mt[threadseqid][chrid]->exonend=pnode->exonend;
            curnode_mt[threadseqid][chrid]->cdsstart=pnode->cdsstart;
            curnode_mt[threadseqid][chrid]->cdsend=pnode->cdsend;
            curnode_mt[threadseqid][chrid]->phase=pnode->phase;
            curnode_mt[threadseqid][chrid]->seqtype=pnode->seqtype;
            curnode_mt[threadseqid][chrid]->pNext=NULL;
         }
     }
     
     maxseglen[threadseqid]=maxseglen[0]; 
     
     for(pidname_mt[0]=idnamehead_mt[0];pidname_mt[0]!=NULL;pidname_mt[0]=pidname_mt[0]->pNext)
     {
         appendidname_mt(threadseqid,pidname_mt[0]->id,pidname_mt[0]->gene,pidname_mt[0]->name);
     }
  }  
  fprintf(stderr,"Done.\n");

  threadid = (pthread_t *)malloc(sizeof(pthread_t)*threadcount);
  threadseqid=-1;
  
  resetsnpdata();
  
  fprintf(stderr,"Load Genome ... ");
  load_genome_seq(argv[5]);
  fprintf(stderr,"Done.\n");
  
  //printf("chr\tsnppos\tAncestralCodon\tCodon\tAncestralAA\tAA\tSnpInCodonPos\tgene\tmid\texonid\tstrand\tExonStart\tExonEnd\tcdsStart\tcdsEnd\tphase\n");
  while(ReadLine(snpfile,&linedata)>0)
  {
    ParseString(linedata,1,"\t",tmpstr);
    chrid=getchrindex(tmpstr);
    ParseString(linedata,2,"\t",tmpstr);
    snppos=atol(tmpstr);
    ParseString(linedata,4,"\t",snpbase);
    ParseString(linedata,3,"\t",snpbase1);
    ParseString(linedata,4,"\t",snpbase2);

    if(chrid<0)   
      continue;

    snpcount++;
    if(snpcount==0)
    {
      threadseqid++;
      resetthreaddata(threadseqid);
    }
    
    if(snpcount<SNPSIZE)
    {
      snpdata[threadseqid][snpcount].chrid = chrid;
      snpdata[threadseqid][snpcount].snppos = snppos;
      strcpy(snpdata[threadseqid][snpcount].snpbase1,snpbase1);
      strcpy(snpdata[threadseqid][snpcount].snpbase2,snpbase2);
    }
    if(snpcount<SNPSIZE-1)
      continue;
    
    if(threadseqid<threadcount)
    {
       pthread_create(&threadid[threadseqid], NULL, findsnpincodon2, (void *)&pthreadid[threadseqid]);
    }
    
    if(threadseqid==threadcount-1)
    {
       for(threadseqid=0;threadseqid<threadcount;threadseqid++)
       {
          pthread_join(threadid[threadseqid],&threadret);
          for(retid=0;retid<SNPRETSIZE;retid++)
          {
             if(threadsdata[threadseqid][retid][0]!=0)
                 printf(threadsdata[threadseqid][retid]);
             else
                 break;
          }
       }
       resetsnpdata();
       threadseqid=-1;
    }
    snpcount=-1;
  }
  
  if(snpcount>=0)
  {
     pthread_create(&threadid[threadseqid], NULL, findsnpincodon2, (void *)&pthreadid[threadseqid]);
  }

  for(;threadseqid>=0;threadseqid--)
  {
     pthread_join(threadid[threadseqid],&threadret);
     for(retid=0;retid<SNPRETSIZE;retid++)
     {
        if(threadsdata[threadseqid][retid][0]!=0)
           printf(threadsdata[threadseqid][retid]);
        else
           break;
     }
  }

  unload_genome_seq();
  
  free(threadid);
  fclose(snpfile);
  free(linedata);
  
  for(threadseqid=0;threadseqid<threadcount;threadseqid++)
  {
     deletenodes_mt(threadseqid);
     deleteidname_mt(threadseqid);
  }
  
  return 1;
}

