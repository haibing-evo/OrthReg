#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "xiestring.h"

#ifdef _MAC
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

void StringRight(char * str,long rlen,char *retstr)
{
	long linelen;
	
	if(rlen==0)
	{
		 retstr[0] = '\0';
		 return;
  }
  
  linelen = strlen(str);
	if(linelen<rlen)
		rlen = linelen;
		
	memmove(retstr,str+linelen-rlen,rlen);
	retstr[rlen] = '\0';
}

void StringLeft(char * str,long llen,char *retstr)
{
	long linelen;
	
	if(llen==0)
	{
		 retstr[0] = '\0';
		 return;
  }
  
  linelen = strlen(str);
  if(linelen<llen)
	   llen = linelen;

	memmove(retstr,str,llen);
	retstr[llen] = '\0';
}

void StringMid(char *str,long start,long length,char *retstr)
{
	long linelen;
	
	if(length==0)
	{
		 retstr[0]='\0';
		 return;
	}
	
	linelen = strlen(str);
	if((linelen-start)<length)
		length = linelen-start;
	
	memmove(retstr,str+start,length);
	retstr[length] = '\0';
}

long StringFind(char *str,char *substr,long start)
{
	char *pDest;
	long tmplen=strlen(str)-start;
	if(tmplen<=0) return -10;

	pDest = strstr(str+start,substr);
	if(pDest!=NULL) 
		return pDest-str;
	else
	    return -1;
}

void ParseString(char *srcstr,int column,char *separator,char *retstr)
{
	long srclen;
	long position,position1,position2;
	
	srclen = strlen(srcstr);
	if(srclen==0||column<=0||separator[0]=='\0')
	{
		retstr[0]='\0';
		return;
	}
	
	if(column==1)
	{
		position = StringFind(srcstr,separator,0);
		if(position>=0)
			 StringLeft(srcstr,position,retstr);
		else
		{
			 memmove(retstr,srcstr,strlen(srcstr));
			 retstr[strlen(srcstr)]=0;
			 return;
		}
	}
	else
	{
		position1 = 0;
		position2 = 0;
		position1 = StringFind(srcstr,separator,0);
		if(position1>=0) 
		   position1+=strlen(separator);
		column--;
		while(column>1&&position1>=0)
		{
			 column--;
			 position1 = StringFind(srcstr,separator,position1);
			 if(position1>=0)
		      position1+=strlen(separator);
		}
		if(position1==-1)
		{
			 retstr[0]='\0';
			 return;
		}
		position2 = StringFind(srcstr,separator,position1);
		if(position2==-1)
		{
			 StringRight(srcstr,srclen-position1,retstr);
			 return;
		}
		StringMid(srcstr,position1,position2-position1,retstr);
	}
	return;
}

long StringColumns(char *srcstr,char *separator)
{
	long position=0;
	long colcount=1;

	if(srcstr[0]==0)
		 return 0;
		 
	position = StringFind(srcstr,separator,0);
	while(position>=0)
	{
		colcount++;
		position = StringFind(srcstr,separator,position+1);
	}
  return colcount;
}

void MultiParseString(char *srcstr,int column1,char *separator1, int column2, char *separator2, char *retstr)
{
  char *tmpdata;
  tmpdata = (char*)malloc(strlen(srcstr)+1);
  ParseString(srcstr,column1,separator1,tmpdata);
  ParseString(tmpdata,column2,separator2,retstr);
  free(tmpdata);
}
