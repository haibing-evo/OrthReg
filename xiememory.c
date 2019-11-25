#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "xiememory.h"

#ifdef _MAC
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

long ReadMemLine(char *pMem,char **linestr,long curpos)
{
	int pos1,pos2;
	pos1 = curpos;
	*linestr = (char*)realloc(*linestr,1);
	(*linestr)[0]=0;
	
	if(pos1>=0)
	if(pMem[pos1]=='\0')
		 return pos1;
		 
	pos1 = curpos+1;
	pos2 = curpos+1;

	while(pMem[pos2]!='\n'&&pMem[pos2]!='\0')
	{
	   pos2++;
	}
	*linestr = (char*)realloc(*linestr,pos2-pos1+2);
	memcpy(*linestr,pMem+pos1,pos2-pos1+1);
	(*linestr)[pos2-pos1+1]=0;

	return pos2;
}

