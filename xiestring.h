#ifndef _xiestring

#define _xiestring

#ifdef __cplusplus
extern "C" {
#endif

void StringRight(char * str,long rlen,char *retstr);
void StringLeft(char * str,long llen,char *retstr);
void StringMid(char *str,long start,long length,char *retstr);
long StringFind(char *str,char *substr,long start);
void ParseString(char *srcstr,int column,char *separator,char *retstr);
void MultiParseString(char *srcstr,int column1,char *separator1, int column2, char *separator2, char *retstr);
long StringColumns(char *srcstr,char *separator);

#ifdef __cplusplus
}
#endif

#endif

