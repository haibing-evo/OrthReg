#ifndef _xiefile

#define _xiefile
#include <stdio.h>
#include "zlib.h"
#include "zconf.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gzfile
{
	gzFile zHandle;
	char readbuf[16385];
	long linepos;
	long buffersize;
} gzline;

long ReadLine(FILE * file,char **linestr);
int IsEndOfFile(FILE *file);
void gz_open_linemode(char *filename,gzline *gz_link);
void gz_close_linemode(gzline *gz_link);
long gz_readline(gzline *gz_link,char **retstr);

#ifdef __cplusplus
}
#endif

#endif


