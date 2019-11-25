/*******************************************************/
/*      Interface for File Operations: Read and Write  */
/*      Date   : 2006.1.22                             */
/*      Author : Haibing Xie                           */
/*      Address: Kunming Insititue of Zoology,CAS.     */
/*      Email  : xiehb@mail.kiz.ac.cn                  */
/*******************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "xiefile.h"

#ifdef _MAC
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

//Read a line from file terminated with \r\n or \n
long ReadLine(FILE *file,char **linestr)
{
    int tmppos=0,tmplen=0;
    char *tmpstr,tmpread[1025];
    int finished = 0;
    
    if (IsEndOfFile(file))
        return -1;
  
    tmpstr = (char*)malloc(1); 
    tmpstr[0] = '\0';
    *linestr = (char*)realloc(*linestr,1);
    (*linestr)[0] = '\0';

    while(!finished&&!IsEndOfFile(file))
    {
        tmpread[0] = '\0';
        tmplen = fread(tmpread,sizeof(char),1024,file);
        tmpread[tmplen] = '\0';

        for(tmppos=0;tmppos<tmplen;tmppos++)
        {
          if((tmpread[tmppos]=='\r'&&tmpread[tmppos+1]=='\n')||tmpread[tmppos]=='\n')
          { 
              if(tmpread[tmppos]=='\r'&&tmpread[tmppos+1]=='\n')
                 fseeko(file,tmppos-tmplen+2,SEEK_CUR);
              else
                 if(tmpread[tmppos]=='\n')
                   fseeko(file,tmppos-tmplen+1,SEEK_CUR);
    
              finished = 1;
              break;
           }
        }
        tmpread[tmppos] = '\0';
        
        tmplen = strlen(*linestr);
        tmpstr = (char*)realloc(tmpstr,tmplen+1);
        memcpy(tmpstr,*linestr,tmplen);
        tmpstr[tmplen] = '\0';
    
        *linestr = (char*)realloc(*linestr,tmplen+tmppos+1);
        memcpy(*linestr,tmpstr,tmplen);
        memcpy(*linestr+tmplen,tmpread,tmppos);
        (*linestr)[tmplen+tmppos] = '\0';
    }

    free(tmpstr);
    return tmplen+tmppos;
}

int IsEndOfFile(FILE *file)
{
    if(fgetc(file)==EOF)
        return 1;
     else
     {
         fseeko(file,-1,SEEK_CUR);
         return 0;
    }
}

void gz_open_linemode(char *filename,gzline *gz_link)
{
    gz_link->zHandle = gzopen(filename, "rb");
    memset(gz_link->readbuf,0,sizeof(gz_link->readbuf));
    gz_link->linepos=-1;
    gz_link->buffersize=-1;
}
void gz_close_linemode(gzline *gz_link)
{
    gzclose(gz_link->zHandle);
}

long gz_readline(gzline *gz_link,char **retstr)
{
    char *linestr,*savedata;
    int datasize;
   
    linestr = NULL;
    savedata = NULL;
    
    if(gz_link->buffersize<0)
    {
       memset(gz_link->readbuf,0,sizeof(gz_link->readbuf));
       gz_link->buffersize = gzread(gz_link->zHandle, gz_link->readbuf, sizeof(gz_link->readbuf)-1);
       gz_link->linepos=-1;
       if (gz_link->buffersize < 0) return -1;
       if (gz_link->buffersize == 0) return -100;
    }
    
    while(gz_link->buffersize>0)
    {
      if(gz_link->linepos==gz_link->buffersize)
      {
         memset(gz_link->readbuf,0,sizeof(gz_link->readbuf));
         gz_link->buffersize = gzread(gz_link->zHandle, gz_link->readbuf, sizeof(gz_link->readbuf)-1);
         gz_link->linepos=-1;
         if (gz_link->buffersize < 0) return -1;
         if (gz_link->buffersize == 0) return -100;
      }
      
      gz_link->linepos = ReadMemLine(gz_link->readbuf,&linestr,gz_link->linepos);
      if(gz_link->linepos<gz_link->buffersize)
      {
         if(savedata!=NULL)
         {
             *retstr = (char*)realloc(*retstr,strlen(savedata)+strlen(linestr)+1);
             memcpy(*retstr,savedata,strlen(savedata));
             memcpy((*retstr)+strlen(savedata),linestr,strlen(linestr)+1);
         }
         else
         {
             *retstr = (char*)realloc(*retstr,strlen(linestr)+1);
             memcpy(*retstr,linestr,strlen(linestr)+1);
         }
         free(savedata);
         free(linestr);
         return strlen(*retstr);
      }
      
      if(gz_link->linepos==gz_link->buffersize)
      {
         if(savedata==NULL)
         {
            savedata = (char*)realloc(savedata,strlen(linestr)+1);
            memcpy(savedata,linestr,strlen(linestr)+1);
         }
         else
         {
         	datasize = strlen(savedata);
            savedata = (char*)realloc(savedata,datasize+strlen(linestr)+1);
            memcpy(savedata+datasize,linestr,strlen(linestr)+1);
         }
         
         //if(gz_link->buffersize<sizeof(gz_link->buffersize)-1)
         if(gz_link->buffersize<sizeof(gz_link->readbuf)-1)
         {
            /*
            if(savedata!=NULL)
            {
                *retstr = (char*)realloc(*retstr,strlen(savedata)+strlen(linestr)+1);
                memcpy(*retstr,savedata,strlen(savedata));
                memcpy((*retstr)+strlen(savedata),linestr,strlen(linestr)+1);
            }
            else
            {
                *retstr = (char*)realloc(*retstr,strlen(linestr)+1);
                memcpy(*retstr,linestr,strlen(linestr)+1);
            }
            */
            *retstr = (char*)realloc(*retstr,strlen(savedata)+1);
            memcpy(*retstr,savedata,strlen(savedata)+1);
            free(savedata);
            free(linestr);
            return strlen(*retstr);
         }
      }
      
    }
    return -1;
}

