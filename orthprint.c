#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#ifdef _MAC
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

struct orthinfo
{
  char chr[6];
  long orthsite;
  char nt;
  char strand;
  float encode;
};

int main(int argc, char *argv[])
{
  if(argc!=3)
  {
     printf("Usage:\n");
     printf("%s orthfile position\n",argv[0]);
     return 0;
  }
  struct orthinfo siteinfo;
  long startpos;
  FILE * pfile;
  startpos = atof(argv[2]);
  if(startpos<=0)
  {
     printf("position should be > 0\n");
     return 0;
  }
  pfile = fopen(argv[1],"rb");
  if(pfile==NULL)
  {
     printf("orthfile not found\n");
     return 0;
  }
  fseek(pfile,sizeof(struct orthinfo)*(startpos-1),SEEK_SET);
  fread(&siteinfo,sizeof(struct orthinfo),1,pfile);
  
  if(siteinfo.orthsite!=0)
    printf("%s\t%c\t%ld\t%c\t%f\n",siteinfo.chr,siteinfo.strand,siteinfo.orthsite,siteinfo.nt,siteinfo.encode);
  else
    printf("-\t-\t-\t-\t-\n");
  
  fclose(pfile);
  return 1;  
}