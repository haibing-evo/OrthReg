/*******************************************************/
/*      Interface for Sequence Operations              */
/*      Date   : 2006.1.22                             */
/*      Author : Haibing Xie                           */
/*      Address: Kunming Institute of Zoology,CAS.     */
/*      Email  : xiehb@mail.kiz.ac.cn                  */
/*******************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "xiesequence.h"

#ifdef _MAC
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

//OK
long GetChrLength(char *filename,long LineLen,int CarriageLen)
{
    //Open File
   FILE *file;
   char ReadBuffer[101];
   int Carriage=0;
   long filelength=0,seqlength=0;
   int isfirstread = 1;
   int Count=0;
   long linecount=0;
   int i=0,tmplen=0,j=0;
    
   if((file = fopen(filename,"rb"))==NULL)
     return -1;

   while (Carriage==0)
   {
      fread(ReadBuffer,sizeof(char),100,file);
      ReadBuffer[100]='\0';
      if ((ReadBuffer[0]!='>')&&isfirstread) break;
      isfirstread = 0;
      Count = 0;
      while (Count<100 && Carriage==0)
      {   
        if (ReadBuffer[Count]=='\r'||ReadBuffer[Count]=='\n')
        {   
          Count++;
          if(ReadBuffer[Count] == '\n')
            Carriage = Count+1;
          else
            Carriage = Count;
        }
        Count++;
      }
   }


   fseeko(file,0,SEEK_END);
   filelength = ftello(file);

   seqlength = filelength-Carriage;
   if(LineLen==0) LineLen=50;
   if(CarriageLen==0) CarriageLen=1;
  
   linecount = seqlength/(LineLen+CarriageLen);
   fseeko(file,Carriage+linecount*(LineLen+CarriageLen),SEEK_SET);
   fread(ReadBuffer,sizeof(char),filelength-Carriage-linecount*(LineLen+CarriageLen),file);
   tmplen = filelength-Carriage-linecount*(LineLen+CarriageLen);

   for(;i<tmplen;i++)
   {
     if(!(ReadBuffer[i]=='\r'||ReadBuffer[i]=='\n'))
     {
       ReadBuffer[j] = ReadBuffer[i];
       j++;
     }
   }
   ReadBuffer[j]='\0';
   seqlength = strlen(ReadBuffer);
   fclose(file);
   return seqlength+linecount*LineLen;
}

void ReverseSequence(char *seq)
{
  long seqlen = strlen(seq);
  long seqhalflen = seqlen/2;
  char base;
  long i=0;
  for(i=0;i<seqhalflen;i++)
  {
     base = seq[i];
    seq[i] = seq[seqlen-1-i];
    seq[seqlen-1-i] = base;
  }
} 

void ComplementSequence(char *seq)
{
    long seqlength,position;
    position = 0;
    seqlength = strlen(seq);
    
    for(position=0;position<seqlength;position++)
  {
      switch(seq[position])
     {
       case 'A':
         seq[position] = 'T';
                 break;
       case 'T':
         seq[position] = 'A';
                 break;
       case 'G':
         seq[position] = 'C';
                 break;
       case 'C':
         seq[position] = 'G';
                 break;
       case 'a':
         seq[position] = 't';
         break;
       case 't':
              seq[position] = 'a';
         break;
       case 'g':
         seq[position] = 'c';
            break;
       case 'c':
         seq[position] = 'g';
         break;
       case 'Y':
         seq[position] = 'R';
         break;
       case 'R':
         seq[position] = 'Y';
         break;
       case 'W':
         seq[position] = 'S';
         break;
       case 'S':
         seq[position] = 'W';
         break;
       case 'M':
         seq[position] = 'K';
         break;
       case 'K':
         seq[position] = 'M';
         break;
       case 'y':
         seq[position] = 'r';
         break;
       case 'r':
         seq[position] = 'y';
         break;
       case 'w':
         seq[position] = 's';
         break;
       case 's':
         seq[position] = 'w';
         break;
       case 'm':
         seq[position] = 'k';
         break;
       case 'k':
         seq[position] = 'm';
         break;
       default:
         //seq[position] = seq[position];
         break;
       }
     }
}

//OK
void SeqUpper(char *seq)
{
    long seqlen,pos;
    seqlen = strlen(seq);
    for(pos=0;pos<seqlen;pos++)
    {
        seq[pos]=(char)toupper(seq[pos]);
    }
}

//OK
void SeqLower(char *seq)
{
    long seqlen,pos;
    seqlen = strlen(seq);
    for(pos=0;pos<seqlen;pos++)
    {
        seq[pos]=(char)tolower(seq[pos]);
    }
}

//OK
long ReadFasta(char *filename,char *strand,long OffSet,long ReadLen,long LineLen,int CarriageLen,
      char *PartData,long *readstart,long *readend,int *errorflag)
{
   FILE *file;
   long Count,Carriage,RelativePos,TotalLen;
   char *ReadBuffer,*SaveBuffer;
   int isfirstread=1;
   long i,j;

   *errorflag = 0; 
   file = fopen(filename,"rb");
   if(file==NULL)
   {
       *errorflag = 1; 
       return -1;
   }    

   ReadBuffer = (char*)malloc(1024);
   Carriage = 0;
   while (Carriage==0)
   {
    fread(ReadBuffer,sizeof(char),1023,file);
    ReadBuffer[1024]='\0';
    if ((ReadBuffer[0]!='>')&&isfirstread) break;
    isfirstread = 0;
    Count = 0;
    while (Count<1023 && Carriage==0)
    {   
     if (ReadBuffer[Count]=='\r'||ReadBuffer[Count]=='\n')
     {   
      Count++;
      if(ReadBuffer[Count] == '\n')
       Carriage = Count+1;
      else
       Carriage = Count;
     }
     Count++;
    }
  }
  if (OffSet<0) 
  {
   ReadLen+=(OffSet-1);
   OffSet = 1;
  }
  
  RelativePos = ((long)((OffSet-1)/LineLen))*CarriageLen;
  fseeko(file,RelativePos+OffSet+Carriage-1,SEEK_SET);
    
  TotalLen = ReadLen+(((long)((OffSet+ReadLen-1)/LineLen))-((long)((OffSet-1)/LineLen)))*CarriageLen;
  SaveBuffer = (char*)malloc(TotalLen+1);
  SaveBuffer[TotalLen] = '\0';
      
     TotalLen = fread(SaveBuffer,sizeof(char),TotalLen,file);
  for(i=0,j=0;i<TotalLen;i++)
  {
   if(!((SaveBuffer[i]=='\r')||(SaveBuffer[i]=='\n')))
   {
    PartData[j] = SaveBuffer[i];
    j++;
   }
  }
  PartData[j]= '\0';

  if(!strcmp(strand,"-"))
  {
   ReverseSequence(PartData);
   ComplementSequence(PartData);
  }
  fclose(file);
  *readstart = OffSet;
  *readend   = j+OffSet-1;
   
  free(ReadBuffer);
  free(SaveBuffer);
  return (*readend-*readstart+1);
} 

char f_codon(char *codon)
{
  if (codon[0] == 'T' || codon[0] == 'U' )
  {
         if (codon[1] == 'T' || codon[1] == 'U' )
    {
           if (codon[2] == 'A' ||codon[2] =='G' )
              //Leu
            return 'L';
           else
              //Phe  
        if (codon[2] == 'C' ||codon[2] =='T' ||codon[2] =='U'  )
                  return 'F' ; 
    }
        if (codon[1] == 'C')
            //Ser
            return 'S';   
        if (codon[1] == 'A')
    {
            if (codon[2] == 'A' || codon[2] == 'G' )
                //Stop
                return '.' ;
            else
               //Tyr
           if (codon[2] == 'C' ||codon[2] =='T' ||codon[2] =='U'  )
                    return 'Y' ; 
    }
 
        if ( codon[1] == 'G')
    {
         if (codon[2] == 'U' || codon[2]== 'T' || codon[2]== 'C')
                 //Cys
         return 'C';
         else
       {
           if (codon[2] =='A')
               //Stop
             return '.';
               else
             //Trp
             if (codon[2] =='G')
                 return 'W';
       }
    }
  }

  if (codon[0]=='C')
  {
    if ( codon[1] == 'T' || codon[1] == 'U' )
        //Leu
      return 'L';
    if ( codon[1] == 'C')
      //Pro
      return 'P';
        if ( codon[1] == 'A')
    {  
      if ( codon[2] == 'A' || codon[2] == 'G' )
           //Gln
          return 'Q';
        else
         //His
         if (codon[2] == 'U' || codon[2]== 'T' || codon[2]== 'C')
           return 'H';
    }
        if ( codon[1] == 'G')
       //Arg
       return 'R';
  }

  if ( codon[0] == 'A')
  {
    if ( codon[1] == 'T'|| codon[1] == 'U')
    {
      if ( codon[2] == 'G')
        //Met
        return 'M';
      else
        //Ile
        if (codon[2] == 'U' || codon[2]== 'T' || codon[2]== 'C' ||codon[2]=='A')
          return 'I';
    }
    if ( codon[1] == 'C')
       //Thr
      return 'T';
    if ( codon[1] == 'A')
    {
      if ( codon[2] == 'A'|| codon[2] == 'G')
         //Lys
         return 'K';
      else
         //Asn
         if (codon[2] == 'U' || codon[2]== 'T' || codon[2]== 'C')
            return 'N';
    }

    if ( codon[1] == 'G')
    {
      if ( codon[2] == 'A'|| codon[2] == 'G')
        //Arg
        return 'R';
      else
        //Ser
        if (codon[2] == 'U' || codon[2]== 'T' || codon[2]== 'C')
           return 'S';
    }

  }

  if ( codon[0] == 'G')
  {
    if ( codon[1] == 'T' ||  codon[1] == 'U')
      //Val
      return 'V';
    if ( codon[1] == 'C')
      //Ala
      return 'A';
    if ( codon[1] == 'A')
    {
      if ( codon[2] == 'A'|| codon[2] == 'G')
        //Glu
        return 'E';
      else
        //Asp
        if (codon[2] == 'U' || codon[2]== 'T' || codon[2]== 'C')
           return 'D';
    }
    if ( codon[1] == 'G')
       //Gly
       return 'G';
  }

  //When Error
     return '#';
}
