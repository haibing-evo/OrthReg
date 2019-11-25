#ifndef _xiesequence

#define _xiesequence

#ifdef __cplusplus
extern "C" {
#endif

long GetChrLength(char *filename,long LineLen,int CarriageLen);
void ReverseSequence(char *seq);
void ComplementSequence(char *seq);
void SeqUpper(char *seq);
void SeqLower(char *seq);
long ReadFasta(char *filename,char *strand,long OffSet,long ReadLen,long LineLen,int CarriageLen,char *PartData,long *readstart,long *readend,int *errorflag);
char f_codon(char *codon);

#ifdef __cplusplus
}
#endif

#endif

