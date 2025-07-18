#ifndef DIRUTILS_H
#define DIRUTILS_H

#define MAXDIRLENGTH 50

typedef struct _sortstruct
{
	int *time;
	char *nomefile;
} sortstruct;


int scandir_pattern(const char *, char ***,const char *);

//void freeList(char **list,int num_files);
void freeListFiles(char **list,int num_files);

int sortFunction(const void *node1,const void *node2);

void sortList(char **list,int num,int (*sortFunction)(const void *,const void *),char *prefix);

void sortImproved(char **list,int num,int (*estraiOrdine)(char *));

void getDirAndPattern(char *dirSlashPattern,char *dir,char *pattern);

void completeName(char *complete,char *dir,char *file);

#endif
