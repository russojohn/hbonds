#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include "dirutils.h"

#define MAXLENGTH 100
#define NUMFILES 10000


#define MAXDIRDEPTH 20




/* scan directory searching for file with particular pattern */
int scandir_pattern(const char *dir, char ***namelist,const char *pattern)
{
	DIR *dp;
	dp=opendir(dir);
	int size=NUMFILES;
	int num_files=0;
	struct dirent *ep;
	int i;
	
	if (dp==NULL)
	{
		printf("Error: Could not open directory '%s'\n",dir);
		exit(1);
	}
	
	(*namelist)=(char **)malloc(size*sizeof(char*));
	for (i=0;i<size;i++)
		(*namelist)[i]=(char*)malloc(MAXLENGTH*sizeof(char));
	
	
	while (  (ep=readdir(dp)) !=NULL )
	{
		/* compare file with pattern */
		if (strncmp(pattern,ep->d_name,strlen(pattern))==0)
		{
			/* match */
			strcpy((*namelist)[num_files++],ep->d_name);
			if (num_files==size)
			{
				size*=2;
				(*namelist)=realloc(*namelist,size*sizeof(char*));
				for (i=(size/2);i<size;i++)
					(*namelist)[i]=(char*)malloc(MAXLENGTH*sizeof(char));
			}
			
		}
	}
	
	closedir(dp);
	
	// resize the namelist
	for (i=num_files;i<size;i++)
		free((*namelist)[i]);
	(*namelist)=realloc(*namelist,num_files*sizeof(char*));
	
	return num_files;
	
}


void freeListFiles(char **list,int num_files)
{
	int i;
	
	for (i=0;i<num_files;i++)
	{
		free(list[i]);
	}
	free(list);
}

int sortFunction(const void *node1,const void *node2)
{
	int t1,t2;
	
	t1=*(((sortstruct*)node1)->time);
	t2=*(((sortstruct*)node2)->time);
	
	if (t1<t2) return -1;
	if (t1>t2) return 1;
	else return 0;
	
}

void sortList(char **list,int num,int (*sortFunction)(const void *,const void *),char *prefix)
{
	int i;
	
	int tempi[num];
	sortstruct s[num];
	
	for (i=0;i<num;i++)
	{
		tempi[i]=atoi(list[i]+strlen(prefix));
		
		s[i].time=tempi+i;
		s[i].nomefile=list[i];
	}
	
	qsort(s,num,sizeof(sortstruct),sortFunction);
	
	for (i=0;i<num;i++)
	{
		list[i]=s[i].nomefile;
	}
	
}

void sortImproved(char **list,int num,int (*estraiOrdine)(char *))
{
	int i;
	
	int tempi[num];
	sortstruct s[num];
	
	for (i=0;i<num;i++)
	{
		tempi[i]=estraiOrdine(list[i]);
		
		s[i].time=tempi+i;
		s[i].nomefile=list[i];
	}
	
	qsort(s,num,sizeof(sortstruct),sortFunction);
	
	for (i=0;i<num;i++)
	{
		list[i]=s[i].nomefile;
	}
	
}

void getDirAndPattern(char *dirSlashPattern,char *dir,char *pattern)
{
        char (*directories)[MAXDIRLENGTH];
        
        directories=calloc(MAXDIRDEPTH,MAXDIRLENGTH*sizeof(char));
        
        char copy[MAXDIRDEPTH*MAXDIRLENGTH+MAXDIRDEPTH];
        
        strcpy(copy,dirSlashPattern);
        
        char *token=strtok(copy,"/");
        int pos=0;
        while (token!=NULL)
        {
                strcpy(directories[pos++],token);
                token=strtok(NULL,"/");
        }
        
        if (pos>MAXDIRDEPTH)
        {
                printf("Error: increare MAXDIRDEPTH\n");
                exit(1);
        }
        
        // pattern è l'ultimo elemento di directories
        strcpy(pattern,directories[pos-1]);
        
        // ricostruiamo il ramo delle directory
        int i;
        strcpy(dir,".");
        for (i=0;i<pos-1;i++)
        {
                strcat(dir,"/");
                strcat(dir,directories[i]);
        }
        
        free(directories);
}



void completeName(char *complete,char *dir,char *file)
{
	strcpy(complete,dir);
	strcat(complete,"/");
	strcat(complete,file);
}
