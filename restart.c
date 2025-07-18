#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>


#include "log.h"
#include "vector.h"
#include "utilities.h"
#include "secure_search.h"
#include "restart.h"

#define NULLKEYWORD "none"

#define MaxFiles 11
#define MaxRestartModules 10

//                            VARIABLES                          //
static char Restart_pos[100];
//static char Restart_pos_read[100];

// copy restarts
static char Restart_copy_prefix[100];
static int Copy=0;


static int NFiles=0;                         // number of files to be restarted
static char* FileNames[MaxFiles];                   // name of files to be restarted
static FILE* Files[MaxFiles];                       // files opened during simulation and to be restarted

void (*saveArray[MaxRestartModules])(unsigned long long int);
unsigned long long int (*getArray[MaxRestartModules])(void);
static int NRObjects=0;

static int Restart;
static int Get_positions=0;

static int Save=1;   // if Save is 0 no restart modules will be written


void restartConstructor(FILE *pfile)
{
	
	SearchTable *s=searchNew();
	int *found_file=searchTryString("Restart_file",Restart_pos,s);
	//int *found_filestart=searchTryString("Restart_file_start",Restart_pos_read,s);
	searchInt("Restart",&Restart,s);
	int *found_copy=searchTryString("Restart_copy_prefix",Restart_copy_prefix,s);
	
	
	searchFile(pfile,s);
	
	// settiamo il file di restart
	if ( (!*found_file) || (!strcmp(Restart_pos,"none")) )
		Save=0;
	else
	{
		Save=1;
		// controlliamo se sono abilitate le copie
		if ( (*found_copy) && (strcmp(Restart_copy_prefix,"none")!=0) )
		{
			Copy=1;
		}
	}
	
	/* chiave disabilitata
	if (Restart%2==1)
	{
		if (!*found_filestart)
		{
			logPrint("Error: cannot read file %s\n",Restart_pos_read);
			exit(1);
		}
	}
	*/
	
	
	searchFree(s);
	
	if (Restart==1)
	{
		Get_positions=1;
	}
	
}


void addRestartModule(void (*s)(unsigned long long int),unsigned long long int (*g)(void))
{
	saveArray[NRObjects]=s;
	getArray[NRObjects++]=g;
}

void closeRestartableFiles()
{
	int i;
	
	for (i=0;i<NFiles;i++)
		fclose(Files[i]);
	
}

FILE* openRestartableFile(char filename[])
{
	
	if (!strcmp(filename,NULLKEYWORD))
	{
		return NULL;
	}
	else
	{
		FileNames[NFiles]=filename;
		
		if (Restart==1)
		{
			Files[NFiles]=fopen(filename,"r+");
			if (Files[NFiles]==NULL)
			{
				printf("Error: missing file %s\n",filename);
				exit(1);
			}
		}
		else if (Restart==0)
			Files[NFiles]=fopen(filename,"a");
		else
			Files[NFiles]=fopen(filename,"w");
		
		return Files[NFiles++];
	}
}


FILE* openFile(char filename[],char mod[])
// Files which are not scheduled for restart
{
	if (!strncmp(filename,NULLKEYWORD,4))
	{
		return NULL;
	}
	else
	{
		FILE *pfile;
		
		pfile=fopen(filename,mod);
		
		if (pfile==NULL)
		{
			printf("Error: missing file %s\n",filename);
			exit(1);
		}
		
		return pfile;
		
	}
}


void closeFile(FILE* pfile)
{
	if (pfile==NULL) return;
	else fclose(pfile);
}


void saveRestart(unsigned long long int step)
{
	int i;
	
	// forcing file writes before writing restart informations
	for (i=0;i<NFiles;i++)
	{
		fflush(Files[i]);
	}
	
	if (!Save)
	{
		return;
	}
	
	restartCopy(Restart_pos);
	
	FILE *file_pos;
	char name[200];
	sprintf(name,"%s%lld",Restart_pos,step);
	
	//file_pos=fopen(name,"wb");
	file_pos=fopen(Restart_pos,"wb");
	
	if (file_pos==NULL)
	{
		printf("Error: cannot write restart files\n");
		exit(1);
	}
	
	// scriviamo il tempo nel file posizioni
	fwrite(&step,sizeof(unsigned long long int),1,file_pos);
	
	// file positions save
	fpos_t position;
	
	memset(&position,0,sizeof(position));
	
	for (i=0;i<NFiles;i++)
	{
		if (!fgetpos(Files[i],&position))
			fwrite(&position,sizeof(fpos_t),1,file_pos);
		else
		{
			printf("Error: cannot save current state in %s\n",FileNames[i]);
			exit(1);
		}
	}
	fclose(file_pos);
	
	// saving restart objects
	for (i=0;i<NRObjects;i++)
	{
		saveArray[i](step);
	}
	
	
}

void getRestart(unsigned long long int *step)
{
	int i;
	size_t status;
	
	if ( (Restart==1) || (Restart==3) )
	{
		if (!Save)
		{
			printf("Error: asking to restart without restart file\n");
			exit(1);
		}
		
		unsigned long long int file_step=0;
		
		if (Get_positions)
		{
			FILE *file_pos;
			
			//file_pos=fopen(Restart_pos_read,"rb");
			file_pos=fopen(Restart_pos,"rb");
			
			if (file_pos==NULL)
			{
				printf("Error: missing restart file '%s'\n",Restart_pos);
				exit(1);
			}
			
			// leggiamo il tempo nel file posizioni
			status=fread(&file_step,sizeof(unsigned long long int),1,file_pos);
			file_step++;
			
			// file position recovery
			fpos_t position;
			
			
			memset(&position,0,sizeof(position));
			
			for (i=0;i<NFiles;i++)
			{
				status=fread(&position,sizeof(fpos_t),1,file_pos);
				fsetpos(Files[i],&position);
				
			}
			
			fclose(file_pos);
		}
		
		unsigned long long int steps[NRObjects];
		
		
		// get restart objects
		for (i=0;i<NRObjects;i++)
		{
			steps[i]=getArray[i]();
		}
		/*
		// check
		for (i=1;i<NRObjects;i++)
		{
			if (steps[1]!=steps[0])
			{
				printf("Error: restarting from different time files");
				exit(1);
			}
		}
		
		*/
		if (Restart==1)
			*step=file_step;
		else
			*step=1;
		
	}
}

int restartCopy(char *restart_name)
{
	if (Copy==0)
		return 1;
	
	FILE *fp1;
	FILE *fp2;
	
	if ((fp1 = fopen(restart_name, "rb")) == 0)
	{
		logPrint("Warning: cannot open file %s for reading\n",restart_name);
		return 1;
	}
	
	char copy_name[100];
	sprintf(copy_name,"%s%s",Restart_copy_prefix,restart_name);
	
	if ((fp2 = fopen(copy_name, "wb")) == 0)
	{
		logPrint("Warning: cannot open file %s for writing\n", copy_name);
		return 1;
	}
	
	 fcopy(fp1, fp2);
	 
	 fclose(fp1);
	 fclose(fp2);
	 
	return 0;
}


void flushRestartableFiles()
{
	int i;
	
	for (i=0;i<NFiles;i++)
	{
		fflush(Files[i]);
	}
}
