#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "log.h"

/*
TODO:

*) introdurre una variabile per il controllo del flush bufferizzato

*) introdurre un puntatore alla variabile timestep in modo tale da poter avere delle funzioni
di log temporizzate. Qualcosa del tipo

logPrint( ...,500);  -> stampa al passo 500

*/

FILE *LogFile=NULL;
int Open=0;
int IsStdout=0;

int logStart(char *filename,char* mode)
{
	if (Open==1)
		return 1;
	
	LogFile=fopen(filename,mode);
	Open=1;
	return 0;
}

int logStartStdout()
{
	if (Open==1)
		return 1;
	
	IsStdout=1;
	Open=1;
	LogFile=stdout;
	return 0;
}

int logStartStderr()
{
	if (Open==1)
		return 1;
	
	IsStdout=1;
	Open=1;
	LogFile=stderr;
	return 0;
}

int logPrint(char *formato,...)
{
	if (Open==0)
	{
		return 1;
	}
	
	va_list args;
	va_start(args,formato);
	vfprintf(LogFile,formato,args);
	va_end(args);
	return 0;
}

int logFlush()
{
	if (Open)
	{
		fflush(LogFile);
		return 0;
	}
	return 1;
}

int logClose()
{
	if ((Open==1) && (IsStdout==0))
	{
		fclose(LogFile);
		Open=0;
		return 0;
	}
	else if (IsStdout==1)
		return 0;
	else
		return 1;
}

