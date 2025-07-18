// FIL files module  //

#ifndef RESTART_H
#define RESTART_H

/*
# 'Restart' possibilies:
	#     0 - continue simulation from 'Initial_conditions_file'
	#     1 - continue simulation from restart files
	#     2 - start new simulation from 'Initial_conditions_file'
	#     3 - start new simulation from restart files
*/

FILE* openFile(char[],char[]);

void closeFile(FILE*);

void restartConstructor(FILE *Configname);

void addRestartModule(void (*s)(unsigned long long int),unsigned long long int (*g)(void));

void closeRestartableFiles();

FILE* openRestartableFile(char filename[]);

void saveRestart(unsigned long long int step);

void getRestart(unsigned long long int *step);

int restartCopy(char *restart_name);

void flushRestartableFiles();

#endif
