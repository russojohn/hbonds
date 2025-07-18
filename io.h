#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#define MAX_LINE_LENGTH 2000

int getNumberParticles(char *_input_name);
void getBoxLength(char *_input_name, double Box[]);
void getHeader(char *_input_name,steps *_step,int *_numparticles,double Box[]);
void readPositions(char *_input_name,vector *pos,steps *time,int *numparticles,double Box[],double INOBox[]);
void readPositionsKeepRealSpace(char *_input_name,vector *pos,steps *time,int *numparticles,double NOBox[],double INOBox[]);
void readPositionsPutZeroAtZero(char *_input_name,vector *pos,steps *time,int *numparticles,double NOBox[],double INOBox[]);
void savePositions(char *output_name,vector *pos,steps time,int numparticles,double Box[]);


// cancellare se ci sono conflitti
int getLine(char *line,FILE *pfile);

#endif

