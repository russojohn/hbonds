#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/types.h>
#include <time.h>

#include "vector.h"
#include "secure_search.h"
#include "global_definitions.h"
#include "restart.h"
#include "random.h"


static const gsl_rng_type * Gsl_T;
static gsl_rng * Gsl_random;

static char Random_restart_file[100];
//static char Random_restart_file_read[100];

gsl_rng* randomConstructor(FILE *config_file)
{
	// build gsl variables
	gsl_rng_env_setup();
	Gsl_T = gsl_rng_rand48;
	Gsl_random = gsl_rng_alloc(Gsl_T);
	
	unsigned long int seed;
	
	SearchTable *s=searchNew();
	int *found_seed=searchTryULInt("Seed",&seed,s);
	searchString("Restart_random",Random_restart_file,s);
	//searchTryString("Restart_random_start",Random_restart_file_read,s);
	searchFile(config_file,s);
	
	if (!*found_seed)
	{
		time_t t1;
		(void)time(&t1);
		seed=(unsigned long int)t1;
	}
	gsl_rng_set(Gsl_random,seed);
	
	searchFree(s);
	
	return Gsl_random;
}

gsl_rng* randomConstructorInteractive(int seed)
{
	// build gsl variables
	gsl_rng_env_setup();
	const gsl_rng_type *Gsl_T = gsl_rng_rand48;
	gsl_rng *gsl_random = gsl_rng_alloc(Gsl_T);
	
	if (seed==-1)
	{
		time_t t1;
		(void)time(&t1);
		seed=(unsigned long int)t1;
	}
	gsl_rng_set(gsl_random,seed);
	
	return gsl_random;
}


void randomFree()
{
	gsl_rng_free(Gsl_random);
}

gsl_rng* randomStructure()
{
	return Gsl_random;
}


void saveRandom(steps step)
{
	if (strcmp(Random_restart_file,"none")!=0)
	{
		restartCopy(Random_restart_file);
		
		char name[200];
		sprintf(name,"%s%lld",Random_restart_file,step);
		
		//FILE *pfile=fopen(name,"wb");
		FILE *pfile=fopen(Random_restart_file,"wb");
		
		if (pfile==NULL)
		{
			printf("Error: cannot write restart file '%s'\n",Random_restart_file);
			exit(1);
		}
		
		fwrite(&step,sizeof(steps),1,pfile);
		gsl_rng_fwrite(pfile,Gsl_random);
		fclose(pfile);
	}
}

steps getRandom()
{
	steps step;
	size_t status;
	//FILE *pfile=fopen(Random_restart_file_read,"rb");
	FILE *pfile=fopen(Random_restart_file,"rb");
	
	if (pfile==NULL)
	{
		printf("Error: cannot read restart file '%s'\n",Random_restart_file);
		exit(1);
	}
	
	status=fread(&step,sizeof(steps),1,pfile);
	// nessun problema di memory leak anche se Gsl_random è già stato aperto
	gsl_rng_fread(pfile,Gsl_random);
	fclose(pfile);
	
	return step;
}

