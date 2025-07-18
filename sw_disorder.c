#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_erf.h>


#include "vector.h"
#include "global_definitions.h"
#include "log.h"
#include "secure_search.h"
#include "interaction_map.h"
#include "smart_allocator.h"
#include "bilista.h"
#include "events.h"
#include "restart.h"
#include "random.h"
#include "utilities.h"
#include "sw_disorder.h"


#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

// wrapper per eventi
void swVirialSystem(args *argv,int argc);
void swEnergySystem(args *argv,int argc);
void swEnergySystemHamiltonian(args *argv,int argc);
void swRejectParticleMove(args *argv,int argc);
void swRejectSwapMove(args *argv,int argc);
void swRejectSystemMove(args *argv,int argc);
void swRejectParticleMoveHamiltonian(args *argv,int argc);
void swRejectSystemMoveHamiltonian(args *argv,int argc);
void swEnergyParticleOld(args *argv,int argc);
void swEnergyParticleNew(args *argv,int argc);
void swEnergySwapOld(args *argv,int argc);
void swEnergySwapNew(args *argv,int argc);
void swEnergyParticleOldHamoltonian(args *argv,int argc);
void swEnergyParticleNewHamoltonian(args *argv,int argc);
void swSaveState(args *argv,int argc);
void swSaveStateHamiltonian(args *argv,int argc);
void swLoadState(args *argv,int argc);
void swLoadStateHamiltonian(args *argv,int argc);
void swRemoveParticle(args *argv,int argc);
void swRemoveAccept(args *argv,int argc);
void swRemoveReject(args *argv,int argc);
void swAddParticle(args *argv,int argc);
void swAddReject(args *argv,int argc);

// parametri del modulo mc
static args *ParametersMc;
static int NumParametersMc;

// parametri privati del modulo
static sw *Sws_1;
static sw *Sws_2;

static double Potential_1;
static double Potential_2;

static double Potential_2body;
static double Potential_3body;

static FILE *Hamiltonian_file;

static int Model;

// hamiltonian integration
static double Hlambda;

// angular disorder
static int Disorder;
static double Disorder_amplitude;

static gsl_rng * Gsl_random;


static inline double intpow (const double x, const int i) {
	if (i < 0) abort();
	if (i == 0) return 1.;
	if (i == 1) return x;
	else return x * intpow(x, i - 1);
}

sw* sw_pointer()
{
  return Sws_1;
}

sw* createSW(int max_colloids,int *ncolloids,double lambda)
{
	sw *sws=malloc(sizeof(sw));

	Matrix2D(sws->TwoBody,max_colloids,10,double);
	Matrix2D(sws->TwoBody_copy,max_colloids,10,double);
	Matrix2D(sws->TwoBody_savestate,max_colloids,10,double);
	sws->WhoCopied=calloc(max_colloids,sizeof(int));
	sws->List_modified=bilistaGet(max_colloids);
	sws->Num_copies=0;
	sws->NColloids=ncolloids;

	sws->A=7.049556277;
	sws->B=0.6022245584;
	sws->p=4;
	sws->q=0;
	sws->Gamma=1.2;
	sws->costeta0=calloc(max_colloids,sizeof(double));
	int i;
	for (i=0;i<max_colloids;i++)
		sws->costeta0[i]=-1./3.;

	sws->lambda=lambda;
	sws->a=1.8;

	return sws;
}

void freeSW(sw *sws)
{
	Free2D(sws->TwoBody);
	Free2D(sws->TwoBody_copy);
	Free2D(sws->TwoBody_savestate);
	bilistaFree(sws->List_modified);
	free(sws->WhoCopied);
	free(sws->costeta0);
}

void swConstructor(FILE *config_file,vector *pos,int max_colloids,int *ncolloids,double box[],double *cutoff,interactionmap *interactionList,double *potential,double *virial,vector *force)
{
	char hamiltonian_name[100];

	// leggiamo se acqua o silica
	Model=1;
	SearchTable *s=searchNew();
	int *found_model=searchTryInt("Model",&Model,s);
	int *found_Hlambda=searchTryDouble("Hamiltonian_lambda",&Hlambda,s);
	int *found_file=searchTryString("Hamiltonian_file",hamiltonian_name,s);
	int *found_disorder=searchTryInt("Disorder",&Disorder,s);
	int *found_disordermax=searchTryDouble("Disorder_amplitude",&Disorder_amplitude,s);

	// 0 Silicon
	// 1 Water
	// 2 Hamiltonian integration H=H_lambda*H_1+(1-H_lambda)*H_2
	double lambda_1;
	double lambda_2;
	double lambda;

	int *found_lambda_1=searchTryDouble("SW_lambda_1",&lambda_1,s);
	int *found_lambda_2=searchTryDouble("SW_lambda_2",&lambda_2,s);
	int *found_lambda=searchTryDouble("SW_lambda",&lambda,s);

	searchFile(config_file,s);


	if (!*found_model)
	{
		logPrint("\n\nWarning: Keyword 'Model' not found - reverting to Water\n\n");
		Model=1;
	}


	// eventi
	int status=0;
	status+=eventAddEvent(ENERGY_SYSTEM);
	status+=eventAddEvent(ENERGY_PARTICLE_OLD);
	status+=eventAddEvent(ENERGY_PARTICLE_NEW);
	status+=eventAddEvent(ENERGY_SWAP_OLD);
	status+=eventAddEvent(ENERGY_SWAP_NEW);
	status+=eventAddEvent(PARTICLEMOVE_REJECTED);
	status+=eventAddEvent(SYSTEMMOVE_REJECTED);
	status+=eventAddEvent(TRAJECTORYMOVE_REJECTED);
	status+=eventAddEvent(TRAJECTORYMOVE_BEGIN);
	status+=eventAddEvent(PAIRMOVE_REJECTED);
	status+=eventAddEvent(REMOVE_PARTICLE);
	status+=eventAddEvent(REMOVE_ACCEPT);
	status+=eventAddEvent(REMOVE_REJECT);
	status+=eventAddEvent(ADD_PARTICLE);
	status+=eventAddEvent(ADD_ACCEPT);
	status+=eventAddEvent(ADD_REJECT);
	status+=eventAddEvent(COMPUTE_VIRIAL);

	/*if (status!=7)
	{
		logPrint("Error in definitions of events in module sw\n");
		exit(1);
	}
	*/

	Sws_1=NULL;
	Sws_2=NULL;


	NumParametersMc=11;
	ParametersMc=calloc(NumParametersMc,sizeof(args));

	ParametersMc[2].argv=(void*)&Hlambda;
	ParametersMc[3].argv=(void*)pos;
	ParametersMc[4].argv=(void*)ncolloids;
	ParametersMc[5].argv=(void*)box;
	ParametersMc[6].argv=(void*)interactionList;
	ParametersMc[7].argv=(void*)cutoff;
	ParametersMc[8].argv=(void*)potential;
	ParametersMc[9].argv=(void*)virial;
	ParametersMc[10].argv=(void*)force;

	eventAddAction(REMOVE_ACCEPT,&swRemoveAccept,ParametersMc,NumParametersMc);
	eventAddAction(REMOVE_REJECT,&swRemoveReject,ParametersMc,NumParametersMc);
	eventAddAction(ADD_REJECT,&swAddReject,ParametersMc,NumParametersMc);

	if (Model==0)
	{
		// silicon
		logPrint("\nModel: Silicon\n");
		Sws_1=createSW(max_colloids,ncolloids,21.0);
		Sws_2=NULL;
		ParametersMc[0].argv=(void*)Sws_1;
		ParametersMc[1].argv=(void*)Sws_2;


		eventAddAction(ENERGY_SYSTEM,&swEnergySystem,ParametersMc,NumParametersMc);
		eventAddAction(ENERGY_PARTICLE_OLD,&swEnergyParticleOld,ParametersMc,NumParametersMc);
		eventAddAction(ENERGY_PARTICLE_NEW,&swEnergyParticleNew,ParametersMc,NumParametersMc);

		eventAddAction(ENERGY_SWAP_OLD,&swEnergySwapOld,ParametersMc,NumParametersMc);
		eventAddAction(ENERGY_SWAP_NEW,&swEnergySwapNew,ParametersMc,NumParametersMc);

		eventAddAction(PARTICLEMOVE_REJECTED,&swRejectParticleMove,ParametersMc,NumParametersMc);
		eventAddAction(SYSTEMMOVE_REJECTED,&swRejectSystemMove,ParametersMc,NumParametersMc);
		eventAddAction(TRAJECTORYMOVE_REJECTED,&swLoadState,ParametersMc,NumParametersMc);
		eventAddAction(TRAJECTORYMOVE_BEGIN,&swSaveState,ParametersMc,NumParametersMc);
		eventAddAction(PAIRMOVE_REJECTED,&swRejectSwapMove,ParametersMc,NumParametersMc);

		eventAddIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__energy,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__energy,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__virial,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__virial,sizeof(double));

		eventAddIO(ENERGY_SWAP_OLD,ENERGY_SWAP_OLD__energy,sizeof(double));
		eventAddIO(ENERGY_SWAP_NEW,ENERGY_SWAP_NEW__energy,sizeof(double));
		eventAddIO(ENERGY_SWAP_OLD,ENERGY_SWAP_OLD__virial,sizeof(double));
		eventAddIO(ENERGY_SWAP_NEW,ENERGY_SWAP_NEW__virial,sizeof(double));

		eventAddAction(REMOVE_PARTICLE,&swRemoveParticle,ParametersMc,NumParametersMc);
		eventAddAction(ADD_PARTICLE,&swAddParticle,ParametersMc,NumParametersMc);
	}
	else if (Model==1)
	{
		// water
		logPrint("\nModel: Water\n");
		Sws_1=createSW(max_colloids,ncolloids,23.15);
		Sws_2=NULL;

		ParametersMc[0].argv=(void*)Sws_1;
		ParametersMc[1].argv=(void*)Sws_2;


		eventAddAction(ENERGY_SYSTEM,&swEnergySystem,ParametersMc,NumParametersMc);
		eventAddAction(ENERGY_PARTICLE_OLD,&swEnergyParticleOld,ParametersMc,NumParametersMc);
		eventAddAction(ENERGY_PARTICLE_NEW,&swEnergyParticleNew,ParametersMc,NumParametersMc);

		eventAddAction(ENERGY_SWAP_OLD,&swEnergySwapOld,ParametersMc,NumParametersMc);
		eventAddAction(ENERGY_SWAP_NEW,&swEnergySwapNew,ParametersMc,NumParametersMc);

		eventAddAction(PARTICLEMOVE_REJECTED,&swRejectParticleMove,ParametersMc,NumParametersMc);
		eventAddAction(SYSTEMMOVE_REJECTED,&swRejectSystemMove,ParametersMc,NumParametersMc);
		eventAddAction(TRAJECTORYMOVE_REJECTED,&swLoadState,ParametersMc,NumParametersMc);
		eventAddAction(TRAJECTORYMOVE_BEGIN,&swSaveState,ParametersMc,NumParametersMc);
		eventAddAction(PAIRMOVE_REJECTED,&swRejectSwapMove,ParametersMc,NumParametersMc);


		eventAddIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__energy,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__energy,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__virial,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__virial,sizeof(double));

		eventAddIO(ENERGY_SWAP_OLD,ENERGY_SWAP_OLD__energy,sizeof(double));
		eventAddIO(ENERGY_SWAP_NEW,ENERGY_SWAP_NEW__energy,sizeof(double));
		eventAddIO(ENERGY_SWAP_OLD,ENERGY_SWAP_OLD__virial,sizeof(double));
		eventAddIO(ENERGY_SWAP_NEW,ENERGY_SWAP_NEW__virial,sizeof(double));

		eventAddAction(REMOVE_PARTICLE,&swRemoveParticle,ParametersMc,NumParametersMc);
		eventAddAction(ADD_PARTICLE,&swAddParticle,ParametersMc,NumParametersMc);
	}
	else if (Model==2)
	{
		logPrint("\nModel: Hamiltonian integration\n");

		if ( (!*found_lambda_1) || (!*found_lambda_2) || (!*found_Hlambda) )
		{
			logPrint("\n\nError: Missing hamiltonian integration parameters\n\n");
			exit(1);
		}

		if (!*found_file)
			strcpy(hamiltonian_name,"energy_hamiltonian.dat");
		Hamiltonian_file=openRestartableFile(hamiltonian_name);

		Sws_1=createSW(max_colloids,ncolloids,lambda_1);
		Sws_2=createSW(max_colloids,ncolloids,lambda_2);

		ParametersMc[0].argv=(void*)Sws_1;
		ParametersMc[1].argv=(void*)Sws_2;

		eventAddNotImplemented(ENERGY_SWAP_OLD,"Hamiltonian Integration does not support swap moves");
		eventAddNotImplemented(ENERGY_SWAP_NEW,"Hamiltonian Integration does not support swap moves");

		eventAddAction(ENERGY_SYSTEM,&swEnergySystemHamiltonian,ParametersMc,NumParametersMc);
		eventAddAction(ENERGY_PARTICLE_OLD,&swEnergyParticleOldHamoltonian,ParametersMc,NumParametersMc);
		eventAddAction(ENERGY_PARTICLE_NEW,&swEnergyParticleNewHamoltonian,ParametersMc,NumParametersMc);

		eventAddAction(PARTICLEMOVE_REJECTED,&swRejectParticleMoveHamiltonian,ParametersMc,NumParametersMc);
		eventAddAction(SYSTEMMOVE_REJECTED,&swRejectSystemMoveHamiltonian,ParametersMc,NumParametersMc);
		eventAddAction(TRAJECTORYMOVE_REJECTED,&swLoadStateHamiltonian,ParametersMc,NumParametersMc);
		eventAddAction(TRAJECTORYMOVE_BEGIN,&swSaveStateHamiltonian,ParametersMc,NumParametersMc);

		eventAddIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__energy,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__energy,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__virial,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__virial,sizeof(double));


		eventAddNotImplemented(REMOVE_PARTICLE,"Hamiltonian Integration does not support gran canonical ensemble");
		eventAddNotImplemented(ADD_PARTICLE,"Hamiltonian Integration does not support gran canonical ensemble");
	}
	else if (Model==3)
	{
		// water
		logPrint("\nModel: Custom SW Potential with lambda = %lf\n",lambda);

		if (!*found_lambda)
		{
			logPrint("\n\nError: Missing 'SW_lambda' key\n\n");
			exit(1);
		}

		Sws_1=createSW(max_colloids,ncolloids,lambda);
		Sws_2=NULL;

		ParametersMc[0].argv=(void*)Sws_1;
		ParametersMc[1].argv=(void*)Sws_2;


		eventAddAction(ENERGY_SYSTEM,&swEnergySystem,ParametersMc,NumParametersMc);
		eventAddAction(ENERGY_PARTICLE_OLD,&swEnergyParticleOld,ParametersMc,NumParametersMc);
		eventAddAction(ENERGY_PARTICLE_NEW,&swEnergyParticleNew,ParametersMc,NumParametersMc);

		eventAddAction(ENERGY_SWAP_OLD,&swEnergySwapOld,ParametersMc,NumParametersMc);
		eventAddAction(ENERGY_SWAP_NEW,&swEnergySwapNew,ParametersMc,NumParametersMc);


		eventAddAction(PARTICLEMOVE_REJECTED,&swRejectParticleMove,ParametersMc,NumParametersMc);
		eventAddAction(SYSTEMMOVE_REJECTED,&swRejectSystemMove,ParametersMc,NumParametersMc);
		eventAddAction(TRAJECTORYMOVE_REJECTED,&swLoadState,ParametersMc,NumParametersMc);
		eventAddAction(TRAJECTORYMOVE_BEGIN,&swSaveState,ParametersMc,NumParametersMc);
		eventAddAction(PAIRMOVE_REJECTED,&swRejectSwapMove,ParametersMc,NumParametersMc);

		eventAddIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__energy,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__energy,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__virial,sizeof(double));
		eventAddIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__virial,sizeof(double));

		eventAddIO(ENERGY_SWAP_OLD,ENERGY_SWAP_OLD__energy,sizeof(double));
		eventAddIO(ENERGY_SWAP_NEW,ENERGY_SWAP_NEW__energy,sizeof(double));
		eventAddIO(ENERGY_SWAP_OLD,ENERGY_SWAP_OLD__virial,sizeof(double));
		eventAddIO(ENERGY_SWAP_NEW,ENERGY_SWAP_NEW__virial,sizeof(double));

		eventAddAction(REMOVE_PARTICLE,&swRemoveParticle,ParametersMc,NumParametersMc);
		eventAddAction(ADD_PARTICLE,&swAddParticle,ParametersMc,NumParametersMc);
	}
	else
	{
		logPrint("\n\Error: Keyword 'Model' unrecognized key\n\n");
		exit(1);
	}

	// introduciamo una funzione ad hoc per il calcolo del viriale
	// che essendo piu' pesante nel caso di interazioni a tre corpi
	// permettiamo di calcolare separatamente
	eventAddAction(COMPUTE_VIRIAL,&swVirialSystem,ParametersMc,NumParametersMc);


	if ((*found_disorder) && (Disorder==1))
	{
		if (!*found_disordermax)
		{
			logPrint("\n\Error: missin Keyword 'Disorder_amplitude'\n\n");
			exit(1);
		}

		Gsl_random=randomStructure();

		int i;
		for (i=0;i<max_colloids;i++)
		{
			double dteta=Disorder_amplitude*(2.*gsl_rng_uniform(Gsl_random)-1.)*M_PI/180.;

			Sws_1->costeta0[i]=cos(acos(-1./3.)+dteta);

			if (Sws_2!=NULL)
			{
				Sws_2->costeta0[i]=cos(acos(-1./3.)+dteta);
			}
		}

		logPrint("\nDisorder amplitude: %lf\n",Disorder_amplitude);

	}



	searchFree(s);


}

void swFree()
{
	freeSW(Sws_1);
	free(Sws_1);

	if (Model==2)
	{
		freeSW(Sws_2);
		free(Sws_2);
	}
	free(ParametersMc);
}


void swVirialSystem(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	//sw *sws_2=(sw*)argv[1].argv;
	//double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	double *box=(double*)argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	double *potential=(double*)argv[8].argv;
	double *virial=(double*)argv[9].argv;
	vector *force=(vector*)argv[10].argv;

	//double old_potential=*potential;

	sw_system_virial(sws_1,pos,im,box,cutoff,potential,virial,force);


// 	if (fabs(old_potential-*potential)>0.00001)
// 	{
// 		logPrint("Checks error: potential energy is %lf instead of %lf after virial calculation\n",*potential,old_potential);
// 		exit(1);
// 	}

}

void swEnergySystem(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	//sw *sws_2=(sw*)argv[1].argv;
	//double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	double *box=(double*)argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	double *potential=(double*)argv[8].argv;
	double *virial=(double*)argv[9].argv;

	*potential=0.;
	sw_system(sws_1,pos,im,box,cutoff,potential,virial);
}

void swEnergySystemHamiltonian(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	sw *sws_2=(sw*)argv[1].argv;
	double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	double *box=(double*)argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	double *potential=(double*)argv[8].argv;
	double *virial=(double*)argv[9].argv;

	Potential_1=0.;
	Potential_2=0.;

	sw_system(sws_1,pos,im,box,cutoff,&Potential_1,virial);
	sw_system(sws_2,pos,im,box,cutoff,&Potential_2,virial);

	*potential=(1.-hlambda)*Potential_1+hlambda*Potential_2;
}

void swRejectParticleMove(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	/*sw *sws_2=(sw*)argv[1].argv;
	double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	int ncolloids=*((int*)argv[4].argv);
	double box[];
	box=(double[])argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	double *potential=(double*)argv[8].argv;
	double *virial=(double*)argv[9].argv;*/

	sw_particle_rejectmove(sws_1);
}

void swRejectSwapMove(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	/*sw *sws_2=(sw*)argv[1].argv;
	double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	int ncolloids=*((int*)argv[4].argv);
	double box[];
	box=(double[])argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	double *potential=(double*)argv[8].argv;
	double *virial=(double*)argv[9].argv;*/

	sw_particle_rejectmove(sws_1);
}

void swRejectSystemMove(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	/*sw *sws_2=(sw*)argv[1].argv;
	double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	int ncolloids=*((int*)argv[4].argv);
	double box[];
	box=(double[])argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	double *potential=(double*)argv[8].argv;
	double *virial=(double*)argv[9].argv;*/

	sw_system_rejectmove(sws_1);
}

void swRejectParticleMoveHamiltonian(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	sw *sws_2=(sw*)argv[1].argv;
	/*double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	int ncolloids=*((int*)argv[4].argv);
	double box[];
	box=(double[])argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	double *potential=(double*)argv[8].argv;
	double *virial=(double*)argv[9].argv;*/

	sw_particle_rejectmove(sws_1);
	sw_particle_rejectmove(sws_2);
}

void swRejectSystemMoveHamiltonian(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	sw *sws_2=(sw*)argv[1].argv;
	/*double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	int ncolloids=*((int*)argv[4].argv);
	double box[];
	box=(double[])argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	double *potential=(double*)argv[8].argv;
	double *virial=(double*)argv[9].argv;*/

	sw_system_rejectmove(sws_1);
	sw_system_rejectmove(sws_2);
}

void swRemoveParticle(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	//sw *sws_2=(sw*)argv[1].argv;
	//double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	double *box=(double*)argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	//double *potential=(double*)argv[8].argv;
	//double *virial=(double*)argv[9].argv;

	double potential,virial;

	sw_particle_remove(sws_1,pos,im,box,cutoff,&potential,&virial);

	eventSetIO(REMOVE_PARTICLE,REMOVE_PARTICLE__energy,(void*)&potential,sizeof(double));
	eventSetIO(REMOVE_PARTICLE,REMOVE_PARTICLE__virial,(void*)&virial,sizeof(double));
	//eventSetIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__virial,(void*)&virial,sizeof(double));

#ifdef DEBUG
	if (status!=0)
	{
		printf("Error: setting output with wrong type\n");
		exit(1);
	}
#endif

}

void swRemoveAccept(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	//sw *sws_2=(sw*)argv[1].argv;
	//double hlambda=*((double*)argv[2].argv);
	//vector *pos=(vector*)argv[3].argv;
	int ncolloids=*((int*)argv[4].argv);
	//double *box=(double*)argv[5].argv;
	//interactionmap *im=(interactionmap *)argv[6].argv;
	//double cutoff=*((double*)argv[7].argv);
	//double *potential=(double*)argv[8].argv;
	//double *virial=(double*)argv[9].argv;

	int movedparticle=*((int*)eventGetIO(REMOVE_PARTICLE,REMOVE_PARTICLE__particle));
	sw_remove_acceptmove(sws_1,movedparticle,ncolloids);

}

void swRemoveReject(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	//sw *sws_2=(sw*)argv[1].argv;
	//double hlambda=*((double*)argv[2].argv);
	//vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	//double *box=(double*)argv[5].argv;
	//interactionmap *im=(interactionmap *)argv[6].argv;
	//double cutoff=*((double*)argv[7].argv);
	//double *potential=(double*)argv[8].argv;
	//double *virial=(double*)argv[9].argv;

	sw_remove_rejectmove(sws_1);
}


void swAddParticle(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	//sw *sws_2=(sw*)argv[1].argv;
	//double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	double *box=(double*)argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	//double *potential=(double*)argv[8].argv;
	//double *virial=(double*)argv[9].argv;

	double potential,virial;

	sw_particle_add(sws_1,pos,im,box,cutoff,&potential,&virial);

	eventSetIO(ADD_PARTICLE,ADD_PARTICLE__energy,(void*)&potential,sizeof(double));
	eventSetIO(ADD_PARTICLE,ADD_PARTICLE__virial,(void*)&virial,sizeof(double));
}

void swAddReject(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	//sw *sws_2=(sw*)argv[1].argv;
	//double hlambda=*((double*)argv[2].argv);
	//vector *pos=(vector*)argv[3].argv;
	int ncolloids=*((int*)argv[4].argv);
	//double *box=(double*)argv[5].argv;
	//interactionmap *im=(interactionmap *)argv[6].argv;
	//double cutoff=*((double*)argv[7].argv);
	//double *potential=(double*)argv[8].argv;
	//double *virial=(double*)argv[9].argv;

	sw_add_rejectmove(sws_1,ncolloids);
}


void swEnergyParticleOld(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	//sw *sws_2=(sw*)argv[1].argv;
	//double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	double *box=(double*)argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	//double *potential=(double*)argv[8].argv;
	//double *virial=(double*)argv[9].argv;

	double potential,virial;

	sw_particle_oldconfiguration(sws_1,pos,im,box,cutoff,&potential,&virial);

	eventSetIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__energy,(void*)&potential,sizeof(double));
	eventSetIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__virial,(void*)&virial,sizeof(double));

#ifdef DEBUG
	if (status!=0)
	{
		printf("Error: setting output with wrong type\n");
		exit(1);
	}
#endif

}

void swEnergyParticleNew(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	//sw *sws_2=(sw*)argv[1].argv;
	//double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	double *box=(double*)argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	//double *potential=(double*)argv[8].argv;
	//double *virial=(double*)argv[9].argv;

	double potential,virial;

	sw_particle_newconfiguration(sws_1,pos,im,box,cutoff,&potential,&virial);

	eventSetIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__energy,(void*)&potential,sizeof(double));
	eventSetIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__virial,(void*)&virial,sizeof(double));

#ifdef DEBUG
	if (status!=0)
	{
		printf("Error: setting output with wrong type\n");
		exit(1);
	}
#endif


}

void swEnergySwapOld(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	//sw *sws_2=(sw*)argv[1].argv;
	//double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	double *box=(double*)argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	//double *potential=(double*)argv[8].argv;
	//double *virial=(double*)argv[9].argv;

	double potential,virial;

	sw_swap_oldconfiguration(sws_1,pos,im,box,cutoff,&potential,&virial);

	eventSetIO(ENERGY_SWAP_OLD,ENERGY_SWAP_OLD__energy,(void*)&potential,sizeof(double));
	eventSetIO(ENERGY_SWAP_OLD,ENERGY_SWAP_OLD__virial,(void*)&virial,sizeof(double));

#ifdef DEBUG
	if (status!=0)
	{
		printf("Error: setting output with wrong type\n");
		exit(1);
	}
#endif

}

void swEnergySwapNew(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	//sw *sws_2=(sw*)argv[1].argv;
	//double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	double *box=(double*)argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	//double *potential=(double*)argv[8].argv;
	//double *virial=(double*)argv[9].argv;

	double potential,virial;

	sw_swap_newconfiguration(sws_1,pos,im,box,cutoff,&potential,&virial);

	eventSetIO(ENERGY_SWAP_NEW,ENERGY_SWAP_NEW__energy,(void*)&potential,sizeof(double));
	eventSetIO(ENERGY_SWAP_NEW,ENERGY_SWAP_NEW__virial,(void*)&virial,sizeof(double));

#ifdef DEBUG
	if (status!=0)
	{
		printf("Error: setting output with wrong type\n");
		exit(1);
	}
#endif


}


void swEnergyParticleOldHamoltonian(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	sw *sws_2=(sw*)argv[1].argv;
	double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	double *box=(double*)argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	//double *potential=(double*)argv[8].argv;
	//double *virial=(double*)argv[9].argv;

	double potential_1=0.;
	double potential_2=0.;
	double virial=0.;

	sw_particle_oldconfiguration(sws_1,pos,im,box,cutoff,&potential_1,&virial);
	sw_particle_oldconfiguration(sws_2,pos,im,box,cutoff,&potential_2,&virial);

	double potential=(1.-hlambda)*potential_1+hlambda*potential_2;

	eventSetIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__energy,(void*)&potential,sizeof(double));
	eventSetIO(ENERGY_PARTICLE_OLD,ENERGY_PARTICLE_OLD__virial,(void*)&virial,sizeof(double));

#ifdef DEBUG
	if (status!=0)
	{
		printf("Error: setting output with wrong type\n");
		exit(1);
	}
#endif
}

void swEnergyParticleNewHamoltonian(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	sw *sws_2=(sw*)argv[1].argv;
	double hlambda=*((double*)argv[2].argv);
	vector *pos=(vector*)argv[3].argv;
	//int ncolloids=*((int*)argv[4].argv);
	double *box=(double*)argv[5].argv;
	interactionmap *im=(interactionmap *)argv[6].argv;
	double cutoff=*((double*)argv[7].argv);
	//double *potential=(double*)argv[8].argv;
	//double *virial=(double*)argv[9].argv;

	double potential_1=0.;
	double potential_2=0.;
	double virial=0.;

	sw_particle_newconfiguration(sws_1,pos,im,box,cutoff,&potential_1,&virial);
	sw_particle_newconfiguration(sws_2,pos,im,box,cutoff,&potential_2,&virial);

	double potential=(1.-hlambda)*potential_1+hlambda*potential_2;

	eventSetIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__energy,(void*)&potential,sizeof(double));
	eventSetIO(ENERGY_PARTICLE_NEW,ENERGY_PARTICLE_NEW__virial,(void*)&virial,sizeof(double));

#ifdef DEBUG
	if (status!=0)
	{
		printf("Error: setting output with wrong type\n");
		exit(1);
	}
#endif
}

void swSaveState(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;

	sw_save_state(sws_1);
}

void swSaveStateHamiltonian(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	sw *sws_2=(sw*)argv[1].argv;

	sw_save_state(sws_1);
	sw_save_state(sws_2);
}

void swLoadState(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;

	sw_load_state(sws_1);
}

void swLoadStateHamiltonian(args *argv,int argc)
{
	// lettura degli argomenti
	sw *sws_1=(sw*)argv[0].argv;
	sw *sws_2=(sw*)argv[1].argv;

	sw_load_state(sws_1);
	sw_load_state(sws_2);
}

void swPrintHamiltonian(steps t)
{
	if (Model==2)
	{
		//printf("Eccomi\n");
		fprintf(Hamiltonian_file,"%lld %lf %lf\n",t,Potential_1,Potential_2);
	}
}

void sw_check(vector *pos,interactionmap *ime,double Box[],double cutoff,double *potential,double *virial,double *twobody,double *threebody)
{
	sw *sws=Sws_1;

	int i,j,k;

	*twobody=0.;
	*threebody=0.;
	*virial=0.;


	for (i=0;i<ime->num;i++)
	{
		//int particlei=ime->who[i];

		for (j=0;j<ime->howmany[i];j++)
		{
			//int particlej=ime->with[i][j];

			// parte a due corpi
			double distij=sqrt(ime->rij2[i][j]);
			vector rij=ime->rij[i][j];

			double v2=sws->A*(sws->B/pow(distij,sws->p)-1.)*exp(1./(distij-sws->a));
			(*twobody)+=0.5*v2;


			// parte a tre corpi
			for (k=j+1;k<ime->howmany[i];k++)
			{
				//int particlek=ime->with[i][k];

				double distik=sqrt(ime->rij2[i][k]);
				vector rik=ime->rij[i][k];

				double costeta=(rij.x*rik.x+rij.y*rik.y+rij.z*rik.z)/(distij*distik);

				(*threebody)+=sws->lambda*exp(sws->Gamma/(distij-sws->a)+sws->Gamma/(distik-sws->a))*SQR(costeta-sws->costeta0[i]);
			}
		}
	}

	*potential=*twobody+*threebody;
}

void sw_check_force(sw *const sws,vector *pos,interactionmap *ime,double Box[],double cutoff,double *potential,double *virial,vector *force)
{
	int i,j,k;

	*potential=0.;
	*virial=0.;

	memset(force,0,ime->num*sizeof(vector));

	for (i=0;i<ime->num;i++)
	{
		int particlei=ime->who[i];

		for (j=0;j<ime->howmany[i];j++)
		{
			int particlej=ime->with[i][j];

			// parte a due corpi
			double distij=sqrt(ime->rij2[i][j]);
			vector rij=ime->rij[i][j];

			double distij4i=1./pow(distij,4);
			double distij5i=1./pow(distij,5);

			double v2=sws->A*(sws->B/pow(distij,sws->p)-1.)*exp(1./(distij-sws->a));
			(*potential)+=0.5*v2;

			double force_module=sws->A*exp(1./(distij-sws->a))*(4.*sws->B*distij5i+(sws->B*distij4i-1.)/(SQR(distij-sws->a)))/distij;

			force[particlei].x+=force_module*rij.x;
			force[particlei].y+=force_module*rij.y;
			force[particlei].z+=force_module*rij.z;


			//force[particlej].x-=force_module*rij.x;
			//force[particlej].y-=force_module*rij.y;
			//force[particlej].z-=force_module*rij.z;


			// parte a tre corpi
			for (k=j+1;k<ime->howmany[i];k++)
			{
				int particlek=ime->with[i][k];

				double distik=sqrt(ime->rij2[i][k]);
				vector rik=ime->rij[i][k];

				double costeta=(rij.x*rik.x+rij.y*rik.y+rij.z*rik.z)/(distij*distik);

				double pot3=sws->lambda*exp(sws->Gamma/(distij-sws->a)+sws->Gamma/(distik-sws->a))*SQR(costeta+1./3.);

				(*potential)+=pot3;

				double prefactor=2.*sws->lambda*exp(sws->Gamma/(distij-sws->a)+sws->Gamma/(distik-sws->a))*(costeta+1./3.);
				// contributi alla forza della particlej

				force_module=sws->Gamma*pot3/(distij*SQR(distij-sws->a));
				force_module+=prefactor*(costeta/SQR(distij));

				force[particlei].x+=force_module*rij.x;
				force[particlei].y+=force_module*rij.y;
				force[particlei].z+=force_module*rij.z;
				force[particlej].x-=force_module*rij.x;
				force[particlej].y-=force_module*rij.y;
				force[particlej].z-=force_module*rij.z;

				// contributi alla forza della particlek

				force_module=sws->Gamma*pot3/(distik*SQR(distik-sws->a));
				force_module+=prefactor*(costeta/SQR(distik));

				force[particlei].x+=force_module*rik.x;
				force[particlei].y+=force_module*rik.y;
				force[particlei].z+=force_module*rik.z;
				force[particlek].x-=force_module*rik.x;
				force[particlek].y-=force_module*rik.y;
				force[particlek].z-=force_module*rik.z;


				// termine incrociato
				force_module=-prefactor/(distij*distik);

				force[particlei].x+=force_module*rij.x;
				force[particlei].y+=force_module*rij.y;
				force[particlei].z+=force_module*rij.z;
				force[particlei].x+=force_module*rik.x;
				force[particlei].y+=force_module*rik.y;
				force[particlei].z+=force_module*rik.z;

				force[particlek].x-=force_module*rij.x;
				force[particlek].y-=force_module*rij.y;
				force[particlek].z-=force_module*rij.z;
				force[particlej].x-=force_module*rik.x;
				force[particlej].y-=force_module*rik.y;
				force[particlej].z-=force_module*rik.z;
			}

		}
	}
}

void sw_system(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial)
{
	int i,j;


	Potential_2body=0.;
	Potential_3body=0.;

	*virial=0.;

	memcpy(sws->TwoBody_copy[0],sws->TwoBody[0],(*sws->NColloids)*10*sizeof(double));
	memset(sws->TwoBody[0],0,(*sws->NColloids)*10*sizeof(double));


	// two body potential
	for (i=0;i<im->num;i++)
	{
		int particle1=im->who[i];

		for (j=0;j<im->howmany[i];j++)
		{

			int particle2=im->with[i][j];

			double dist=sqrt(im->rij2[i][j]);
			double inorm=1./dist;


			// vdist is the vector going from j to i
			vector vdist=im->rij[i][j];

			double gij=exp(sws->Gamma/(dist-sws->a));

			double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
			double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[particle2]));

			//Potential_2body+=v2+Uc;
			Potential_2body+=v2;
			Potential_3body+=Uc;

			sws->TwoBody[particle1][0]+=gij;
			sws->TwoBody[particle2][0]+=gij;

			sws->TwoBody[particle1][1]+=gij*vdist.x*inorm;
			sws->TwoBody[particle1][2]+=gij*vdist.y*inorm;
			sws->TwoBody[particle1][3]+=gij*vdist.z*inorm;

			sws->TwoBody[particle2][1]-=gij*vdist.x*inorm;
			sws->TwoBody[particle2][2]-=gij*vdist.y*inorm;
			sws->TwoBody[particle2][3]-=gij*vdist.z*inorm;

			double inorm2=SQR(inorm);

			sws->TwoBody[particle1][4]+=gij*SQR(vdist.x)*inorm2;
			sws->TwoBody[particle1][5]+=gij*SQR(vdist.y)*inorm2;
			sws->TwoBody[particle1][6]+=gij*SQR(vdist.z)*inorm2;
			sws->TwoBody[particle1][7]+=gij*vdist.x*vdist.y*inorm2;
			sws->TwoBody[particle1][8]+=gij*vdist.x*vdist.z*inorm2;
			sws->TwoBody[particle1][9]+=gij*vdist.y*vdist.z*inorm2;

			sws->TwoBody[particle2][4]+=gij*SQR(vdist.x)*inorm2;
			sws->TwoBody[particle2][5]+=gij*SQR(vdist.y)*inorm2;
			sws->TwoBody[particle2][6]+=gij*SQR(vdist.z)*inorm2;
			sws->TwoBody[particle2][7]+=gij*vdist.x*vdist.y*inorm2;
			sws->TwoBody[particle2][8]+=gij*vdist.x*vdist.z*inorm2;
			sws->TwoBody[particle2][9]+=gij*vdist.y*vdist.z*inorm2;
		}



	}

	for (i=0;i<im->num;i++)
	{
		int particle1=im->who[i];

		Potential_3body+=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));

	}

	*potential=Potential_2body+Potential_3body;

}

void swGetTwoThreeBodyParts(double *potential2body,double *potential3body)
{
	// questa funzione ha senso solo dopo la chiamata a sw_system
	*potential2body=Potential_2body;
	*potential3body=Potential_3body;
}

void sw_system_force(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial,vector *force)
{
	int i,j;
	double alpha=1./3.;

	*potential=0.;
	*virial=0.;
	memset(force,0,im->num*sizeof(vector));

	memcpy(sws->TwoBody_copy[0],sws->TwoBody[0],(*sws->NColloids)*10*sizeof(double));
	memset(sws->TwoBody[0],0,(*sws->NColloids)*10*sizeof(double));


	// two body potential
	for (i=0;i<im->num;i++)
	{
		int particle1=im->who[i];

		for (j=0;j<im->howmany[i];j++)
		{

			int particle2=im->with[i][j];

			double dist=sqrt(im->rij2[i][j]);
			double inorm=1./dist;
			double distij4i=1./pow(dist,4);
			double distij5i=1./pow(dist,5);


			// vdist is the vector going from j to i
			vector vdist=im->rij[i][j];
			double gij=exp(sws->Gamma/(dist-sws->a));

			double v2=sws->A*(sws->B/pow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
			double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[particle2]));

			(*potential)+=v2+Uc;

			// forza dalla parte a due corpi del potenziale
			double force_module=sws->A*exp(1./(dist-sws->a))*(4.*sws->B*distij5i+(sws->B*distij4i-1.)/(SQR(dist-sws->a)))/dist;

			force_module+=-2*sws->lambda*(16./9.)*(gij/dist)*sws->Gamma*exp(sws->Gamma/(dist-sws->a))/((dist-sws->a)*(dist-sws->a));

			force[particle1].x+=force_module*vdist.x;
			force[particle1].y+=force_module*vdist.y;
			force[particle1].z+=force_module*vdist.z;
			force[particle2].x-=force_module*vdist.x;
			force[particle2].y-=force_module*vdist.y;
			force[particle2].z-=force_module*vdist.z;



			sws->TwoBody[particle1][0]+=gij;
			sws->TwoBody[particle2][0]+=gij;

			sws->TwoBody[particle1][1]+=gij*vdist.x*inorm;
			sws->TwoBody[particle1][2]+=gij*vdist.y*inorm;
			sws->TwoBody[particle1][3]+=gij*vdist.z*inorm;

			sws->TwoBody[particle2][1]-=gij*vdist.x*inorm;
			sws->TwoBody[particle2][2]-=gij*vdist.y*inorm;
			sws->TwoBody[particle2][3]-=gij*vdist.z*inorm;

			double inorm2=SQR(inorm);

			sws->TwoBody[particle1][4]+=gij*SQR(vdist.x)*inorm2;
			sws->TwoBody[particle1][5]+=gij*SQR(vdist.y)*inorm2;
			sws->TwoBody[particle1][6]+=gij*SQR(vdist.z)*inorm2;
			sws->TwoBody[particle1][7]+=gij*vdist.x*vdist.y*inorm2;
			sws->TwoBody[particle1][8]+=gij*vdist.x*vdist.z*inorm2;
			sws->TwoBody[particle1][9]+=gij*vdist.y*vdist.z*inorm2;

			sws->TwoBody[particle2][4]+=gij*SQR(vdist.x)*inorm2;
			sws->TwoBody[particle2][5]+=gij*SQR(vdist.y)*inorm2;
			sws->TwoBody[particle2][6]+=gij*SQR(vdist.z)*inorm2;
			sws->TwoBody[particle2][7]+=gij*vdist.x*vdist.y*inorm2;
			sws->TwoBody[particle2][8]+=gij*vdist.x*vdist.z*inorm2;
			sws->TwoBody[particle2][9]+=gij*vdist.y*vdist.z*inorm2;
		}



	}

	for (i=0;i<im->num;i++)
	{
		int particle1=im->who[i];

		(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));


		for (j=0;j<im->howmany[i];j++)
		{
			int particle2=im->with[i][j];

			double dist=sqrt(im->rij2[i][j]);
			//double inorm=1./dist;

			// vdist is the vector going from j to i
			vector vdist=im->rij[i][j];
			double gij=exp(sws->Gamma/(dist-sws->a));

			double dgdr=-sws->Gamma*exp(sws->Gamma/(dist-sws->a))/((dist-sws->a)*(dist-sws->a));

			double cpc=sws->lambda*alpha*alpha*dgdr*(sws->TwoBody[particle1][0]+sws->TwoBody[particle2][0])/dist;

			cpc+=2.*sws->lambda*alpha*(dgdr-gij/dist)*(vdist.x*(sws->TwoBody[particle1][1]-sws->TwoBody[particle2][1])+vdist.y*(sws->TwoBody[particle1][2]-sws->TwoBody[particle2][2])+vdist.z*(sws->TwoBody[particle1][3]-sws->TwoBody[particle2][3]))/(dist*dist);

			double new_a=sws->TwoBody[particle1][4]+sws->TwoBody[particle2][4];
			double new_d=sws->TwoBody[particle1][5]+sws->TwoBody[particle2][5];
			double new_f=sws->TwoBody[particle1][6]+sws->TwoBody[particle2][6];
			double new_b=sws->TwoBody[particle1][7]+sws->TwoBody[particle2][7];
			double new_c=sws->TwoBody[particle1][8]+sws->TwoBody[particle2][8];
			double new_e=sws->TwoBody[particle1][9]+sws->TwoBody[particle2][9];

			double tensor_product=vdist.x*(new_a*vdist.x+new_b*vdist.y+new_c*vdist.z)+vdist.y*(new_b*vdist.x+new_d*vdist.y+new_e*vdist.z)+vdist.z*(new_c*vdist.x+new_e*vdist.y+new_f*vdist.z);

			cpc+=sws->lambda*(dgdr-2.*gij/dist)*tensor_product/(dist*dist*dist);

			vector first_component,second_component,third_component;

			first_component.x=cpc*vdist.x;
			first_component.y=cpc*vdist.y;
			first_component.z=cpc*vdist.z;

			second_component.x=2.*sws->lambda*alpha*gij*(sws->TwoBody[particle1][1]-sws->TwoBody[particle2][1])/dist;
			second_component.y=2.*sws->lambda*alpha*gij*(sws->TwoBody[particle1][2]-sws->TwoBody[particle2][2])/dist;
			second_component.z=2.*sws->lambda*alpha*gij*(sws->TwoBody[particle1][3]-sws->TwoBody[particle2][3])/dist;

			third_component.x=2.*sws->lambda*gij*(new_a*vdist.x+new_b*vdist.y+new_c*vdist.z)/(dist*dist);
			third_component.y=2.*sws->lambda*gij*(new_b*vdist.x+new_d*vdist.y+new_e*vdist.z)/(dist*dist);
			third_component.z=2.*sws->lambda*gij*(new_c*vdist.x+new_e*vdist.y+new_f*vdist.z)/(dist*dist);

			force[particle1].x-=first_component.x+second_component.x+third_component.x;
			force[particle1].y-=first_component.y+second_component.y+third_component.y;
			force[particle1].z-=first_component.z+second_component.z+third_component.z;

			force[particle2].x+=first_component.x+second_component.x+third_component.x;
			force[particle2].y+=first_component.y+second_component.y+third_component.y;
			force[particle2].z+=first_component.z+second_component.z+third_component.z;
		}

	}

}


void sw_system_virial(sw *const sws,vector *pos,interactionmap *im,double box[],double cutoff,double *potential,double *virial,vector *force)
{
	int i,j;
	double alpha=1./3.;

	*potential=0.;
	*virial=0.;

	// warning sws->TwoBody must be already set

	for (i=0;i<im->num;i++)
	{
		int particle1=im->who[i];

		// parte a tre corpi del potenziale
		(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));

		// portiamo particle1 nella cella centrale
		vector p1_rint;
		p1_rint.x=pos[particle1].x-rint(pos[particle1].x);
		p1_rint.y=pos[particle1].y-rint(pos[particle1].y);
		p1_rint.z=pos[particle1].z-rint(pos[particle1].z);

		// coordinate ortogonali
		vector p1;
		p1.x=box[0]*p1_rint.x+box[1]*p1_rint.y+box[2]*p1_rint.z;
		p1.y=box[3]*p1_rint.y+box[4]*p1_rint.z;
		p1.z=box[5]*p1_rint.z;

		for (j=0;j<im->howmany[i];j++)
		{
			int particle2=im->with[i][j];

			vector p2_pbc;
			p2_pbc.x=pos[particle2].x;
			p2_pbc.y=pos[particle2].y;
			p2_pbc.z=pos[particle2].z;

			pbcNearestImage(&p2_pbc,&p1_rint);

			vector p2;
			p2.x=box[0]*p2_pbc.x+box[1]*p2_pbc.y+box[2]*p2_pbc.z;
			p2.y=box[3]*p2_pbc.y+box[4]*p2_pbc.z;
			p2.z=box[5]*p2_pbc.z;

			double dist=sqrt(im->rij2[i][j]);
			double distij4i=1./pow(dist,4);
			double distij5i=1./pow(dist,5);

			// vdist is the vector going from j to i
			vector vdist=im->rij[i][j];
			double gij=exp(sws->Gamma/(dist-sws->a));


			// PARTE A DUE CORPI
			double v2=sws->A*(sws->B/pow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
			double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[particle2]));

			// parte a due corpi del potenziale
			(*potential)+=v2+Uc;

			double force_module=sws->A*exp(1./(dist-sws->a))*(4.*sws->B*distij5i+(sws->B*distij4i-1.)/(SQR(dist-sws->a)))/dist;
			force_module+=-2*sws->lambda*(16./9.)*(gij/dist)*sws->Gamma*exp(sws->Gamma/(dist-sws->a))/((dist-sws->a)*(dist-sws->a));

			vector fg;
			fg.x=force_module*vdist.x;
			fg.y=force_module*vdist.y;
			fg.z=force_module*vdist.z;

			(*virial)+=fg.x*p1.x+fg.y*p1.y+fg.z*p1.z;
			(*virial)-=fg.x*p2.x+fg.y*p2.y+fg.z*p2.z;



			// PARTE A TRE CORPI

			double dgdr=-sws->Gamma*exp(sws->Gamma/(dist-sws->a))/((dist-sws->a)*(dist-sws->a));

			double cpc=sws->lambda*alpha*alpha*dgdr*(sws->TwoBody[particle1][0]+sws->TwoBody[particle2][0])/dist;

			cpc+=2.*sws->lambda*alpha*(dgdr-gij/dist)*(vdist.x*(sws->TwoBody[particle1][1]-sws->TwoBody[particle2][1])+vdist.y*(sws->TwoBody[particle1][2]-sws->TwoBody[particle2][2])+vdist.z*(sws->TwoBody[particle1][3]-sws->TwoBody[particle2][3]))/(dist*dist);

			double new_a=sws->TwoBody[particle1][4]+sws->TwoBody[particle2][4];
			double new_d=sws->TwoBody[particle1][5]+sws->TwoBody[particle2][5];
			double new_f=sws->TwoBody[particle1][6]+sws->TwoBody[particle2][6];
			double new_b=sws->TwoBody[particle1][7]+sws->TwoBody[particle2][7];
			double new_c=sws->TwoBody[particle1][8]+sws->TwoBody[particle2][8];
			double new_e=sws->TwoBody[particle1][9]+sws->TwoBody[particle2][9];

			double tensor_product=vdist.x*(new_a*vdist.x+new_b*vdist.y+new_c*vdist.z)+vdist.y*(new_b*vdist.x+new_d*vdist.y+new_e*vdist.z)+vdist.z*(new_c*vdist.x+new_e*vdist.y+new_f*vdist.z);

			cpc+=sws->lambda*(dgdr-2.*gij/dist)*tensor_product/(dist*dist*dist);

			vector first_component,second_component,third_component;

			first_component.x=cpc*vdist.x;
			first_component.y=cpc*vdist.y;
			first_component.z=cpc*vdist.z;

			second_component.x=2.*sws->lambda*alpha*gij*(sws->TwoBody[particle1][1]-sws->TwoBody[particle2][1])/dist;
			second_component.y=2.*sws->lambda*alpha*gij*(sws->TwoBody[particle1][2]-sws->TwoBody[particle2][2])/dist;
			second_component.z=2.*sws->lambda*alpha*gij*(sws->TwoBody[particle1][3]-sws->TwoBody[particle2][3])/dist;

			third_component.x=2.*sws->lambda*gij*(new_a*vdist.x+new_b*vdist.y+new_c*vdist.z)/(dist*dist);
			third_component.y=2.*sws->lambda*gij*(new_b*vdist.x+new_d*vdist.y+new_e*vdist.z)/(dist*dist);
			third_component.z=2.*sws->lambda*gij*(new_c*vdist.x+new_e*vdist.y+new_f*vdist.z)/(dist*dist);

			fg.x=first_component.x+second_component.x+third_component.x;
			fg.y=first_component.y+second_component.y+third_component.y;
			fg.z=first_component.z+second_component.z+third_component.z;

			(*virial)-=fg.x*p1.x+fg.y*p1.y+fg.z*p1.z;
			(*virial)+=fg.x*p2.x+fg.y*p2.y+fg.z*p2.z;
		}

	}

}


double sw_particle_ghost_cheat(interactionmap *im)
{
	sw *sws=Sws_1;

	int j;

	double potential=0.;


	int particle1=im->who[0];

	// tensore
	double twobody_particle1[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	double twobody_particle2[10];

	for (j=0;j<im->howmany[0];j++)
	{
		int particle2=im->with[0][j];


		double dist=sqrt(im->rij2[0][j]);
		double inorm=1./dist;
		double inorm2=SQR(inorm);

		// vdist is the vector going from j to i
		vector vdist=im->rij[0][j];

		double gij=exp(sws->Gamma/(dist-sws->a));

		double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
		double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[particle2]));

		// constributo a due corpi di particle1
		potential+=v2+Uc;


		// sottraiamo il vecchio tensore dei vicinis
		potential-=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9]));
		sws->Num_copies++;

		// incrementiamo il tensore della particella1
		twobody_particle1[0]+=gij;
		twobody_particle1[1]+=gij*vdist.x*inorm;
		twobody_particle1[2]+=gij*vdist.y*inorm;
		twobody_particle1[3]+=gij*vdist.z*inorm;
		twobody_particle1[4]+=gij*SQR(vdist.x)*inorm2;
		twobody_particle1[5]+=gij*SQR(vdist.y)*inorm2;
		twobody_particle1[6]+=gij*SQR(vdist.z)*inorm2;
		twobody_particle1[7]+=gij*vdist.x*vdist.y*inorm2;
		twobody_particle1[8]+=gij*vdist.x*vdist.z*inorm2;
		twobody_particle1[9]+=gij*vdist.y*vdist.z*inorm2;

		// aggiorniamo il tensore dei vicini
		memcpy(twobody_particle2,sws->TwoBody[particle2],10*sizeof(double));


		twobody_particle2[0]+=gij;
		twobody_particle2[1]-=gij*vdist.x*inorm;
		twobody_particle2[2]-=gij*vdist.y*inorm;
		twobody_particle2[3]-=gij*vdist.z*inorm;
		twobody_particle2[4]+=gij*SQR(vdist.x)*inorm2;
		twobody_particle2[5]+=gij*SQR(vdist.y)*inorm2;
		twobody_particle2[6]+=gij*SQR(vdist.z)*inorm2;
		twobody_particle2[7]+=gij*vdist.x*vdist.y*inorm2;
		twobody_particle2[8]+=gij*vdist.x*vdist.z*inorm2;
		twobody_particle2[9]+=gij*vdist.y*vdist.z*inorm2;

		//nuovo contributo
		potential+=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(twobody_particle2[0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(twobody_particle2[1])+SQR(twobody_particle2[2])+SQR(twobody_particle2[3]))+0.5*sws->lambda*(SQR(twobody_particle2[4])+SQR(twobody_particle2[5])+SQR(twobody_particle2[6])+2.*SQR(twobody_particle2[7])+2.*SQR(twobody_particle2[8])+2.*SQR(twobody_particle2[9]));

	}

	// nuovo contributo particella1
	potential+=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(twobody_particle1[0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(twobody_particle1[1])+SQR(twobody_particle1[2])+SQR(twobody_particle1[3]))+0.5*sws->lambda*(SQR(twobody_particle1[4])+SQR(twobody_particle1[5])+SQR(twobody_particle1[6])+2.*SQR(twobody_particle1[7])+2.*SQR(twobody_particle1[8])+2.*SQR(twobody_particle1[9]));

	return potential;
}

void sw_particle_add(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *pot,double *virial)
{
	int j;

	double potential=0.;

	int particle1=im->who[0];

	sws->Num_copies=0;


	for (j=0;j<im->howmany[0];j++)
	{
		int particle2=im->with[0][j];


		double dist=sqrt(im->rij2[0][j]);
		double inorm=1./dist;
		double inorm2=SQR(inorm);

		// vdist is the vector going from j to i
		vector vdist=im->rij[0][j];

		double gij=exp(sws->Gamma/(dist-sws->a));

		double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
		double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[particle2]));

		// constributo a due corpi di particle1
		potential+=v2+Uc;


		// sottraiamo il vecchio tensore dei vicini
		potential-=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9]));

		// incrementiamo il tensore della particella1
		sws->TwoBody[particle1][0]+=gij;
		sws->TwoBody[particle1][1]+=gij*vdist.x*inorm;
		sws->TwoBody[particle1][2]+=gij*vdist.y*inorm;
		sws->TwoBody[particle1][3]+=gij*vdist.z*inorm;
		sws->TwoBody[particle1][4]+=gij*SQR(vdist.x)*inorm2;
		sws->TwoBody[particle1][5]+=gij*SQR(vdist.y)*inorm2;
		sws->TwoBody[particle1][6]+=gij*SQR(vdist.z)*inorm2;
		sws->TwoBody[particle1][7]+=gij*vdist.x*vdist.y*inorm2;
		sws->TwoBody[particle1][8]+=gij*vdist.x*vdist.z*inorm2;
		sws->TwoBody[particle1][9]+=gij*vdist.y*vdist.z*inorm2;


		memcpy(sws->TwoBody_copy[sws->Num_copies],sws->TwoBody[particle2],10*sizeof(double));
		sws->WhoCopied[sws->Num_copies]=particle2;
		sws->Num_copies++;


		sws->TwoBody[particle2][0]+=gij;
		sws->TwoBody[particle2][1]-=gij*vdist.x*inorm;
		sws->TwoBody[particle2][2]-=gij*vdist.y*inorm;
		sws->TwoBody[particle2][3]-=gij*vdist.z*inorm;
		sws->TwoBody[particle2][4]+=gij*SQR(vdist.x)*inorm2;
		sws->TwoBody[particle2][5]+=gij*SQR(vdist.y)*inorm2;
		sws->TwoBody[particle2][6]+=gij*SQR(vdist.z)*inorm2;
		sws->TwoBody[particle2][7]+=gij*vdist.x*vdist.y*inorm2;
		sws->TwoBody[particle2][8]+=gij*vdist.x*vdist.z*inorm2;
		sws->TwoBody[particle2][9]+=gij*vdist.y*vdist.z*inorm2;

		//nuovo contributo
		potential+=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9]));

	}

	// nuovo contributo particella1
	potential+=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));

	*pot=potential;
}

void sw_particle_remove(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial)
{
	// l'idea e' di fare energia con particle 1 - energia senza particle 2

	int j;

	*potential=0.;
	*virial=0.;

	sws->Num_copies=0;

	int particle1=im->who[0];

	// si puo' evitare di fare la copia della particella centrale
	// se la mossa viene accettata il tensore viene resettato
	// altrimenti resta uguale
	//memcpy(sws->TwoBody_copy[sws->Num_copies],sws->TwoBody[particle1],10*sizeof(double));
	//sws->WhoCopied[sws->Num_copies]=particle1;

	// parte tensoriale di particle1
	(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));


	for (j=0;j<im->howmany[0];j++)
	{
		int particle2=im->with[0][j];

		double dist=sqrt(im->rij2[0][j]);
		double inorm=1./dist;


		// vdist is the vector going from j to i
		vector vdist=im->rij[0][j];
		double gij=exp(sws->Gamma/(dist-sws->a));

		double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
		double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[particle2]));

		// parte scalare di particle1-particle2
		(*potential)+=v2+Uc;

		memcpy(sws->TwoBody_copy[sws->Num_copies],sws->TwoBody[particle2],10*sizeof(double));
		sws->WhoCopied[sws->Num_copies]=particle2;
		sws->Num_copies++;

		// parte tensoriale di particle2
		(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9]));


		/* IMPORTANTE */
		// Invertiamo il segno di gij in modo tale da sottrarre il contributo della particella particle1
		// a tutte le matrici
		gij*=-1.;

		sws->TwoBody[particle2][0]+=gij;

		sws->TwoBody[particle2][1]-=gij*vdist.x*inorm;
		sws->TwoBody[particle2][2]-=gij*vdist.y*inorm;
		sws->TwoBody[particle2][3]-=gij*vdist.z*inorm;

		double inorm2=SQR(inorm);

		sws->TwoBody[particle2][4]+=gij*SQR(vdist.x)*inorm2;
		sws->TwoBody[particle2][5]+=gij*SQR(vdist.y)*inorm2;
		sws->TwoBody[particle2][6]+=gij*SQR(vdist.z)*inorm2;
		sws->TwoBody[particle2][7]+=gij*vdist.x*vdist.y*inorm2;
		sws->TwoBody[particle2][8]+=gij*vdist.x*vdist.z*inorm2;
		sws->TwoBody[particle2][9]+=gij*vdist.y*vdist.z*inorm2;

		// sottraiamo la parte tensoriale di particle 2 in assenza di particle 1
		// parte tensoriale di particle2
		(*potential)-=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9]));
	}

	// da resettare se la mossa viene accettata
	//memset(sws->TwoBody[particle1],0,10*sizeof(double));
}


void sw_particle_oldconfiguration(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial)
{
	int j;

	*potential=0.;
	*virial=0.;

	sws->Num_copies=0;
	//bilistaReset(sws->List_modified,sws->NColloids); // dovrebbe essere sempre vuota a questo punto

	int particle1=im->who[0];

	memcpy(sws->TwoBody_copy[sws->Num_copies],sws->TwoBody[particle1],10*sizeof(double));
	sws->WhoCopied[sws->Num_copies]=particle1;
	(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));
	sws->Num_copies++;


	for (j=0;j<im->howmany[0];j++)
	{
		int particle2=im->with[0][j];

		double dist=sqrt(im->rij2[0][j]);
		double inorm=1./dist;


		// vdist is the vector going from j to i
		vector vdist=im->rij[0][j];
		double gij=exp(sws->Gamma/(dist-sws->a));

		double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
		double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[particle2]));

		(*potential)+=v2+Uc;

		memcpy(sws->TwoBody_copy[sws->Num_copies],sws->TwoBody[particle2],10*sizeof(double));
		sws->WhoCopied[sws->Num_copies]=particle2;
		bilistaInsert(sws->List_modified,particle2);
		(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9]));
		sws->Num_copies++;

		/* IMPORTANTE */
		// Invertiamo il segno di gij in modo tale da sottrarre il contributo della particella particle1
		// a tutte le matrici
		gij*=-1.;

		sws->TwoBody[particle2][0]+=gij;

		sws->TwoBody[particle2][1]-=gij*vdist.x*inorm;
		sws->TwoBody[particle2][2]-=gij*vdist.y*inorm;
		sws->TwoBody[particle2][3]-=gij*vdist.z*inorm;

		double inorm2=SQR(inorm);

		sws->TwoBody[particle2][4]+=gij*SQR(vdist.x)*inorm2;
		sws->TwoBody[particle2][5]+=gij*SQR(vdist.y)*inorm2;
		sws->TwoBody[particle2][6]+=gij*SQR(vdist.z)*inorm2;
		sws->TwoBody[particle2][7]+=gij*vdist.x*vdist.y*inorm2;
		sws->TwoBody[particle2][8]+=gij*vdist.x*vdist.z*inorm2;
		sws->TwoBody[particle2][9]+=gij*vdist.y*vdist.z*inorm2;
	}

	memset(sws->TwoBody[particle1],0,10*sizeof(double));
}

void sw_particle_newconfiguration(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial)
{
	int j;

	*potential=0.;
	*virial=0.;


	// two body potential
	int particle1=im->who[0];

	for (j=0;j<im->howmany[0];j++)
	{

		int particle2=im->with[0][j];

		double dist=sqrt(im->rij2[0][j]);
		double inorm=1./dist;


		// vdist is the vector going from j to i
		vector vdist=im->rij[0][j];

		double gij=exp(sws->Gamma/(dist-sws->a));

		double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
		double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[particle2]));

		(*potential)+=v2+Uc;

		double inorm2=SQR(inorm);

		if (bilistaIsIn(sws->List_modified,particle2)==0)
		{
			// bisogna sottrarre il vecchio contributo
			(*potential)-=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9]));

			// aggiungiamolo alla lista
			memcpy(sws->TwoBody_copy[sws->Num_copies],sws->TwoBody[particle2],10*sizeof(double));
			sws->WhoCopied[sws->Num_copies]=particle2;
			sws->Num_copies++;
		}
		else
		{
			bilistaRemove(sws->List_modified,particle2);
		}

		sws->TwoBody[particle1][0]+=gij;
		sws->TwoBody[particle1][1]+=gij*vdist.x*inorm;
		sws->TwoBody[particle1][2]+=gij*vdist.y*inorm;
		sws->TwoBody[particle1][3]+=gij*vdist.z*inorm;
		sws->TwoBody[particle1][4]+=gij*SQR(vdist.x)*inorm2;
		sws->TwoBody[particle1][5]+=gij*SQR(vdist.y)*inorm2;
		sws->TwoBody[particle1][6]+=gij*SQR(vdist.z)*inorm2;
		sws->TwoBody[particle1][7]+=gij*vdist.x*vdist.y*inorm2;
		sws->TwoBody[particle1][8]+=gij*vdist.x*vdist.z*inorm2;
		sws->TwoBody[particle1][9]+=gij*vdist.y*vdist.z*inorm2;


		sws->TwoBody[particle2][0]+=gij;
		sws->TwoBody[particle2][1]-=gij*vdist.x*inorm;
		sws->TwoBody[particle2][2]-=gij*vdist.y*inorm;
		sws->TwoBody[particle2][3]-=gij*vdist.z*inorm;
		sws->TwoBody[particle2][4]+=gij*SQR(vdist.x)*inorm2;
		sws->TwoBody[particle2][5]+=gij*SQR(vdist.y)*inorm2;
		sws->TwoBody[particle2][6]+=gij*SQR(vdist.z)*inorm2;
		sws->TwoBody[particle2][7]+=gij*vdist.x*vdist.y*inorm2;
		sws->TwoBody[particle2][8]+=gij*vdist.x*vdist.z*inorm2;
		sws->TwoBody[particle2][9]+=gij*vdist.y*vdist.z*inorm2;

		//nuovo contributo
		(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9]));
	}

	// nuovo contributo particella1
	(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));


	// ora in bilista sws->List_modified sono rimaste le particelle che dopo la mossa della particella non interagiscono piu' con essa
	int particle2;
	while (bilistaPop(sws->List_modified,&particle2))
	{
		(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9]));
	}

}

void sw_swap_oldconfiguration(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial)
{
	int j;

	*potential=0.;
	*virial=0.;

	sws->Num_copies=0;
	//bilistaReset(sws->List_modified,sws->NColloids); // dovrebbe essere sempre vuota a questo punto

	int particle1=im->who[0];
	int particle2=im->who[1];

	// parte tensoriale di particle1
	memcpy(sws->TwoBody_copy[sws->Num_copies],sws->TwoBody[particle1],10*sizeof(double));
	sws->WhoCopied[sws->Num_copies]=particle1;
	(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));
	sws->Num_copies++;

	// parte tensoriale di particle2
	memcpy(sws->TwoBody_copy[sws->Num_copies],sws->TwoBody[particle2],10*sizeof(double));
	sws->WhoCopied[sws->Num_copies]=particle2;
	(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9]));
	sws->Num_copies++;


	for (j=0;j<im->howmany[0];j++)
	{
		int neighbour=im->with[0][j];

		// assumiamo che la parte 2 body tra particle1 e particle2 rimane invariata dopo la swap move
		if (neighbour!=particle2)
		{
			double dist=sqrt(im->rij2[0][j]);
			double inorm=1./dist;


			// vdist is the vector going from j to i
			vector vdist=im->rij[0][j];
			double gij=exp(sws->Gamma/(dist-sws->a));

			double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
			double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[neighbour]));

			(*potential)+=v2+Uc;

			memcpy(sws->TwoBody_copy[sws->Num_copies],sws->TwoBody[neighbour],10*sizeof(double));
			sws->WhoCopied[sws->Num_copies]=neighbour;
			bilistaInsert(sws->List_modified,neighbour);
			(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[neighbour])*SQR(sws->TwoBody[neighbour][0])+sws->lambda*(-sws->costeta0[neighbour])*(SQR(sws->TwoBody[neighbour][1])+SQR(sws->TwoBody[neighbour][2])+SQR(sws->TwoBody[neighbour][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[neighbour][4])+SQR(sws->TwoBody[neighbour][5])+SQR(sws->TwoBody[neighbour][6])+2.*SQR(sws->TwoBody[neighbour][7])+2.*SQR(sws->TwoBody[neighbour][8])+2.*SQR(sws->TwoBody[neighbour][9]));
			sws->Num_copies++;

			/* IMPORTANTE */
			// Invertiamo il segno di gij in modo tale da sottrarre il contributo della particella particle1
			// a tutte le matrici
			gij*=-1.;

			sws->TwoBody[neighbour][0]+=gij;

			sws->TwoBody[neighbour][1]-=gij*vdist.x*inorm;
			sws->TwoBody[neighbour][2]-=gij*vdist.y*inorm;
			sws->TwoBody[neighbour][3]-=gij*vdist.z*inorm;

			double inorm2=SQR(inorm);

			sws->TwoBody[neighbour][4]+=gij*SQR(vdist.x)*inorm2;
			sws->TwoBody[neighbour][5]+=gij*SQR(vdist.y)*inorm2;
			sws->TwoBody[neighbour][6]+=gij*SQR(vdist.z)*inorm2;
			sws->TwoBody[neighbour][7]+=gij*vdist.x*vdist.y*inorm2;
			sws->TwoBody[neighbour][8]+=gij*vdist.x*vdist.z*inorm2;
			sws->TwoBody[neighbour][9]+=gij*vdist.y*vdist.z*inorm2;
		}
	}

	for (j=0;j<im->howmany[1];j++)
	{
		int neighbour=im->with[1][j];

		// assumiamo che la parte 2 body tra particle1 e particle2 rimane invariata dopo la swap move
		if (neighbour!=particle1)
		{
			double dist=sqrt(im->rij2[1][j]);
			double inorm=1./dist;


			// vdist is the vector going from j to i
			vector vdist=im->rij[1][j];
			double gij=exp(sws->Gamma/(dist-sws->a));

			double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
			double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle2])+SQR(1.-sws->costeta0[neighbour]));

			(*potential)+=v2+Uc;

			// se neighbour era anche vicino di particle1 non va ricalcolata la parte a tre corpi
			if (bilistaIsIn(sws->List_modified,neighbour)==0)
			{
				memcpy(sws->TwoBody_copy[sws->Num_copies],sws->TwoBody[neighbour],10*sizeof(double));
				sws->WhoCopied[sws->Num_copies]=neighbour;
				bilistaInsert(sws->List_modified,neighbour);
				(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[neighbour])*SQR(sws->TwoBody[neighbour][0])+sws->lambda*(-sws->costeta0[neighbour])*(SQR(sws->TwoBody[neighbour][1])+SQR(sws->TwoBody[neighbour][2])+SQR(sws->TwoBody[neighbour][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[neighbour][4])+SQR(sws->TwoBody[neighbour][5])+SQR(sws->TwoBody[neighbour][6])+2.*SQR(sws->TwoBody[neighbour][7])+2.*SQR(sws->TwoBody[neighbour][8])+2.*SQR(sws->TwoBody[neighbour][9]));
				sws->Num_copies++;
			}

			/* IMPORTANTE */
			// Invertiamo il segno di gij in modo tale da sottrarre il contributo della particella particle1
			// a tutte le matrici
			gij*=-1.;

			sws->TwoBody[neighbour][0]+=gij;

			sws->TwoBody[neighbour][1]-=gij*vdist.x*inorm;
			sws->TwoBody[neighbour][2]-=gij*vdist.y*inorm;
			sws->TwoBody[neighbour][3]-=gij*vdist.z*inorm;

			double inorm2=SQR(inorm);

			sws->TwoBody[neighbour][4]+=gij*SQR(vdist.x)*inorm2;
			sws->TwoBody[neighbour][5]+=gij*SQR(vdist.y)*inorm2;
			sws->TwoBody[neighbour][6]+=gij*SQR(vdist.z)*inorm2;
			sws->TwoBody[neighbour][7]+=gij*vdist.x*vdist.y*inorm2;
			sws->TwoBody[neighbour][8]+=gij*vdist.x*vdist.z*inorm2;
			sws->TwoBody[neighbour][9]+=gij*vdist.y*vdist.z*inorm2;
		}
	}

	memset(sws->TwoBody[particle1],0,10*sizeof(double));
	memset(sws->TwoBody[particle2],0,10*sizeof(double));
}

void sw_swap_newconfiguration(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial)
{
	int j;

	*potential=0.;
	*virial=0.;


	// two body potential
	int particle1=im->who[0];
	int particle2=im->who[1];

	for (j=0;j<im->howmany[0];j++)
	{

		int neighbour=im->with[0][j];



		double dist=sqrt(im->rij2[0][j]);
		double inorm=1./dist;


		// vdist is the vector going from j to i
		vector vdist=im->rij[0][j];

		double gij=exp(sws->Gamma/(dist-sws->a));

		double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
		double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[neighbour]));

		// la parte 2-body tra particle1 e particle2 non e' stata considerata nella old move
		if (neighbour!=particle2)
			(*potential)+=v2+Uc;

		double inorm2=SQR(inorm);


		sws->TwoBody[particle1][0]+=gij;
		sws->TwoBody[particle1][1]+=gij*vdist.x*inorm;
		sws->TwoBody[particle1][2]+=gij*vdist.y*inorm;
		sws->TwoBody[particle1][3]+=gij*vdist.z*inorm;
		sws->TwoBody[particle1][4]+=gij*SQR(vdist.x)*inorm2;
		sws->TwoBody[particle1][5]+=gij*SQR(vdist.y)*inorm2;
		sws->TwoBody[particle1][6]+=gij*SQR(vdist.z)*inorm2;
		sws->TwoBody[particle1][7]+=gij*vdist.x*vdist.y*inorm2;
		sws->TwoBody[particle1][8]+=gij*vdist.x*vdist.z*inorm2;
		sws->TwoBody[particle1][9]+=gij*vdist.y*vdist.z*inorm2;


		sws->TwoBody[neighbour][0]+=gij;
		sws->TwoBody[neighbour][1]-=gij*vdist.x*inorm;
		sws->TwoBody[neighbour][2]-=gij*vdist.y*inorm;
		sws->TwoBody[neighbour][3]-=gij*vdist.z*inorm;
		sws->TwoBody[neighbour][4]+=gij*SQR(vdist.x)*inorm2;
		sws->TwoBody[neighbour][5]+=gij*SQR(vdist.y)*inorm2;
		sws->TwoBody[neighbour][6]+=gij*SQR(vdist.z)*inorm2;
		sws->TwoBody[neighbour][7]+=gij*vdist.x*vdist.y*inorm2;
		sws->TwoBody[neighbour][8]+=gij*vdist.x*vdist.z*inorm2;
		sws->TwoBody[neighbour][9]+=gij*vdist.y*vdist.z*inorm2;
	}

	for (j=0;j<im->howmany[1];j++)
	{

		int neighbour=im->with[1][j];

		// se il vicino e' particle1 la parte a tre corpi e' gia' stata aggiornata al loop precedente
		if (neighbour!=particle1)
		{

			double dist=sqrt(im->rij2[1][j]);
			double inorm=1./dist;


			// vdist is the vector going from j to i
			vector vdist=im->rij[1][j];

			double gij=exp(sws->Gamma/(dist-sws->a));

			double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
			double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle2])+SQR(1.-sws->costeta0[neighbour]));

			(*potential)+=v2+Uc;

			double inorm2=SQR(inorm);


			sws->TwoBody[particle2][0]+=gij;
			sws->TwoBody[particle2][1]+=gij*vdist.x*inorm;
			sws->TwoBody[particle2][2]+=gij*vdist.y*inorm;
			sws->TwoBody[particle2][3]+=gij*vdist.z*inorm;
			sws->TwoBody[particle2][4]+=gij*SQR(vdist.x)*inorm2;
			sws->TwoBody[particle2][5]+=gij*SQR(vdist.y)*inorm2;
			sws->TwoBody[particle2][6]+=gij*SQR(vdist.z)*inorm2;
			sws->TwoBody[particle2][7]+=gij*vdist.x*vdist.y*inorm2;
			sws->TwoBody[particle2][8]+=gij*vdist.x*vdist.z*inorm2;
			sws->TwoBody[particle2][9]+=gij*vdist.y*vdist.z*inorm2;


			sws->TwoBody[neighbour][0]+=gij;
			sws->TwoBody[neighbour][1]-=gij*vdist.x*inorm;
			sws->TwoBody[neighbour][2]-=gij*vdist.y*inorm;
			sws->TwoBody[neighbour][3]-=gij*vdist.z*inorm;
			sws->TwoBody[neighbour][4]+=gij*SQR(vdist.x)*inorm2;
			sws->TwoBody[neighbour][5]+=gij*SQR(vdist.y)*inorm2;
			sws->TwoBody[neighbour][6]+=gij*SQR(vdist.z)*inorm2;
			sws->TwoBody[neighbour][7]+=gij*vdist.x*vdist.y*inorm2;
			sws->TwoBody[neighbour][8]+=gij*vdist.x*vdist.z*inorm2;
			sws->TwoBody[neighbour][9]+=gij*vdist.y*vdist.z*inorm2;
		}
	}


	// nuovo contributo particella1
	(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));

	// nuovo contributo particella2
	(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9]));


	// ora in bilista sws->List_modified sono rimasti i vicini di entrambe le particelle (senza ripetizioni)
	int neighbour;
	while (bilistaPop(sws->List_modified,&neighbour))
	{
		(*potential)+=0.5*sws->lambda*SQR(sws->costeta0[neighbour])*SQR(sws->TwoBody[neighbour][0])+sws->lambda*(-sws->costeta0[neighbour])*(SQR(sws->TwoBody[neighbour][1])+SQR(sws->TwoBody[neighbour][2])+SQR(sws->TwoBody[neighbour][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[neighbour][4])+SQR(sws->TwoBody[neighbour][5])+SQR(sws->TwoBody[neighbour][6])+2.*SQR(sws->TwoBody[neighbour][7])+2.*SQR(sws->TwoBody[neighbour][8])+2.*SQR(sws->TwoBody[neighbour][9]));
	}

}


void sw_remove_rejectmove(sw *const sws)
{
	int i;

	for (i=0;i<sws->Num_copies;i++)
	{
		int particle=sws->WhoCopied[i];

		memcpy(sws->TwoBody[particle],sws->TwoBody_copy[i],10*sizeof(double));
	}
}

void sw_remove_acceptmove(sw *const sws,int particle,int ncolloids)
{

	if (particle<ncolloids-1)
	{
		memcpy(sws->TwoBody[particle],sws->TwoBody_copy[ncolloids-1],10*sizeof(double));
	}
	memset(sws->TwoBody[ncolloids-1],0,10*sizeof(double));
}

void sw_add_rejectmove(sw *const sws,int ncolloids)
{
	int i;

	for (i=0;i<sws->Num_copies;i++)
	{
		int particle=sws->WhoCopied[i];

		memcpy(sws->TwoBody[particle],sws->TwoBody_copy[i],10*sizeof(double));
	}

	memset(sws->TwoBody[ncolloids],0,10*sizeof(double));
}



void sw_particle_rejectmove(sw *const sws)
{
	int i;

	for (i=0;i<sws->Num_copies;i++)
	{
		int particle=sws->WhoCopied[i];

		memcpy(sws->TwoBody[particle],sws->TwoBody_copy[i],10*sizeof(double));
	}
}

void sw_system_rejectmove(sw *const sws)
{
	double **buffer;

	buffer=sws->TwoBody;
	sws->TwoBody=sws->TwoBody_copy;
	sws->TwoBody_copy=buffer;
}

void sw_save_state(sw *const sws)
{
	memcpy(sws->TwoBody_savestate[0],sws->TwoBody[0],(*sws->NColloids)*10*sizeof(double));
}

void sw_load_state(sw *const sws)
{
	memcpy(sws->TwoBody[0],sws->TwoBody_savestate[0],(*sws->NColloids)*10*sizeof(double));
}

void getSW(FILE *pfile)
{
	int status;

	if (*Sws_1->NColloids>0)
		status=fread(Sws_1->TwoBody[0],sizeof(double),(*Sws_1->NColloids)*10,pfile);

	if (Model==2)
		status=fread(Sws_2->TwoBody[0],sizeof(double),(*Sws_2->NColloids)*10,pfile);
}

void saveSW(FILE *pfile)
{
	int status;

	if (*Sws_1->NColloids>0)
		status=fwrite(Sws_1->TwoBody[0],sizeof(double),(*Sws_1->NColloids)*10,pfile);

	if (Model==2)
		status=fwrite(Sws_2->TwoBody[0],sizeof(double),(*Sws_2->NColloids)*10,pfile);
}


void bonds_energy(sw *sws,interactionmap *im,double *energy)
{
	// ATTENZIONE: RICHIEDE CHE sws->TwoBody SIA CALCOLATO

	double twobody_buffer[10];


	// SCORRIAMO SU TUTTE LE COPPIE DI VICINI
	int i;
	for (i=0;i<im->num;i++)
	{
		int particle1=im->who[i];

		int j;
		for (j=0;j<im->howmany[i];j++)
		{
			int particle2=im->with[i][j];

			// parte a due corpi
			double dist=sqrt(im->rij2[i][j]);
			double inorm=1./dist;
			double inorm2=SQR(inorm);
			vector vdist=im->rij[i][j];
			double gij=exp(sws->Gamma/(dist-sws->a));
			double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));

			energy[particle1]+=v2;
			energy[particle2]+=v2;

			double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[particle2]));

			energy[particle1]+=0.5*Uc;
			energy[particle2]+=0.5*Uc;

			// sommiamo la parte a tre corpi di i e j

			energy[particle1]+=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9])));
			energy[particle2]+=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9])));
			energy[particle1]+=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9])));
			energy[particle2]+=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9])));

			// sottraiamo al tensore di i il contributo di j
			memcpy(twobody_buffer,sws->TwoBody[particle1],10*sizeof(double));

			// vdist is the vector going from j to i
			twobody_buffer[0]-=gij;
			twobody_buffer[1]-=gij*vdist.x*inorm;
			twobody_buffer[2]-=gij*vdist.y*inorm;
			twobody_buffer[3]-=gij*vdist.z*inorm;
			twobody_buffer[4]-=gij*SQR(vdist.x)*inorm2;
			twobody_buffer[5]-=gij*SQR(vdist.y)*inorm2;
			twobody_buffer[6]-=gij*SQR(vdist.z)*inorm2;
			twobody_buffer[7]-=gij*vdist.x*vdist.y*inorm2;
			twobody_buffer[8]-=gij*vdist.x*vdist.z*inorm2;
			twobody_buffer[9]-=gij*vdist.y*vdist.z*inorm2;

			energy[particle1]-=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(twobody_buffer[0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(twobody_buffer[1])+SQR(twobody_buffer[2])+SQR(twobody_buffer[3]))+0.5*sws->lambda*(SQR(twobody_buffer[4])+SQR(twobody_buffer[5])+SQR(twobody_buffer[6])+2.*SQR(twobody_buffer[7])+2.*SQR(twobody_buffer[8])+2.*SQR(twobody_buffer[9])));
			energy[particle2]-=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(twobody_buffer[0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(twobody_buffer[1])+SQR(twobody_buffer[2])+SQR(twobody_buffer[3]))+0.5*sws->lambda*(SQR(twobody_buffer[4])+SQR(twobody_buffer[5])+SQR(twobody_buffer[6])+2.*SQR(twobody_buffer[7])+2.*SQR(twobody_buffer[8])+2.*SQR(twobody_buffer[9])));

			// sottraiamo al tensore di j il contributo di i
			memcpy(twobody_buffer,sws->TwoBody[particle2],10*sizeof(double));

			twobody_buffer[0]-=gij;
			twobody_buffer[1]+=gij*vdist.x*inorm;
			twobody_buffer[2]+=gij*vdist.y*inorm;
			twobody_buffer[3]+=gij*vdist.z*inorm;
			twobody_buffer[4]-=gij*SQR(vdist.x)*inorm2;
			twobody_buffer[5]-=gij*SQR(vdist.y)*inorm2;
			twobody_buffer[6]-=gij*SQR(vdist.z)*inorm2;
			twobody_buffer[7]-=gij*vdist.x*vdist.y*inorm2;
			twobody_buffer[8]-=gij*vdist.x*vdist.z*inorm2;
			twobody_buffer[9]-=gij*vdist.y*vdist.z*inorm2;

			energy[particle1]-=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(twobody_buffer[0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(twobody_buffer[1])+SQR(twobody_buffer[2])+SQR(twobody_buffer[3]))+0.5*sws->lambda*(SQR(twobody_buffer[4])+SQR(twobody_buffer[5])+SQR(twobody_buffer[6])+2.*SQR(twobody_buffer[7])+2.*SQR(twobody_buffer[8])+2.*SQR(twobody_buffer[9])));
			energy[particle2]-=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(twobody_buffer[0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(twobody_buffer[1])+SQR(twobody_buffer[2])+SQR(twobody_buffer[3]))+0.5*sws->lambda*(SQR(twobody_buffer[4])+SQR(twobody_buffer[5])+SQR(twobody_buffer[6])+2.*SQR(twobody_buffer[7])+2.*SQR(twobody_buffer[8])+2.*SQR(twobody_buffer[9])));

		}

	}
}
