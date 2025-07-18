#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "vector.h"
#include "global_definitions.h"
#include "io.h"
#include "interaction_map.h"
#include "events.h"
#include "cell_list.h"
#include "dirutils.h"
#include "smart_allocator.h"
#include "ordinator.h"
#include "bilista.h"
#include "random.h"
#include "sw_disorder.h"


#define SQR(x) ((x)*(x))
#define SCALAR(v,u) (((v)->x)*((u)->x)+((v)->y)*((u)->y)+((v)->z)*((u)->z))


typedef struct _bonds {
	int num;		// number particles interacting
	int *who;		// particles interacting
	int **with;		// interacting partners
	int *howmany;		// number of interacting partners
	double **energy;
} bonds;

static inline double intpow (const double x, const int i) {
	if (i < 0) abort();
	if (i == 0) return 1.;
	if (i == 1) return x;
	else return x * intpow(x, i - 1);
}

void calculate_potential(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *total_potential,double *twobody,double *threebody);
void check(sw *const sws,vector *pos,interactionmap *ime,double Box[],double cutoff,double *potential,double *virial,double *twobody,double *threebody);
double check_bonds(sw *const sws,interactionmap *ime,bonds *b);
double bonds_energy_local(sw *const sws,interactionmap *im,bonds *b,double cutoff);
bonds* createBondsMap(int max_elements,int max_neighbours);
void freeBonds(bonds *i);
void ImToBonds(bonds *b,interactionmap *im);
void calculate_potential_complete(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *total_potential,double *twobody,double *threebody,double *tottwobody,double *totthreebody);
void BimeFromBim(bonds *bim,bonds *bime);
void BimeFromBim_EnergyOrdered(bonds *bim,bonds *bime);
void insertionSort_energy(int *number,int num_el,double *energy,double energy_el,int *length);

int main(int argc,char *argv[])
{

	if (argc!=3)
	{
		printf("%s [input file] [lambda=23.15 for water]\n",argv[0]);
		exit(1);
	}

	char input[100];
	strcpy(input,argv[1]);
	double lambda=atof(argv[2]);
	double range=2.5;

	// leggiamo la configurazione
	int ncolloids=getNumberParticles(input);

	int *lista_atoms=calloc(ncolloids,sizeof(int));
	int num_atoms_lista=0;

		int i;
		for (i=0;i<ncolloids;i++)
			lista_atoms[num_atoms_lista++]=i;


	vector *pos=calloc(ncolloids,sizeof(vector));
	double *tetrahedral=calloc(ncolloids,sizeof(double));
	steps time;
	double box[6],ibox[6];
	readPositions(input,pos,&time,&ncolloids,box,ibox);

	// celle
	listcell *cells=getList(box,range,ncolloids);
	fullUpdateList(cells,pos,ncolloids,box,range);

	// struttura per il calcolo del potenziale a tre corpi
	sw *Sws=createSW(ncolloids,&ncolloids,lambda);

	// lista delle interazioni
	interactionmap *im=createInteractionMap(ncolloids,ncolloids);
	calculateSystemInteractionMap(cells,im,pos,box,range);

	double potential;

	// two and three body terms for each particle
	double *twobody=calloc(ncolloids,sizeof(double));
	double *threebody=calloc(ncolloids,sizeof(double));

	double totthreebody,tottwobody;

	calculate_potential_complete(Sws,pos,im,box,1.8,&potential,twobody,threebody,&tottwobody,&totthreebody);

	bonds *bim=createBondsMap(ncolloids,ncolloids);
	ImToBonds(bim,im);

	double bonds_potential=bonds_energy_local(Sws,im,bim,1.8);

	bonds *bime=createBondsMap(ncolloids,ncolloids);
	BimeFromBim_EnergyOrdered(bim,bime);

	assert(fabs(bonds_potential-potential)<0.000001);

	printf("%d\n",ncolloids);
	int ilist;
	for (ilist=0;ilist<num_atoms_lista;ilist++)
	{
		int i=lista_atoms[ilist];

		int j,k;

		double ienergy=0.;
		int numneighbours=0;

		//assert(bime->howmany[i]>4);
		int particle=i;
		bime->howmany[particle]=(bime->howmany[particle]<4 ? bime->howmany[particle] : 4 );

		int howmany=bime->howmany[particle];

		for (j=0;j<howmany;j++)
		{
			int particle_j=bime->with[particle][j];

			for (k=j+1;k<howmany;k++)
			{
				int particle_k=bime->with[particle][k];

				vector dist_j,dist_k;
				vector olddist_j,olddist_k;

				olddist_j.x=pos[particle_j].x-pos[particle].x;
				olddist_j.y=pos[particle_j].y-pos[particle].y;
				olddist_j.z=pos[particle_j].z-pos[particle].z;

				olddist_j.x-=rint(olddist_j.x);
				olddist_j.y-=rint(olddist_j.y);
				olddist_j.z-=rint(olddist_j.z);

				dist_j.x=box[0]*olddist_j.x+box[1]*olddist_j.y+box[2]*olddist_j.z;
				dist_j.y=box[3]*olddist_j.y+box[4]*olddist_j.z;
				dist_j.z=box[5]*olddist_j.z;


				double distnorm_j=sqrt(SQR(dist_j.x)+SQR(dist_j.y)+SQR(dist_j.z));

				olddist_k.x=pos[particle_k].x-pos[particle].x;
				olddist_k.y=pos[particle_k].y-pos[particle].y;
				olddist_k.z=pos[particle_k].z-pos[particle].z;

				olddist_k.x-=rint(olddist_k.x);
				olddist_k.y-=rint(olddist_k.y);
				olddist_k.z-=rint(olddist_k.z);

				dist_k.x=box[0]*olddist_k.x+box[1]*olddist_k.y+box[2]*olddist_k.z;
				dist_k.y=box[3]*olddist_k.y+box[4]*olddist_k.z;
				dist_k.z=box[5]*olddist_k.z;

				double distnorm_k=sqrt(SQR(dist_k.x)+SQR(dist_k.y)+SQR(dist_k.z));

				double costeta=(dist_j.x*dist_k.x+dist_j.y*dist_k.y+dist_j.z*dist_k.z)/(distnorm_j*distnorm_k);

				tetrahedral[i]+=SQR(costeta+1./3.);
			}
		}

		tetrahedral[i]=1-9.*tetrahedral[i]/(2.*howmany*(howmany-1.));

		/*
		for (j=0;j<4;j++)
		{

			ienergy+=bime->energy[i][j];
			numneighbours++;
			printf("%d %d\n",i,bime->with[i][j]);

			//printf("%lf\n",bime->energy[i][j]);
		}
		*/

	}

	// attenzione: non aggiorniamo bime->energy perche' non ci serve piu'
	for (i=0;i<ncolloids;i++)
	{
		int j=0;
		while (j<bime->howmany[i])
		//for (j=0;j<bime->howmany[i];j++)
		{
			// looking for non-reciprocated bonds
			int neighbour=bime->with[i][j];




			int k=0;
			while ((bime->with[neighbour][k]!=i) && (++k<bime->howmany[neighbour]));

			if (k==bime->howmany[neighbour])
			{
				// non reciprocated bond

				if (tetrahedral[i]>tetrahedral[neighbour])
				{
					// add bond between neighbour and i
					bime->with[neighbour][bime->howmany[neighbour]++]=i;
				}
				else
				{
					// remove bond between i and neighbour
					int l=j;
					bime->howmany[i]--;
					while (l<bime->howmany[i])
					{
						bime->with[i][l]=bime->with[i][l+1];
						l++;
					}
				}

			}


			j++;
		}
	}

	for (i=0;i<ncolloids;i++)
	{
		int j;
		for (j=0;j<bime->howmany[i];j++)
			printf("%d %d\n",i,bime->with[i][j]);
	}

	freeInteractionMap(im);
	freeList(cells);
	free(pos);
	free(twobody);
	free(threebody);
	freeSW(Sws);
	freeBonds(bim);
	freeBonds(bime);

	return 0;
}


void calculate_potential(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *total_potential,double *twobody,double *threebody)
{
	int i,j;


	double Potential_2body=0.;
	double Potential_3body=0.;


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
			twobody[particle1]+=v2;
			twobody[particle2]+=v2;

			Potential_3body+=Uc;
			threebody[particle1]+=Uc;
			threebody[particle2]+=Uc;

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

		double buffer=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));
		Potential_3body+=buffer;
		threebody[particle1]+=2.*buffer;
	}

	//printf("Total two-body %lf; three-body %lf total %lf\n",Potential_2body,Potential_3body,Potential_2body+Potential_3body);

	*total_potential=Potential_2body+Potential_3body;

	//printf("total_potential %lf\n",*total_potential);
}




void calculate_potential_complete(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *total_potential,double *twobody,double *threebody,double *tottwobody,double *totthreebody)
{
	int i,j;


	double Potential_2body=0.;
	double Potential_3body=0.;


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

			if (dist<cutoff)
			{
				// vdist is the vector going from j to i
				vector vdist=im->rij[i][j];

				double gij=exp(sws->Gamma/(dist-sws->a));

				double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));
				double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[particle2]));

				//Potential_2body+=v2+Uc;
				Potential_2body+=v2;
				twobody[particle1]+=v2;
				twobody[particle2]+=v2;

				Potential_3body+=Uc;
				threebody[particle1]+=Uc;
				threebody[particle2]+=Uc;

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



	}

	for (i=0;i<im->num;i++)
	{
		int particle1=im->who[i];

		double buffer=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));
		Potential_3body+=buffer;
		threebody[particle1]+=2.*buffer;
	}


	*tottwobody=Potential_2body;
	*totthreebody=Potential_3body;
	*total_potential=Potential_2body+Potential_3body;

}


double bonds_energy_local(sw *const sws,interactionmap *im,bonds *b,double cutoff)
{
	// ATTENZIONE: RICHIEDE CHE sws->TwoBody SIA CALCOLATO

	double twobody_buffer[10];

	double total_energy=0.;

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

			if (dist<cutoff)
			{
				double inorm=1./dist;
				double inorm2=SQR(inorm);
				vector vdist=im->rij[i][j];
				double gij=exp(sws->Gamma/(dist-sws->a));
				double v2=sws->A*(sws->B/intpow(dist,sws->p)-1.)*exp(1./(dist-sws->a));

				b->energy[i][j]=v2;

				double Uc=-0.5*sws->lambda*SQR(gij)*(SQR(1.-sws->costeta0[particle1])+SQR(1.-sws->costeta0[particle2]));

				b->energy[i][j]+=0.5*Uc;

				// sommiamo la parte a tre corpi di i e j

				b->energy[i][j]+=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9])));
				b->energy[i][j]+=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(sws->TwoBody[particle2][0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(sws->TwoBody[particle2][1])+SQR(sws->TwoBody[particle2][2])+SQR(sws->TwoBody[particle2][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle2][4])+SQR(sws->TwoBody[particle2][5])+SQR(sws->TwoBody[particle2][6])+2.*SQR(sws->TwoBody[particle2][7])+2.*SQR(sws->TwoBody[particle2][8])+2.*SQR(sws->TwoBody[particle2][9])));

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

				b->energy[i][j]-=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(twobody_buffer[0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(twobody_buffer[1])+SQR(twobody_buffer[2])+SQR(twobody_buffer[3]))+0.5*sws->lambda*(SQR(twobody_buffer[4])+SQR(twobody_buffer[5])+SQR(twobody_buffer[6])+2.*SQR(twobody_buffer[7])+2.*SQR(twobody_buffer[8])+2.*SQR(twobody_buffer[9])));

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

				b->energy[i][j]-=0.5*(0.5*sws->lambda*SQR(sws->costeta0[particle2])*SQR(twobody_buffer[0])+sws->lambda*(-sws->costeta0[particle2])*(SQR(twobody_buffer[1])+SQR(twobody_buffer[2])+SQR(twobody_buffer[3]))+0.5*sws->lambda*(SQR(twobody_buffer[4])+SQR(twobody_buffer[5])+SQR(twobody_buffer[6])+2.*SQR(twobody_buffer[7])+2.*SQR(twobody_buffer[8])+2.*SQR(twobody_buffer[9])));

				total_energy+=b->energy[i][j];
			}
		}

	}

// 	for (i=0;i<im->num;i++)
// 	{
// 		int particle1=im->who[i];
//
// 		double buffer=0.5*sws->lambda*SQR(sws->costeta0[particle1])*SQR(sws->TwoBody[particle1][0])+sws->lambda*(-sws->costeta0[particle1])*(SQR(sws->TwoBody[particle1][1])+SQR(sws->TwoBody[particle1][2])+SQR(sws->TwoBody[particle1][3]))+0.5*sws->lambda*(SQR(sws->TwoBody[particle1][4])+SQR(sws->TwoBody[particle1][5])+SQR(sws->TwoBody[particle1][6])+2.*SQR(sws->TwoBody[particle1][7])+2.*SQR(sws->TwoBody[particle1][8])+2.*SQR(sws->TwoBody[particle1][9]));
//
// 	}
//
	return total_energy;
}


double check_bonds(sw *const sws,interactionmap *ime,bonds *b)
{
	int i,j,k;

	double Potential_2body=0.;
	double Potential_3body=0.;

	for (i=0;i<ime->num;i++)
	{
		int particlei=ime->who[i];

		for (j=0;j<ime->howmany[i];j++)
		{

			int particlej=ime->with[i][j];

			// parte a due corpi
			double distij=sqrt(ime->rij2[i][j]);
			vector rij=ime->rij[i][j];

			double v2=sws->A*(sws->B/pow(distij,sws->p)-1.)*exp(1./(distij-sws->a));
			Potential_2body+=v2;

			b->energy[particlei][j]+=v2;


			// dobbiamo cercare dove particlej interagisce con particlei
			int pos=0;
			while (b->with[particlej][pos]!=particlei)
				pos++;




			// parte a tre corpi
			for (k=0;k<ime->howmany[i];k++)
			{
				if (k!=j)
				{
					//int particlek=ime->with[i][k];

					double distik=sqrt(ime->rij2[i][k]);
					vector rik=ime->rij[i][k];

					double costeta=(rij.x*rik.x+rij.y*rik.y+rij.z*rik.z)/(distij*distik);

					double buffer=sws->lambda*exp(sws->Gamma/(distij-sws->a)+sws->Gamma/(distik-sws->a))*SQR(costeta-sws->costeta0[i]);
					Potential_3body+=buffer;
					//threebody[i]+=buffer;
					b->energy[particlei][j]+=0.5*buffer;
					b->energy[particlej][pos]+=0.5*buffer;
				}
			}


		}
	}

	//return Potential_2body+Potential_3body;
	return 0.5*(Potential_2body+Potential_3body);

	//printf("check completo; total %lf 2-body %lf 3-body %lf\n",*potential/(double)ime->num,Potential_2body,Potential_3body);
}


void check(sw *const sws,vector *pos,interactionmap *ime,double Box[],double cutoff,double *potential,double *virial,double *twobody,double *threebody)
{
	int i,j,k;

	double Potential_2body=0.;
	double Potential_3body=0.;

	//*twobody=0.;
	//*threebody=0.;

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
			Potential_2body+=0.5*v2;
			twobody[i]+=v2;

			// parte a tre corpi
			for (k=j+1;k<ime->howmany[i];k++)
			{
				//int particlek=ime->with[i][k];

				double distik=sqrt(ime->rij2[i][k]);
				vector rik=ime->rij[i][k];

				double costeta=(rij.x*rik.x+rij.y*rik.y+rij.z*rik.z)/(distij*distik);

				double buffer=sws->lambda*exp(sws->Gamma/(distij-sws->a)+sws->Gamma/(distik-sws->a))*SQR(costeta-sws->costeta0[i]);
				Potential_3body+=buffer;
				threebody[i]+=buffer;
			}
		}
	}

	*potential=Potential_2body+Potential_3body;

	//printf("check completo; total %lf 2-body %lf 3-body %lf\n",*potential/(double)ime->num,Potential_2body,Potential_3body);
}


bonds* createBondsMap(int max_elements,int max_neighbours)
{
	bonds*i=malloc(sizeof(bonds));;

	i->num=0;
	i->who=calloc(max_elements,sizeof(int));
	i->howmany=calloc(max_elements,sizeof(int));

	// la creazione degli array bidimensionali con un massimo numero di vicini
	// e' delicata. Dobbiamo lasciare dello spazio in coda a questi array nel
	// caso ci sia uno sforamento nel numero di vicini dell'ultima particella.
	// lasciamo in fondo max_elements spazi vuoti
	Matrix2DSafe(i->with,max_elements,max_neighbours,max_elements,int);
	Matrix2DSafe(i->energy,max_elements,max_neighbours,max_elements,double);

	return i;
}


void freeBonds(bonds *i)
{
	free(i->who);
	Free2D(i->with);
	Free2D(i->energy);
	free(i->howmany);
	free(i);
}


void ImToBonds(bonds *b,interactionmap *im)
{
	int i;

	b->num=im->num;

	for (i=0;i<b->num;i++)
	{
		b->who[i]=im->who[i];
		b->howmany[i]=im->howmany[i];

		int j;
		for (j=0;j<b->howmany[i];j++)
		{
			b->with[i][j]=im->with[i][j];
			b->energy[i][j]=0.;
		}
	}
}

void BimeFromBim(bonds *bim,bonds *bime)
{
	bime->num=bim->num;

	int i,index,j;
	for (i=0;i<bim->num;i++)
	{
		bime->who[i]=bim->who[i];

		for (index=0;index<bim->howmany[i];index++)
		{
			j=bim->with[i][index];
			double energy=bim->energy[i][index];

			bime->energy[i][bime->howmany[i]]=energy;
			bime->energy[j][bime->howmany[j]]=energy;

			bime->with[i][bime->howmany[i]++]=j;
			bime->with[j][bime->howmany[j]++]=i;
		}
	}
}

void BimeFromBim_EnergyOrdered(bonds *bim,bonds *bime)
{
	bime->num=bim->num;

	int i,index,j;
	for (i=0;i<bim->num;i++)
	{
		bime->who[i]=bim->who[i];

		for (index=0;index<bim->howmany[i];index++)
		{
			j=bim->with[i][index];
			double energy=bim->energy[i][index];

			insertionSort_energy(bime->with[i],j,bime->energy[i],energy,bime->howmany+i);
			insertionSort_energy(bime->with[j],i,bime->energy[j],energy,bime->howmany+j);
		}
	}
}

void insertionSort_energy(int *number,int num_el,double *energy,double energy_el,int *length)
{
	int l=*length;
	number[l]=num_el;
	energy[l]=energy_el;
	while ((l>0) && (energy[l]<energy[l-1]))
	{
		double buffer_energy;
		buffer_energy=energy[l];
		energy[l]=energy[l-1];
		energy[l-1]=buffer_energy;

		int buffer_num;
		buffer_num=number[l];
		number[l]=number[l-1];
		number[l-1]=buffer_num;

		l--;
	}
	(*length)++;
}
