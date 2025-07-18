#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "vector.h"
#include "interaction_map.h"
#include "secure_search.h"
#include "log.h"
#include "smart_allocator.h"
#include "events.h"
#include "cell_list.h"


#define SQR(x) ((x)*(x))
#define MAX_NEIGHBOURS 250

// variabili private del modulo
static listcell *ListMc;
static int Max_neighbours;

static double Cellsize;

// variabili del modulo mc
static args *ParametersMc;
static int NumParametersMc;

// lista degli eventi
void celllistGetSystemInteractionMap(args *argv,int argc);
void celllistGetParticleInteractionMap(args *argv,int argc);
void celllistGetPairInteractionMap(args *argv,int argc);
void celllistRemoveParticle(args *argv,int argc);
void celllistAddParticle(args *argv,int argc);

static int module(int n,int mo)
{
	n=n%mo;
	
	while (n<0)
	{
		n+=mo;
	}
	
	return n;
}


listcell* getList(double Box[],double cutoff,int num_particles)
{
	listcell *l=malloc(sizeof(listcell));
	
	double volume=Box[0]*Box[3]*Box[5];
	l->NumberCells_x=(int)(volume/(sqrt(Box[3]*Box[3]*Box[5]*Box[5]+Box[5]*Box[5]*Box[1]*Box[1]+Box[3]*Box[3]*Box[2]*Box[2]+Box[4]*Box[4]*Box[1]*Box[1])*cutoff));
	l->NumberCells_y=(int)(volume/(Box[0]*sqrt(Box[5]*Box[5]+Box[4]*Box[4])*cutoff));
	l->NumberCells_z=(int)(volume/(Box[0]*Box[3]*cutoff));
	
	
	int *maxncells;
	while ((l->NumberCells_x*l->NumberCells_y*l->NumberCells_z>27) && (l->NumberCells_x*l->NumberCells_y*l->NumberCells_z>num_particles))
	{
		maxncells=(l->NumberCells_x>l->NumberCells_y ? &(l->NumberCells_x) : &(l->NumberCells_y) );
		maxncells=( *maxncells > l->NumberCells_z ? maxncells : &(l->NumberCells_z));
		
		(*maxncells)--;
	}
	
	int *mincells;
	mincells=(l->NumberCells_x<l->NumberCells_y ? &(l->NumberCells_x) : &(l->NumberCells_y) );
	mincells=( *mincells < l->NumberCells_z ? mincells : &(l->NumberCells_z));
	
	if (*mincells<3)
	{
		logPrint("Error: not enough cells\n");
		exit(1);
	}
	
	// le coordinate sono periodiche fra [0,1], quindi i lati della cella nello
	// spazio s sono sempre 1
	l->CellSize_x=1./(double)l->NumberCells_x;
	l->CellSize_y=1./(double)l->NumberCells_y;
	l->CellSize_z=1./(double)l->NumberCells_z;
	
	
	l->HoC=(int*)calloc(l->NumberCells_x*l->NumberCells_y*l->NumberCells_z,sizeof(int));
	l->LinkedList=(int*)calloc(num_particles,sizeof(int));
	
	int i;
	int nnn=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
	for (i=0;i<nnn;i++)
		(l->HoC)[i]=-1;
	
	l->MyCell=calloc(num_particles,sizeof(int));
	
	return l;
}


void copyList(listcell *dst,listcell *src,int ncolloids)
{
	dst->NumberCells_x=src->NumberCells_x;
	dst->NumberCells_y=src->NumberCells_y;
	dst->NumberCells_z=src->NumberCells_z;
	dst->CellSize_x=src->CellSize_x;
	dst->CellSize_y=src->CellSize_y;
	dst->CellSize_z=src->CellSize_z;
	
	int size=src->NumberCells_x*src->NumberCells_y*src->NumberCells_z;
	
	memcpy(dst->HoC,src->HoC,size*sizeof(int));
	memcpy(dst->LinkedList,src->LinkedList,ncolloids*sizeof(int));
	memcpy(dst->MyCell,src->MyCell,ncolloids*sizeof(int));
	
}


void celllistConstructor(FILE *config_file,vector *pos,int max_number_colloids,int *ncolloids,double box[],double *cutoff,interactionmap *interactionList,int *movedparticle)
{
	
	SearchTable *s=searchNew();
	searchDouble("Cell_size",&Cellsize,s);
	int *found_maxneighbours=searchTryInt("Max_neighbours",&Max_neighbours,s);
	searchFile(config_file,s);
	
	
	if (!*found_maxneighbours)
	{
		Max_neighbours=MAX_NEIGHBOURS;
	}
	
	if (Cellsize<*cutoff)
	{
		logPrint("Warning: Cell size too small - was %lf; expanding to %lf\n",Cellsize,*cutoff);
		Cellsize=*cutoff+0.0001;
	}
	
	ListMc=getList(box,Cellsize,max_number_colloids);
	
	NumParametersMc=8;
	ParametersMc=calloc(NumParametersMc,sizeof(args));
	
	ParametersMc[0].argv=(void*)ListMc;
	ParametersMc[1].argv=(void*)pos;
	ParametersMc[2].argv=(void*)ncolloids;
	ParametersMc[3].argv=(void*)box;
	ParametersMc[4].argv=(void*)interactionList;
	ParametersMc[5].argv=(void*)&Cellsize;
	ParametersMc[6].argv=(void*)cutoff;
	ParametersMc[7].argv=(void*)movedparticle;
	
	// eventi
	int status=0;
	status+=eventAddEvent(INTERACTIONMAP_SYSTEM);
	status+=eventAddEvent(INTERACTIONMAP_PARTICLE);
	status+=eventAddEvent(INTERACTIONMAP_PAIR);
	status+=eventAddEvent(TRAJECTORYMOVE_REJECTED);
	status+=eventAddEvent(SYSTEMMOVE_REJECTED);
	status+=eventAddEvent(REMOVE_ACCEPT);
	status+=eventAddEvent(ADD_ACCEPT);
	
	
	
	if (status!=6)
	{
		logPrint("Warning: events %d and %d already defined\n",INTERACTIONMAP_SYSTEM,INTERACTIONMAP_PARTICLE);
	}
	
	eventAddAction(INTERACTIONMAP_SYSTEM,&celllistGetSystemInteractionMap,ParametersMc,NumParametersMc);
	eventAddAction(INTERACTIONMAP_PARTICLE,&celllistGetParticleInteractionMap,ParametersMc,NumParametersMc);
	eventAddAction(INTERACTIONMAP_PAIR,&celllistGetPairInteractionMap,ParametersMc,NumParametersMc);
	eventAddAction(TRAJECTORYMOVE_REJECTED,&celllistGetSystemInteractionMap,ParametersMc,NumParametersMc);
	eventAddAction(SYSTEMMOVE_REJECTED,&celllistGetSystemInteractionMap,ParametersMc,NumParametersMc);
	eventAddAction(REMOVE_ACCEPT,&celllistRemoveParticle,ParametersMc,NumParametersMc);
	eventAddAction(ADD_ACCEPT,&celllistAddParticle,ParametersMc,NumParametersMc);
	
	eventAddIO(INTERACTIONMAP_PAIR,INTERACTIONMAP_PAIR__p1,sizeof(int));
	eventAddIO(INTERACTIONMAP_PAIR,INTERACTIONMAP_PAIR__p2,sizeof(int));
	
	logPrint("Cell size: %lf\n",Cellsize);
	
	searchFree(s);
}

void celllistFree()
{
	freeList(ListMc);
	free(ParametersMc);
}


void celllistGetSystemInteractionMap(args *argv,int argc)
{
	// lettura degli argomenti
	listcell *list=(listcell*)argv[0].argv;
	vector *pos=(vector*)argv[1].argv;
	int ncolloids=*((int*)argv[2].argv);
	double *box=(double*)argv[3].argv;
	interactionmap *im=(interactionmap *)argv[4].argv;
	double cellsize=*((double*)argv[5].argv);
	double cutoff=*((double*)argv[6].argv);
	
	fullUpdateList(list,pos,ncolloids,box,cellsize);
	int neighbours=calculateSystemInteractionMap(list,im,pos,box,cutoff);
	
	if (neighbours>Max_neighbours)
	{
		logPrint("Error: bad estimate of MAX_NEIGHBOURS (%d), since %d neighbours were found\n",MAX_NEIGHBOURS,neighbours);
		exit(1);
	}
}

void celllistUpdateList()
{
	args *argv=ParametersMc;
	//int argc=NumParametersMc;
	// lettura degli argomenti
	listcell *list=(listcell*)argv[0].argv;
	vector *pos=(vector*)argv[1].argv;
	int ncolloids=*((int*)argv[2].argv);
	double *box=(double*)argv[3].argv;
	//interactionmap *im=(interactionmap *)argv[4].argv;
	double cellsize=*((double*)argv[5].argv);
	//double cutoff=*((double*)argv[6].argv);
	
	fullUpdateList(list,pos,ncolloids,box,cellsize);
}

void celllistGetParticleInteractionMap(args *argv,int argc)
{
	// lettura degli argomenti
	listcell *list=(listcell*)argv[0].argv;
	vector *pos=(vector*)argv[1].argv;
	//int ncolloids=*((int*)argv[2].argv);
	double *box=(double*)argv[3].argv;
	interactionmap *im=(interactionmap *)argv[4].argv;
	//double cellsize=*((double*)argv[5].argv);
	double cutoff=*((double*)argv[6].argv);
	int movedparticle=*((int*)argv[7].argv);
	
	int neighbours=getParticleInteractionMap(list,pos,movedparticle,im,box,cutoff);
	
	if (neighbours>Max_neighbours)
	{
		logPrint("Error: bad estimate of MAX_NEIGHBOURS (%d), since %d neighbours were found\n",MAX_NEIGHBOURS,neighbours);
		exit(1);
	}
	
}

void celllistRemoveParticle(args *argv,int argc)
{
	// lettura degli argomenti
	listcell *list=(listcell*)argv[0].argv;
	//vector *pos=(vector*)argv[1].argv;
	int ncolloids=*((int*)argv[2].argv);
	//double *box=(double*)argv[3].argv;
	//interactionmap *im=(interactionmap *)argv[4].argv;
	//double cellsize=*((double*)argv[5].argv);
	//double cutoff=*((double*)argv[6].argv);
	//int movedparticle=*((int*)argv[7].argv);
	
	int movedparticle=*((int*)eventGetIO(REMOVE_PARTICLE,REMOVE_PARTICLE__particle));
	
	if (movedparticle==ncolloids-1)
	{
		removeFromList(list,movedparticle);
	}
	else
	{
		changeIdentityInList(list,ncolloids-1,movedparticle);
	}
	
}

void celllistAddParticle(args *argv,int argc)
{
	// lettura degli argomenti
	listcell *list=(listcell*)argv[0].argv;
	vector *pos=(vector*)argv[1].argv;
	int ncolloids=*((int*)argv[2].argv);
	//double *box=(double*)argv[3].argv;
	//interactionmap *im=(interactionmap *)argv[4].argv;
	//double cellsize=*((double*)argv[5].argv);
	//double cutoff=*((double*)argv[6].argv);
	//int movedparticle=*((int*)argv[7].argv);
	
	addToList(list,pos+ncolloids,ncolloids);
}

void celllistGetPairInteractionMap(args *argv,int argc)
{
	// lettura degli argomenti
	listcell *list=(listcell*)argv[0].argv;
	vector *pos=(vector*)argv[1].argv;
	//int ncolloids=*((int*)argv[2].argv);
	double *box=(double*)argv[3].argv;
	interactionmap *im=(interactionmap *)argv[4].argv;
	//double cellsize=*((double*)argv[5].argv);
	double cutoff=*((double*)argv[6].argv);
	//int movedparticle=*((int*)argv[7].argv);
	
	int particle1=*((int*)eventGetIO(INTERACTIONMAP_PAIR,INTERACTIONMAP_PAIR__p1));
	int particle2=*((int*)eventGetIO(INTERACTIONMAP_PAIR,INTERACTIONMAP_PAIR__p2));
	
	im->num=0;
	int neighbours=getParticleInteractionMap_ImPos(list,pos,particle1,im,0,box,cutoff);
	
	if (neighbours>Max_neighbours)
	{
		logPrint("Error: bad estimate of MAX_NEIGHBOURS (%d), since %d neighbours were found\n",MAX_NEIGHBOURS,neighbours);
		exit(1);
	}
	
	neighbours=getParticleInteractionMap_ImPos(list,pos,particle2,im,1,box,cutoff);
	
	if (neighbours>Max_neighbours)
	{
		logPrint("Error: bad estimate of MAX_NEIGHBOURS (%d), since %d neighbours were found\n",MAX_NEIGHBOURS,neighbours);
		exit(1);
	}
}

void freeList(listcell *l)
{
	free(l->HoC);
	free(l->LinkedList);
	free(l->MyCell);
	free(l);
}

void resetList(listcell *l)
{
	int i;
	int nnn=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
	
	for (i=0;i<nnn;i++)
		(l->HoC)[i]=-1;
}

void updateList(listcell *l,const vector *pos,int num)
{
	int i;
	int ncell;                       // cell number
	int posx,posy,posz;              // cell coordinates
	int nnn;                         // total number of cells
	
	
	nnn=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
	
	// HoC initialization
	for (i=0;i<nnn;i++)
		(l->HoC)[i]=-1;
	
	// colloids loop
	for (i=0;i<num;i++)
	{
		posx=(int)floor(pos[i].x/l->CellSize_x);
		posy=(int)floor(pos[i].y/l->CellSize_y);
		posz=(int)floor(pos[i].z/l->CellSize_z);
		
		posx=module(posx,l->NumberCells_x);
		posy=module(posy,l->NumberCells_y);
		posz=module(posz,l->NumberCells_z);
		
		ncell=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);
		
		(l->LinkedList)[i]=(l->HoC)[ncell];
		(l->HoC)[ncell]=i;
		(l->MyCell)[i]=ncell;
	}
	
}

void fullUpdateList(listcell *l,const vector *pos,int num,double box[],double cutoff)
{
	int i;
	int ncell;                       // cell number
	int posx,posy,posz;              // cell coordinates
	
	int old=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
	
	double volume=box[0]*box[3]*box[5];
	l->NumberCells_x=(int)(volume/(sqrt(box[3]*box[3]*box[5]*box[5]+box[5]*box[5]*box[1]*box[1]+box[3]*box[3]*box[2]*box[2]+box[4]*box[4]*box[1]*box[1])*cutoff));
	l->NumberCells_y=(int)(volume/(box[0]*sqrt(box[5]*box[5]+box[4]*box[4])*cutoff));
	l->NumberCells_z=(int)(volume/(box[0]*box[3]*cutoff));
	
	int *maxncells;
	while ((l->NumberCells_x*l->NumberCells_y*l->NumberCells_z>27) && (l->NumberCells_x*l->NumberCells_y*l->NumberCells_z>num))
	{
		maxncells=(l->NumberCells_x>l->NumberCells_y ? &(l->NumberCells_x) : &(l->NumberCells_y) );
		maxncells=( *maxncells > l->NumberCells_z ? maxncells : &(l->NumberCells_z));
		
		(*maxncells)--;
	}
	
	int *mincells;
	mincells=(l->NumberCells_x<l->NumberCells_y ? &(l->NumberCells_x) : &(l->NumberCells_y) );
	mincells=( *mincells < l->NumberCells_z ? mincells : &(l->NumberCells_z));
	
	if (*mincells<3)
	{
		logPrint("Error: not enough cells\n");
		exit(1);
	}
	
	
	int nnn=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
	
	if ( (nnn)>(old) )
		l->HoC=realloc(l->HoC,nnn*sizeof(int));
	
	
	// le coordinate sono periodiche fra [0,1], quindi i lati della cella nello
	// spazio s sono sempre 1
	l->CellSize_x=1./(double)l->NumberCells_x;
	l->CellSize_y=1./(double)l->NumberCells_y;
	l->CellSize_z=1./(double)l->NumberCells_z;
		
	// HoC initialization
	for (i=0;i<nnn;i++)
		(l->HoC)[i]=-1;
	
	// colloids loop
	for (i=0;i<num;i++)
	{
		
		posx=(int)floor(pos[i].x/l->CellSize_x);
		posy=(int)floor(pos[i].y/l->CellSize_y);
		posz=(int)floor(pos[i].z/l->CellSize_z);
		
		posx=module(posx,l->NumberCells_x);
		posy=module(posy,l->NumberCells_y);
		posz=module(posz,l->NumberCells_z);
		
		ncell=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);
		
		(l->LinkedList)[i]=(l->HoC)[ncell];
		(l->HoC)[ncell]=i;
		(l->MyCell)[i]=ncell;
	}
}

int getParticleInteractionMap_ImPos(listcell *l,vector *pos,int label,interactionmap *im,int im_pos,double Box[],double cutoff)
{
	int j;
	int nx,ny,nz,nn;               // Number of cells on surface and in box
	int neighbour;                 // Neighbour cell
	int ix,iy,iz;                  // present cell coordinates
	int dx,dy,dz;
	int xv[3];                     // xlist values
	int yv[3];
	int zv[3];
	
	nx=l->NumberCells_x;
	ny=l->NumberCells_y;
	nz=l->NumberCells_z;
	nn=l->NumberCells_x*l->NumberCells_y;
	
	//pos+=label;
	
	// reconstruct cartesian coordinates
	ix=(int)floor(pos[label].x/l->CellSize_x);
	iy=(int)floor(pos[label].y/l->CellSize_y);
	iz=(int)floor(pos[label].z/l->CellSize_z);
	
	
	// the next step shoudn't be necessary
	xv[0]=module(ix,nx);
	yv[0]=module(iy,ny);
	zv[0]=module(iz,nz);
	
	
	// cell number
	//int ncell=xv[0]+yv[0]*l->NumberCells_x+zv[0]*(l->NumberCells_x*l->NumberCells_y);
	
	
	xv[1]=module(ix+1,nx);
	yv[1]=module(iy+1,ny);
	zv[1]=module(iz+1,nz);
	
	
	xv[2]=module(ix-1,nx);
	yv[2]=module(iy-1,ny);
	zv[2]=module(iz-1,nz);
	
	
	// interactions with particles in the same cell
	// and in neighbouring cells
	
	double cutoff2=cutoff*cutoff;
	
	int ncounter=0;
	
	vector dist;
	double dist2;
	
	for (dx=0;dx<3;dx++)
	{
		for (dy=0;dy<3;dy++)
		{
			for (dz=0;dz<3;dz++)
			{
				neighbour=xv[dx]+yv[dy]*nx+zv[dz]*nn;
				
				j=(l->HoC)[neighbour];
				
				while (j!=-1)
				{
					if (j!=label)
					{
						vector olddist;
						
						olddist.x=pos[label].x-pos[j].x;
						olddist.y=pos[label].y-pos[j].y;
						olddist.z=pos[label].z-pos[j].z;
						
						olddist.x-=rint(olddist.x);
						olddist.y-=rint(olddist.y);
						olddist.z-=rint(olddist.z);
						
						dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
						dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
						dist.z=Box[5]*olddist.z;
						
						dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
						
						if (dist2<cutoff2)
						{
							im->rij2[im_pos][ncounter]=dist2;
							im->rij[im_pos][ncounter]=dist;
							im->with[im_pos][ncounter++]=j;
						}
					}
					j=(l->LinkedList)[j];
				}
			}
		}
	}
	
	im->num++;
	im->who[im_pos]=label;
	im->howmany[im_pos]=ncounter;
	
	return ncounter;
}

int getParticleInteractionMap(listcell *l,vector *pos,int label,interactionmap *im,double Box[],double cutoff)
{
	int j;
	int nx,ny,nz,nn;               // Number of cells on surface and in box
	int neighbour;                 // Neighbour cell
	int ix,iy,iz;                  // present cell coordinates
	int dx,dy,dz;
	int xv[3];                     // xlist values
	int yv[3];
	int zv[3];
	
	nx=l->NumberCells_x;
	ny=l->NumberCells_y;
	nz=l->NumberCells_z;
	nn=l->NumberCells_x*l->NumberCells_y;
	
	//pos+=label;
	
	// reconstruct cartesian coordinates
	ix=(int)floor(pos[label].x/l->CellSize_x);
	iy=(int)floor(pos[label].y/l->CellSize_y);
	iz=(int)floor(pos[label].z/l->CellSize_z);
	
	
	// the next step shoudn't be necessary
	xv[0]=module(ix,nx);
	yv[0]=module(iy,ny);
	zv[0]=module(iz,nz);
	
	
	// cell number
	//int ncell=xv[0]+yv[0]*l->NumberCells_x+zv[0]*(l->NumberCells_x*l->NumberCells_y);
	
	
	xv[1]=module(ix+1,nx);
	yv[1]=module(iy+1,ny);
	zv[1]=module(iz+1,nz);
	
	
	xv[2]=module(ix-1,nx);
	yv[2]=module(iy-1,ny);
	zv[2]=module(iz-1,nz);
	
	
	// interactions with particles in the same cell
	// and in neighbouring cells
	
	double cutoff2=cutoff*cutoff;
	
	int ncounter=0;
	
	vector dist;
	double dist2;
	
	for (dx=0;dx<3;dx++)
	{
		for (dy=0;dy<3;dy++)
		{
			for (dz=0;dz<3;dz++)
			{
				neighbour=xv[dx]+yv[dy]*nx+zv[dz]*nn;
				
				j=(l->HoC)[neighbour];
				
				while (j!=-1)
				{
					if (j!=label)
					{
						vector olddist;
						
						olddist.x=pos[label].x-pos[j].x;
						olddist.y=pos[label].y-pos[j].y;
						olddist.z=pos[label].z-pos[j].z;
						
						olddist.x-=rint(olddist.x);
						olddist.y-=rint(olddist.y);
						olddist.z-=rint(olddist.z);
						
						dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
						dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
						dist.z=Box[5]*olddist.z;
						
						dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
						
						if (dist2<cutoff2)
						{
							im->rij2[0][ncounter]=dist2;
							im->rij[0][ncounter]=dist;
							im->with[0][ncounter++]=j;
						}
					}
					j=(l->LinkedList)[j];
				}
			}
		}
	}
	
	im->num=1;
	im->who[0]=label;
	im->howmany[0]=ncounter;
	
	return ncounter;
}


int celllistCountParticles(listcell *l,int cell_label)
{
	int nparticles=0;
	
	int j=(l->HoC)[cell_label];
	
	while (j!=-1)
	{
		nparticles++;
		
		j=(l->LinkedList)[j];
	}
	
	return nparticles;
}


int celllistGetCellNumber(listcell *l,vector *pos)
{
	int posx=(int)floor(pos->x/l->CellSize_x);
	int posy=(int)floor(pos->y/l->CellSize_y);
	int posz=(int)floor(pos->z/l->CellSize_z);
	
	posx=module(posx,l->NumberCells_x);
	posy=module(posy,l->NumberCells_y);
	posz=module(posz,l->NumberCells_z);
	
	return posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);
}

void changeCell(listcell *l,const vector *oldpos,const vector *newpos,int num)
{
	int ncell_old,ncell_new;         // cell number
	int posx,posy,posz;              // cell coordinates
	
	
	posx=(int)floor(oldpos->x/l->CellSize_x);
	posy=(int)floor(oldpos->y/l->CellSize_y);
	posz=(int)floor(oldpos->z/l->CellSize_z);
	
	posx=module(posx,l->NumberCells_x);
	posy=module(posy,l->NumberCells_y);
	posz=module(posz,l->NumberCells_z);
	
	ncell_old=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);
	
	posx=(int)floor(newpos->x/l->CellSize_x);
	posy=(int)floor(newpos->y/l->CellSize_y);
	posz=(int)floor(newpos->z/l->CellSize_z);
	
	posx=module(posx,l->NumberCells_x);
	posy=module(posy,l->NumberCells_y);
	posz=module(posz,l->NumberCells_z);
	
	ncell_new=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);
	
	if (ncell_old==ncell_new)
		return;
	
	int *old_item,*new_item;
	
	// delete the old position
	int j=(l->HoC)[ncell_old];
	
	old_item=l->HoC+ncell_old;
	
	while (j!=-1)
	{
		new_item=l->LinkedList+j;
		
		if (j==num)
		{
			*old_item=*new_item;
			break;
		}
		
		j=(l->LinkedList)[j];
		
		old_item=new_item;
	}
	
	// you should never be here because there must be a cancellation
	// in the list. This means that there are some problems in the list
	assert(j!=-1);
	
	// add the new position
	(l->LinkedList)[num]=(l->HoC)[ncell_new];
	(l->HoC)[ncell_new]=num;
	(l->MyCell)[num]=ncell_new;
}

void calculateExtendedInteractionMapWithCutoffDistance(listcell *l,interactionmap *im,interactionmap *ime,vector *pos,double Box[],double cutoff)
{
	
	int i;
	int nn,nnn;                      // Number of cells on surface and in box
	int particle1,particle2;         // interacting particles
	int neighbour;                   // Neighbour cell
	int ix,iy,iz,alpha;              // present cell coordinates
	int jx,jy,jz;                    // next cell coordinates
	int kx,ky,kz;                    // previous cell coordinates
	
	vector dist,olddist;
	double dist2;
	
	nn=l->NumberCells_x*l->NumberCells_y;
	nnn=l->NumberCells_z*nn;
	
	
	int number_particles=0;
	
	double cutoff2=SQR(cutoff);
	
	im->num_bonds=0;
	
	// scan cells
	for (i=0;i<nnn;i++)
	{
		// reconstruct cartesian coordinates
		iz=i/nn;
		alpha=i%nn;
		ix=alpha%l->NumberCells_x;
		iy=alpha/l->NumberCells_x;
		
		jx=(ix+1)%l->NumberCells_x;
		jy=(iy+1)%l->NumberCells_y;
		jz=(iz+1)%l->NumberCells_z;
		
		kx=module(ix-1,l->NumberCells_x);
		ky=module(iy-1,l->NumberCells_y);
		kz=module(iz-1,l->NumberCells_z);
		
		
		particle1=(l->HoC)[i];
		
		while (particle1!=-1)
		{
			//interactions with particles in the same cell
			// and in neighbouring cells
			
			im->who[particle1]=particle1;
			ime->who[particle1]=particle1;
			int howmany_particle1=0;
			number_particles++;
			
			// 0 same cell
			particle2=(l->LinkedList)[particle1];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
					
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 1 cell y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 2 cell x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 3 cell x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 4 cell x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 5 cell z+1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 6 cell z+1
			neighbour=(ix)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 7 cell z+1 x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 8 cell z+1 x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 9 cell z+1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 10 cell z-1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 11 cell x+1 y+1 z-1
			neighbour=(jx)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 12 cell x+1 z-1
			neighbour=(jx)+(iy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 13 cell z-1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					im->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			im->howmany[particle1]=howmany_particle1;
			
			particle1=(l->LinkedList)[particle1];
		}
	}
	
	im->num=number_particles;
	
	ime->num_bonds=im->num_bonds;
	
	assert(number_particles==im->size);
}

static void insertionSort_dist2(int *number,int num_el,double *dist2,double dist2_el,int *length)
{
	int l=*length;
	number[l]=num_el;
	dist2[l]=dist2_el;
	while ((l>0) && (dist2[l]<dist2[l-1]))
	{
		double buffer_dist2;
		buffer_dist2=dist2[l];
		dist2[l]=dist2[l-1];
		dist2[l-1]=buffer_dist2;
		
		int buffer_num;
		buffer_num=number[l];
		number[l]=number[l-1];
		number[l-1]=buffer_num;
		
		l--;
	}
	(*length)++;
}

void calculateInteractionMapWithCutoffDistanceOrdered(listcell *l,interactionmap *ime,vector *pos,double Box[],double cutoff)
{
	
	
	vector olddist,dist;
	double dist2;
	
	int i;
	int nn,nnn;                      // Number of cells on surface and in box
	int particle1,particle2;         // interacting particles
	int neighbour;                   // Neighbour cell
	int ix,iy,iz,alpha;              // present cell coordinates
	int jx,jy,jz;                    // next cell coordinates
	int kx,ky,kz;                    // previous cell coordinates
	
	
	nn=l->NumberCells_x*l->NumberCells_y;
	nnn=l->NumberCells_z*nn;
	
	
	int number_particles=0;
	
	double cutoff2=SQR(cutoff);
	
	ime->num_bonds=0;
	
	
	// abbiamo bisogno della mappa delle distanze
	// usiamo ime->rij2
	/*
	int ncolloids=ime->size;
	double **dist2map;
	Matrix2D(dist2map,ncolloids,ncolloids,double);
	*/
	
	
	// scan cells
	for (i=0;i<nnn;i++)
	{
		// reconstruct cartesian coordinates
		iz=i/nn;
		alpha=i%nn;
		ix=alpha%l->NumberCells_x;
		iy=alpha/l->NumberCells_x;
		
		jx=(ix+1)%l->NumberCells_x;
		jy=(iy+1)%l->NumberCells_y;
		jz=(iz+1)%l->NumberCells_z;
		
		kx=module(ix-1,l->NumberCells_x);
		ky=module(iy-1,l->NumberCells_y);
		kz=module(iz-1,l->NumberCells_z);
		
		
		particle1=(l->HoC)[i];
		
		while (particle1!=-1)
		{
			//interactions with particles in the same cell
			// and in neighbouring cells
			ime->who[particle1]=particle1;
			number_particles++;
			
			// 0 same cell
			particle2=(l->LinkedList)[particle1];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 1 cell y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 2 cell x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 3 cell x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 4 cell x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 5 cell z+1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 6 cell z+1
			neighbour=(ix)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 7 cell z+1 x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 8 cell z+1 x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 9 cell z+1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 10 cell z-1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 11 cell x+1 y+1 z-1
			neighbour=(jx)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 12 cell x+1 z-1
			neighbour=(jx)+(iy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 13 cell z-1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					insertionSort_dist2(ime->with[particle1],particle2,ime->rij2[particle1],dist2,ime->howmany+particle1);
					insertionSort_dist2(ime->with[particle2],particle1,ime->rij2[particle2],dist2,ime->howmany+particle2);
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			particle1=(l->LinkedList)[particle1];
		}
	}
	
	ime->num=number_particles;
	
	//assert(number_particles==ime->size);
	
	//Free2D(dist2map);
}

int calculateSystemInteractionMapExtended(listcell *l,interactionmap *ime,vector *pos,double Box[],double cutoff)
{
	int i;
	int nn,nnn;                      // Number of cells on surface and in box
	int particle1,particle2;         // interacting particles
	int neighbour;                   // Neighbour cell
	int ix,iy,iz,alpha;              // present cell coordinates
	int jx,jy,jz;                    // next cell coordinates
	int kx,ky,kz;                    // previous cell coordinates
	
	
	nn=l->NumberCells_x*l->NumberCells_y;
	nnn=l->NumberCells_z*nn;
	
	vector olddist,dist;
	double dist2;
	
	int number_particles=0;
	
	double cutoff2=SQR(cutoff);
	
	ime->num_bonds=0;
	
	// scan cells
	for (i=0;i<nnn;i++)
	{
		// reconstruct cartesian coordinates
		iz=i/nn;
		alpha=i%nn;
		ix=alpha%l->NumberCells_x;
		iy=alpha/l->NumberCells_x;
		
		jx=(ix+1)%l->NumberCells_x;
		jy=(iy+1)%l->NumberCells_y;
		jz=(iz+1)%l->NumberCells_z;
		
		kx=module(ix-1,l->NumberCells_x);
		ky=module(iy-1,l->NumberCells_y);
		kz=module(iz-1,l->NumberCells_z);
		
		
		particle1=(l->HoC)[i];
		
		while (particle1!=-1)
		{
			//interactions with particles in the same cell
			// and in neighbouring cells
			ime->who[particle1]=particle1;
			number_particles++;
			
			
			// 0 same cell
			particle2=(l->LinkedList)[particle1];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					
					
					ime->num_bonds++;
					
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 1 cell y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 2 cell x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 3 cell x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 4 cell x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 5 cell z+1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 6 cell z+1
			neighbour=(ix)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 7 cell z+1 x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 8 cell z+1 x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 9 cell z+1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 10 cell z-1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 11 cell x+1 y+1 z-1
			neighbour=(jx)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 12 cell x+1 z-1
			neighbour=(jx)+(iy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 13 cell z-1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					int howmany1=ime->howmany[particle1];
					int howmany2=ime->howmany[particle2];
					
// 					ime->rij[particle2][howmany2].x=dist.x;
// 					ime->rij[particle2][howmany2].y=dist.y;
// 					ime->rij[particle2][howmany2].z=dist.z;
					ime->rij2[particle2][howmany2]=dist2;
					
// 					ime->rij[particle1][howmany1].x=-dist.x;
// 					ime->rij[particle1][howmany1].y=-dist.y;
// 					ime->rij[particle1][howmany1].z=-dist.z;
					ime->rij2[particle1][howmany1]=dist2;
					
					ime->with[particle1][ime->howmany[particle1]++]=particle2;
					ime->with[particle2][ime->howmany[particle2]++]=particle1;
					
					ime->num_bonds++;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			particle1=(l->LinkedList)[particle1];
		}
	}
	
	ime->num=number_particles;
	
	assert(number_particles==ime->size);
	
	int max_neighbours=0;
	for (i=0;i<number_particles;i++)
	{
		if (ime->howmany[i]>max_neighbours)
			max_neighbours=ime->howmany[i];
	}
	
	return max_neighbours;
}

int calculateSystemInteractionMap(listcell *l,interactionmap *im,vector *pos,double Box[],double cutoff)
{
	
	int i;
	int nn,nnn;                      // Number of cells on surface and in box
	int particle1,particle2;         // interacting particles
	int neighbour;                   // Neighbour cell
	int ix,iy,iz,alpha;              // present cell coordinates
	int jx,jy,jz;                    // next cell coordinates
	int kx,ky,kz;                    // previous cell coordinates
	
	vector olddist,dist;
	double dist2;
	
	nn=l->NumberCells_x*l->NumberCells_y;
	nnn=l->NumberCells_z*nn;
	
	
	int number_particles=0;
	
	double cutoff2=SQR(cutoff);
	
	im->num_bonds=0;
	
	// scan cells
	for (i=0;i<nnn;i++)
	{
		// reconstruct cartesian coordinates
		iz=i/nn;
		alpha=i%nn;
		ix=alpha%l->NumberCells_x;
		iy=alpha/l->NumberCells_x;
		
		jx=(ix+1)%l->NumberCells_x;
		jy=(iy+1)%l->NumberCells_y;
		jz=(iz+1)%l->NumberCells_z;
		
		kx=module(ix-1,l->NumberCells_x);
		ky=module(iy-1,l->NumberCells_y);
		kz=module(iz-1,l->NumberCells_z);
		
		
		particle1=(l->HoC)[i];
		
		while (particle1!=-1)
		{
			//interactions with particles in the same cell
			// and in neighbouring cells
			im->who[particle1]=particle1;
			number_particles++;
			int howmany_particle1=0;
			
			// 0 same cell
			particle2=(l->LinkedList)[particle1];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 1 cell y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 2 cell x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 3 cell x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 4 cell x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 5 cell z+1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 6 cell z+1
			neighbour=(ix)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 7 cell z+1 x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 8 cell z+1 x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 9 cell z+1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 10 cell z-1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 11 cell x+1 y+1 z-1
			neighbour=(jx)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 12 cell x+1 z-1
			neighbour=(jx)+(iy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			// 13 cell z-1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->rij2[particle1][howmany_particle1]=dist2;
					im->rij[particle1][howmany_particle1]=dist;
					im->with[particle1][howmany_particle1++]=particle2;
				}
				particle2=(l->LinkedList)[particle2];
			}
			
			im->howmany[particle1]=howmany_particle1;
			particle1=(l->LinkedList)[particle1];
		}
	}
	
	im->num=number_particles;
	
	
	//assert(number_particles==im->size);
	
	int max_neighbours=0;
	for (i=0;i<number_particles;i++)
	{
		if (im->howmany[i]>max_neighbours)
			max_neighbours=im->howmany[i];
	}
	
	return max_neighbours;
}

void addToList(listcell *l,const vector *pos,int num)
{
	int ncell;                       // cell number
	int posx,posy,posz;              // cell coordinates
	
	
	posx=(int)floor(pos->x/l->CellSize_x);
	posy=(int)floor(pos->y/l->CellSize_y);
	posz=(int)floor(pos->z/l->CellSize_z);
	
	posx=module(posx,l->NumberCells_x);
	posy=module(posy,l->NumberCells_y);
	posz=module(posz,l->NumberCells_z);
	
	ncell=posx+posy*l->NumberCells_x+posz*(l->NumberCells_x*l->NumberCells_y);
	
	(l->LinkedList)[num]=(l->HoC)[ncell];
	(l->HoC)[ncell]=num;
	(l->MyCell)[num]=ncell;
}


void removeFromList(listcell *l,int num)
{
	int *current_item,*next_item;
	int ncell_old=l->MyCell[num];
	
	int item=(l->HoC)[ncell_old];
	current_item=l->HoC+ncell_old;
	
	while (item!=num)
	{
		current_item=l->LinkedList+item;
		item=(l->LinkedList)[item];
	}
	
#ifdef DEBUG
	assert(item!=-1);
#endif
	
	next_item=l->LinkedList+item;
	*current_item=*next_item;
	
// 	l->MyCell[num]=-1;
}


void changeIdentityInList(listcell *l,int oldnum,int newnum)
{
	int cell=l->MyCell[oldnum];
	
	// a particle oldnum changes identity in newnum
	
	removeFromList(l,oldnum);
	removeFromList(l,newnum);
	
	// ora la newnum deve prendere il posto della oldnum
	(l->LinkedList)[newnum]=(l->HoC)[cell];
	(l->HoC)[cell]=newnum;
	(l->MyCell)[newnum]=cell;
}


void calculateLinksWithinCutoff(listcell *l,vector *pos,double Box[],double cutoff,interactionmap *im,int *num_links)
{
	
	int i;
	int nn,nnn;                      // Number of cells on surface and in box
	int particle1,particle2;         // interacting particles
	int neighbour;                   // Neighbour cell
	int ix,iy,iz,alpha;              // present cell coordinates
	int jx,jy,jz;                    // next cell coordinates
	int kx,ky,kz;                    // previous cell coordinates
	
	
	nn=l->NumberCells_x*l->NumberCells_y;
	nnn=l->NumberCells_z*nn;
	
	vector olddist,dist;
	double dist2;
	
	int number_particles=0;
	
	double cutoff2=SQR(cutoff);
	
	// scan cells
	for (i=0;i<nnn;i++)
	{
		// reconstruct cartesian coordinates
		iz=i/nn;
		alpha=i%nn;
		ix=alpha%l->NumberCells_x;
		iy=alpha/l->NumberCells_x;
		
		jx=(ix+1)%l->NumberCells_x;
		jy=(iy+1)%l->NumberCells_y;
		jz=(iz+1)%l->NumberCells_z;
		
		kx=module(ix-1,l->NumberCells_x);
		ky=module(iy-1,l->NumberCells_y);
		kz=module(iz-1,l->NumberCells_z);
		
		particle1=(l->HoC)[i];
		
		while (particle1!=-1)
		{
			//interactions with particles in the same cell
			// and in neighbouring cells
			im->who[particle1]=particle1;
			number_particles++;
			int howmany_particle1=0;
			
			// 0 same cell
			particle2=(l->LinkedList)[particle1];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 1 cell y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 2 cell x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 3 cell x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 4 cell x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(iz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 5 cell z+1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 6 cell z+1
			neighbour=(ix)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 7 cell z+1 x+1 y+1
			neighbour=(jx)+(jy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 8 cell z+1 x+1
			neighbour=(jx)+(iy)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 9 cell z+1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(jz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 10 cell z-1 y+1
			neighbour=(ix)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 11 cell x+1 y+1 z-1
			neighbour=(jx)+(jy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 12 cell x+1 z-1
			neighbour=(jx)+(iy)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			// 13 cell z-1 x+1 y-1
			neighbour=(jx)+(ky)*l->NumberCells_x+(kz)*nn;
			particle2=(l->HoC)[neighbour];
			while (particle2!=-1)
			{
				olddist.x=pos[particle1].x-pos[particle2].x;
				olddist.y=pos[particle1].y-pos[particle2].y;
				olddist.z=pos[particle1].z-pos[particle2].z;
				
				olddist.x-=rint(olddist.x);
				olddist.y-=rint(olddist.y);
				olddist.z-=rint(olddist.z);
				
				dist.x=Box[0]*olddist.x+Box[1]*olddist.y+Box[2]*olddist.z;
				dist.y=Box[3]*olddist.y+Box[4]*olddist.z;
				dist.z=Box[5]*olddist.z;
				
				dist2=SQR(dist.x)+SQR(dist.y)+SQR(dist.z);
				
				if (dist2<cutoff2)
				{
					im->with[particle1][howmany_particle1++]=particle2;
					num_links[particle1]++;
					num_links[particle2]++;
				}
				
				particle2=(l->LinkedList)[particle2];
			}
			
			im->howmany[particle1]=howmany_particle1;
			
			particle1=(l->LinkedList)[particle1];
		}
	}
	
	im->num=number_particles;
}

vector getCellSize(listcell *l)
{
	vector size;
	size.x=l->CellSize_x;
	size.y=l->CellSize_y;
	size.z=l->CellSize_z;
	return size;
}

listcell *selfList()
{
	return ListMc;
}

