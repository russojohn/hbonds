#ifndef SW_H
#define SW_H

typedef struct _sw {
	int *NColloids;

	double **TwoBody;
	double **TwoBody_copy;
	double **TwoBody_savestate;

	int *WhoCopied;
	int Num_copies;

	bilista *List_modified;

	// parameters
	double A;
	double B;
	int p;
	int q;
	double Gamma;
	double *costeta0;
	double lambda;
	double a;
} sw;

sw* createSW(int max_colloids,int *ncolloids,double lambda);
void freeSW(sw *sws);

void swConstructor(FILE *config_file,vector *pos,int max_colloids,int *ncolloids,double box[],double *cutoff,interactionmap *interactionList,double *potential,double *virial,vector *force);
void swFree();
void sw_system(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial);
void sw_particle_oldconfiguration(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial);
void sw_particle_newconfiguration(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial);
void sw_swap_oldconfiguration(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial);
void sw_swap_newconfiguration(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial);
void sw_particle_rejectmove(sw *const sws);
void sw_system_rejectmove(sw *const sws);
double sw_particle_ghost_cheat(interactionmap *im);

void sw_remove_rejectmove(sw *const sws);
void sw_remove_acceptmove(sw *const sws,int particle,int ncolloids);
void sw_add_rejectmove(sw *const sws,int ncolloids);
void sw_particle_add(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *pot,double *virial);
void sw_particle_remove(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial);

void sw_check(vector *pos,interactionmap *ime,double Box[],double cutoff,double *potential,double *virial,double *twobody,double *threebody);

void sw_check_force(sw *const sws,vector *pos,interactionmap *ime,double Box[],double cutoff,double *potential,double *virial,vector *force);
void sw_system_force(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial,vector *force);
void sw_system_virial(sw *const sws,vector *pos,interactionmap *im,double Box[],double cutoff,double *potential,double *virial,vector *force);

void sw_save_state(sw *const sws);
void sw_load_state(sw *const sws);

void getSW(FILE *pfile);
void saveSW(FILE *pfile);

void swPrintHamiltonian(steps t);

void swGetTwoThreeBodyParts(double *potential2body,double *potential3body);

void bonds_energy(sw *sws,interactionmap *im,double *energy);

sw* sw_pointer();

#endif
