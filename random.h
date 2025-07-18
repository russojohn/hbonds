#ifndef RANDOM_H
#define RANDOM_H

gsl_rng* randomConstructor(FILE *config_file);

gsl_rng* randomConstructorInteractive(int seed);

void randomFree();

gsl_rng* randomStructure();

void saveRandom(steps step);

steps getRandom();


#endif

