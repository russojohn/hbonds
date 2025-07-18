## Edit here to the gsl library location

GSL_INCLUDE = /opt/homebrew/Cellar/gsl/2.7.1/include

GSL_LIB = /opt/homebrew/Cellar/gsl/2.7.1/lib

#############


OPTIMIZATION_FLAG = -O3 -ffast-math -Wno-unused-variable

SRC = .

OBJBONDSMAPTOPOLOGIC = bonds_map_topologic.o io.o ordinator.o cell_list.o bilista.o sw_disorder.o log.o secure_search.o interaction_map.o random.o restart.o utilities.o events.o cluster.o dirutils.o

# ----- COMPILER -----#
CC = gcc -Wall

#CC = icc

# --- PREPROCESSOR FLAGS --- #
CPPFLAGS = $(INCLUDE_DIRS) -I$(SRC)

# -- CC FLAGS   ------#
CFLAGS = $(OPTIMIZATION_FLAG) $(ARCH) -I$(GSL_INCLUDE)

# --- LINKER FLAGS  -------#
LFLAGS = $(LIBRARY_DIRS) $(STATIC_FLAG) -L$(GSL_LIB)


# ---------  EXECUTABLE -----------#



EXEBONDSMAPTOPOLOGIC = ./bonds_map_topologic

all: bondsmap

bondsmap: $(OBJBONDSMAPTOPOLOGIC)
	$(CC) $(OBJBONDSMAPTOPOLOGIC) $(LFLAGS) -lgsl -lgslcblas -lm -o $(EXEBONDSMAPTOPOLOGIC)



# -------- SUFFIX RULES ------#
.c.o: local
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

io.o: $(SRC)/io.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

cluster.o: $(SRC)/cluster.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

ordinator.o: $(SRC)/ordinator.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

cell_list.o: $(SRC)/cell_list.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

verlet_list.o: $(SRC)/verlet_list.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

bilista.o: $(SRC)/bilista.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

order_parameters.o: $(SRC)/order_parameters.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

secure_search.o: $(SRC)/secure_search.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

log.o: $(SRC)/log.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

interaction_map.o: $(SRC)/interaction_map.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

random.o: $(SRC)/random.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

restart.o: $(SRC)/restart.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

utilities.o: $(SRC)/utilities.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

events.o: $(SRC)/events.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?

sw_disorder.o: $(SRC)/sw_disorder.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $?




clean:
	rm -f $(EXEBONDSMAPTOPOLOGIC) *.o
