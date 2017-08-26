CC=g++

ODIR=obj

IDIR_RECON=../include/
INCLUDES=-I$(IDIR_RECON)

LOCAL_INCLUDES=recon_structs.h

LIBS=-lm -lyaml-cpp -lctbb_read

CFLAGS=-g -DHAVE_CONFIG_H -std=gnu++11

_OBJ=main.o parse_config.o setup.o generate_system_matrix.o generate_system_matrix_ffs.o initialize_recon_volume.o icd_iteration.o rotate_slices.o penalties.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ));

all: ../icd

../icd: object_dir $(OBJ) 
	(cd obj && $(CC) -g -fopenmp -o ../../icd *.o $(LIBS))

object_dir:
	mkdir -p obj

obj/main.o: main.cpp recon_structs.h	
	$(CC) $(CFLAGS) -c -o $@ $<

obj/parse_config.o: parse_config.cpp parse_config.h recon_structs.h
	$(CC) $(CFLAGS) -c -o $@ $<

obj/setup.o: setup.cpp setup.h recon_structs.h parse_config.h
	$(CC) $(CFLAGS) -c -o $@ $<

obj/generate_system_matrix.o: generate_system_matrix.cpp generate_system_matrix.h recon_structs.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -fopenmp -o $@ $<

obj/generate_system_matrix_ffs.o: generate_system_matrix_ffs.cpp generate_system_matrix_ffs.h recon_structs.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -fopenmp -o $@ $<

obj/initialize_recon_volume.o: initialize_recon_volume.cpp initialize_recon_volume.h recon_structs.h
	$(CC) $(CFLAGS) -c -o $@ $<

obj/icd_iteration.o: icd_iteration.cpp icd_iteration.h recon_structs.h
	$(CC) $(CFLAGS) -c -fopenmp -o $@ $<

obj/rotate_slices.o: rotate_slices.cpp rotate_slices.h recon_structs.h interp.h
	$(CC) $(CFLAGS) -c -o $@ $<

obj/penalties.o: penalties.cpp penalties.h
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: all test clean object_dir

clean:
	rm obj/*.o; rm -rf ../*~; rm -f ../icd

