CC      = gcc
INCLUDE = -I../../../lmip
CCFLAGS = -Wall -pedantic -O2
LIBS    = -lm
TARGETS = image_ced_LP

all:
	make clean
	make $(TARGETS)

.c:
	$(CC) $(INCLUDE) $(CCFLAGS) $@.c $(LIBS) -o $@

clean:
	rm -f $(TARGETS)

cleanbkp:
	rm -f *.c~
