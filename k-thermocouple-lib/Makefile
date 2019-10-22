CC=gcc
OC=objcopy
ODIR=./

NAME = k-c
STATIC_LIB=lib$(NAME).a

CFLAGS += -DDEBUG=1
all: 	libk-c clean


libk-c: 
	$(CC) thermocouple.c -c 
	ar -cr $(STATIC_LIB) thermocouple.o



clean:
	rm -f *.o





#libk-c.o: thermocouple.o
#	$(OC) -N main thermocouple.o 

