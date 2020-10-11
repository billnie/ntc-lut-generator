CC=gcc
OC=objcopy
ODIR=thermistor1.0

CFLAGS += -DDEBUG=1
all: 	lut rtot ttor coeff rtd clean


lut:	libcoeff.o librtot.o thermistor.h lut.c
	$(CC) lut.c -o lut libcoeff.o librtot.o -lm

rtd: rtd.c
	$(CC) rtd.c -o rtd

ttor :  $(ODIR)/ttor.c
	$(CC) $(ODIR)/ttor.c -o ttor -lm

rtot :  $(ODIR)/rtot.c
	$(CC) $(ODIR)/rtot.c -o rtot -lm	

coeff :  $(ODIR)/coeff.c
	$(CC) $(ODIR)/coeff.c -o coeff -lm	

clean:
	rm -f $(ODIR)/coeff.o $(ODIR)/rtot.o $(ODIR)/ttor.o libcoeff.o librtot.o rtd.o


$(ODIR)/coeff.o: $(ODIR)/coeff.c
	$(CC) -c $(ODIR)/coeff.c -o $(ODIR)/coeff.o -lm

$(ODIR)/rtot.o: $(ODIR)/rtot.c
	$(CC) -c $(ODIR)/rtot.c -o $(ODIR)/rtot.o -lm

libcoeff.o: $(ODIR)/coeff.o
	$(OC) -N main $(ODIR)/coeff.o libcoeff.o

librtot.o: $(ODIR)/rtot.o
	$(OC) -N main $(ODIR)/rtot.o librtot.o

