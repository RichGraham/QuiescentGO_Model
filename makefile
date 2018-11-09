# simple make file (Migrate.c)
# September 2004
#HOME= /Users/rsgraham
#NR=$(HOME)/recipes_c-ansi
#GSL=$(HOME)/GNU_Library


VPATH = ./: $(NR)/misc/ :  $(NR)/recipes/


SOURCES= main.cpp ran1.c 


OBJECTS= main.o ran1.o	

PRODUCT=a.out

CC=g++


CPPFLAGS= -O3 #-I$(NR)/include -O3 #-I$(GSL)/include #-L$(GSL)/lib -O2 #-fast -B 
CFLAGS= -O3 #-I$(NR)/include -O3  #-I$(GSL)/include #-L$(GSL)/lib -O2 #-fast -B 



all: $(PRODUCT)
$(PRODUCT): $(OBJECTS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $(PRODUCT) $(OBJECTS) -lm 

.c.o:	
	$(CC) $(CFLAGS)  $(CPPFLAGS)-c $<

.cpp.o:	
	g++ $(CFLAGS)  $(CPPFLAGS)-c $< 

clean:
	rm -f *.o
