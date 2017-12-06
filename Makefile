#Makefile

# CXX = icc
# CXXFLAGS = -Wall -parallel -lm -fPIC

CC = gcc
all:a-estr

a-estr:graph.o bmp.o ford_fulkerson.o queue_stack.o str_drv.o a_estr.o
	$(CC) graph.o bmp.o ford_fulkerson.o queue_stack.o str_drv.o a_estr.o -o a-estr -Wall -lm

.c.o:
	$(CC) -c $< -Wall -lm
remove:
	rm -f *.o
clean:
	rm -f *.o *~ a-expan a-estr
