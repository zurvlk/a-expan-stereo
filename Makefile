#Makefile

# CXX = icc
# CXXFLAGS = -Wall -parallel -lm -fPIC

CC = gcc
all:a-expan

a-expan:graph.o bmp.o ford_fulkerson.o queue_stack.o a_edrv.o a-expan.o
	$(CC) graph.o bmp.o ford_fulkerson.o queue_stack.o a_edrv.o a-expan.o -o a-expan -Wall -lm

a-estr:graph.o bmp.o ford_fulkerson.o queue_stack.o str_drv.o a_estr.o
	$(CC) graph.o bmp.o ford_fulkerson.o queue_stack.o str_drv.o a_estr.o -o a-estr -Wall -lm

.c.o:
	$(CC) -c $< -Wall -lm
remove:
	rm -f *.o
clean:
	rm -f *.o *~ a-expan a-estr
