MF=	Makefile

CC=	mpicc
CFLAGS= -std=c99 -O3

LFLAGS=	-lm

EXE=	imagempi

SRC= \
	imagempi.c \
	pgmio.c \
	arralloc.c

INC=\
	pgmio.h \
	arralloc.h

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(INC)

$(OBJ):	$(MF)

clean:
	rm -f $(OBJ) $(EXE) core
