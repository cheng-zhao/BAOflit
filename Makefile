CC = gcc
FC = gfortran
CFLAGS = -std=c99 -O3 -Wall
FFLAGS = -ffree-line-length-none -O2

MNDIR = MultiNest
CLIBS = -lm -lnest3 -L$(MNDIR)
FLIBS = -llapack

# Settings for LAPACK
LAPACK_DIR =
ifneq ($(LAPACK_DIR),)
  FLIBS += -L$(LAPACK_DIR)
endif

INCL = -Ilib -Imath -Isrc
SRCS = $(wildcard */*.c)
OBJS = $(SRCS:.c=.o)
EXEC = BAOflit

export FC FFLAGS FLIBS

all: mn baoflit

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $*.c $(CLIBS) $(INCL)

mn:
	make -C $(MNDIR)

baoflit: $(OBJS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS) $(CLIBS) $(FLIBS)

clean:
	rm $(EXEC) $(OBJS)

cleanmn:
	make -C $(MNDIR) clean

cleanall: cleanmn clean

