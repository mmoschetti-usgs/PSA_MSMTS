INST_DIR = $(HOME)/bin
BIN  = sac2ascii_c

# need to include the math library when the math.h header file is used
cflags = -lm
CFLAGS = $(DBG) $(cflags)

CC = gcc

DBG =
FOBJS = sac2ascii_c.c read_sac64.c 

$(BIN) : $(FOBJS)
	$(CC) $(CFLAGS) $(FOBJS) -o $(BIN)

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core.*  
