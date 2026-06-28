EXE		= dtmf2num
CFLAGS	+= -O2 -s
PREFIX	= /usr/local
BINDIR	= $(PREFIX)/bin
LIBS	= -lm

all:
	$(CC) $(CFLAGS) -o $(EXE) dtmf2num.c $(LIBS)

install:
	install -m 755 -d $(BINDIR)
	install -m 755 $(EXE) $(BINDIR)/$(EXE)

.PHONY:
	install
