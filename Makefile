CFLAGS=-O2 -Wunused

logft: logft.o
	cc -o $@ $< -lsndfile -lpng -lm

logft.o: logft.c frmhdr.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f logft logft.o
