CFLAGS=-O2 -Wunused

logft: logft.o
	cc -o $@ $< -lsndfile -lpng -lm

logft.o: logft.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f logft logft.o
