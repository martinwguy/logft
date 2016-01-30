CFLAGS=-O2 -Wunused

logft: logft.o
	cc -o $@ $< -lm

logft.o: logft.c frmhdr.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f logft logft.o
