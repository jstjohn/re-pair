CC=g++

SRC=re-pair.cpp
CFLAGS=-O3
re-pair: $(SRC)
	$(CC) $(CFLAGS) -o $@ $<

PREFIX=$(HOME)/bin

install: re-pair
	cp $< $(PREFIX)

