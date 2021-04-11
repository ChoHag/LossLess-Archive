CTANGLE?= ctangle
CWEAVE?=  cweave
PDFTEX?=  pdftex
CFLAGS+=  -Wall -Wpedantic -Wextra

lossless: lossless.c
	$(CC) $(CFLAGS) -o lossless lossless.c

lltest: lossless.c
	$(CC) $(CFLAGS) -DLL_TEST -o lltest lossless.c

lossless.pdf: lossless.tex
	$(PDFTEX) lossless.tex

all: lossless lossless.pdf

lossless.c: lossless.w
	$(CTANGLE) lossless.w

lossless.tex: lossless.w
	$(CWEAVE) lossless.w

# test depends on both binaries to guard against testing one thing and
# shipping another.
test: lossless lltest
	PATH=..:$$PATH $(MAKE) -C t

clean:
	rm -f *.core *.idx *.log *.scn *.toc *.o
	rm -f lossless.c lossless lltest lossless.tex lossless.pdf
	$(MAKE) -C t clean
