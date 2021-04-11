CTANGLE?= ctangle
CWEAVE?=  cweave
PDFTEX?=  pdftex
CFLAGS+=  -Wall -Wpedantic -Wextra

lossless: lossless.c
	$(CC) $(CFLAGS) -o lossless lossless.c

# If lltest doesn't depend on lossless it can result in an old llbuild
# of lltest not being replaced under some conditions.
lltest: lossless.c lossless
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
	rm -f core *.core *.idx *.log *.scn *.toc *.o
	rm -f lossless.c lossless lltest lossless.tex lossless.pdf
	$(MAKE) -C t clean
