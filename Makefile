lossless: lossless.c
	cc $(CFLAGS) -Wall -Wpedantic -o lossless lossless.c

lltest: lossless.c
	cc $(CFLAGS) -DLL_TEST -Wall -Wpedantic -o lltest lossless.c

lossless.c: lossless.w
	ctangle lossless.w

lossless.pdf: lossless.tex
	pdftex lossless.tex

lossless.tex: lossless.w
	cweave lossless.w

clean:
	rm -f *.core *.idx *.log *.scn *.toc *.o
	rm -f lossless.c lossless lltest lossless.tex lossless.pdf
