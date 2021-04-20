CTANGLE?= ctangle
CWEAVE?=  cweave
PDFTEX?=  pdftex
CFLAGS+=  -Wall -Wpedantic -Wextra

SOURCES:=       lossless.c repl.c
OBJECTS:=       lossless.o repl.o
TESTS:=         \
	t/allocator.t \
	t/eval.t \
	t/exception.t \
	t/if.t \
	t/lambda.t \
	t/pair.t \
	t/sanity.t \
	t/vov.t
TEST_SOURCES:= \
	t/allocator.c \
	t/eval.c \
	t/exception.c \
	t/if.c \
	t/lambda.c \
	t/pair.c \
	t/sanity.c \
	t/vov.c

all: lossless lossless.pdf

full: test all

lossless.pdf: lossless.tex
	$(PDFTEX) lossless.tex

lossless.tex: lossless.w
	$(CWEAVE) lossless.w

lossless: $(OBJECTS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o lossless $(OBJECTS)
	strip lossless

$(OBJECTS): $(SOURCES)

$(SOURCES) $(TEST_SOURCES): lossless.w
	$(CTANGLE) lossless.w

test: $(TESTS)
	prove -vr -e '' t

$(TESTS): t $(TEST_SOURCES)

t:
	mkdir -p t

.SUFFIXES: .t

.c.t:
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -I. -o $@ $<

clean:
	rm -f core *.core *.idx *.log *.scn *.toc *.o
	rm -f lossless *.o
	rm -f repl.c lossless.c lossless.tex lossless.pdf
	rm -f t/*.c t/*.o t/*.t
