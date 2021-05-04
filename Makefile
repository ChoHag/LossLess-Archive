CTANGLE?= ctangle
CWEAVE?=  cweave
PDFTEX?=  pdftex
CFLAGS+=  -Wall -Wpedantic -Wextra

SOURCES:=       lossless.c repl.c
OBJECTS:=       lossless.o repl.o
TESTS:=         \
	t/cell-heap.t \
	t/environments.t \
	t/eval.t \
	t/exception.t \
	t/gc-mark.t \
	t/gc-sweep.t \
	t/gc-vector.t \
	t/if.t \
	t/lambda.t \
	t/pair.t \
	t/sanity.t \
	t/vector-heap.t \
	t/vov.t
TEST_SOURCES:= \
	t/cell-heap.c \
	t/environments.c \
	t/eval.c \
	t/exception.c \
	t/gc-mark.c \
	t/gc-sweep.c \
	t/gc-vector.c \
	t/if.c \
	t/lambda.c \
	t/pair.c \
	t/sanity.c \
	t/vector-heap.c \
	t/vov.c

ALLOC_TESTS:= t/cell-heap.t t/vector-heap.t

TEST_OBJECTS:= llalloc.o lltest.o

OTHER_SOURCES:= \
	lossless.h \
	t/llalloc.c \
	t/llt.h \
	t/lltest.c

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

$(SOURCES) $(TEST_SOURCES) $(OTHER_SOURCES): t lossless.w
	$(CTANGLE) lossless.w

test: $(TESTS)
	prove -vr -e '' t

$(TESTS): $(TEST_OBJECTS)

$(TEST_OBJECTS): $(TEST_SOURCES) $(OTHER_SOURCES)

t:
	mkdir -p t

.SUFFIXES: .t

.c.t:
	if echo " $(ALLOC_TESTS) " | grep -qF " $@ "; then                     \
		$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -I. llalloc.o -o $@ $<; \
	else                                                                   \
		$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -I. lltest.o -o $@ $<;  \
	fi

clean:
	rm -f core *.core *.idx *.log *.scn *.toc *.o
	rm -f lossless *.o
	rm -f repl.c lossless.c lossless.tex lossless.pdf
	rm -fr t
