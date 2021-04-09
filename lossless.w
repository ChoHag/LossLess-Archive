% Sorry human but first we need to define a bunch of TeX macros for
% pretty-printing some common things. Just ignore this bit.

\def\LL/{\.{LossLess}}
\def\point#1.{\yskip\indent#1.\quad\ignorespaces}
\def\*{$\bullet$}
\def\title{\LL/ Programming Environment}
\def\qo{$\langle\!\langle$}
\def\qc{$\rangle\!\rangle$}
\def\qclosebra/{\qo\.{]}\qc}
\def\qclosepar/{\qo\.{)}\qc}
\def\qdot/{\qo\.{.}\qc}
\def\qopenbra/{\qo\.{[}\qc}
\def\qopenpar/{\qo\.{(}\qc}
\def\qquasi/{\qo\.{`}\qc}
\def\qquote/{\qo\.{'}\qc}
\def\qunquote/{\qo\.{,}\qc}
\def\qunsplice/{\qo\.{,@@}\qc}
\def\qcomment/{\qo\.{;}\qc}
\def\qstring/{\qo\.{"}\qc}
\def\qspecial/{\qo\.{\#}\qc}
\def\qsymbol/{\qo\.{|}\qc}
\def\qsymboldots/{\qo\.{|...|}\qc}

% Ignore this bit as well because we're only working around C's and
% TeX's inability to fully work with each other.

@s cell int
@s error return

% And now for your unusual programming...

@** Introduction. \LL/ is a programming language and environment
similar to scheme. This document describes the implementation of a
\LL/ runtime written in \CEE/ and \LL/ itself will be described
elsewhere. In unambiguous cases \LL/ may be used to refer specifically
to the implementation.

% TODO: URLs
This code started off its life as \.{s9fes} by Nils M. Holm
(\.{http://t3x.org/s9fes/index.html}). After a few iterations including
being briefly ported to perl this rather different code is the result,
although at its core it follows the same design.

The structure is of a virtual machine with a single accumulator
register and a stack. There is a single entry point to the
VM---|interpret|---called after parsed source code has been put
into the accumulator, where the result will also be left.

@c
@<System headers@>@/
@h
@<Global constants@>@/
@<Type definitions@>@/
@<Function declarations@>@/
@<Global variables@>@/

@ @<Global initialisation@>=
/* There is no code here; this section exists to give it a name */

@ \LL/ has few external dependencies, primarily |stdio| and
|stdlib|, plus some obvious memory mangling functions from the
\CEE/ library there's no point in duplicating.

@<System headers@>=
#include <ctype.h>
#include <limits.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for |memset| */
#include <sys/types.h>

@ The |boolean| and |predicate| {\bf \CEE/} types are used to
distinguish between |boolean|-returning functions reporting \CEE/
truth (0 or 1) or |predicate|-returning functions reporting \LL/
truth (|FALSE| or |TRUE|). Otherwise-untyped \CEE/ macros always
report \CEE/ truth.

@d bfalse 0
@d btrue 1
@<Type definitions@>=
typedef int32_t cell;
typedef int boolean;
typedef cell predicate;

@** Error Handling. When the VM has finished then again when it
begins interpretation it establishes two jump buffers. To understand
jump buffers it's necessary to understand how \CEE/'s stack works
and we have enough stacks already.

The main thing to know is that whenever \CEE/ code calls a function
it grows its own stack with the caller's return address. When
|setjmp| is called the position in this stack is saved. When jumping
back to that position with |longjmp|, anything which has been added
to {\it\CEE/'s} stack since the corresponding call to |setjmp| is
discarded which has the effect of returning to exactly the point
in the program where the corresponding |setjmp| was called, this
time with a non-zero return value (the value that was given as an
argument to |longjmp|; this facility is not used by \LL/ for
anything and it always sends 1).

The other thing that you don't need to know is that sometimes \CEE/
compilers can make the previous paragraph a tissue of lies.

@d ERR_UNIMPLEMENTED "unimplemented"
@<Global var...@>=
jmp_buf Goto_Begin;
jmp_buf Goto_Error;

@ @<Function dec...@>=
void handle_error (char *, cell, cell) __dead;
void warn (char *, cell);

@ When the initialisation has finished \LL/ sets the |Goto_Begin|
jump buffer to return to when all else fails.

@<Initialise error handling@>=
setjmp(Goto_Begin);

@ To support user error handling a second jump buffer |Goto_Error|
is established immediately after {\it beginning} run-time computation.
Because there is no support yet for exceptions this ``handler''
will never be entered---it exists here as a placeholder and
demonstration of how \CEE/ will handle them when \LL/ does.

@<Set up run-time error handling@>=
if (setjmp(Goto_Error)) {
        Ip = -1; /* call the handler, wherever that is */
        if (Ip < 0)
                longjmp(Goto_Begin, 1);
}

@ When an error occurs if a handler has been established (how?)
then control passes immediately there, never to return. Otherwise
the error is displayed and control returns to the beginning.

Raised errors may either be a \CEE/-`string'\footnote{$^1$}{\CEE/ does not
have strings, it has pointers to memory buffers that probably contain
ASCII and might also happen to have a |NULL| in them somewhere.} when raised
by an internal process or a |symbol| when raised at run-time.

@d error(x,d) handle_error((x), NIL, (d))
@c
void
handle_error(char *message,
             cell  id,
             cell  detail)
{
        int len;

        if (!null_p(id)) {
                message = symbol_store(id);
                len = symbol_length(id);
        } else
                len = strlen(message);
        if (0) { /* handled */
                vms_push(detail);
                if (null_p(id))
                        id = sym(message);
                Acc = atom(id, detail, FORMAT_EXCEPTION);
                vms_clear();
                longjmp(Goto_Error, 1);
        }
        printf("UNHANDLED ERROR: ");
        for (; len--; message++)
                putchar(*message);
        putchar(':');
        putchar(' ');
        write_form(detail, 0);
        printf("\n");
        longjmp(Goto_Begin, 1);
}

@ Run-time errors are raised by the |OP_ERROR| opcode which passes
control to |handle_error| (and never returns).

@<Opcode implementations@>=
case OP_ERROR:@/
        handle_error(NULL, Acc, rts_pop(1));
        break; /* superfluous */

@ We additionally define |warn| here because where else is it going
to go?

@c
void
warn (char *message,
      cell  detail)
{
        printf("WARNING: %s: ", message);
        write_form(detail, 0);
        printf("\n");
}

@** Memory Management. The most commonly used data type in lisp-like
languages is the |pair|, also called a ``cons cell'' for histerical
raisins, which is a datum consisting of two equally-sized halves.
For reasons that don't bear thinking about they are called the |car|
for the ``first'' half and the |cdr| for the ``second'' half. In
this code \AM\ document, |cell| refers to each half of a |pair|.
|cell| is not used to refer to a whole cons cell in order to avoid
confusion.

A |pair| in \LL/ is stored in 2 equally-sized areas of memory.  On
64-bit x86 implementations, which are all I'm considering at the
moment, each half is 32 bits wide. Each pair additionally has an 8
bit tag (1 byte) associated with it, stored in a third array.

Internally a |pair| is represented by an offset into these memory
areas. Negative numbers are therefore available for a few global
constants.

The |pair|'s tag is treated as a bitfield. The garbage collector
uses two bits (|TAG_MARK| and |TAG_STATE|). The other 6 bits are
used to identify what data is stored in the |cell|s.

@d NIL         -1 /* Not |NULL|, but not |nil_p|/{\it nil?} either */
@d FALSE       -2 /* Yes, */
@d TRUE        -3 /* really. */
@d END_OF_FILE -4 /* stdio has |EOF| */
@d VOID        -5
@d UNDEFINED   -6
@#
@d TAG_NONE   0x00
@d TAG_MARK   0x80 /* GC mark bit */
@d TAG_STATE  0x40 /* GC state bit */
@d TAG_ACARP  0x20 /* CAR is a pair */
@d TAG_ACDRP  0x10 /* CDR is a pair */
@d TAG_FORMAT 0x3f /* Mask lower 6 bits */
@#
@d HEAP_SEGMENT 0x8000
@<Global var...@>=
cell *CAR = NULL;
cell *CDR = NULL;
char *TAG = NULL;

cell  Cells_Free = NIL;
int   Cells_Poolsize = 0;
int   Cells_Segment = HEAP_SEGMENT;

@ The pool is spread across |CAR|, |CDR| and |TAG| and starts off
with a size of zero |cell|s, growing by |Cells_Segment| |cell|s
each time it's enlarged. When the heap is enlarged newly allocated
memory is set to zero and the segment size set to half of the total
pool size.

@d ERR_OOM "out-of-memory"
@c
void
new_cells_segment(void)
{
        cell *new_car, *new_cdr;
        char *new_tag;
        new_car = reallocarray(CAR, Cells_Poolsize + Cells_Segment,
                sizeof (cell));
        new_cdr = reallocarray(CAR, Cells_Poolsize + Cells_Segment,
                sizeof (cell));
        new_tag = reallocarray(CAR, Cells_Poolsize + Cells_Segment,
                sizeof (char));
        if (new_car == NULL || new_cdr == NULL || new_tag == NULL) {
                free(new_car);
                free(new_cdr);
                free(new_tag);
                error(ERR_OOM, NIL);
        }
        bzero(((char *) new_car) + Cells_Poolsize,
                Cells_Segment * sizeof (cell));
        bzero(((char *) new_cdr) + Cells_Poolsize,
                Cells_Segment * sizeof (cell));
        bzero(((char *) new_tag) + Cells_Poolsize,
                Cells_Segment * sizeof (char));
        CAR = new_car;
        CDR = new_cdr;
        TAG = new_tag;
        Cells_Poolsize += Cells_Segment;
        Cells_Segment = Cells_Poolsize / 2;
}

@ Preprocessor directives provide precidates to interrogate a
|pair|'s tag and find out what it is.

Although not all of these \\{cXr} macros are used they are all
defined here for completeness (and it's easier than working out
which ones really are needed).

@d special_p(p)   ((p) < 0)
@d boolean_p(p)   ((p) == FALSE or (p) == TRUE)
@d eof_p(p)       ((p) == END_OF_FILE)
@d false_p(p)     ((p) == FALSE)
@d null_p(p)      ((p) == NIL)
@d true_p(p)      ((p) == TRUE)
@d void_p(p)      ((p) == VOID)
@d undefined_p(p) ((p) == UNDEFINED)
@#
@d mark_p(p)  (!special_p(p) && (TAG[(p)] & TAG_MARK))
@d state_p(p) (!special_p(p) && (TAG[(p)] & TAG_STATE))
@d acar_p(p)  (!special_p(p) && (TAG[(p)] & TAG_ACARP))
@d acdr_p(p)  (!special_p(p) && (TAG[(p)] & TAG_ACDRP))
@d mark_clear(p)  (TAG[(p)] &= ~TAG_MARK)
@d mark_set(p)    (TAG[(p)] |=  TAG_MARK)
@d state_clear(p) (TAG[(p)] &= ~TAG_STATE)
@d state_set(p)   (TAG[(p)] |=  TAG_STATE)
@d format(p)      (TAG[(p)] & TAG_FORMAT)
@#
@d tag(p) (TAG[(p)])
@#
@d car(p)    (            CAR[(p)]   )
@d cdr(p)    (            CDR[(p)]   )
@d caar(p)   (        CAR[CAR[(p)]]  )
@d cadr(p)   (        CAR[CDR[(p)]]  )
@d cdar(p)   (        CDR[CAR[(p)]]  )
@d cddr(p)   (        CDR[CDR[(p)]]  )
@d caaar(p)  (    CAR[CAR[CAR[(p)]]] )
@d caadr(p)  (    CAR[CAR[CDR[(p)]]] )
@d cadar(p)  (    CAR[CDR[CAR[(p)]]] )
@d caddr(p)  (    CAR[CDR[CDR[(p)]]] )
@d cdaar(p)  (    CDR[CAR[CAR[(p)]]] )
@d cdadr(p)  (    CDR[CAR[CDR[(p)]]] )
@d cddar(p)  (    CDR[CDR[CAR[(p)]]] )
@d cdddr(p)  (    CDR[CDR[CDR[(p)]]] )
@d caaaar(p) (CAR[CAR[CAR[CAR[(p)]]]])
@d caaadr(p) (CAR[CAR[CAR[CDR[(p)]]]])
@d caadar(p) (CAR[CAR[CDR[CAR[(p)]]]])
@d caaddr(p) (CAR[CAR[CDR[CDR[(p)]]]])
@d cadaar(p) (CAR[CDR[CAR[CAR[(p)]]]])
@d cadadr(p) (CAR[CDR[CAR[CDR[(p)]]]])
@d caddar(p) (CAR[CDR[CDR[CAR[(p)]]]])
@d cadddr(p) (CAR[CDR[CDR[CDR[(p)]]]])
@d cdaaar(p) (CDR[CAR[CAR[CAR[(p)]]]])
@d cdaadr(p) (CDR[CAR[CAR[CDR[(p)]]]])
@d cdadar(p) (CDR[CAR[CDR[CAR[(p)]]]])
@d cdaddr(p) (CDR[CAR[CDR[CDR[(p)]]]])
@d cddaar(p) (CDR[CDR[CAR[CAR[(p)]]]])
@d cddadr(p) (CDR[CDR[CAR[CDR[(p)]]]])
@d cdddar(p) (CDR[CDR[CDR[CAR[(p)]]]])
@d cddddr(p) (CDR[CDR[CDR[CDR[(p)]]]])
@c

@ Both |atom|s and cons cells are stored in |pair|s. The lower 6
bits of the tag define the format of data stored in that |pair|.
The |atom|s are grouped into three types depending on whether both
|cell|s point to another |pair|, whether only the |cdr| does, or
whether both |cell|s are opaque. From this we obtain the core data
types.

@d FORMAT_CONS        (TAG_ACARP | TAG_ACDRP | 0x00)
@d FORMAT_APPLICATIVE (TAG_ACARP | TAG_ACDRP | 0x01)
@d FORMAT_OPERATIVE   (TAG_ACARP | TAG_ACDRP | 0x02)
@d FORMAT_SYNTAX      (TAG_ACARP | TAG_ACDRP | 0x03)
@d FORMAT_ENVIRONMENT (TAG_ACARP | TAG_ACDRP | 0x04)
@d FORMAT_EXCEPTION   (TAG_ACARP | TAG_ACDRP | 0x05)
@d FORMAT_INTEGER     (TAG_ACDRP | 0x00) /* value   : next/|NIL| */
@d FORMAT_SYMBOL      (TAG_NONE | 0x00) /* length   : offset */
@d FORMAT_VECTOR      (TAG_NONE | 0x01) /* gc-index : offset */
@d FORMAT_COMPILER    (TAG_NONE | 0x02) /* offset   : |NIL| */
@#
@d atom_p(p)          (!special_p(p) && ((tag(p) & TAG_FORMAT) != (TAG_ACARP | TAG_ACDRP)))
@d pair_p(p)          (!special_p(p) && ((tag(p) & TAG_FORMAT) == (TAG_ACARP | TAG_ACDRP)))
@d applicative_p(p)   (!special_p(p) && ((tag(p) & TAG_FORMAT) == FORMAT_APPLICATIVE))
@d compiler_p(p)      (!special_p(p) && ((tag(p) & TAG_FORMAT) == FORMAT_COMPILER))
@d environment_p(p)   (!special_p(p) && ((tag(p) & TAG_FORMAT) == FORMAT_ENVIRONMENT))
@d integer_p(p)       (!special_p(p) && ((tag(p) & TAG_FORMAT) == FORMAT_INTEGER))
@d operative_p(p)     (!special_p(p) && ((tag(p) & TAG_FORMAT) == FORMAT_OPERATIVE))
@d symbol_p(p)        (!special_p(p) && ((tag(p) & TAG_FORMAT) == FORMAT_SYMBOL))
@d syntax_p(p)        (!special_p(p) && ((tag(p) & TAG_FORMAT) == FORMAT_SYNTAX))
@d vector_p(p)        (!special_p(p) && ((tag(p) & TAG_FORMAT) == FORMAT_VECTOR))
@c

@ Allocating a new |pair| may require garbage collection to be
performed. If the data being put into either half of the new |pair|
is itself a |pair| it may be discarded by the collector. To avoid
this happening the data are saved into preallocated temporary storage
while a new |pair| is being located.

@<Global var...@>=
cell Tmp_CAR = NIL;
cell Tmp_CDR = NIL;

@ @<Protect...@>=
&Tmp_CAR, &Tmp_CDR,

@ @<Function dec...@>=
cell atom (cell, cell, char);

@ @d cons(a, d) atom((a), (d), FORMAT_CONS)
@c
cell
atom (cell ncar,
      cell ncdr,
      char ntag)
{
        cell r;
        if (null_p(Cells_Free)) {
                if (ntag & TAG_ACARP)
                        Tmp_CAR = ncar;
                if (ntag & TAG_ACDRP)
                        Tmp_CDR = ncdr;
                if (gc() <= (Cells_Poolsize / 2)) {
                        new_cells_segment();
                        gc();
                }
                Tmp_CAR = Tmp_CDR = NIL;
        }
        r = Cells_Free;
        Cells_Free = cdr(Cells_Free);
        car(r) = ncar;
        cdr(r) = ncdr;
        tag(r) = ntag;
        return r;
}

@* Vectors. A |vector| stores a contiguous sequence of |cell|s,
each referring to a |pair| on the heap. Unlike |pair|s |vector|s
are compacted during garbage collection to avoid fragmentation.

Storage is largely the same as |cell|s except for how the free
pointer is maintained: an index into the next unused |cell| in |VECTOR|.

@<Global var...@>=
cell *VECTOR = NULL;
int   Vectors_Free = 0;
int   Vectors_Poolsize = 0;
int   Vectors_Segment = HEAP_SEGMENT;

@ @c
void
new_vector_segment(void)
{
        cell *new_vector;
        new_vector = reallocarray(VECTOR, Vectors_Poolsize + Vectors_Segment,
                sizeof (cell));
        if (new_vector == NULL)
                error(ERR_OOM, NIL);
        bzero(((char *) new_vector) + Vectors_Poolsize,
                Vectors_Segment * sizeof (cell));
        VECTOR = new_vector;
        Vectors_Poolsize += Vectors_Segment;
        Vectors_Segment = Vectors_Poolsize / 2;
        return;
}

@ When a |pair| holds a |vector| its tag is |FORMAT_VECTOR|,
the |car| is used by the garbage collecter and the |cdr| is an
index into |VECTOR|.

Each |vector| contains 2 additional pieces of metadata (which are
{\bf above} the index), the length of the |vector| and a reference
back to the |pair| holding the |vector|.

A |vector| of length 0 is treated as a global constant akin to
|NIL| but it must be stored in a variable and created during
initialisation.

@d VECTOR_CELL 0
@d VECTOR_SIZE 1
@d VECTOR_HEAD 2
@#
@d vector_realsize(s) ((s) + VECTOR_HEAD)
@#
@d vector_cell(v)   (VECTOR[vector_offset(v) - (VECTOR_HEAD - VECTOR_CELL)])
@d vector_index(v)  (car(v))
@d vector_length(v) (VECTOR[vector_offset(v) - (VECTOR_HEAD - VECTOR_SIZE)])
@d vector_offset(v) (cdr(v))
@d vector_ref(v,o)  (VECTOR[vector_offset(v) + (o)])
@<Global var...@>=
cell Zero_Vector = NIL;

@ @<Global init...@>=
Zero_Vector = vector_new_imp(0, 0, 0);

@ Separate storage means separate garbage collection and a different
allocator. |vector_new_imp|, again, is broadly similar to |atom|
without the need for preallocated storage.

@ @c
cell
vector_new_imp (int  size,
                int  fill_p,
                cell fill)
{
        int wsize, off, i;
        cell r;
        wsize = vector_realsize(size);
        if (Vectors_Free + wsize >= Vectors_Poolsize) {
                gc_vectors();
                while (Vectors_Free + wsize >= (Vectors_Poolsize - (Vectors_Poolsize / 2))) {
                        new_vector_segment();
                        gc_vectors(); /* Is this really necessary? */
                }
        }
        r = atom(NIL, NIL, FORMAT_VECTOR);
        off = Vectors_Free;
        Vectors_Free += wsize;
        vector_offset(r) = off + VECTOR_HEAD; /* must be first */
        vector_length(r) = size;
        vector_cell(r) = r;
        vector_index(r) = 0;
        if (fill_p)
                for (i = VECTOR_HEAD; i <= size + (VECTOR_HEAD - 1); i++)@/
                        vector_ref(r, i) = fill;
        return r;
}

@ @c
cell
vector_new (int  size,
            cell fill)
{
        if (size == 0)
                return Zero_Vector;
        return vector_new_imp(size, 1, fill);
}

@ |vector_new_list| turns a |list| of |pair|s into a |vector|.

@c
cell
vector_new_list (cell list,
                 int  len)
{
        cell r;
        int i;
        r = vector_new(len, 0);
        for (i = 0; i < len; i++) {
                vector_ref(r, i) = car(list);
                list = cdr(list);
        }
        return r;
}

@ Although a little early in the narrative |vector_sub| is defined
here because it's the only other function substantially dealing
with |vector| data.

@c
cell
vector_sub (cell src,
            int srcfrom,
            int srcto,
            int dstfrom,
            int dstto,
            cell fill)
{
        cell dst;
        int copy, i;
        copy = srcto - srcfrom;
        if (dstto < 0)
                dstto = dstfrom + copy;
        dst = vector_new_imp(dstto, 0, 0);
        for (i = 0; i < dstfrom; i++)
                vector_ref(dst, i) = fill;
        for (i = srcfrom; i < srcto; i++)@/
                vector_ref(dst, (dstfrom - srcfrom) + i) = vector_ref(src, i);
        for (i = dstfrom + copy; i < dstto; i++)
                vector_ref(dst, i) = fill;
        return dst;
}

@* Garbage Collection. The garbage collector is a straightforward
mark and sweep collector. |mark| is called for every entry in |ROOTS|
to recursively set the mark bit on every reachable |pair|, then
the whole pool is scanned and any |pair|s which aren't marked are
added to the free list.

% TODO: break this algorithm down into explained pieces.

@<Function dec...@>=
int gc (void);
int gc_vectors (void);

@ |ROOTS| is a |NULL|-terminated \CEE/ array of objects to protect
from collection. I can't think of any better way of declaring it
but hard-coding it right here.

@c
cell *ROOTS[] = { @<Protected Globals@>@t, @> NULL };

@ @c
void
mark(cell next)
{
        cell parent, prev;
        int i;
        parent = prev = NIL;
        while (1) {
                if (!(special_p(next) || mark_p(next))) {
                        if (vector_p(next)) {         /* S0 $\to$ S.1 */
                                mark_set(next);
                                vector_cell(next) = next;
                                if (vector_length(next) > 0) {
                                        state_set(next);
                                        vector_index(next) = 0;
                                        prev = vector_ref(next, 0);
                                        vector_ref(next, 0) = parent;
                                        parent = next;
                                        next = prev;
                                }
                        } else if (!acar_p(next)
                                   && acdr_p(next)) { /* S0 $\to$ S2 */
                                prev = cdr(next);
                                cdr(next) = parent;
                                parent = next;
                                next = prev;
                                mark_set(parent);
                        } else if (acar_p(next)) {    /* S0 $\to$ S1 */
                                prev = car(next);
                                car(next) = parent;
                                mark_set(next);
                                parent = next;
                                next = prev;
                                state_set(parent);
                        } else {                      /* S0 $\to$ S1 */
                                mark_set(next);
                        }
                } else {
                        if (null_p(parent))
                                break;
                        if (vector_p(parent)) {       /* S.1 $\to$ S.1/done */
                                i = vector_index(parent);
                                if ((i + 1) < vector_length(parent)) {
                                        prev = vector_ref(parent, i + 1);
                                        vector_ref(parent, i + 1) = vector_ref(parent, i);
                                        vector_ref(parent, i) = next;
                                        next = prev;
                                        vector_index(parent) = i + 1;
                                } else {              /* S.1 $\to$ done */
                                        state_clear(parent);
                                        prev = parent;
                                        parent = vector_ref(prev, i);
                                        vector_ref(prev, i) = next;
                                        next = prev;
                                }
                        } else if (state_p(parent)) { /* S1 $\to$ S2 */
                                prev = cdr(parent);
                                cdr(parent) = car(parent);
                                car(parent) = next;
                                state_clear(parent);
                                next = prev;
                        } else if (acdr_p(parent)) {  /* S2 $\to$ done */
                                prev = parent;
                                parent = cdr(prev);
                                cdr(prev) = next;
                                next = prev;
                        } else {
                                error(ERR_UNIMPLEMENTED, NIL);
                        }
                }
        }
}

int
gc (void)
{
        int count, sk, i;
        if (!null_p(RTS)) {
                sk = vector_length(RTS);
                vector_length(RTS) = RTSp + 1;
        }
        for (i = 0; ROOTS[i]; i++)
                mark(*ROOTS[i]);
        for (i = SCHAR_MIN; i <= SCHAR_MAX; i++) {
                mark(Small_Int[(unsigned char) i]);
        }
        if (!null_p(RTS))
                vector_length(RTS) = sk;
        Cells_Free = NIL;
        count = 0;
        for (i = 0; i < Cells_Poolsize; i++) {
                if (!mark_p(i)) {
                        cdr(i) = Cells_Free;
                        Cells_Free = i;
                        count++;
                } else {
                        mark_clear(i);
                }
        }
        return count;
}

@ |vector| garbage collection works by using the |pair|s garbage
collector to scan |ROOTS| and determine which vectors are really
in use then removes any which aren't from |VECTORS|, decrementing
|Vectors_Free| if it can.

@c
int
gc_vectors (void)
{
        int to, from, d, i, r;
        @<Unmark all vectors@>@/
        gc();
        from = to = 0;
        while (from < Vectors_Free) {
                d = vector_realsize(VECTOR[from + VECTOR_SIZE]);
                if (!null_p(VECTOR[from + VECTOR_CELL])) {
                        if (to != from) {
                                for (i = 0; i < d; i++)
                                        VECTOR[to + i] = VECTOR[from + i];
                                vector_offset(VECTOR[to + VECTOR_CELL])
                                        = to + VECTOR_HEAD;
                        }
                        to += d;
                }
                from += d;
        }
        r = Vectors_Free - to;
        Vectors_Free = to;
        return r;
}

@ To ``unmark'' a |vector|, all the links in |VECTOR| back to the
cell which refers to it (|vector_cell|) are set to |NIL|. |gc| will
re-set the link in any vectors that it can reach.

@<Unmark all vectors@>=
i = 0;
while (i < Vectors_Free) {
        VECTOR[i + VECTOR_CELL] = NIL;
        i += vector_realsize(VECTOR[i + VECTOR_SIZE]);
}

@* Objects. Although not \.{object}s per se, the first \.{object}s
which will be defined are three stacks. We could define the run-time
stack later because it's not used until the virtual machine is
implemented but the implementations mirror each other and the
internal VM stack is required before real objects can be defined.
Also the runtime stack uses the VM stack in its implementation.

The compiler stack is included here because it's identical to the
VM stack.

The VM stack is a pointer to the head of a |list|. This means
that accessing the top few elements of the stack---especially pushing
and popping a single \.{object}---is effectively free but accessing an
arbitrary part of the stack requires an expensive walk over each
item in turn.

On the other hand the run-time stack is stored in a |vector| with
a pointer |RTSp| to the current head of the stack, which is -1 if
the stack is empty.

This has the obvious disadvantage that its storage space is finite
and occasionally the whole stack will need to be copied into a new,
larger |vector| (and conversely it may waste space or require
occasional trimming). On the other hand random access to any part
of the stack has the same (negligable) cost.

When it's not ambiguous ``stack'' in this document refers to the
run-time stack; the VM stack is an implementation detail. In fact
the run-time stack is also an implementation detail but the VM stack
is an implementation detail of that implementation detail; do you
like recursion yet?.

The main interface to each stack is its |push|/|pop|/|ref|/|clear|
functions. There are some additional handlers for the run-time stack.

@d ERR_UNDERFLOW "underflow"
@d ERR_OVERFLOW "overflow"
@d CHECK_UNDERFLOW(s) if (null_p(s)) error(ERR_UNDERFLOW, VOID)
@d RTS_UNDERFLOW(p) if ((p) < -1) error(ERR_UNDERFLOW, RTS)
@d RTS_OVERFLOW(p) if ((p) > RTSp) error(ERR_OVERFLOW, RTS)
@<Global var...@>=
cell CTS = NIL;
cell RTS = NIL;
cell VMS = NIL;
int  RTS_Size = 0;
int  RTSp = -1;

@ @<Protected...@>=
&CTS, &RTS, &VMS,

@ @<Function dec...@>=
cell vms_pop (void);
void vms_push (cell);


@ The VM and compiler stacks |VMS| and |CTS| are built on |list|s.

@d vms_clear() ((void) vms_pop())
@c
cell
vms_pop (void)
{
        cell r;
        CHECK_UNDERFLOW(VMS);
        r = car(VMS);
        VMS = cdr(VMS);
        return r;
}

void
vms_push (cell item)
{@+
        VMS = cons(item, VMS);@+
}

cell
vms_ref (void)
{
        CHECK_UNDERFLOW(VMS);
        return car(VMS);
}

void
vms_set (cell item)
{
        CHECK_UNDERFLOW(VMS);
        car(VMS) = item;
}

@ |CTS| is treated identically to |VMS|. Using the \CEE/ preprocessor
for this would be unnecessarily inelegant so instead here is a
delicious bowl of pasta.

@d cts_clear() ((void) cts_pop())
@d cts_reset() CTS = NIL
@c
cell
cts_pop ()
{
        cell r;
        CHECK_UNDERFLOW(CTS);
        r = car(CTS);
        CTS = cdr(CTS);
        return r;
}

void
cts_push (cell item)
{@+
        CTS = cons(item, CTS);@+
}

cell
cts_ref (void)
{
        CHECK_UNDERFLOW(CTS);
        return car(CTS);
}

void
cts_set (cell item)
{
        CHECK_UNDERFLOW(CTS);
        car(CTS) = item;
}

@ Being built on a |vector| the run-time stack needs to increase
its size when it's full. Functions can call |rts_prepare| to ensure
that the stack is big enough for their needs.

@d RTS_SEGMENT 0x1000
@c
void
rts_prepare (int need)
{
        int b, s;
        if (RTSp + need >= RTS_Size) {
                b = RTS_SEGMENT * ((need + RTS_SEGMENT) / RTS_SEGMENT);
                s = RTS_Size + b;
                RTS = vector_sub(RTS, 0, RTS_Size, 0, s, UNDEFINED);
                RTS_Size = s;
        }
}

@ Otherwise, the run-time stack has the same interface but a different
implementation.

@d rts_clear(c) ((void) rts_pop(c))
@c
cell
rts_pop (int count)
{
        RTS_UNDERFLOW(RTSp - count);
        RTSp -= count;
        return vector_ref(RTS, RTSp + 1);
}

void
rts_push (cell o)
{
        vms_push(o);
        rts_prepare(1);
        vector_ref(RTS, ++RTSp) = vms_pop();
}

cell
rts_ref (int d)
{
        RTS_UNDERFLOW(RTSp - d);
        RTS_OVERFLOW(RTSp - d);
        return vector_ref(RTS, RTSp - d);
}

cell
rts_ref_abs (int d)
{
        RTS_UNDERFLOW(d);
        RTS_OVERFLOW(d);
        return vector_ref(RTS, d);
}

void
rts_set (int  d,
         cell v)
{
        RTS_UNDERFLOW(RTSp - d);
        RTS_OVERFLOW(RTSp - d);
        vector_ref(RTS, RTSp - d) = v;
}

void
rts_set_abs (int  d,
             cell v)
{
        RTS_UNDERFLOW(d);
        RTS_OVERFLOW(d);
        vector_ref(RTS, d) = v;
}

@*1 Symbols. With the basics in place, the first thing to define
is |symbol|s; they're not needed yet but everything becomes easier
with them extant and they depend on nothing but themselves since
they are themselves.

|symbol|s are never garbage collected. This was not a conscious
decision it just doesn't seem like it matters. Instead, every
|symbol| once created is immediately added to the |Symbol_Table|
|list|. When a reference to a |symbol| is requested, the object
in this |list| is returned.

Eventually this should implement a hash table but I'm not making
one of those this morning.

Owing to the nasty c-to-perl-to-c route that I've taken, combined
with plans for vector/byte storage, the storage backing symbols is
going to be hairy without explanation (for now it's a mini duplicate
of vector storage).

@d sym(s) symbol((s), 1)
@d symbol_length car
@d symbol_offset cdr
@d symbol_store(s) (SYMBOL + symbol_offset(s))
@<Global var...@>=
cell  Symbol_Table = NIL;
char *SYMBOL = NULL;
int   Symbol_Free = 0;
int   Symbol_Poolsize = 0;

@ @<Protected...@>=
&Symbol_Table,

@ @<Function dec...@>=
cell symbol (char *, int);

@ @c
void
symbol_expand (void)
{
        char *new;
        new = realloc(SYMBOL, Symbol_Poolsize + HEAP_SEGMENT);
        if (new == NULL)
                error(ERR_OOM, NIL);
        Symbol_Poolsize += HEAP_SEGMENT;
        SYMBOL = new;
}

@ A |symbol| can ``steal'' storage from |SYMBOL| which results
in an \.{object} which can be mostly treated like a normal |symbol|,
used to compare a potentially new |symbol| with those currently
stored in |Symbol_Table|. This is the closest that |symbol|s get
to being garbage collected.

@c
cell
symbol_steal (char *cstr)
{
        cell r;
        int len;
        len = strlen(cstr);
        while (Symbol_Free + len > Symbol_Poolsize)
                symbol_expand();
        r = atom(len, Symbol_Free, FORMAT_SYMBOL);
        memcpy(SYMBOL + Symbol_Free, cstr, len);
        /* |Symbol_Free| is {\bf not} incremented here */
        return r;
}

@ Temporary |symbol|s compare byte-by-byte with existing |symbol|s.
This is not efficient at all.

@c
boolean
symbol_same_p (cell maybe,
               cell match)
{
        char *pmaybe, *pmatch;
        int i, len;
        len = symbol_length(match);
        if (symbol_length(maybe) != len)
                return bfalse;
        pmaybe = symbol_store(maybe);
        pmatch = symbol_store(match);
        if (maybe == match) /* This shouldn't happen */
                return btrue;
        for (i = 0; i < len; i++) {
                if (pmaybe[i] != pmatch[i])
                        return bfalse;
        }
        return btrue;
}

@ @c
void
symbol_reify (cell s)
{
        Symbol_Free += symbol_length(s);
        Symbol_Table = cons(s, Symbol_Table);
}

cell
symbol (char *cstr,
        int   permanent_p)
{
        cell st, s;
        s = symbol_steal(cstr);
        st = Symbol_Table;
        while (!null_p(st)) {
                if (symbol_same_p(s, car(st)))
                        return car(st);
                st = cdr(st);
        }
        if (permanent_p)
                symbol_reify(s);
        return s;
}

@*1 Numbers. The only numbers supported by this early implementation
of \LL/ are signed integers that fit in a single |cell| (ie.
32-bit integers).

The 256 numbers closest to 0 (ie. from |-0x80| to |+0x7f|) are
preallocated during initialisation. If you live in a parallel
universe where the char type isn't 8 bits then adjust those numbers
accordingly.

% TODO: need to verify that char casting works as I expect.

@d fixint_p(p) (integer_p(p) && null_p(int_next(p)))
@d int_value(p) ((int) (car(p)))
@d int_next cdr
@<Global var...@>=
cell Small_Int[UCHAR_MAX + 1];

@ Even though the |Small_Int| objects are about to be created, in
order to create objects garbage collection will happen and assume
that |Small_Int| has already been initialised and attempt to protect
data which don't exist from collection. This is a silly solution
but I'm leaving it alone until I have a better memory model.

@<Pre-init...@>=
for (i = 0; i < 256; i++)
        Small_Int[i] = NIL;

@ @<Global init...@>=
for (i = SCHAR_MIN; i <= SCHAR_MAX; i++)@/
        Small_Int[(unsigned char) i] = int_new_imp(i, NIL);

@ As with |vector|s, |int_new| checks whether it should return
an |object| from |Small_Int| or build a new one.

@c
cell
int_new_imp (int  value,
             cell next)
{
        if (!null_p(next))
                error(ERR_UNIMPLEMENTED, NIL);
        return atom((cell) value, next, FORMAT_INTEGER);
}

cell
int_new (int value)
{
        if (value >= SCHAR_MIN && value <= SCHAR_MAX)@/
                return Small_Int[(unsigned char) value];
        return int_new_imp(value, NIL);
}

@*1 Pairs \AM\ Lists. Of course |pair|s---and so by definition
|list|s---have already been implemented but so far only enough
to implement core features. Here we define handlers for operations
specifically on |list| \.{object}s.

First to count it's length a |list| is simply walked from head
to tail. It is not considered an error if the |list| is improper
(or not a |list| at all). To indicate this case the returned
length is negated.

@c
int
list_length (cell l)
{
        int c = 0;
        if (null_p(l))
                return 0;
        for (; pair_p(l); l = cdr(l))
                c++;
        if (!null_p(l))
                c = -(c + 1);
        return c;
}

@ A |list| is either NIL or a pair with one restriction, that its
|cdr| must itself be a |list|. The size of the |list| is also counted
to avoid walking it twice but nothing uses that (yet?).

@c
predicate
list_p (cell o,
        predicate improper_p,
        cell *sum)
{
        int c = 0;
        if (null_p(o)) {
                if (sum != NULL)
                        *sum = int_new(0);
                return TRUE;
        }
        while (pair_p(o)) {
                o = cdr(o);
                c++;
        }
        if (sum != NULL)
                *sum = int_new(c);
        if (null_p(o))
                return TRUE;
        if (sum != NULL)
                *sum = int_new(-(c + 1));
        return improper_p;
}

@ A proper |list| can be reversed simply into a new |list|.

@d ERR_IMPROPER_LIST "improper-list"
@c
cell
list_reverse (cell  l,
              cell *improper,
              cell *sum)
{
        cell saved, r;
        int c;
        saved = l;
        c = 0;
        vms_push(NIL);
        while (!null_p(l)) {
                if (!pair_p(l)) {
                        r = vms_pop();
                        if (improper != NULL) {
                                *improper = l;
                                if (sum != NULL)
                                        *sum = c;
                                return r;
                        } else
                                error(ERR_IMPROPER_LIST, saved);
                }
                vms_set(cons(car(l), vms_ref()));
                l = cdr(l);
                c++;
        }
        if (sum != NULL)
                *sum = int_new(c);
        return vms_pop();
}

@ Reversing a |list| in-place means maintaining a link to the
previous |pair| (or |NIL|) and replacing each |pair|'s |cdr|.
The new head |pair| is returned, or |FALSE| if the |list| turned
out to be improper.

@c
cell
list_reverse_m (cell    l,
                boolean error_p)
{
        cell m, t, saved;
        saved = l;
        m = NIL;
        while (!null_p(l)) {
                if (!pair_p(l)) {
                        if (!error_p) /* TODO: repair? */
                                return FALSE;
                        error(ERR_IMPROPER_LIST, saved);
                }
                t = cdr(l);
                cdr(l) = m;
                m = l;
                l = t;
        }
        return m;
}

@*1 Environments. In order to associate a value with a |symbol|
(a variable) they are paired together in an |environment|.

Like an onion or an ogre\footnote{$^1$}{Or a cake.}, an |environment|
has layers. The top layer is both the current layer and the current
|environment|. The bottom layer is the root |environment| |Root|.

An |environment| is stored in an |atom| with the |car| pointing
to the previous layer (or |NIL| in the root |environment|).

The |cdr| is a |list| of association |pair|s representing the
variables in that layer. An association |pair| is a proper |list|
with two items: an identifier, in this case a |symbol|, and a
value.

|environment|-handling functions and macros are generally named
``env''.

@d ERR_BOUND "already-bound"
@d ERR_UNBOUND "unbound"
@#
@d env_empty() atom(NIL, NIL, FORMAT_ENVIRONMENT)
@d env_extend(e) atom((e), NIL, FORMAT_ENVIRONMENT)
@d env_layer cdr
@d env_parent car
@d env_root_p(e) (environment_p(e) && null_p(car(e))
@ Searching through an |environment| starts at its top layer and
walks along each |pair|. If it encounters a |pair| who's
|symbol| matches, the value is returned. If not then the search
repeats layer by layer until the |environment| is exhausted and
|UNDEFINED| is returned.

|env_search| does not raise an error if a |symbol| isn't found.
This means that |UNDEFINED| is the only value which cannot be stored
in a variable as there is no way to distinguish its return from
this function.

@c
cell
env_search (cell haystack,
            cell needle)
{
        cell n;
        for (; !null_p(haystack); haystack = env_parent(haystack))
                for (n = env_layer(haystack); !null_p(n); n = cdr(n))
                        if (caar(n) == needle)
                                return cadar(n);
        return UNDEFINED;
}

@ To set a variable's value the |environment|'s top layer is first
searched to see if the |symbol| is already bound. An |error| is
raised if the symbol is bound (when running on behalf of {\it
define!\/}) or not bound (when running on behalf of {\it set!\/}).

@c
void
env_set (cell e,
         cell name,
         cell value,
         boolean new_p)
{
        cell ass, t;
        ass = cons(name, cons(value, NIL));
        if (new_p) { @<Mutate if unbound@> }
        else { @<Mutate if bound@> }
}

@ Updating an already-bound variable means removing the existing
binding from the |environment| and inserting the new binding. During
the walk over the layer |t| is one pair ahead of the pair being
considered so that when |name| is found |t|'s |cdr| can be changed,
snipping the old binding out, so the first pair is checked specially.
@<Mutate if bound@>=
if (null_p(env_layer(e)))
        error(ERR_UNBOUND, name);
if (caar(env_layer(e)) == name) {
        env_layer(e) = cons(ass, cdr(env_layer(e)));
        return;
}
for (t = env_layer(e); !null_p(cdr(t)); t = cdr(t)) {
        if (caadr(t) == name) {
                cdr(t) = cddr(t);
                env_layer(e) = cons(ass, env_layer(e));
                return;
        }
}
error(ERR_UNBOUND, name);

@ The case is simpler if the |name| must {\bf not} be bound already
as the new binding can be prepended to the layer after searching with
no need for special cases.
@<Mutate if unbound@>=
for (t = env_layer(e); !null_p(t); t = cdr(t))
        if (caar(t) == name)
                error(ERR_BOUND, name);
env_layer(e) = cons(ass, env_layer(e));

@ Values are passed to functions on the stack. |env_lift_stack|
moves these values from the stack into an |environment|.

@c

cell
env_lift_stack (cell e,
                int nargs,
                cell formals)
{
        cell p, name, value, ass;
        vms_push(env_extend(e));
        p = NIL; /* prepare a new layer */
        vms_push(p);
        while (nargs--) {
                if (pair_p(formals)) {
                        name = car(formals);
                        formals = cdr(formals);
                } else {@+
                        name = formals;@+
                }
                value = rts_pop(1);
                if (!null_p(name)) {
                        ass = cons(name, cons(value, NIL));
                        vms_set((p = cons(ass, p)));
                }
        }
        vms_pop();
        cdr(vms_ref()) = p; /* place the new layer in the extended
                               environment */
        return vms_pop();
}

@*1 Closures \AM\ Compilers. Finally we have data structures to
save run-time state: |closure|s. The way the compiler and virtual
machine work to get |closure| objects built is described below---here
is only a description of their backing stores.

\LL/ has two types of |closure|, |applicative| and |operative|.
They store the same data in identical containers; the difference
is in how they're used.

The data required to define a |closure| are a program \AM\ the
|environment| to run it in. A |closure| in \LL/ also contains
the formals given in the |lambda| or |vov| expression that
was used to define it.

Program code in \LL/ is stored as compiled bytecode in a
|vector| with an instruction pointer indicating the entry point
(0 is not implied). The |closure|s then look like this:

$$(APPLICATIVE\ \langle formals \rangle\ \langle environment \rangle\
\langle code \rangle\ \langle pointer \rangle)$$
$$(OPERATIVE\ \langle formals \rangle\ \langle environment \rangle\
\langle code \rangle\ \langle pointer \rangle)$$

However the |environment|, code and pointer are never referred to
directly until the closure is unpicked by |OP_APPLY|/|OP_APPLY_TAIL|.
Instead the objects effectively look like this:

$$(A\vert O\ \langle formals\rangle\ .\ \langle opaque\_closure\rangle)$$

@d applicative_closure cdr
@d applicative_formals car
@d applicative_new(f,e,p,i)
        closure_new_imp(FORMAT_APPLICATIVE, (f), (e), (p), (i))
@d operative_closure cdr
@d operative_formals car
@d operative_new(f,e,p,i)
        closure_new_imp(FORMAT_OPERATIVE, (f), (e), (p), (i))
@c
cell
closure_new_imp (char ntag,
                 cell formals,
                 cell env,
                 cell prog,
                 cell ip)
{
        cell r;
        r = cons(int_new(ip), NIL);
        r = cons(prog, r);
        r = cons(env, r);
        return atom(formals, r, ntag);
}

@ Other than closures, and required in order to make them, the
evaluator uses ``|compiler|'' \.{object}s that compile \LL/ source
code to VM bytecode. Each |compiler| is described in the structure
|primitive|, containing the native function pointer to it.

% Apparently cannot parse function pointer typedefs
%s native void
@d compiler_cname(c) COMPILER[car(c)].name
@d compiler_fn(c) COMPILER[car(c)].fn
@<Type def...@>=
typedef void @[@] (*native) (cell, cell, boolean);
typedef struct {
        char *name;
        native fn;
} primitive;

@ The contents of |COMPILER| are populated by the \CEE/ compiler
of this source. During initialisation |Root| then becomes the root
environment filled with an association |pair| for each one.

@<Global var...@>=
@q This is littered with @@+ so it appears on one line. The final @@+@>
@q should be @@& but that causes a literal @@& to appear in the output@>
primitive COMPILER[] =@+ {@+
        @<List of opcode primitives@>@t,@>@+
        { NULL, NULL }@+ }@+;

@** Virtual Machine. This implementation of \LL/ compiles user
source code to an internal bytecode representation which is then
executed sequentially by a virtual machine (VM).

Additionally to the myriad stacks already mentioned, the VM maintains
(global!) state primarily in 6 registers. Two of these are simple
flags (booleans) which indicate whether interpretation should
continue.

\point 1. |Running| is a flag raised (1) when the VM begins and
lowered by user code to indicate that it should halt cleanly. This
flag is checked on the beginning of each iteration of the VM's main
loop.

\point 2. |Interrupt| is normally lowered (0) and is raised in
response to external events such as a unix signal. Long-running
operations---especially those which could potentially run
unbounded---check frequently for the state of this flag and abort
and return immediately when it's raised.

The other four registers represent the computation.

\point 3. |Acc| is the accumulator. Opcodes generally read and/or
write this register to do their work. This is where the final result
of computation will be found.

\point 4. |Env| holds the current |environment|. Changing this
is the key to implementing |closure|s.

\point 5. |Prog| is the compiled bytecode of the currently running
computation, a |vector| of VM opcodes with their in-line arguments.

\point 6. |Ip| is the instruction pointer. This is an |int|, not a
|cell| and must be boxed to be used outside of the VM.

|Root| and |Prog_Main| are also defined here which hold, respectively,
the root |environment| and the virtual machine's starting program.

@<Global var...@>=
boolean Interrupt = 0;
boolean Running   = 0;
cell    Acc       = NIL;
cell    Env       = NIL;
cell    Prog      = NIL;
cell    Prog_Main = NIL;
cell    Root      = NIL;
int     Ip        = 0;

@ @<Protected...@>=
&Acc, &Env, &Prog, &Prog_Main, &Root,

@ The virtual machine is initialised in two stages. First |vm_init|
performs global initialisation; this must be called exactly once.

@c
void
vm_init (void)
{
        cell t;
        int i;
        primitive *n;
        @<Pre-initialise |Small_Int|@>@;
        @<Global init...@>@;
        Prog_Main = compile_main();
        i = 0;
        Root = atom(NIL, NIL, FORMAT_ENVIRONMENT);
        for (n = COMPILER + i; n->fn != NULL; n = COMPILER + (++i)) {
                t = atom(i, NIL, FORMAT_COMPILER);
                t = cons(t, NIL);
                t = cons(sym(n->name), t);
                env_layer(Root) = cons(t, env_layer(Root));
        }
        Env = Root;
}

@ To complete initialisation |vm_clear| is called to ready the VM
for running another program. |vm_clear| does {\bf not} reset |Env|
(or |Acc|), which is what allows state to be maintained between
instructions in the \.{REPL}, for example.

@c
void
vm_clear (void)
{
        Prog = Prog_Main;
        Running = Interrupt = Ip = 0;
}

@* Frames. The VM enters a |closure|---aka. calls a function---by
appending a |frame| header to the stack. A |frame| consists of any
work-in-progress items on the stack followed by a fixed-size header.
A |frame|'s header captures the state of computation at the time
it's created which is what lets another subroutine run and then
return. The |frame| header contains 4 objects: $\ll$|Ip| |Prog|
|Env| |Fp|$\gg$.

|Fp| is a quasi-register which points into the stack to the current
|frame|'s header. It's saved when entering a |frame| and its value
set to that of the stack pointer |RTSp|. |RTSp| is restored to the
saved value when returning from a |frame|.

@d FRAME_HEAD 4
@d frame_ip(f)         rts_ref_abs((f) + 1)
@d frame_prog(f)       rts_ref_abs((f) + 2)
@d frame_env(f)        rts_ref_abs((f) + 3)
@d frame_fp(f)         rts_ref_abs((f) + 4)
@d frame_set_ip(f,v)   rts_set_abs((f) + 1, (v));
@d frame_set_prog(f,v) rts_set_abs((f) + 2, (v));
@d frame_set_env(f,v)  rts_set_abs((f) + 3, (v));
@d frame_set_fp(f,v)   rts_set_abs((f) + 4, (v));
@<Global var...@>=
int Fp = -1;

@ Creating a |frame| is pushing the header items onto the stack.
Entering it is changing the VM's registers that are now safe. This
is done in two stages for some reason.

@c
void
frame_push (int ipdelta)
{
        rts_push(int_new(Ip + ipdelta));
        rts_push(Prog);
        rts_push(Env);
        rts_push(int_new(Fp));
}

void
frame_enter (cell e,
             cell p,
             cell i)
{
        Env = e;
        Prog = p;
        Ip = i;
        Fp = RTSp - FRAME_HEAD;
}

@ Leaving a |frame| means restoring the registers that were saved
in it by |frame_push| and then returning |RTSp| and |Fp| to their
previous values; |Fp| from the header and |RTSp| as the current
|Fp| minus the |frame| header in case there were previously any
in-progress items on top of the stack.

@c
void
frame_leave (void)
{
        int prev;
        Ip = int_value(frame_ip(Fp));
        Prog = frame_prog(Fp);
        Env = frame_env(Fp);
        prev = int_value(frame_fp(Fp));
        rts_clear(FRAME_HEAD);
        Fp = prev;
}

@*1 Tail Recursion. {\bf TODO}

This is a straight copy of what I wrote in perl which hasn't been
used there. Looks about right. Might work.

@c
void
frame_consume (void)
{
        int src, dst, i;
        src = Fp;
        dst = int_value(frame_fp(src));
        /* Copy the parts of the old frame header that are needed */
        frame_set_prog(src, frame_prog(dst));
        frame_set_ip(src,   frame_ip(dst));
        frame_set_fp(src,   frame_fp(dst));
        /* Move the active frame over the top of the previous one */
        for (i = 1; i <= FRAME_HEAD; i++)@/
                rts_set_abs(dst + i, rts_ref_abs(src + i));
        rts_clear(src - dst);
        Fp -= src - dst;
}

@* Interpreter. The workhorse of the virtual machine is |interpret|.
After being reset with |vm_clear|, parsed (but not compiled) source
code is put into |Acc| and the VM can be started by calling
|interpret|.

@d ERR_INTERRUPTED "interrupted"
@c
void
interpret (void)
{
        int ins;
        cell tmp; /* {\bf not} saved in |ROOTS| */
        @<Set up run-time...@>@;
        Running = 1;
        RTSp = -1;
        while (Running && !Interrupt) {
                ins = int_value(vector_ref(Prog, Ip));
                switch (ins) {
                        @<Opcode implementations@>@;
                }
        }
        if (Interrupt)
                error(ERR_INTERRUPTED, NIL);
}

@** I/O. Before embarking on the meat of the interpreter a final
detour to describe routines to parse a string (or stream) of source
code into s-expressions, and because it's useful to see what's being
done routines to write them back again.

These routines use \CEE/'s |stdio| for now to get a simple
implementation finished.

@* Reader (or Parser). The s-expression reader is an ad-hoc LALR
parser; a single byte is read to determine which type of form to
parse. Bytes are then read one at a time to validate the syntax and
create the appropriate \.{object}.

The reading routines call into themselves recursively (for which
it cheats and relies on \CEE/'s stack). To prevent it running out
of control |Read_Level| records the recursion depth and |read_form|
aborts if it exceeds |READER_MAX_DEPTH|.

The compiler's rather than the VM's stack is used for temporary
storage so that error handling doesn't need to clean it up. This is
safe provided the reader and compiler are never used simultaneously.

The parser often needs the byte that was used to determine which
kind of form to parse (the one that was ``looked ahead'' at).
|Putback| is a small buffer to contain this byte. In fact this
buffer can hold {\it two} bytes to accomodate lisp's {\it
unquote-splicing} operator \qunsplice/.
% TODO: make unquote-splicing indexable

@d ERR_RECURSION "recursion"
@d ERR_UNEXPECTED "unexpected"
@d WARN_AMBIGUOUS_SYMBOL "ambiguous"
@#
@d READER_MAX_DEPTH 1024 /* gotta pick something */
@d READ_SPECIAL       -10
@d READ_DOT           -10 /* \qdot/ */
@d READ_CLOSE_BRACKET -11 /* \qclosebra/ */
@d READ_CLOSE_PAREN   -12 /* \qclosepar/ */
@#
@d SYNTAX_DOTTED   "dotted"           /* \qdot/ */
@d SYNTAX_QUOTE    "quote"            /* \qquote/ */
@d SYNTAX_QUASI    "quasiquote"       /* \qquasi/ */
@d SYNTAX_UNQUOTE  "unquote"          /* \qunquote/ */
@d SYNTAX_UNSPLICE "unquote-splicing" /* \qunsplice/ */
@<Global var...@>=
char Putback[2] = { '\0', '\0' };
int Read_Level = 0;
cell Sym_ERR_UNEXPECTED = NIL;
cell Sym_SYNTAX_DOTTED = NIL;
cell Sym_SYNTAX_QUASI = NIL;
cell Sym_SYNTAX_QUOTE = NIL;
cell Sym_SYNTAX_UNQUOTE = NIL;
cell Sym_SYNTAX_UNSPLICE = NIL;

@ @<Global init...@>=
Sym_ERR_UNEXPECTED = sym(ERR_UNEXPECTED);
Sym_SYNTAX_DOTTED = sym(SYNTAX_DOTTED);
Sym_SYNTAX_QUASI = sym(SYNTAX_QUASI);
Sym_SYNTAX_QUOTE = sym(SYNTAX_QUOTE);
Sym_SYNTAX_UNQUOTE = sym(SYNTAX_UNQUOTE);
Sym_SYNTAX_UNSPLICE = sym(SYNTAX_UNSPLICE);

@ @<Function dec...@>=
cell read_symbol (void);
cell read_form (void);
cell read_list (cell);
cell read_number (void);
cell read_symbol (void);
void unread_byte (char);

@ @c
int
read_byte (void)
{
        int r;
        if ((r = Putback[0]) != '\0') {
                Putback[0] = Putback[1];
                Putback[1] = '\0';
                return r;
        }
        return getchar();
}

void
unread_byte(char c)
{
        @q assert(Putback[1] == '\0')@>
        Putback[1] = Putback[0];
        Putback[0] = c;
}

@ Even this primitive parser should support primitive comments.

@c
int
useful_byte (void)
{
        int c;
        while (!Interrupt) {
                c = read_byte();
                switch (c) {
                case ' ':
                case '\n':
                case '\r':
                case '\t':
                        continue;
                case ';':
                        c = read_byte();
                        while (c != '\n' && !Interrupt) {
                                /* read up to but not beyond the next
                                   newline */
                                c = read_byte();
                                if (c == EOF)
                                        return c;
                        }
                        break; /* go around again */
                default:
                        return c; /* includes |EOF| (which |!=|
                                     |END_OF_FILE|) */
                }
        }
        return EOF;
}

@ The public entry point to the reader is |read_sexp|. This simply
resets the reader's global state and calls |read_form|.

@c
cell
read_sexp (void)
{
        cts_clear();
        Read_Level = 0;
        Putback[0] = Putback[1] = '\0';
        return read_form();
}

@ |read_form| reads a single (in most cases) byte which it uses to
determine which parser function to dispatch to. The parser function
will then return a complete s-expression (or raise an error).

@c
cell
read_form (void)
{
        cell r;
        int c, n;
        if (Interrupt)
                return VOID;
        if (Read_Level > READER_MAX_DEPTH)
                error(ERR_RECURSION, NIL);
        c = useful_byte();
        switch (c) { @<Reader forms@> }
        error(ERR_UNEXPECTED, NIL);
}

@ Here are the different bytes which |read_form| can understand,
starting with the non-byte value |EOF| which is an error if
the reader is part-way through parsing an expression.

@<Reader forms@>=
case EOF:
        if (!Read_Level)
                return END_OF_FILE;
        else
                error(ERR_ARITY_SYNTAX, NIL);

@ Lists and vectors are read in exactly the same way, differentiating
by being told to expect the appropriate delimiter.

@ @<Reader...@>=
case '(':
        return read_list(READ_CLOSE_PAREN);
case '[':
        return read_list(READ_CLOSE_BRACKET);

case ')':
case ']':
        /* If |Read_Level > 0| then |read_form| was called by
           |read_list|, otherwise |read_sexp| */
        if (!Read_Level)
                error(ERR_ARITY_SYNTAX, NIL);
        else
                return c == ')' ? READ_CLOSE_PAREN : READ_CLOSE_BRACKET;

@ A lone dot can only appear in a |list| and only before precisely
one more expression. This is verified later by |read_list|.

@<Reader...@>=
case '.':
        if (!Read_Level)
                error(ERR_ARITY_SYNTAX, NIL);
        c = useful_byte();
        if (c == EOF)
                error(ERR_ARITY_SYNTAX, NIL);
        unread_byte(c);
        return READ_DOT;

@ Special forms and strings aren't supported yet.

@<Reader...@>=
case '"':
case '#':
case '|':
        error(ERR_UNIMPLEMENTED, NIL);

@ In addition to the main syntactic characters, three other characters
commonly have special meaning in lisps: \qquote/, \qquasi/ and
\qunquote/. \qunquote/ can also appear as \qunsplice/. Primarily these
are for working with the macro expander.

In \LL/ this syntax is unnecessary thanks to its first-class
operatives but it's helpful so it's been retained. To differentiate
between having parsed the syntactic form of these operators (eg.
\qo\.{'foo}\qc\ or \qo\.{'(bar baz)}\qc) and their symbolic form
(eg. \qo\.{(quote . foo)}\qc\ or \qo\.{(quote bar baz)}\qc) an
otherwise ordinary |pair| with the operative's |symbol| in the |car|
is created with the tag |FORMAT_SYNTAX|. These |syntax| \.{object}s
are treated specially by the compiler and the writer.

@<Reader...@>=
case '\'':
case '`':
        n = useful_byte();
        if (n == EOF)
                error(ERR_ARITY_SYNTAX, NIL);
        unread_byte(n);
        if (n == ')' || n == ']')
                error(ERR_ARITY_SYNTAX, NIL);
        r = sym(c == '`' ? SYNTAX_QUASI : SYNTAX_QUOTE);
        return atom(r, read_form(), FORMAT_SYNTAX);

case ',':
        c = read_byte();
        if (c == EOF)
                error(ERR_ARITY_SYNTAX, NIL);
        if (c == ')' || c == ']')
                error(ERR_ARITY_SYNTAX, NIL);
        if (c == '@@') {
                r = sym(SYNTAX_UNSPLICE);
                return atom(r, read_form(), FORMAT_SYNTAX);
        } else {
                unread_byte(c);
                r = sym(SYNTAX_UNQUOTE);
                return atom(r, read_form(), FORMAT_SYNTAX);
        }

@ Anything else is a number or a symbol (and this byte is part of
it) provided it's ASCII.

@<Reader...@>=
default:
        if (!isprint(c))
                error(ERR_ARITY_SYNTAX, NIL);
        unread_byte(c);
        if (isdigit(c))
                return read_number();
        else
                return read_symbol();

@ |read_list| sequentially reads complete forms until it encounters
the closing delimiter \qclosepar/ or \qclosebra/.

A pointer to the head of the |list| is saved and another pointer
to its tail, |write|, is updated and used to insert the next object
after it's been read, avoiding the need to reverse the |list| at
the end.

@c
cell
read_list (cell delimiter)
{
        cell write, next, r;
        int count = 0;
        Read_Level++;
        write = cons(NIL, NIL);
        cts_push(write);
        while (1) {
                if (Interrupt) {
                        cts_pop();
                        Read_Level--;
                        return VOID;
                }
                next = read_form();
                if (special_p(next)) {
                        /* These must return or terminate unless |n| is
                           a `real' |special| */
                        @<Handle terminable `forms' during |list| construction@>@;
                }
                count++;
                cdr(write) = cons(NIL, NIL);
                write = cdr(write);
                car(write) = next;
        }
        Read_Level--;
        r = cdr(cts_pop());
        if (delimiter == READ_CLOSE_BRACKET)@/
                return vector_new_list(r, count);
        return count ? r : NIL;
}

@ |read_form| is expected to return an s-expression or raise an
error if the input is invalid. In order to recognise when a closing
parenthesis/bracket is read 3 `special' special forms are defined,
|READ_CLOSE_PAREN|, |READ_CLOSE_BRACKET| and |READ_DOT|. Although
these look and act like the other global constants they don't exist
outside of the parser.

@<Handle terminable `forms'...@>=
if (eof_p(next))
        error(ERR_ARITY_SYNTAX, NIL);
else if (next == READ_CLOSE_BRACKET
                || next == READ_CLOSE_PAREN) {
        if (next != delimiter)
                error(ERR_ARITY_SYNTAX, NIL);
        break;
} else if (next == READ_DOT) { @<Read dotted pair@> }

@ Encountering a \qdot/ requires more special care than it deserves,
made worse because if a |list| is dotted, a |syntax| \.{object}
is created instead of a normal s-expression so that the style in
which it's written out will be in the same that was read in.

@<Read dotted pair@>=
if (count < 1 || delimiter != READ_CLOSE_PAREN)
        /* There must be at least one item already and we must be
           parsing a \.{list}. */
        error(ERR_ARITY_SYNTAX, NIL);
next = read_form();
if (special_p(next) && next <= READ_SPECIAL)
        /* Check that the next `form' isn't one of \qdot/, \qclosepar/
           or \qclosebra/ */
        error(ERR_ARITY_SYNTAX, NIL);
cdr(write) = atom(sym(SYNTAX_DOTTED), next, FORMAT_SYNTAX);
next = read_form();
if (next != delimiter)@t\hbox to 42mm{}@>
        /* Check that the next `form' is really the closing delimiter */
        error(ERR_ARITY_SYNTAX, NIL);
break;

@ If it's not a |list| or a |vector| (or a |string| (\qstring/),
|special form| (\qspecial/), |raw symbol| (\qsymbol/) or |comment|) then
the form being read is an |atom|. If the atom starts with a numeric
digit then control proceeds directly to |read_number| otherwise
|read_symbol| reads enough to determine whether the atom is a number
beginning with $\pm$ or a valid or invalid symbol.

@d CHAR_TERMINATE "()[]\"; \t\r\n"
@d terminable_p(c) strchr(CHAR_TERMINATE, (c))
@c
cell
read_number (void)
{
        char buf[12] = {0}; /* $2^{32}$ is 10 digits, also $\pm$ and |NULL| */
        int c, i;
        long r;
        i = 0;
        while (1) {
                c = read_byte();
                if (c == EOF)
                        /* TODO: If |Read_Level| is 0 is this an error? */
                        error(ERR_ARITY_SYNTAX, NIL);
                if (i == 0 && (c == '-' || c == '+'))
                        buf[i++] = c;
                else if (isdigit(c))
                        buf[i++] = c;
                else if (!terminable_p(c))
                        error(ERR_ARITY_SYNTAX, NIL);
                else {
                        unread_byte(c);
                        break;
                }
                if (i > 11)
                        error(ERR_UNIMPLEMENTED, NIL);
        }
        r = atol(buf);
        if (r > INT_MAX || r < INT_MIN)
                error(ERR_UNIMPLEMENTED, NIL);
        return int_new(r);
}

@ Although \LL/ specifices (read: would specify) that there are no
restrictions on the value of a |symbol|'s label, memory permitting,
an artificial limit is being placed on the length of |symbol|s
of 16KB\footnote{$^1$}{640KB was deemed to be far more than enough
for anyone's needs.}.

That said, there are no restrictions on the value of a |symbol|'s
label, memory permitting. There are limits on what can be {\it
parsed} as a |symbol| in source code. The limits on plain |symbol|s
are primarily to avoid things that look vaguely like numbers to the
human eye being parsed as |symbol|s when the programmer thinks
they should be parsed as a number. This helps to avoid mistakes
like `\.{3..14159}' and harder to spot human errors being silently
ignored.

\point \* A |symbol| must not begin with a numeric digit or a
syntactic character (comments (\qcomment/), whitespace and everything
recognised by |read_form|).

\point \* The syntactic characters \qopenpar/, \qclosepar/, \qopenbra/,
\qclosebra/, \qcomment/ and \qstring/ cannot appear anywhere in the
|symbol|.

nb. This means that the following otherwise syntactic characters
{\it are} permitted in a symbol provided they do not occupy the
first byte: \qdot/, \qunquote/, \qquote/, \qquasi/, \qspecial/ and
\qsymbol/. You probably shouldn't do that lightly though.

\point \* If the first character of a |symbol| is \qo\.{-}\qc\
or \qo\.{+}\qc\ then it cannot be followed a numeric digit.
\.{++$\langle digit\rangle$} is valid.

\point \* A \qo\.{-}\qc\ character or \qo\.{+}\qc\ followed by a
\qdot/ is a valid if strange |symbol| but a warning should probably
be emitted by the parser if it finds that.

@d CHUNK_SIZE 0x80
@d READSYM_EOF_P if (c == EOF) error(ERR_ARITY_SYNTAX, NIL)
@c /* TODO: If |Read_Level| is 0 is this an error? */
cell
read_symbol (void)
{
        cell r;
        char *buf, *nbuf;
        int c, i, s;
        c = read_byte();
        READSYM_EOF_P;
        buf = malloc(CHUNK_SIZE);
        if (!buf)
                error(ERR_OOM, NIL);
        s = CHUNK_SIZE;
        @<Read the first two bytes to check for a number@>@;
        while (1) { @<Read bytes until an invalid or terminating character@> }
        buf[i] = '\0'; /* |NULL|-terminate the \CEE/-`string' */
        r = sym(buf);
        free(buf);
        return r;
}

@ Reading the first two bytes of a symbol is done specially to
detect numbers beginning with $\pm$. The first byte---which has
already been read to check for |EOF|---is put into |buf| then if
it matches $\pm$ the next byte is also read and also put into |buf|.

If that second byte is a digit then we're actually reading a number
so put the bytes that were read so far into |Putback| and go to
|read_number|, which will read them again. If the second byte is
\qdot/ then the |symbol| is valid but possibly a typo, so emit a
warning and carry on.

@<Read the first two...@>=
buf[0] = c;
i = 1;
if (c == '-' || c == '+') {
        c = read_byte();
        READSYM_EOF_P;
        buf[1] = c;
        i++;
        if (isdigit(buf[1])) {
                /* This is a number! */
                unread_byte(buf[1]);
                unread_byte(buf[0]);
                free(buf);
                return read_number();
        } else if (buf[1] == '.')
                warn(WARN_AMBIGUOUS_SYMBOL, NIL);
        else if (!isprint(c))
                error(ERR_ARITY_SYNTAX, NIL);
}

@ After the first two bytes we're definitely reading a |symbol|
so anything goes except non-printable characters (which are an
error) or syntactic terminators which indicate the end of the
|symbol|.

@<Read bytes until...@>=
c = read_byte();
READSYM_EOF_P;
if (terminable_p(c)) {
        unread_byte(c);
        break;
}
if (!isprint(c))
        error(ERR_ARITY_SYNTAX, NIL);
buf[i++] = c;
if (i == s) {
        /* Enlarge |buf| if it's now full (this will also allow the
           |NULL|-terminator to fit) */
        nbuf = realloc(buf, s *= 2);
        if (nbuf == NULL) {
                free(buf);
                error(ERR_OOM, NIL);
        }
        buf = nbuf;
}

@* Writer. Although not an essential part of the language itself,
the ability to display an s-expression to the user/programmer is
obviously invaluable.

It is expected that this will (very!) shortly be changed to return
a |string| representing the s-expression which can be passed on
to an output routine but for the time being \LL/ has no support for
|string|s or output routines so the expression is written directly
to |stdout|.

@d WRITER_MAX_DEPTH 1024 /* gotta pick something */
@<Function dec...@>=
void write_form (cell, int);

@*1 Opaque Objects. |applicative|s, |compiler|s and |operative|s
don't have much to say.

@c
boolean
write_applicative(cell sexp,
                  int depth __unused)
{
        if (!applicative_p(sexp))
                return bfalse;
        printf("#<applicative ...>");
        return btrue;
}

boolean
write_compiler(cell sexp,
               int depth __unused)
{
        if (!compiler_p(sexp))
                return bfalse;
        printf("#<compiler-%s>", compiler_cname(sexp));
        return btrue;
}

boolean
write_operative(cell sexp,
                int depth __unused)
{
        if (!operative_p(sexp))
                return bfalse;
        printf("#<operative ...>");
        return btrue;
}

@*1 As-Is Objects. |integer|s and |symbol|s print themselves.

@c
boolean
write_integer(cell sexp,
              int depth __unused)
{
        if (!integer_p(sexp))
                return bfalse;
        printf("%d", int_value(sexp));
        return btrue;
}

boolean
write_symbol(cell sexp,
             int depth __unused)
{
        int i;
        if (!symbol_p(sexp))
                return bfalse;
        for (i = 0; i < symbol_length(sexp); i++)
                putchar(symbol_store(sexp)[i]);
        return btrue;
}

@*1 Secret Objects. The hidden |syntax| object prints its syntactic
form and then itself.

@c
boolean
write_syntax(cell sexp,
             int depth)
{
        if (!syntax_p(sexp))
                return bfalse;
        else if (car(sexp) == sym(SYNTAX_DOTTED))   printf(". ");
        else if (car(sexp) == sym(SYNTAX_QUASI))    printf("`");
        else if (car(sexp) == sym(SYNTAX_QUOTE))    printf("'");
        else if (car(sexp) == sym(SYNTAX_UNQUOTE))  printf(",");
        else if (car(sexp) == sym(SYNTAX_UNSPLICE)) printf(",@@");
        write_form(cdr(sexp), depth + 1);
        return btrue;
}

@*1 Environment Objects. An |environment| prints its own layer
and then the layers above it.

@c
boolean
write_environment(cell sexp,
                  int depth)
{
        if (!environment_p(sexp))
                return bfalse;
        printf("#<environment ");
        write_form(env_layer(sexp), depth + 1);
        if (!null_p(env_parent(sexp))) {
                printf(" ON ");
                write_form(env_parent(sexp), depth + 1);
                printf(">");
        } else
                printf(" ROOT>");
        return btrue;
}

@*1 Sequential Objects. The routines for a |list| and |vector| are
more or less the same -- write each item in turn with whitespace
after each form but the last, with the appropriate delimiters.
|list|s also need to deal with being improper.

@c
boolean
write_list(cell sexp,
           int depth)
{
        if (!pair_p(sexp))
                return bfalse;
        printf("(");
        while (pair_p(sexp)) {
                write_form(car(sexp), depth + 1);
                if (pair_p(cdr(sexp)) || syntax_p(cdr(sexp)))
                        printf(" ");
                else if (!null_p(cdr(sexp))
                         && !pair_p(cdr(sexp))
                         && !syntax_p(cdr(sexp)))
                        printf(" . ");
                sexp = cdr(sexp);
        }
        if (!null_p(sexp))
                write_form(sexp, depth + 1);
        printf(")");
        return btrue;
}

boolean
write_vector(cell sexp,
             int depth)
{
        int i;
        if (!vector_p(sexp))
                return bfalse;
        printf("[");
        for (i = 0; i < vector_length(sexp); i++) {
                write_form(vector_ref(sexp, i), depth + 1);
                if (i + 1 < vector_length(sexp))
                        printf(" ");
        }
        printf("]");
        return btrue;
}

@ |write_form| simply calls each writer in turn, stopping after the
first one returning (\CEE/'s) true.

@c
void
write_form (cell sexp,
            int  depth)
{
        if (Interrupt) {
                if (!depth)
                        printf("... \n");
                return;
        }
        if (depth > WRITER_MAX_DEPTH)
                error(ERR_RECURSION, NIL);
        if (undefined_p(sexp))
                printf("#><"); /* nothing should ever print this */
        else if (eof_p(sexp))
                printf("#<eof>");
        else if (false_p(sexp))
                printf("#f");
        else if (null_p(sexp))
                printf("()");
        else if (true_p(sexp))
                printf("#t");
        else if (void_p(sexp))
                printf("#<>");
        else if (write_applicative(sexp, depth)) /* NOP */@+;
        else if (write_compiler(sexp, depth))    /* NOP */@+;
        else if (write_environment(sexp, depth)) /* NOP */@+;
        else if (write_integer(sexp, depth))     /* NOP */@+;
        else if (write_list(sexp, depth))        /* NOP */@+;
        else if (write_operative(sexp, depth))   /* NOP */@+;
        else if (write_symbol(sexp, depth))      /* NOP */@+;
        else if (write_syntax(sexp, depth))      /* NOP */@+;
        else if (write_vector(sexp, depth))      /* NOP */@+;
        else printf("#<wtf?>");                  /* impossibru! */
}

@** Opcodes. With the core infrastucture out of the way we can
finally turn our attention to the virtual machine implementation,
or the implementation of the opcodes that the compiler will turn
\LL/ code into.

The opcodes that the virtual machine can perform must be declared
before anything can be said about them. They take the form of an
|enum|, this one unnamed. This list is sorted alphabetically for
want of anything else.

Also defined here are |fetch| and |skip| which |opcode| implementations
will use to obtain their argument(s) from |Prog| or advance |Ip|,
respectively.

@d skip(d) Ip += (d)
@d fetch(d) vector_ref(Prog, Ip + (d))
@<Global constants@>=
enum {
        OP_APPLY,
        OP_APPLY_TAIL,
        OP_CAR,
        OP_CDR, /* 3 */
        OP_COMPILE,
        OP_CONS,
        OP_CYCLE,
        OP_ENVIRONMENT_P, /* 7 */
        OP_ENV_MUTATE_M,
        OP_ENV_QUOTE,
        OP_ENV_ROOT,
        OP_ENV_SET_ROOT_M, /* 11 */
        OP_ERROR,
        OP_HALT,
        OP_JUMP,
        OP_JUMP_FALSE, /* 15 */
        OP_JUMP_TRUE,
        OP_LAMBDA,
        OP_LIST_P,
        OP_LIST_REVERSE, /* 19 */
        OP_LIST_REVERSE_M,
        OP_LOOKUP,
        OP_NIL,
        OP_NOOP, /* 23 */
        OP_NULL_P,
        OP_PAIR_P,
        OP_PEEK,
        OP_POP, /* 27 */
        OP_PUSH,
        OP_QUOTE,
        OP_RETURN,
        OP_RUN, /* 31 */
        OP_RUN_THERE,
        OP_SET_CAR_M,
        OP_SET_CDR_M,
        OP_SNOC, /* 35 */
        OP_SWAP,
        OP_SYNTAX,
        OP_VOV,
        OPCODE_MAX
};

@* Basic Flow Control. The most basic opcodes that the virtual machine
needs are those which control whether to operate and where.

@ @<Opcode imp...@>=
case OP_HALT:@/
        Running = 0;
        break;
case OP_JUMP:@/
        Ip = int_value(fetch(1));
        break;
case OP_JUMP_FALSE:@/
        if (false_p(Acc))
                Ip = int_value(fetch(1));
        else
                skip(2);
        break;
case OP_JUMP_TRUE:@/
        if (true_p(Acc))
                Ip = int_value(fetch(1));
        else
                skip(2);
        break;
case OP_NOOP:@/
        skip(1);
        break;

@ |OP_QUOTE| isn't really flow control but I don't know where else
to put it.

@<Opcode imp...@>=
case OP_QUOTE:@/
        Acc = fetch(1);
        skip(2);
        break;

@* Pairs \AM\ Lists. |OP_CAR|, |OP_CDR|, |OP_NULL_P| and |OP_PAIR_P|
are self explanatory.

@<Opcode imp...@>=
case OP_CAR:@/
        Acc = car(Acc);
        skip(1);
        break;
case OP_CDR:@/
        Acc = cdr(Acc);
        skip(1);
        break;
case OP_NULL_P:@/
        Acc = null_p(Acc) ? TRUE : FALSE;
        skip(1);
        break;
case OP_PAIR_P:@/
        Acc = pair_p(Acc) ? TRUE : FALSE;
        skip(1);
        break;

@ |OP_CONS| consumes one stack item (for the |cdr|) and puts the
new pair in |Acc|. |OP_SNOC| does the opposite, pushing |Acc|'s
|cdr| to the stack and leaving its |car| in |Acc|.
@<Opcode imp...@>=
case OP_CONS:@/
        Acc = cons(Acc, rts_pop(1));
        skip(1);
        break;
case OP_SNOC:@/
        rts_push(cdr(Acc));
        Acc = car(Acc);
        skip(1);
        break;

@ Cons cell mutators clear take an item from the stack and clear |Acc|.
@<Opcode imp...@>=
case OP_SET_CAR_M:@/
        car(rts_pop(1)) = Acc;
        Acc = VOID;
        skip(1);
        break;
case OP_SET_CDR_M:@/
        car(rts_pop(1)) = Acc;
        Acc = VOID;
        skip(1);
        break;

@* Other Objects. There is not much to say about these.

@<Opcode imp...@>=
case OP_LIST_P:@/
        if (!false_p(fetch(2)))
                error(ERR_UNIMPLEMENTED, NIL);
        Acc = list_p(Acc, fetch(1), NULL);
        skip(3);
        break;
case OP_LIST_REVERSE:@/
        if (!true_p(fetch(1)) || !false_p(fetch(2)))
                error(ERR_UNIMPLEMENTED, NIL);
        Acc = list_reverse(Acc, NULL, NULL);
        skip(3);
        break;
case OP_LIST_REVERSE_M:@/
        Acc = list_reverse_m(Acc, btrue);
        skip(1);
        break;
case OP_SYNTAX:@/
        Acc = atom(fetch(1), Acc, FORMAT_SYNTAX);
        skip(2);
        break;

@* Stack. |OP_PUSH| and |OP_POP| push the \.{object} in |Acc| onto
the stack, or remove the top stack \.{object} into |Acc|, respectively.
|OP_PEEK| is |OP_POP| without removing the item from the stack.

|OP_SWAP| swaps the \.{object} in |Acc| with the \.{object} on top
of the stack.

|OP_CYCLE| swaps the top two stack items with each other.

|OP_NIL| pushes a |NIL| straight onto the stack without the need
to quote it first.

@<Opcode imp...@>=
case OP_CYCLE:@/
        tmp = rts_ref(0);
        rts_set(0, rts_ref(1));
        rts_set(1, tmp);
        skip(1);
        break;
case OP_PEEK:@/
        Acc = rts_ref(0);
        skip(1);
        break;
case OP_POP:@/
        Acc = rts_pop(1);
        skip(1);
        break;
case OP_PUSH:@/
        rts_push(Acc);
        skip(1);
        break;
case OP_SWAP:@/
        tmp = Acc;
        Acc = rts_ref(0);
        rts_set(0, tmp);
        skip(1);
        break;
case OP_NIL:@/
        rts_push(NIL);
        skip(1);
        break;

@* Environments. Get or mutate |environment| objects. |OP_ENV_SET_ROOT_M|
isn't used yet.

@<Opcode imp...@>=
case OP_ENVIRONMENT_P:@/
        Acc = environment_p(Acc) ? TRUE : FALSE;
        skip(1);
        break;
case OP_ENV_MUTATE_M:@/
        env_set(rts_pop(1), fetch(1), Acc, true_p(fetch(2)));
        Acc = VOID;
        skip(3);
        break;
case OP_ENV_QUOTE:@/
        Acc = Env;
        skip(1);
        break;
case OP_ENV_ROOT:@/
        Acc = Root;
        skip(1);
        break;
case OP_ENV_SET_ROOT_M:@/
        Root = Acc; /* |Root| is `lost'! */
        skip(1);
        break;

@ To look up the value of a variable in an |environment| we use
|OP_LOOKUP| which calls the (recursive) |env_search|, interpreting
the |UNDEFINED| it might return.

@<Opcode imp...@>=
case OP_LOOKUP:@/
        vms_push(Acc);
        Acc = env_search(Env, vms_ref());
        if (undefined_p(Acc)) {
                Acc = vms_pop();
                error(ERR_UNBOUND, Acc);
        }
        vms_pop();
        skip(1);
        break;

@* Closures. A |closure| is the combination of code to interpret
and an |environment| to interpret it in. Usually a closure has
arguments---making it useful---although in some cases a closure may
work with global state or be idempotent.

In order to apply the arguments (if any) to the |closure| it must
be entered by one of the opcodes |OP_APPLY| or |OP_APPLY_TAIL|.
|OP_APPLY_TAIL| works identically to |OP_APPLY| and then consumes
the stack frame which |OP_APPLY| created, allowing for {\it proper
tail recursion} with further support from the compiler.

@<Opcode imp...@>=
case OP_APPLY:@/
        @<Enter a |closure|@>@;
        break;
case OP_APPLY_TAIL:@/
        @<Enter a |closure|@>@;
        frame_consume();
        break;
case OP_RETURN:@/
        frame_leave();
        break;

@ Whether in tail position or not, entering a |closure| is the same.

% This should be a function.
@<Enter a |closure|@>=
{
        cell e, i, p;
        tmp = fetch(2);
        e = env_lift_stack(cadr(tmp), int_value(fetch(1)), car(tmp));
        p = caddr(tmp);
        i = int_value(cadddr(tmp));
        frame_push(3);
        frame_enter(e, p, i);
}

@ Creating a closure in the first place follows an identical procedure
whether it's an applicative or an operative but creates a different
type of \.{object} in each case.

@<Opcode imp...@>=
case OP_LAMBDA:@/ /* The |applicative| */
        Acc = applicative_new(rts_pop(1), Env, Prog, int_value(fetch(1)));
        skip(2);
        break;
case OP_VOV:@/ /* The |operative| */
        Acc = operative_new(rts_pop(1), Env, Prog, int_value(fetch(1)));
        skip(2);
        break;

@* Compiler. The compiler needs to instruct the interpreter to
compile more code and then run it, so these |opcode|s do that.
|OP_COMPILE| compiles an s-expression into \LL/ bytecode.

@ @<Opcode imp...@>=
case OP_COMPILE:@/
        Acc = compile(Acc);
        skip(1);
        break;

@ |OP_RUN| interprets the bytecode in |Acc| in the current
|environment|; the VM's live state is saved into a new stack
|frame| then that |frame| is entered by executing the bytecode
in |Acc|, starting at instruction 0.

@<Opcode imp...@>=
case OP_RUN:@/
        frame_push(1);
        frame_enter(Env, Acc, 0);
        break;

@ |OP_RUN_THERE| is like |OP_RUN| except that the |environment|
to interpret the bytecode in is taken from the stack rather than
staying in the active |environment|.

@<Opcode imp...@>=
case OP_RUN_THERE:@/
        vms_push(rts_pop(1));
        frame_push(1);
        frame_enter(vms_pop(), Acc, 0);
        break;

@** Compiler. Speaking of the compiler, we can now turn our attention
to writing it. The compiler is not advanced in any way but it is a
little unusual. Due to the nature of first-class operatives, how
to compile any expression can't be known until the combinator has
been evaluated (read: compiled and then interpreted) in order to
distinguish an applicative from an operative so that it knows whether
to evaluate the arguments in the expression. I don't know if this
qualifies it for a {\it Just-In-Time} compiler; I think {\it
Finally-Able-To} is more suitable.

The compiler uses a small set of \CEE/ macros which grow and fill
|Compilation|---a |vector| holding the compilation in-progress.

@d ERR_COMPILE_DIRTY "compiler"
@d ERR_UNCOMBINABLE "uncombinable"
@d COMPILATION_SEGMENT 0x80
@#
@d emitop(o)  emit(int_new(o))
@d emitq(o)   do@+ { @+emitop(OP_QUOTE);@+ emit(o);@+ } while (0) /* \CEE/... */
@d patch(i,v) (vector_ref(Compilation, (i)) = (v))
@d undot(p)   ((syntax_p(p) && car(p) == Sym_Dotted) ? cdr(p) : (p))
@<Global var...@>=
int Here = 0;
cell Compilation = NIL;
cell Sym_Dotted = UNDEFINED;

@ @<Function dec...@>=
cell compile (cell);
cell compile_main (void);
void compile_car (cell, cell, boolean);
void compile_cdr (cell, cell, boolean);
void compile_conditional (cell, cell, boolean);
void compile_cons (cell, cell, boolean);
void compile_define_m (cell, cell, boolean);
void compile_env_current (cell, cell, boolean);
void compile_env_root (cell, cell, boolean);
void compile_error (cell, cell, boolean);
void compile_eval (cell, cell, boolean);
void compile_expression (cell, boolean);
void compile_lambda (cell, cell, boolean);
void compile_null_p (cell, cell, boolean);
void compile_pair_p (cell, cell, boolean);
void compile_quasiquote (cell, cell, boolean);
void compile_quote (cell, cell, boolean);
void compile_set_car_m (cell, cell, boolean);
void compile_set_cdr_m (cell, cell, boolean);
void compile_set_m (cell, cell, boolean);
void compile_vov (cell, cell, boolean);

@ @<Protected...@>=
&Compilation,

@ @<Global init...@>=
Sym_Dotted = sym(SYNTAX_DOTTED);

@ @c
void
emit (cell bc)
{
        int l;
        l = vector_length(Compilation);
        if (Here >= l)@/
                Compilation = vector_sub(Compilation,
                        0, l,@|
                        0, l + COMPILATION_SEGMENT,@|
                        OP_HALT);
        vector_ref(Compilation, Here++) = bc;
}

@ While compiling it frequently occurs that the value to emit isn't
known at the time it's being emitted. The most common and obvious
example of this is a forward jump who's address must immediately
follow the opcode but the address won't be known until more compilation
has been performed.

To make this work |comefrom| emits a |NIL| as a placeholder and
returns its offset, which can later be passed in the first argument
of |patch| to replace the |NIL| with the desired address etc.

@c
int
comefrom (void)
{
        emit(NIL);
        return Here - 1;
}

@ Compilation begins by preparing |Compilation| and |CTS| then
recursively walks the tree in |source| dispatching to individual
compilation routines to emit the appropriate bytecode.

@c
cell
compile (cell source)
{
        cell r;
        vms_push(source);
        Compilation = vector_new(COMPILATION_SEGMENT, int_new(OP_HALT));
        Here = 0;
        cts_reset();
        compile_expression(source, 1);
        emitop(OP_RETURN);
        r = vector_sub(Compilation, 0, Here, 0, Here, VOID);
        Compilation = NIL;
        vms_clear();
        if (!null_p(CTS))
                error(ERR_COMPILE_DIRTY, source);
        return r;
}

@ |compile_main| is used during initialisation to build the bytecode
$\ll$|OP_COMPILE|\ |OP_RUN|\ |OP_HALT|$\gg$ which is the program
installed initially into the virtual machine.

@c
cell
compile_main (void)
{
        cell r;
        r = vector_new_imp(3, 0, 0);
        vector_ref(r, 0) = int_new(OP_COMPILE);
        vector_ref(r, 1) = int_new(OP_RUN);
        vector_ref(r, 2) = int_new(OP_HALT);
        return r;
}

@ The first job of the compiler is to figure out what type of
expression it's compiling, chiefly whether it's a |list| to combine
or an |atom| which is itself.

@c
void
compile_expression (cell sexp,
                    int  tail_p)
{
        if (!pair_p(sexp) && !syntax_p(sexp)) {
                @<Compile an atom@>
        } else {
                @<Compile a combiner@>
        }
}

@ The only |atom| which doesn't evaluate to itself is a |symbol|.
A |symbol| being evaluated references a variable which must be
looked up in the active environment.

@<Compile an atom@>=
if (symbol_p(sexp)) {
        emitq(sexp);
        emitop(OP_LOOKUP);
} else {@+
        emitq(sexp);@+
}

@ Combining a |list| requires more work. This is also where
operatives obtain the property of being first-class objects by
delaying compilation of all but the first expression in the |list|
until after that compiled bytecode has been interpreted.

@<Compile a combiner@>=
cell args, combiner;
combiner = car(sexp);
args = undot(cdr(sexp));
@<Search |Root| for syntactic combiners@>@;
if (compiler_p(combiner)) {
        @<Compile native combiner@>
} else if (applicative_p(combiner)) {
        @<Compile applicative combiner@>
} else if (operative_p(combiner)) {
        @<Compile operative combiner@>
} else if (symbol_p(combiner) || pair_p(combiner)) {
        @<Compile unknown combiner@>
} else {@+
        error(ERR_UNCOMBINABLE, combiner);@+
}

@ If the combiner (|sexp|'s |car|) is a |syntax| object then it
represents the result of parsing (for example) \qo\.{'(expression)}\qc\
into \qo\.{(quote expression)}\qc\ and it must always mean the {\it
real} |quote| operator, so |syntax| combiners are always looked for
directly (and only) in |Root|.

@<Search |Root|...@>=
if (syntax_p(sexp)) {
        cell c;
        c = env_search(Root, combiner);
        if (undefined_p(c))
                error(ERR_UNBOUND, combiner); /* should never happen */
        combiner = c;
}

@ A native compiler is simple; look up its address in |COMPILER|
and go there. The individual native compilers are defined below.

@<Compile native...@>=
compiler_fn(combiner)(combiner, args, tail_p);

@ If the compiler doesn't know whether |combiner| is applicative
or operative then that must be determined before |args| can be
considered.

@<Compile unknown...@>=
emitq(args);
emitop(OP_PUSH);    /* save |args| onto the stack */
compile_expression(combiner, 0);
                    /* evaluate the combiner, leaving it in |Acc| */
emitop(OP_CONS);    /* rebuild |sexp| with the evaluated combiner */
emitop(OP_COMPILE); /* continue compiling |sexp| */
emitop(OP_RUN);     /* run that code in the same |environment| */

@* Function Bodies. Nearly everything has arguments to process and
it's nearly always done in the same way. |arity| and |arity_next|
work in concert to help the compiler implementations check how many
arguments there are (but not their value or type) and raise any
errors encountered.

|arity| pushes the minimum required arguments onto the compiler
stack (in reverse) and returns a pointer to the rest of the argument
list.

@d ERR_ARITY_EXTRA   "extra"
@d ERR_ARITY_MISSING "missing"
@d ERR_ARITY_SYNTAX  "syntax"
@d arity_error(e,c,a) error((e), cons((c), (a)))
@c
cell
arity (cell    op,
       cell    args,
       int     min,
       boolean more_p)
{
        cell a = args;
        int i = 0;
        for (; i < min; i++) {
                if (null_p(a)) {
                        if (compiler_p(op) || operative_p(op))
                                arity_error(ERR_ARITY_SYNTAX, op, args);
                        else
                                arity_error(ERR_ARITY_MISSING, op, args);
                }
                if (!pair_p(a))
                        arity_error(ERR_ARITY_SYNTAX, op, args);
                cts_push(car(a));
                a = cdr(a);
        }
        if (min && !more_p && !null_p(a)) {
                if (pair_p(a))
                        arity_error(ERR_ARITY_EXTRA, op, args);
                else
                        arity_error(ERR_ARITY_SYNTAX, op, args);
        }
        return a;
}

@ |arity_next|, given the remainder of the arguments that were
returned from |arity|, checks whether another one is present and
whether it's allowed to be, then returns a value suitable for another
call to |arity_next|.

@c
cell
arity_next (cell    op,
            cell    args,
            cell    more,
            boolean required_p,
            boolean last_p)
{
        if (null_p(more)) {
                if (required_p)
                        arity_error(ERR_ARITY_MISSING, op, args);
                else {
                        cts_push(UNDEFINED);
                        return NIL;
                }
        } else if (!pair_p(more))@/
                arity_error(ERR_ARITY_SYNTAX, op, args);
        else if (last_p && !null_p(cdr(more))) {
                if (operative_p(op) && pair_p(cdr(more)))
                        arity_error(ERR_ARITY_EXTRA, op, args);
                else
                        arity_error(ERR_ARITY_SYNTAX, op, args);
        }
        cts_push(car(more));
        return cdr(more);
}

@ |closure| bodies, and the contents of a |begin| expression, are
compiled by simply walking the list and recursing into |compile_expression|
for everything on it. When compiling the last item in the list the
|tail_p| flag is raised so that the expression can use |OP_APPLY_TAIL|
if appropriate, making tail recursion proper.

@c
void
compile_list (cell    op,
              cell    sexp,
              boolean tail_p)
{
        boolean t;
        cell body, next, this;
        body = undot(sexp);
        t = null_p(body);
        if (t) {
                emitq(VOID);
                return;
        }
        while (!t) {
                if (!pair_p(body))
                        arity_error(ERR_ARITY_SYNTAX, op, sexp);
                this = car(body);
                next = undot(cdr(body));
                t = null_p(next);
                compile_expression(this, t && tail_p);
                body = next;
        }
}

@* Closures (Applicatives \AM\ Operatives). The first thing to
understand is that at their core |applicative|s and |operative|s
work in largely the same way and have the same internal representation:

\point \* The static |environment| which will expand into a local
|environment| when entering the |closure|. This is where the variables
that were ``closed over'' are stored.

\point \* The program which the |closure| will perform, as compiled
bytecode and a starting instruction pointer.

\point \* A list of formals naming any arguments which will be
passed to the |closure|, so that they can be put into the newly-extended
|environment|.

Entering a |closure| means extracting these saved values and restoring
them to the virtual machine's registers, |Env|, |Prog| \AM\ |Ip|.

A |closure| can (usually does) have arguments and it's how they're
handled that differentiates an |applicative| from an |operative|.

@ The main type of |closure| everyone is familiar with already even
if they don't know it is a function\footnote{$^1$}{The word
``function'' is horribly misused everywhere and this trend will
continue without my getting in its way.} or |applicative|.

An |applicative| is created in response to evaluating a |lambda|
expression. The bytecode which does this evaluating is created by
|compile_lambda|.

@c
void
compile_lambda (cell op,
                cell args,
                boolean tail_p)
{
        cell body, in, formals, f;
        int begin_address, comefrom_end;
        body = arity(op, args, 1, 1);
        body = undot(body);
        formals = cts_pop();
        formals = undot(formals);
        if (!symbol_p(formals)) { @<Process lambda formals@> }
        emitq(formals); /* push |formals| onto the stack */
        emitop(OP_PUSH);
        emitop(OP_LAMBDA); /* create the |applicative| */
        begin_address = comefrom(); /* start address; argument to
                                       |OP_LAMBDA| */
        emitop(OP_JUMP); /* jump over the compiled |closure| body */
        comefrom_end = comefrom();
        patch(begin_address, int_new(Here));
        compile_list(op, body, tail_p); /* compile the code that entering
                                           the |closure| will interpret */
        emitop(OP_RETURN); /* returns from the |closure| at run-time */
        patch(comefrom_end, int_new(Here));
}

@ If the |formals| given in the |lambda| expression are not in fact
a single |symbol| then it must be a list of |symbol|s which is
verified here. At the same time if the list is a dotted pair then
the |syntax| wrapper is removed.

@<Process lambda formals@>=
cts_push(f = cons(NIL, NIL));
in = formals;
while (pair_p(in)) {
        if (!symbol_p(car(in)) && !null_p(car(in)))
                arity_error(ERR_ARITY_SYNTAX, op, args);
        cdr(f) = cons(car(in), NIL);
        f = cdr(f);
        in = undot(cdr(in));
}
if (!null_p(in)) {
        if (!symbol_p(in) && !null_p(in))
                arity_error(ERR_ARITY_SYNTAX, op, args);
        cdr(f) = in;
}
formals = cdr(cts_pop());

@ To enter this |closure| at run-time---aka. to call the function
returned by |lambda|---the arguments it's called with must be
evaluated (after being arity checked) then |OP_APPLY| or |OP_APPLY_TAIL|
enters the |closure|, consuming a stack |frame| in the latter case.

The arguments and the formals saved in the |applicative| are walked
together and saved in |direct|. If the formals list ends in a dotted
|pair| then the remainder of the arguments are saved in |collect|.

When |collect| and |direct| have been prepared, being a copy of the
unevaluated arguments in reverse order, they are walked again to
emit the opcodes which will evaluate each argument and put the
results onto the stack.

@<Compile applicative...@>=
cell collect, direct, formals, a;
int nargs = 0;

formals = applicative_formals(combiner);
cts_push(direct = NIL);
a = undot(args);
@<Look for required arguments@>@;
@<Look for optional arguments@>@;
if (pair_p(a))
        arity_error(ERR_ARITY_EXTRA, combiner, args);
else if (!null_p(a))
        arity_error(ERR_ARITY_SYNTAX, combiner, args);
@<Evaluate optional arguments into a |list|@>@;
@<Evaluate required arguments onto the stack@>@;
cts_clear();
emitop(tail_p ? OP_APPLY_TAIL : OP_APPLY);
emit(int_new(nargs));
emit(combiner);

@ It's a syntax error if the arguments are not a proper list,
otherwise there is nothing much to say about this.

@<Look for required...@>=
while (pair_p(formals)) {
        if (!pair_p(a)) {
                if (null_p(a))
                        arity_error(ERR_ARITY_SYNTAX, combiner, args);
                else
                        arity_error(ERR_ARITY_SYNTAX, combiner, args);
        }
        direct = cons(car(a), direct);
        cts_set(direct);
        a = undot(cdr(a));
        formals = cdr(formals);
        nargs++;
}

@ If the |applicative| formals indicate that it can be called with
a varying number of arguments then that counts as one more argument
which will be a list of whatever arguments remain.

@<Look for optional...@>=
if (symbol_p(formals)) {
        nargs++;
        cts_push(collect = NIL);
        while (pair_p(a)) {
                collect = cons(car(a), collect);
                cts_set(collect);
                a = undot(cdr(a));
        }
}

@ To perform the evaluation, each argument in the (now reversed)
list |direct| is compiled followed by an |OP_PUSH| to save the
result on the stack.

@<Evaluate required...@>=
while (!null_p(direct)) {
        compile_expression(car(direct), 0);
        emitop(OP_PUSH);
        direct = cdr(direct);
}

@ If the |applicative| expects a varying number of arguments then
the (also reversed) list in |collect| is compiled in the same way
but before |OP_PUSH|, |OP_CONS| removes the growing list from the
stack and prepends the new result to it and it's this |list| which
is pushed.

@<Evaluate optional...@>=
if (symbol_p(formals)) {
        emitop(OP_NIL);
        while (!null_p(collect)) {
                compile_expression(car(collect), 0);
                emitop(OP_CONS);
                emitop(OP_PUSH);
                collect = cdr(collect);
        }
        cts_clear();
}

@ Analogous to |compile_lambda| for |applicative|s is |compile_vov|
for |operative|s. An |operative| |closure| is a simpler than an
|applicative| because the arguments are not evaluated. Instead
|compile_vov| needs to handle |vov|'s very different way of specifying
its formals.

Resembling |let| rather than |lambda|, |vov|'s formals specify what
run-time detail the |operative| needs: The unevaluated arguments,
the active environment and/or (unimplemented) a continuation
delimiter. To do this each entry in the formals list is an association
pair with the |symbol|ic name for that detail associated with another
symbol specifying what: {\it vov/arguments}, {\it vov/environment}
or {\it vov/continuation}. Because no-one wants RSI these have the
abbreviations {\it vov/args}, {\it vov/env} and {\it vov/cont}.

@<Global var...@>=
cell Sym_vov_args      = UNDEFINED;
cell Sym_vov_args_long = UNDEFINED;
cell Sym_vov_cont      = UNDEFINED;
cell Sym_vov_cont_long = UNDEFINED;
cell Sym_vov_env       = UNDEFINED;
cell Sym_vov_env_long  = UNDEFINED;

@ @<Global init...@>=
Sym_vov_args      = sym("vov/args");
Sym_vov_args_long = sym("vov/arguments");
Sym_vov_cont      = sym("vov/cont");
Sym_vov_cont_long = sym("vov/continuation");
Sym_vov_env       = sym("vov/env");
Sym_vov_env_long  = sym("vov/environment");

@ @c
void
compile_vov (cell op,
             cell args,
             boolean tail_p)
{
        cell body, formals;
        int begin_address, comefrom_end;
        cell a = NIL;
        cell c = NIL;
        cell e = NIL;
        body = arity(op, args, 1, 1);
        body = undot(body);
        formals = cts_pop();
        formals = undot(formals);
        @<Scan operative informals@>@;
        emitop(OP_NIL); /* push formals onto the stack */
        emitq(c);@+ emitop(OP_CONS);@+ emitop(OP_PUSH);
        emitq(e);@+ emitop(OP_CONS);@+ emitop(OP_PUSH);
        emitq(a);@+ emitop(OP_CONS);@+ emitop(OP_PUSH);
        emitop(OP_VOV); /* create the |operative| */
        /* The rest of |compile_vov| is identical to |compile_lambda|: */
        begin_address = comefrom(); /* start address; argument to opcode */
        emitop(OP_JUMP); /* jump over the compiled |closure| body */
        comefrom_end = comefrom();
        patch(begin_address, int_new(Here));
        compile_list(op, body, tail_p); /* compile the code that entering
                                           the |closure| will interpret */
        emitop(OP_RETURN); /* return from the run-time |closure| */
        patch(comefrom_end, int_new(Here)); /* finish building the |closure| */
}

@ To scan the ``informals'' three variables, |a|, |c| and |e| are
prepared with |NIL| representing the |symbol| for the arguments,
|continuation| and |environment| respectively. Each ``informal''
is checked in turn using |arity| and the appropriate placeholder's
|NIL| replaced with the |symbol|.

@<Scan operative informals@>=
cell r, s;
if (!pair_p(formals))
        arity_error(ERR_ARITY_SYNTAX, op, args);
#define CHECK_AND_ASSIGN(v) {                            \
        if (!null_p(v))                                  \
                arity_error(ERR_ARITY_SYNTAX, op, args); \
        (v) = s;                                         \
}
while (pair_p(formals)) {
        arity(op, car(formals), 2, 0);
        r = cts_pop();
        s = cts_pop();
        if (!symbol_p(s))
                arity_error(ERR_ARITY_SYNTAX, op, args);
        else if (r == Sym_vov_args || r == Sym_vov_args_long)
                CHECK_AND_ASSIGN(a)@;
        else if (r == Sym_vov_env || r == Sym_vov_env_long)
                CHECK_AND_ASSIGN(e)@;
        else if (r == Sym_vov_cont || r == Sym_vov_cont_long)
                CHECK_AND_ASSIGN(c)@;
        formals = cdr(formals);
}
if (!null_p(formals))
        arity_error(ERR_ARITY_SYNTAX, op, args);

@ Entering an |operative| involves pushing the 3 desired run-time
properties, or |NIL|, onto the stack as though arguments to an
|applicative| |closure| (remember that the unevaluated run-time
arguments of the |closure| are potentially one of those run-time
properties).

@<Compile operative...@>=
cell a, c, e, f;
f = operative_formals(combiner);
a = !null_p(car(f));@+ f = cdr(f);
e = !null_p(car(f));@+ f = cdr(f);
c = !null_p(car(f));@+ f = cdr(f);

if (c)
        error(ERR_UNIMPLEMENTED, NIL);
else
        emitop(OP_NIL);
if (e) {
        emitop(OP_ENV_QUOTE);
        emitop(OP_PUSH);
} else
        emitop(OP_NIL);
if (a) {
        emitq(args);
        emitop(OP_PUSH);
} else
        emitop(OP_NIL);

emitop(tail_p ? OP_APPLY_TAIL : OP_APPLY);
emit(int_new(3));
emit(combiner);

@* Conditionals (|if|). Although you could define a whole language
with just |lambda| and |vov|\footnote{$^1$}{In fact I think
conditionals can be achieved in both somehow, so you only need one.}
that way lies Church Numerals and other madness, so we will define
the basic conditional, |if|.

@c
void
compile_conditional (cell op,
                     cell args,
                     boolean tail_p)
{
        cell alternate, condition, consequent, more;
        int jump_false, jump_true;
        more = arity(op, args, 2, 1);
        arity_next(op, args, more, 0, 1);
        alternate = cts_pop();
        consequent = cts_pop();
        condition = cts_pop();
        compile_expression(condition, 0);
        emitop(OP_JUMP_FALSE);
        jump_false = comefrom();
        compile_expression(consequent, tail_p);
        emitop(OP_JUMP);
        jump_true = comefrom();
        patch(jump_false, int_new(Here));
        if (undefined_p(alternate))
                emitq(VOID);
        else
                compile_expression(alternate, tail_p);
        patch(jump_true, int_new(Here));
}

@* Run-time Evaluation (|eval|). |eval| must evaluate its 1 or 2
arguments in the current environment, and then enter the environment
described by the second to execute the program in the first.

@c
void
compile_eval (cell op,
              cell args,
              boolean tail_p __unused)
{
        cell more, sexp, eenv;
        int goto_env_p;
        more = arity(op, args, 1, 1);
        sexp = cts_pop();
        arity_next(op, args, more, 0, 1);
        eenv = cts_pop();
        if (undefined_p(eenv)) {
                emitop(OP_ENV_QUOTE);
                emitop(OP_PUSH);
        } else {
                compile_expression(eenv, 0);
                emitop(OP_PUSH);
                emitop(OP_ENVIRONMENT_P);
                emitop(OP_JUMP_TRUE);
                goto_env_p = comefrom();
                emitq(Sym_ERR_UNEXPECTED);
                emitop(OP_ERROR);
                patch(goto_env_p, int_new(Here));
        }
        compile_expression(sexp, 0);
        emitop(OP_COMPILE);
        emitop(OP_RUN_THERE);
}

@* Run-time Errors. |error| expects a symbol an the first position
and an optional form to evaluate in the second.

@c
void
compile_error (cell op,
               cell args,
               boolean tail_p __unused)
{
        cell id, more, value;
        more = arity(op, args, 1, 1);
        arity_next(op, args, more, 0, 1);
        value = cts_pop();
        id = cts_pop();
        if (!symbol_p(id))
                arity_error(ERR_ARITY_SYNTAX, op, args);
        if (undefined_p(value))
                emitq(NIL);
        else
                compile_expression(value, 0);
        emitop(OP_PUSH);
        emitq(id);
        emitop(OP_ERROR);
}

@* Cons Cells. These operators have been written out directly despite
the obvious potential for refactoring into reusable pieces. This
is short-lived until more compiler routines have been written and
the similarity patterns between them become apparent.

Cons cells are defined by the |cons|, |car|, |cdr|, {\it null?} and
{\it pair?} symbols with {\it set-car!} and {\it set-cdr!} providing
for mutation.

@c
void
compile_cons (cell op,
              cell args,
              boolean tail_p __unused)
{ /* pattern 0; |arity==(O,O)| */
        cell ncar, ncdr;
        arity(op, args, 2, 0);
        ncdr = cts_pop();
        ncar = cts_pop();
        compile_expression(ncdr, 0);
        emitop(OP_PUSH);
        compile_expression(ncar, 0);
        emitop(OP_CONS);
}

void
compile_car (cell op,
             cell args,
             boolean tail_p __unused)
{ /* pattern 1; |arity=(OP_PAIR_P)| */
        int comefrom_pair_p;
        arity(op, args, 1, 0);
        compile_expression(cts_pop(), 0);
        emitop(OP_PUSH);
        emitop(OP_PAIR_P);
        emitop(OP_JUMP_TRUE);
        comefrom_pair_p = Here;
        emit(NIL);
        emitq(sym(ERR_UNEXPECTED)); /* TODO */
        emitop(OP_ERROR);
        patch(comefrom_pair_p, int_new(Here));
        emitop(OP_POP);
        emitop(OP_CAR);
}

void
compile_cdr (cell op,
             cell args,
             boolean tail_p __unused)
{
        int comefrom_pair_p;
        arity(op, args, 1, 0);
        compile_expression(cts_pop(), 0);
        emitop(OP_PUSH);
        emitop(OP_PAIR_P);
        emitop(OP_JUMP_TRUE);
        comefrom_pair_p = Here;
        emit(NIL);
        emitq(sym(ERR_UNEXPECTED)); /* TODO */
        emitop(OP_ERROR);
        patch(comefrom_pair_p, int_new(Here));
        emitop(OP_POP);
        emitop(OP_CDR); /* this is the only difference from the above */
}

void
compile_null_p (cell op,
                cell args,
                boolean tail_p __unused)
{ /* pattern 2 = predicate */
        arity(op, args, 1, 0);
        compile_expression(cts_pop(), 0);
        emitop(OP_NULL_P);
}

void
compile_pair_p (cell op,
                cell args,
                boolean tail_p __unused)
{
        arity(op, args, 1, 0);
        compile_expression(cts_pop(), 0);
        emitop(OP_PAIR_P);
}

void
compile_set_car_m (cell op,
                   cell args,
                   boolean tail_p __unused)
{ /* pattern 3 = |arity=(OP_PAIR_P, O)| */
        cell value, object;
        int goto_pair_p;
        arity(op, args, 2, 0);
        value = cts_pop();
        object = cts_pop();
        compile_expression(object, bfalse);
        emitop(OP_PUSH);
        emitop(OP_PAIR_P);
        emitop(OP_JUMP_TRUE);
        goto_pair_p = comefrom();
        emitq(Sym_ERR_UNEXPECTED);
        emitop(OP_ERROR);
        compile_expression(value, bfalse);
        patch(goto_pair_p, int_new(Here));
        emitop(OP_SET_CAR_M);
}

void
compile_set_cdr_m (cell op,
                   cell args,
                   boolean tail_p __unused)
{
        cell value, object;
        int goto_pair_p;
        arity(op, args, 2, 0);
        value = cts_pop();
        object = cts_pop();
        compile_expression(object, bfalse);
        emitop(OP_PUSH);
        emitop(OP_PAIR_P);
        emitop(OP_JUMP_TRUE);
        goto_pair_p = comefrom();
        emitq(Sym_ERR_UNEXPECTED);
        emitop(OP_ERROR);
        compile_expression(value, bfalse);
        patch(goto_pair_p, int_new(Here));
        emitop(OP_SET_CDR_M);
}

@* Environment. The |environment| mutators are the same except for
the flag given to the final opcode.

@c
void
compile_set_m (cell op,
               cell args,
               boolean tail_p __unused)
{ /* pattern 4, |arity = (OP_ENV_P #<> symbol?)| */
        cell env, name, value;
        int goto_env_p;
        arity(op, args, 3, bfalse);
        value = cts_pop();
        name = cts_pop();
        env = cts_pop();
        if (!symbol_p(name))
                error(ERR_ARITY_SYNTAX, NIL);
        compile_expression(env, bfalse);
        emitop(OP_PUSH);
        emitop(OP_ENVIRONMENT_P);
        emitop(OP_JUMP_TRUE);
        goto_env_p = comefrom();
        emitq(Sym_ERR_UNEXPECTED);
        emitop(OP_ERROR);
        patch(goto_env_p, int_new(Here));
        compile_expression(value, bfalse);
        emitop(OP_ENV_MUTATE_M);
        emit(name);
        emit(FALSE);
}

void
compile_define_m (cell op,
                  cell args,
                  boolean tail_p __unused)
{
        cell env, name, value;
        int goto_env_p;
        arity(op, args, 3, bfalse);
        value = cts_pop();
        name = cts_pop();
        env = cts_pop();
        if (!symbol_p(name))
                error(ERR_ARITY_SYNTAX, NIL);
        compile_expression(env, bfalse);
        emitop(OP_PUSH);
        emitop(OP_ENVIRONMENT_P);
        emitop(OP_JUMP_TRUE);
        goto_env_p = comefrom();
        emitq(Sym_ERR_UNEXPECTED);
        emitop(OP_ERROR);
        patch(goto_env_p, int_new(Here));
        compile_expression(value, bfalse);
        emitop(OP_ENV_MUTATE_M);
        emit(name);
        emit(TRUE);
}

void
compile_env_root (cell op,
                  cell args,
                  boolean tail_p __unused)
{ /* pattern 5 = no args */
        arity(op, args, 0, bfalse);
        emitop(OP_ENV_ROOT);
}

void
compile_env_current (cell op,
                     cell args,
                     boolean tail_p __unused)
{
        arity(op, args, 0, bfalse);
        emitop(OP_ENV_QUOTE);
}

@* Quotation \AM\ Quasiquotation. A quoted object is one which is not
evaluated and we have an opcode to do just that, used by many of the
implementations above.

@c
void
compile_quote (cell op __unused,
               cell args,
               boolean tail_p __unused)
{@+
        emitq(args);@+
}

@ Quasiquoting an object is almost, but not quite, entirely different.
The end result is the same however---a run-time object which (almost)
exactly matches the unevaluated source code that it was created from.

A quasiquoted object is converted into its final form by changing any
|unquote| (and {\it unquote-splicing\/}) within it to the result of
evaulating them. This is complicated enough because we're now writing a
compiler within our compiler\footnote{$^1$}{Yo!} but additionally the
quasiquoted object may contain quasiquoted objects, changing the nature
of the inner-|unquote| operators.

@ The compiler for compiling quasiquoted code only calls directly
into the recursive quasicompiler engine (let's call it the
quasicompiler).

@<Function dec...@>=
void compile_quasicompiler (cell, cell, cell, int, boolean);

@ @c
void
compile_quasiquote (cell op,
                    cell args,
                    boolean tail_p __unused)
{ /* pattern Q */
        compile_quasicompiler(op, args, args, 0, bfalse);
}

@ As with any compiler, the first task is to figure out what sort
of expression is being quasicompiled. Atoms are themselves. Otherwise
lists and vectors must be recursively compiled item-by-item, and the
syntactic operators must operate when encountered.

Quasiquoting |vector|s is not supported but I'm not anticipating it
being difficult, just not useful yet.

@c
void
compile_quasicompiler (cell    op,
                       cell    oargs,
                       cell    arg,
                       int     depth,
                       boolean in_list_p)
{
        if (pair_p(arg)) {
                @<Quasiquote a pair/list@>
        } else if (vector_p(arg)) {@+
                error(ERR_UNIMPLEMENTED, NIL);@+
        } else if (syntax_p(arg)) {
                @<Quasiquote syntax@>
        } else {
                emitq(arg);
                if (in_list_p)
                        emitop(OP_CONS);
        }
}

@ Dealing first with the simple case of a list the quasicompiler
reverses the list to find its tail, which may or may not be |NIL|,
and recursively calling |compile_quasicompiler| for every item.

After each item has been quasicompiled it will be combined with the
transformed list being grown on top of the stack.

When quasicompiling the list's tail there is no partial list to
prepend it to so the quasicompiler is entered in atomic mode.
|compile_quasicompiler| can be relied on to handle the tail of a
proper or improper list.

@<Quasiquote a pair/list@>=
cell todo, tail;
tail = NIL;
todo = list_reverse(arg, &tail, NULL);
compile_quasicompiler(op, oargs, tail, depth, bfalse);
for (; !null_p(todo); todo = cdr(todo)) {
        emitop(OP_PUSH); /* Push the list so far */
        compile_quasicompiler(op, oargs, car(todo), depth, btrue);
}
if (in_list_p)
        emitop(OP_CONS);

@ The quote \AM\ unquote syntax is where the quasicompiler starts
to get interesting. |quote|s and |quasiquote|s (and a |dotted| tail)
recurse back into the quasicompiler to emit the transformation of
the quoted object, then re-apply the syntax operator.

|depth| is increased when recursing into a |quasiquote| so that the
compiler knows whether to evaluate an unquote operator.

@<Quasiquote syntax@>=
int d;
if (car(arg) == Sym_SYNTAX_DOTTED@|
        || car(arg) == Sym_SYNTAX_QUOTE@|
        || car(arg) == Sym_SYNTAX_QUASI) {
        d = (car(arg) == Sym_SYNTAX_QUASI) ? 1 : 0;
        compile_quasicompiler(op, oargs, cdr(arg), depth + d, bfalse);
        emitop(OP_SYNTAX);
        emit(car(arg));
        if (in_list_p)
                emitop(OP_CONS);
}

@ |unquote| evaluates the unquoted object. If quasiquote is
quasicompiling an inner quasiquote then the unquoted object isn't
evaluated but compiled at a decreased |depth|. This enables the
correct unquoting-or-not of quasiquoting quasiquoted quotes.

@<Quasiquote syntax@>=
else if (car(arg) == Sym_SYNTAX_UNQUOTE) {
        if (depth > 0) {
                compile_quasicompiler(op, oargs, cdr(arg), depth - 1, bfalse);
                emitop(OP_SYNTAX);
                emit(Sym_SYNTAX_UNQUOTE);
        } else
                compile_expression(cdr(arg), bfalse);
        if (in_list_p)
                emitop(OP_CONS);
}

@ Similarly to |unquote|, {\it unquote-splicing} recurses back into
the quasicompiler at a lower depth when unquoting an inner quasiquote.

@<Quasiquote syntax@>=
else if (car(arg) == Sym_SYNTAX_UNSPLICE) {
        if (depth > 0) {
                compile_quasicompiler(op, oargs, cdr(arg), depth - 1, bfalse);
                emitop(OP_SYNTAX);
                emit(Sym_SYNTAX_UNSPLICE);
                if (in_list_p)
                        emitop(OP_CONS);
        } else { @<Compile unquote-splicing@> }
} else
        error(ERR_UNIMPLEMENTED, NIL);

@*1 Splicing Lists. If not recursing back into the quasicompiler
at a lower depth then we are quasicompiling at the lowest depth and
need to do the work.

When splicing into the tail position of a list we can replace its
|NIL| with the evaluation with minimal further processing. Unfortunately
we don't know until runtime whether we are splicing into the tail
position -- consider constructs like {\tt `(,@@foo ,@@bar)} where {\tt
bar} evaluates to |NIL|.

@<Compile unquote-splicing@>=
int goto_inject_iterate, goto_inject_start, goto_finish;
int goto_list_p, goto_null_p, goto_nnull_p;
if (!in_list_p)
        error(ERR_UNEXPECTED, arg);
emitop(OP_PEEK);
emitop(OP_NULL_P);
emitop(OP_JUMP_TRUE);@+
goto_null_p = comefrom();
emitop(OP_PUSH); /* save |FALSE| */
emitop(OP_JUMP);@+
goto_nnull_p = comefrom();
patch(goto_null_p, int_new(Here));
emitop(OP_SWAP); /* become the tail, save |TRUE| */
patch(goto_nnull_p, int_new(Here));

@ |FALSE| or |TRUE| is now atop the stack indicating whether a new
list is being built otherwise the remainder of the list is left on
the stack. Now we can evaluate and validate the expression.

@<Compile unquote-splicing@>=
compile_expression(cdr(arg), 0);
emitop(OP_PUSH);
emitop(OP_LIST_P);@+
emit(TRUE);@+
emit(FALSE);
emitop(OP_JUMP_TRUE);@+
goto_list_p = comefrom();
emitq(Sym_ERR_UNEXPECTED);
emitop(OP_ERROR);

@ If we have a list we can leave it as-is if we were originally in
the tail position.

@<Compile unquote-splicing@>=
patch(goto_list_p, int_new(Here));
emitop(OP_POP);
emitop(OP_SWAP);
emitop(OP_JUMP_TRUE);@+
goto_finish = comefrom();

@ Splicing a list into the middle of another list is done item-by-item
in reverse. A small efficiency could be gained here by not walking
the list a second time (the first to validate it above) at the cost
of more complex bytecode.

By now the evaluated list to splice in is first on the stack followed
by the partial result.

@<Compile unquote-splicing@>=
emitop(OP_POP);
emitop(OP_LIST_REVERSE);@+
emit(TRUE);@+
emit(FALSE);
@<Walk through the splicing list@>@;

@ @<Walk through the splicing list@>=
emitop(OP_JUMP);@+
goto_inject_start = comefrom();

goto_inject_iterate = Here;
emitop(OP_POP);
emitop(OP_SNOC);
emitop(OP_CYCLE);
emitop(OP_CONS);
emitop(OP_SWAP);

@ If this was the last item (the first of the evaluated list's) or
the evaluation was |NIL| then we're done otherwise we go around
again. This is also where the loop starts to handle the case of
evaluating an empty list.

@<Walk through the splicing list@>=
patch(goto_inject_start, int_new(Here));
emitop(OP_PUSH);
emitop(OP_NULL_P);
emitop(OP_JUMP_FALSE);@+
emit(int_new(goto_inject_iterate));
emitop(OP_POP);
patch(goto_finish, int_new(Here));
emitop(OP_POP);

@** Testing. \LL/ includes (hah!) a comprehensive test suite. It
also includes utilities for self-testing which require a version
of \LL/ with a different entry-point and some internal facilities
exposed. The \CEE/ pre-precessor is used to keep these out of the
\LL/ executable and rename its |main| function in the testing
executable.

@<Function dec...@>=
#ifdef LL_TEST
int test_main (int, char **);
#endif

@ @c
#ifdef LL_TEST
int
main (int    argc,
      char **argv __unused)
{
        vm_init();
        if (argc > 1)
                error(ERR_UNIMPLEMENTED, NIL);
        @<Initialise error...@>@;
        /* Now what? */
        printf("1..1\n");
        printf("not ok 1 There are tests\n");
        return 0;
}
#endif

@** TODO.

@s eval if
@s lambda if
@s vov if
@<List of opcode primitives@>=
        /* Core: */
{ "error", compile_error },
{ "eval", compile_eval },
{ "if", compile_conditional },
{ "lambda", compile_lambda },
{ "vov", compile_vov },
{ "quote", compile_quote },
{ "quasiquote", compile_quasiquote },
        /* Pairs: */
{ "car", compile_car },
{ "cdr", compile_cdr },
{ "cons", compile_cons },
{ "null?", compile_null_p },
{ "pair?", compile_pair_p },
{ "set-car!", compile_set_car_m },
{ "set-cdr!", compile_set_cdr_m },
       /* Mutation: */
{ "current-environment", compile_env_current },
{ "root-environment", compile_env_root },
{ "set!", compile_set_m },
{ "define!", compile_define_m },

@** REPL. The |main| loop is a simple repl.

@c
#ifdef LL_TEST
int
test_main (int    argc,
           char **argv)@;
#else
int
main (int    argc,
      char **argv __unused)
#endif
{
        vm_init();
        if (argc > 1)
                error(ERR_UNIMPLEMENTED, NIL);
        @<Initialise error...@>@;
        while (1) {
                vm_clear();
                printf("> ");
                Acc = read_form();
                if (eof_p(Acc) || Interrupt)
                        break;
                interpret();
                if (!void_p(Acc)) {
                        write_form(Acc, 0);
                        printf("\n");
                }
        }
        if (Interrupt)
                printf("Interrupted");
        return 0;
}

@** Index.
