LLTEST_SCRIPTS := \
	eval.log  \
	if.log    \
	pair.log

test: sanity.log $(LLTEST_SCRIPTS)
	if sh -c "prove -V" >/dev/null 2>&1; then \
	        prove --ext=log -e cat -rv .;     \
	else    ! grep ^not *.log;                \
	fi

clean:
	rm -f *.log

# Does it compile?
sanity.log:
	lltest 0 > sanity.log

$(LLTEST_SCRIPTS): sanity.log

# cons car cdr null? pair? set-car! set-cdr!
pair.log:
	lltest p > pair.log

# eval environments
eval.log:
	lltest e > eval.log

# if
if.log:
	lltest i > if.log
