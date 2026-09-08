all: rs-test

rs-test:
	gcc test/* rs/src/* -O3 -o rs-test -I rs/inc
clean:
	-rm -f ./rs-test.exe
	-rm -f ./rs-test
