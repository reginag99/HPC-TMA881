

all: newton

newton : function_test.c
	gcc -o newton function_test.c -O2 -g -lpthread -lm

run: newton
	./newton -t1 -l1000 1
	./newton -t1 -l1000 2
	./newton -t1 -l1000 5
	./newton -t1 -l1000 7

	./newton -t1 -l1000 5
	./newton -t2 -l1000 5
	./newton -t3 -l1000 5
	./newton -t4 -l1000 5

	./newton -t10 -l1000 7
	./newton -t10 -l2000 7
	./newton -t10 -l3000 7

	make clean

clean:
	rm -f newton

.PHONY: run clean
