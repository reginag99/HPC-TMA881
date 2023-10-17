read_in : read_in.c
	gcc -o $@ $< -O1 -g
ellens : ellens.c
	gcc -o $@ $< -O1 -g

Lisas : Lisas.c
	gcc -o $@ $< -O1 -g

min : min.c
	gcc -o $@ $< -O1 -g

Ass2OpenMP : Ass2OpenMP.c
	gcc -o $@ $< -O2 -g

Idas : Idas.c
	gcc -o $@ $< -O2 -g

Ass2Open_hemma : Ass2Open_hemma.c
	gcc -o $@ $< -O1 -g -lm


ass3_version2: ass3_version2.c
	gcc -o $@ $< -O1 -g -lm

