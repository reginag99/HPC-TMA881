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

Assignment3 : Assignment3.c
	gcc -o $@ $< -O2 -g
