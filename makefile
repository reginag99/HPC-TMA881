read_in : read_in.c
	gcc -o $@ $< -O1 -g
ellens : ellens.c
	gcc -o $@ $< -O1 -g

Lisas : Lisas.c
	gcc -o $@ $< -O1 -g

min : min.c
	gcc -o $@ $< -O1 -g
