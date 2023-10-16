read_in : read_in.c
	gcc -o $@ $< -O1 -g
ellens : ellens.c
	gcc -o $@ $< -O1 -g

Lisas : Lisas.c
	gcc -o $@ $< -O1 -g

min : min.c
	gcc -o $@ $< -O1 -g

Ass2OpenMP : Ass2OpenMP.c
<<<<<<< HEAD
	gcc -o $@ $< -O2 -g

Idas : Idas.c
	gcc -o $@ $< -O2 -g
=======
	gcc -o $@ $< -O1 -g -lm

Ass2Open_hemma : Ass2Open_hemma.c
	gcc -o $@ $< -O1 -g -lm
>>>>>>> f69cdec (nästa klar)
