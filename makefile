all: main
main: main.c functions.o solvers.o
		gcc -g -std=c99 -lSDL -lpthread -Wall -o main main.c functions.o solvers.o
functions.o: functions.c functions.h
			gcc -g -std=c99 -c functions.c
solvers.o: solvers.c functions.h
					gcc -g -std=c99 -c solvers.c
clean:
		rm -f *.o
