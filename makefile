all: main
	./main

main: bresenham_line.o color.o curves.o lighting.o matrix_transform.o matrix.o shapes.o vector_math.o runscript.o
	g++ -o main bresenham_line.o color.o curves.o lighting.o matrix_transform.o matrix.o shapes.o vector_math.o runscript.o

bresenham_line.o: bresenham_line.hpp bresenham_line.cpp
	g++ -c bresenham_line.cpp

color.o: color.hpp color.cpp
	g++ -c color.cpp

curves.o: curves.hpp curves.cpp
	g++ -c curves.cpp

lighting.o: lighting.hpp lighting.cpp
	g++ -c lighting.cpp

matrix_transform.o: matrix_transform.hpp matrix_transform.cpp
	g++ -c matrix_transform.cpp

matrix.o: matrix.hpp matrix.cpp
	g++ -c matrix.cpp

shapes.o: shapes.hpp shapes.cpp
	g++ -c shapes.cpp

vector_math.o: vector_math.hpp vector_math.cpp
	g++ -c vector_math.cpp

runscript.o: runscript.cpp
	g++ -c runscript.cpp

clean:
	rm -rf main *.o