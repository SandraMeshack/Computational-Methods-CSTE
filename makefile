run: compile
	./main

compile:
	g++ *.cpp -o main -std=c++11
docs:
	doxygen ./Doxyfile
