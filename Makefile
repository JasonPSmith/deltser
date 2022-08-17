build: deltser

deltser: deltser.cpp
	g++ -std=c++14 -O3 -pthread deltser.cpp -o deltser

