build: deltser

deltser: deltser.cpp
	c++ -std=c++11 deltser.cpp -o deltser -Ofast -D NDEBUG

