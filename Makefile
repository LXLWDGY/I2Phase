prog: CoordCube.o CubieCube.o I2Phase.o Search.o Util.o
		g++ -g CoordCube.o CubieCube.o I2Phase.o Search.o Util.o -o I2Phase
CoordCube.o : I2Phase.h CoordCube.cpp
		g++ -g -Wall -c CoordCube.cpp
CubieCube.o : I2Phase.h CubieCube.cpp
		g++ -g -Wall -c CubieCube.cpp
I2Phase.o : I2Phase.h I2Phase.cpp
		g++ -g -Wall -c I2Phase.cpp
Search.o : I2Phase.h Search.cpp
		g++ -g -Wall -c Search.cpp
Util.o : I2Phase.h Util.cpp
		g++ -g -Wall -c Util.cpp

