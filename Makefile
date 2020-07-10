sources = common.cpp
headers = common.h params.h

PEP: PEP.cpp $(sources) $(headers)
	g++ -o PEP  PEP.cpp $(sources) -std=c++11 -DPEP -O3 -g -mavx2 -lpthread

LEP: LEP.cpp $(sources) $(headers)
	g++ -o LEP  LEP.cpp $(sources) -std=c++11 -DLEP -O3 -g -mavx2 -lpthread
