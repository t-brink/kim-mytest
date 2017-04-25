#  Copyright (c) 2013,2014,2017 Tobias Brink
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# First set the location of the kim-api-build-config util.
KIMCFG=kim-api-build-config

CXX=$(shell $(KIMCFG) --cxx)
LD=$(shell $(KIMCFG) --ld)

KIMINCLUDES=$(shell $(KIMCFG) --includes)
INCLUDES="-I./nlopt/include/"

CXXFLAGS=$(shell $(KIMCFG) --cxxflags)
LDFLAGS=$(shell $(KIMCFG) --ldflags)
LDLIBS=$(shell $(KIMCFG) --ldlibs)

PWD=$(shell pwd)

all: mytest

clean:
	rm -f mytest *.o

mytest: core.o elastic.o dsl.o mytest.o numer_forces.o poisson_disk_sampling.o
	$(LD) $(LDFLAGS) -L$(PWD)/nlopt/lib/ -Wl,-rpath $(PWD)/nlopt/lib/ -std=c++11 -o mytest mytest.o core.o elastic.o dsl.o numer_forces.o poisson_disk_sampling.o $(LDLIBS) -lnlopt

core.o: core.cpp core.hpp
	$(CXX) $(CXXFLAGS) -std=c++11 $(INCLUDES) $(KIMINCLUDES) -c -o core.o -c core.cpp

elastic.o: elastic.cpp elastic.hpp
	$(CXX) $(CXXFLAGS) -std=c++11 $(INCLUDES) $(KIMINCLUDES) -c -o elastic.o -c elastic.cpp

dsl.o: dsl.cpp dsl.hpp
	$(CXX) $(CXXFLAGS) -std=c++11 $(INCLUDES) $(KIMINCLUDES) -c -o dsl.o -c dsl.cpp

mytest.o: mytest.cpp
	$(CXX) $(CXXFLAGS) -std=c++11 $(INCLUDES) $(KIMINCLUDES) -c -o mytest.o mytest.cpp

numer_forces.o: numer_forces.cpp
	$(CXX) $(CXXFLAGS) -std=c++11 $(INCLUDES) $(KIMINCLUDES) -c -o numer_forces.o numer_forces.cpp

poisson_disk_sampling.o: poisson_disk_sampling.cpp
	$(CXX) $(CXXFLAGS) -std=c++11 $(INCLUDES) $(KIMINCLUDES) -c -o poisson_disk_sampling.o poisson_disk_sampling.cpp
