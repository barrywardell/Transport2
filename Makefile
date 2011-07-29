CXXFLAGS=-Wall -g -I./MosquitoTensor/srcs -I/opt/local/include
LDFLAGS=-L/opt/local/lib -lgsl -lcblas -latlas -lm
CXX = g++

.PHONY: all

all: VanVleck

VanVleck: VanVleck.cc MosquitoTensor/srcs/Tensor.C Schwarzschild.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ VanVleck.cc Schwarzschild.cc MosquitoTensor/srcs/Tensor.C

clean:
	rm -f Test
