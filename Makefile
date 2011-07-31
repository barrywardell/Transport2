CXXFLAGS=-Wall -g -O2 -I./MosquitoTensor/srcs -I/opt/local/include
LDFLAGS=-L/opt/local/lib -lgsl -lcblas -latlas -lm -lhdf5
CXX = g++

.PHONY: all

all: VanVleck

VanVleck: *.cc *.h MosquitoTensor/srcs/*.C MosquitoTensor/srcs/*.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ VanVleck.cc Schwarzschild.cc MosquitoTensor/srcs/Tensor.C MosquitoTensor/srcs/TensorList.C

clean:
	rm -f VanVleck
