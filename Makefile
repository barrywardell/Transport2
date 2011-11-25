CXXFLAGS=-Wall -g -O2 -I./MosquitoTensor/srcs -I/opt/local/include
LDFLAGS=-L/opt/local/lib -lgsl -lcblas -latlas -lm -lhdf5
CXX = g++

.PHONY: all

all: transport

transport: *.cc *.h MosquitoTensor/srcs/*.C MosquitoTensor/srcs/*.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ transport.cc Schwarzschild.cc MosquitoTensor/srcs/Tensor.C MosquitoTensor/srcs/TensorBase.C MosquitoTensor/srcs/IndexedTensor.C MosquitoTensor/srcs/TensorList.C

clean:
	rm -rf transport transport.dSYM
