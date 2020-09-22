host_parasite_dispersal : host_parasite_dispersal.cpp
	g++ -Wall -g -o host_parasite_dispersal host_parasite_dispersal.cpp -lgsl -lgslcblas -lm
	
all:	host_parasite_dispersal

clean:
	rm -f host_parasite_dispersal
