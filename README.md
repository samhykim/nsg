# NSG string graph construction algorithm

###The NSG algorithm consists of two steps:

1. Constructing a string graph with a max of two outgoing and two ingoing edges for each read. The greedy and 'not-so-greedy' read is found. The greedy of read r is defined by the read with the largest overlap with r. The 'not-so-greedy' is defined by the read with the next largest overlap which does not align with the greedy read. 

2. Identifying X and Z structures. We transform X-like structures (nodes in the graph with in-degree = 2 and out-degree = 2) and remove Z structures that satisfy the followering condition:

Assuming edges (u1, v1), (u2, v1), (u2, v2),

	a. w(u1, v1) >= w(u2, v1) 

	b. w(u2, v2) >= w(u2, v1)

where one of the conditions is strictly less. 

##Creating the binary: 

The main file is located in nsg.cpp

The binary can be created by :  g++ -O3 nsg.cpp -o nsg -std=c++11

###Usage: NSG [-flags] [inputFilename] [outputFilename]

valid flags are:

   t : run with t threads (default is max number of cores on your machine )
       must be followed by a number 

   x : don't create output file

   e : do read extension step

example: './nsg -xt 1 bact1.reads' will run in serial, disable output, and use reads from bact1.reads.

