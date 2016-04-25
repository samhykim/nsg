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

The nsg binary is currently compiled on a Linux machine.

###Usage: NSG [-flags] [inputFilename] [outputFilename]

valid flags are:

   t : run with t threads (default is max number of cores on your machine )
       must be followed by a number 

   x : don't create output file

   e : do read extension step

example: './nsg -xt 1 bact1.reads' will run in serial, disable output, and use reads from bact1.reads.

###File Input

The input to the nsg algorithm should be a list of reads, one per line. Alternatively, fasta files with read headers (e.g., >read5) can also be used. The .revreads sample files contain 8000 simulated reads of length 3000, followed by their reverse complements for E. coli K12, R. sphaeroides, and S. aureus.



### File Output Format
The .revxzgraph files contain sample outputs of the nsg algorithm. The format of each line are as follows:

	[read ID] [greedy ID] [greedy overlap] [not-so-greedy ID] [not-so-greedy overlap]

The [greedy overlap] is the length of the overlap the prefix of read with [greedy ID] has with the suffix of read with [read ID].
The same applies to the not-so-greedy read. 

