# ParGol (Parallel Game of Life)

### Conway's Game of Life in parallel for testing threading frameworks
<br>
ParGol runs Conway's Game of Life, using MPI, in a straightforward way that is  
easy to parallelize further with threads. Additionally, It requires tasking and  
overlapping of communication and computation for best performance. The goal of  
ParGol is to create a small, easy-to-understand, but non-trivial, application  
that can be used to test the performance of different threading frameworks.
<br>

### Running the Program
**Usage: ./pargol [input file] [no. generations] [print frequency]**

Example input files are included. Each contains a header with two numbers, size  
and number of repetitions. The first gives the size n of the game board given in  
the input file, where *n* is the number of rows and columns (game board must be  
a square). The full board size will be *nr x nr*, where *r* is the number of  
repetitions. The remainder of the file is the game board (see examples). Note  
that only the asterisk indicates a live cell. So any other character can be  
used for the other cells.

The number of generations indicates the number of iterations of the game to run.

The print frequency, *p*, is optional. By default, only the first and last  
generations are printed. If given, the board is printed every p iterations  
(1 prints all iterations, 2 every other iteration, etc.).

The number of MPI ranks must be a perfect square and must evenly divide the  
game board (size of game board must be divisible by the square root of the  
number of ranks).
<br>

### Parallelizing with Threads
For threading, modify the "next\_gen" function, which consists of three parts,  
computing internal cells, communication, and computing external cells. The only  
dependency is that step 2 must be done before step 3. Thus, you can test tasking  
strategies and overlapping of computation and communication. Note that, with  
a sufficiently large board, the internal cell computations dwarf the external  
computations (quadratic vs linear growth), and thus overlapping the first two  
steps becomes key.
