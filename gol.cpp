#include <cmath>
#include <fstream>
#include <mpi.h>
#include <sstream>
#include <stdio.h>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

// Generic function for computing neighbor indices of a cell
struct grid_nbrs {int tl,t,tr,l,r,bl,b,br;};
grid_nbrs get_grid_nbrs(int idx, int size)
{
    grid_nbrs n;

    int row = idx / size;
    int col = idx % size;
    int row_above = row == 0 ? size-1 : row - 1;
    int col_left  = col == 0 ? size-1 : col - 1;
    int row_below = row == size-1 ? 0 : row + 1;
    int col_right = col == size-1 ? 0 : col + 1;

    n.tl = row_above * size + col_left;
    n.t  = row_above * size + col;
    n.tr = row_above * size + col_right;
    n.l  = row       * size + col_left;
    n.r  = row       * size + col_right;
    n.bl = row_below * size + col_left;
    n.b  = row_below * size + col;
    n.br = row_below * size + col_right;

    return n;
}

class GOL
{
    using world_grid = std::vector<char>;

    public:
    // Constructor for a boring, default world
    GOL(int s, int ws)
    {
        init(s,ws);

        // Populate world with non-adjacent columns that never change
        bool isAlive = false;
        for (int x=1; x<size-1; x++)
        {
            for (int y=1; y<size-1; y++)
            {
                int c = x*size+y;
                old_world[c] = isAlive;
                isAlive = !isAlive;
            }
        }
    }

    // Constructor for a world from a file
    GOL(const char *file_name, int ws)
    {
        // Read and compute size information from file header
        int rep_size, reps;
        std::fstream input(file_name, std::ios_base::in);
        std::string header_string;
        getline(input, header_string);
        std::stringstream header(header_string);
        header >> rep_size >> reps;

        init(rep_size*reps, ws);

        // Populate world from rest of file

        // Rank 0 reads the data and sends each chunk to the appropriate other rank
        if (rank==0)
        {
            std::vector<std::vector<char>> input_buffers(rep_size, std::vector<char>(rep_size*reps));
            int buf_num = 0;
            // Rows of MPI ranks
            for (int wrow=0; wrow < wsize; wrow++)
            {
                // Rows of cells
                for (int row=1; row < size-1; row++)
                {
                    // After input exhausted, repeat input from buffers
                    auto &buf = input_buffers[buf_num % rep_size];
                    if (buf_num < rep_size)
                    {
                        input.read(buf.data(), rep_size);
                        input.ignore(1); // newline
                        // Duplicate the row horizontally
                        for (int c=1; c<reps; c++)
                        {
                            memcpy(&buf[c*rep_size], buf.data(), rep_size*sizeof(char));
                        }
                        chars_to_binary(buf);
                    }
                    buf_num++;

                    // Divide cells among MPI ranks
                    for (int c=0; c<wsize; c++)
                    {
                        int pid=wrow*wsize+c;
                        if (pid==0)
                        {
                            // +1 to avoid external rows and columns meant for communication
                            memcpy(&old_world[row*size+1], &buf[c*(size-2)], (size-2)*sizeof(char));
                        }
                        else
                        {
                            MPI_Send(&buf[c*(size-2)], size-2, MPI::UNSIGNED_CHAR, pid, 0, MPI_COMM_WORLD);
                        }
                    }
                }
            }
            input.close();
        }
        else
        {
            input.close();
            // Receive our chunks (rows) of data from file
            // +1 to avoid external rows and columns meant for communication
            for (int r=1; r<size-1; r++)
            {
                MPI_Recv(&old_world[r*size+1], size-2, MPI::UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

    /*
    ~GOL()
    {
        // This should not be done here because it most likely will be called
        // after MPI_Finalize is called by the main program.
        // TODO: Figure out alternate strategy for freeing the memory, because
        // now we have a memory leak...
        MPI_Type_free(&GRID_COLUMNS);
    }
    */

    int local_size() {return size-2;}
    int world_size() {return wsize;}

    void next_gen()
    {
        int s = size;

        // Compute internal cells
        for (int x=2; x<s-2; x++)
        {
            for (int y=2; y<s-2; y++)
            {
                next_gen_cell(x*s+y);
            }
        }

        grid_nbrs n = get_grid_nbrs(rank,wsize);
        world_grid &w = old_world;

        // Exchange external cells with other ranks

        // Exchange top and bottom rows
        MPI_Sendrecv(&w[s+1],       s-2, MPI::UNSIGNED_CHAR, n.t,  0,
                     &w[s*(s-1)+1], s-2, MPI::UNSIGNED_CHAR, n.b,  0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&w[s*(s-2)+1], s-2, MPI::UNSIGNED_CHAR, n.b,  1,
                     &w[1],         s-2, MPI::UNSIGNED_CHAR, n.t,  1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Exchange diagonal elements
        MPI_Sendrecv(&w[s+1],         1, MPI::UNSIGNED_CHAR, n.tl, 2,
                     &w[s*s-1],       1, MPI::UNSIGNED_CHAR, n.br, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&w[s*(s-2)+s-2], 1, MPI::UNSIGNED_CHAR, n.br, 3,
                     &w[0],           1, MPI::UNSIGNED_CHAR, n.tl, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&w[s+s-2],       1, MPI::UNSIGNED_CHAR, n.tr, 4,
                     &w[s*(s-1)],     1, MPI::UNSIGNED_CHAR, n.bl, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&w[s*(s-2)+1],   1, MPI::UNSIGNED_CHAR, n.bl, 5,
                     &w[s-1],         1, MPI::UNSIGNED_CHAR, n.tr, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       
        // Exchange left and right columns
        MPI_Sendrecv(&w[s+1],         1, GRID_COLUMNS,       n.l,  6,
                     &w[s+s-1],       1, GRID_COLUMNS,       n.r,  6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&w[s+s-2],       1, GRID_COLUMNS,       n.r,  7,
                     &w[s],           1, GRID_COLUMNS,       n.l,  7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Compute external cells

        // Compute top and bottom rows
        for (int y=1; y<s-1; y++) {
            next_gen_cell(s+y);
            next_gen_cell(s*(s-2)+y);
        }

        // Compute left and right columns
        for (int x=1; x<s-1; x++) {
            next_gen_cell(s*x+1);
            next_gen_cell(s*x+s-2);
        }

        // Commit changes by swapping worlds
        old_world.swap(new_world);
    }

    void print()
    {
        if (rank==0)
        {
            std::vector<char> buf(size-2);

            // Rows of MPI ranks
            for (int wrow=0; wrow < wsize; wrow++)
            {
                // Rows of cells
                for (int row=1; row < size-1; row++)
                {
                    // MPI ranks containing cells on current row
                    for (int pid=wrow*wsize; pid < wrow*wsize+wsize; pid++)
                    {
                        if (pid==0)
                        {
                            memcpy(&buf[0], &new_world[row*size+1], size-2);
                        }
                        else
                        {
                            MPI_Recv(&buf[0], size-2, MPI::UNSIGNED_CHAR, pid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }

                        // Individual cells
                        for (int x=0; x < size-2; x++)
                        {
                            if (buf[x]) printf("*");
                            else printf(" ");
                        }
                    }
                    printf("\n");
                }
            }
        }

        else
        {
            for (int row=1; row < size-1; row++)
            {
                MPI_Send(&new_world[row*size+1], size-2, MPI::UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

    // Also print communication cells (primarily for debugging)
    void printall()
    {
        if (rank==0)
        {
            std::vector<char> buf(size-2);

            // Rows of MPI ranks
            for (int wrow=0; wrow < wsize; wrow++)
            {
                // Rows of cells
                for (int row=0; row < size; row++)
                {
                    // MPI ranks containing cells on current row
                    for (int pid=wrow*wsize; pid < wrow*wsize+wsize; pid++)
                    {
                        if (pid==0)
                        {
                            memcpy(&buf[0], &new_world[row*size], size);
                        }
                        else
                        {
                            MPI_Recv(&buf[0], size, MPI::UNSIGNED_CHAR, pid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }

                        // Individual cells
                        for (int x=0; x < size; x++)
                        {
                            if (buf[x]) printf("*");
                            else printf(" ");
                        }
                    }
                    printf("\n");
                }
            }
        }

        else
        {
            for (int row=0; row < size; row++)
            {
                MPI_Send(&new_world[row*size], size, MPI::UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

    private:
    void init(int s, int ws)
    {
        if (s % ws != 0) throw std::invalid_argument("GOL: Population cannot be evenly divided among MPI ranks");
        // +2 for extra outer rows and columns for communicating with neighbors
        size =s / ws + 2;
        wsize = ws;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Type_vector(size-2, 1, size, MPI::UNSIGNED_CHAR, &GRID_COLUMNS);
        MPI_Type_commit(&GRID_COLUMNS);

        old_world.assign(size*size,0);
        new_world.assign(size*size,0);
    }
    // Convert input characters to binary
    void chars_to_binary(std::vector<char>& v)
    {
        for (char& c : v)
        {
            if (c=='*') c=1;
            else c=0;
        }
    }
    void next_gen_cell(int c)
    {
        // Count neighbors
        grid_nbrs n = get_grid_nbrs(c,size);
        int num_nbrs = 0;
        if (old_world[n.tl]) num_nbrs++;
        if (old_world[n.t])  num_nbrs++;
        if (old_world[n.tr]) num_nbrs++;
        if (old_world[n.l])  num_nbrs++;
        if (old_world[n.r])  num_nbrs++;
        if (old_world[n.bl]) num_nbrs++;
        if (old_world[n.b])  num_nbrs++;
        if (old_world[n.br]) num_nbrs++;

        // Apply rules
        switch(num_nbrs)
        {
            case 2:
                new_world[c] = old_world[c];
                break;
            case 3:
                new_world[c] = 1;
                break;
            default:
                new_world[c] = 0;

        }
    }

    int size;  // Local population size
    int rank;  // MPI rank
    int wsize; // Global size ( sqrt(no. of MPI ranks) )
    MPI_Datatype GRID_COLUMNS; // Used for communicating columns
    world_grid old_world;
    world_grid new_world;
};

void print_sep(int length)
{
    for (int x=0; x<length; x++) printf("=");
    printf("\n");
}

int main(int argc, char **argv)
{
    int rank;
    int nranks;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s <input file> <no. generations> <print frequency>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    // Define GOL constants
    const int wsize      = sqrt(nranks);
    const int num_gens   = atoi(argv[2]);
    int print_freq       = num_gens-1; // Print first and last generations only by default
    if (argc > 3) print_freq = atoi(argv[3]);
    if (wsize*wsize != nranks)
    {
        fprintf(stderr, "Error: No. MPI ranks must be a perfect square\n");
        MPI_Finalize();
        exit(1);
    }

    GOL world(argv[1],wsize);
    for (int gen_num = 0; gen_num < num_gens; gen_num++)
    {
        world.next_gen();
        if (gen_num % print_freq == 0)
        {
            world.print();
            if (rank==0) print_sep(world.local_size() * world.world_size());
        }
    }

    // Exit
    MPI_Finalize();
    return 0;
}
