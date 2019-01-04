#include <cmath>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
    using world_grid = std::vector<unsigned char>;

    public:
    // +2 for extra outer rows and columns for communicating with neighbors
    GOL(int s, int ws) :size(s+2), wsize(ws)
    {
        // Gather and create MPI data
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Type_vector(size-2, 1, size, MPI::UNSIGNED_CHAR, &GRID_COLUMNS);
        MPI_Type_commit(&GRID_COLUMNS);

        // Build initial population
        bool isAlive = false;
        old_world.assign(size*size,0);
        new_world.assign(size*size,0);
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
            std::vector<unsigned char> buf(size-2);

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
            std::vector<unsigned char> buf(size-2);

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

    // Define GOL constants
    const int size       = 2;
    const int num_gens   = 100;
    const int wsize      = sqrt(nranks);
    if (wsize*wsize != nranks)
    {
        fprintf(stderr, "Error: No. MPI ranks must be a perfect square\n");
        MPI_Finalize();
        exit(1);
    }

    // Create and print initial world
    GOL world(size,wsize);
    if (rank==0) print_sep(size*wsize);

    // Run main loop
    for (int gen_num = 0; gen_num < num_gens; gen_num++)
    {
        world.next_gen();
        world.print();
        if (rank==0) print_sep(size*wsize);
    }

    // Exit
    MPI_Finalize();
    return 0;
}
