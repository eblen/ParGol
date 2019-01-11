#include <cmath>
#include <fstream>
#include <mpi.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

#include "timestamps.h"

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
    enum task_name {internal, comm, external, NUM_TASKS};

    // Constructor for a boring, default world
    GOL(int s, int ws) :is_valid(false), time_stamps(task_name::NUM_TASKS)
    {
        if (!init(s,ws)) return;

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
        is_valid = true;
    }

    // Constructor for a world from a file
    GOL(const char *file_name, int ws) :is_valid(false), time_stamps(task_name::NUM_TASKS)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // Read and compute size information from file header
        int rep_size, reps;
        std::fstream input(file_name, std::ios_base::in);
        if (!input.is_open())
        {
            if (rank==0) fprintf(stderr, "GOL: Unable to open input file\n");
            return;
        }
        std::string header_string;
        getline(input, header_string);
        std::stringstream header(header_string);
        header >> rep_size >> reps;
        if (rep_size < 1 || reps < 1)
        {
            if (rank==0) fprintf(stderr, "GOL: Invalid header values in input file\n");
            return;
        }
        if (!init(rep_size*reps, ws)) return;

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
                        if (input.eof())
                        {
                            fprintf(stderr, "GOL: Unexpected eol while reading input file\n");
                            return;
                        }
                        char nline = ' ';
                        input.read(&nline, 1); // newline
                        if (nline != '\n')
                        {
                            fprintf(stderr, "GOL: Bad format for input file (missing newline at end of row)\n");
                            return;
                        }
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
                            MPI_Send(&buf[c*(size-2)], size-2, MPI::CHAR, pid, 0, MPI_COMM_WORLD);
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
                MPI_Status status;
                MPI_Recv(&old_world[r*size+1], size-2, MPI::CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                if (status.MPI_TAG == 1) return;
            }
        }
        is_valid = true;
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

    // Should be called before running an instance created from an input file
    bool world_is_valid() const {return is_valid;}
    int local_size() const {return size-2;}
    int world_size() const {return wsize;}
    const TimeStamps& get_time_stamps() const {return time_stamps;}

    void next_gen()
    {
        int s = size;

        // Compute internal cells
        time_stamps.start(task_name::internal);
        for (int x=2; x<s-2; x++)
        {
            for (int y=2; y<s-2; y++)
            {
                next_gen_cell(x*s+y);
            }
        }
        time_stamps.stop(task_name::internal);

        grid_nbrs n = get_grid_nbrs(rank,wsize);
        world_grid &w = old_world;

        // Exchange external cells with other ranks
        time_stamps.start(task_name::comm);

        // Exchange top and bottom rows
        MPI_Sendrecv(&w[s+1],       s-2, MPI::CHAR, n.t,  0,
                     &w[s*(s-1)+1], s-2, MPI::CHAR, n.b,  0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&w[s*(s-2)+1], s-2, MPI::CHAR, n.b,  1,
                     &w[1],         s-2, MPI::CHAR, n.t,  1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Exchange diagonal elements
        MPI_Sendrecv(&w[s+1],         1, MPI::CHAR, n.tl, 2,
                     &w[s*s-1],       1, MPI::CHAR, n.br, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&w[s*(s-2)+s-2], 1, MPI::CHAR, n.br, 3,
                     &w[0],           1, MPI::CHAR, n.tl, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&w[s+s-2],       1, MPI::CHAR, n.tr, 4,
                     &w[s*(s-1)],     1, MPI::CHAR, n.bl, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&w[s*(s-2)+1],   1, MPI::CHAR, n.bl, 5,
                     &w[s-1],         1, MPI::CHAR, n.tr, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       
        // Exchange left and right columns
        MPI_Sendrecv(&w[s+1],         1, GRID_COLUMNS,       n.l,  6,
                     &w[s+s-1],       1, GRID_COLUMNS,       n.r,  6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&w[s+s-2],       1, GRID_COLUMNS,       n.r,  7,
                     &w[s],           1, GRID_COLUMNS,       n.l,  7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        time_stamps.stop(task_name::comm);

        // Compute external cells
        time_stamps.start(task_name::external);

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
        time_stamps.stop(task_name::external);

        // Commit changes by swapping worlds
        old_world.swap(new_world);
    }

    void print(int max_size) const
    {
        if (rank==0)
        {
            std::vector<char> buf(size-2);

            // Rows of MPI ranks
            int rows_printed = 0;
            for (int wrow=0; wrow < wsize; wrow++)
            {
                // Rows of cells
                for (int row=1; row < size-1; row++)
                {
                    int cells_printed = 0;
                    // MPI ranks containing cells on current row
                    for (int pid=wrow*wsize; pid < wrow*wsize+wsize; pid++)
                    {
                        if (pid==0)
                        {
                            memcpy(&buf[0], &new_world[row*size+1], size-2);
                        }
                        else
                        {
                            MPI_Recv(&buf[0], size-2, MPI::CHAR, pid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        }

                        // Individual cells
                        int cells_to_print = std::min(size-2, max_size - cells_printed);
                        for (int x=0; x < cells_to_print; x++)
                        {
                            if (buf[x]) printf("*");
                            else printf(" ");
                        }
                        cells_printed += cells_to_print;
                        if (cells_printed >= max_size) break;
                    }
                    printf("\n");
                    rows_printed++;
                    if (rows_printed >= max_size) return;
                }
            }
        }

        else
        {
            int start_col = (rank % wsize) * (size-2);
            if (start_col >= max_size) return;
            int start_row = (rank / wsize);

            int num_rows = std::min(max_size - start_row*(size-2), size-2);
            for (int row=1; row < num_rows+1; row++)
            {
                MPI_Send(&new_world[row*size+1], size-2, MPI::CHAR, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

    // Also print communication cells (primarily for debugging)
    void printall() const
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
                            MPI_Recv(&buf[0], size, MPI::CHAR, pid, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
                MPI_Send(&new_world[row*size], size, MPI::CHAR, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

    private:
    bool init(int s, int ws)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (s % ws != 0)
        {
            if (rank==0) fprintf(stderr, "GOL: Population cannot be evenly divided among MPI ranks\n");
            return false;
        }
        // +2 for extra outer rows and columns for communicating with neighbors
        size =s / ws + 2;
        wsize = ws;

        MPI_Type_vector(size-2, 1, size, MPI::CHAR, &GRID_COLUMNS);
        MPI_Type_commit(&GRID_COLUMNS);

        old_world.assign(size*size,0);
        new_world.assign(size*size,0);
        return true;
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

    bool is_valid; // whether world was able to be created (useful for worlds created from an input file)
    int size;  // Local population size
    int rank;  // MPI rank
    int wsize; // Global size ( sqrt(no. of MPI ranks) )
    MPI_Datatype GRID_COLUMNS; // Used for communicating columns
    TimeStamps time_stamps;
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
        if (rank==0) fprintf(stderr, "Usage: %s <input file> <no. generations> <print size> <print frequency>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    // Define GOL constants
    const int wsize      = sqrt(nranks);
    const int num_gens   = atoi(argv[2]);
    int print_max_size = 100;
    if (argc > 3) print_max_size = atoi(argv[3]);
    int print_freq       = num_gens-1; // Print first and last generations only by default
    if (argc > 4) print_freq = atoi(argv[4]);
    if (wsize*wsize != nranks)
    {
        if (rank==0) fprintf(stderr, "Error: No. MPI ranks must be a perfect square\n");
        MPI_Finalize();
        exit(1);
    }

    GOL world(argv[1],wsize);
    if (rank==0 && !world.world_is_valid()) MPI_Abort(MPI_COMM_WORLD,1);
    if (world.world_is_valid())
    {
        for (int gen_num = 0; gen_num < num_gens; gen_num++)
        {
            world.next_gen();
            if (gen_num % print_freq == 0)
            {
                world.print(print_max_size);
                if (rank==0) print_sep(std::min(print_max_size, world.local_size() * world.world_size()));
            }
        }
    }

    // Print performance statistics
    if (rank==0)
    {
        const TimeStamps& ts = world.get_time_stamps();
        printf("Time for internal compute: %1.2f ms\n", ts.get_total_time(GOL::task_name::internal) / 1000.0);
        printf("Time for communication:    %1.2f ms\n", ts.get_total_time(GOL::task_name::comm)     / 1000.0);
        printf("Time for external compute: %1.2f ms\n", ts.get_total_time(GOL::task_name::external) / 1000.0);
        printf("Internal and comm overlap: %1.2f ms\n", ts.get_overlap(GOL::task_name::internal, GOL::task_name::comm) / 1000.0);
    }

    // Exit
    MPI_Finalize();
    return 0;
}
