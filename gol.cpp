#include <stdio.h>
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
    GOL(int s) :size(s)
    {
        old_world.assign(size*size,0);
        new_world.assign(size*size,0);
        for (int x=0; x<size; x++)
        {
            for (int y=0; y<size; y++)
            {
                int c = x*size+y;
                old_world[c] = (y==0) || (x % 2 == 0) ? 0 : 1;
            }
        }
    }

    void next_gen()
    {
        for (int c=0; c<size*size; c++) next_gen_cell(c);

        world_grid &tmp = old_world;
        old_world = new_world;
        new_world = tmp;
    }

    void print()
    {
        for (int x=0; x<size; x++)
        {
            for (int y=0; y<size; y++)
            {
                int c = x*size+y;
                if (old_world[c]) printf("*");
                else printf(" ");
            }
            printf("\n");
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

    int size;
    world_grid WORLD1;
    world_grid WORLD2;
    world_grid &old_world = WORLD1;
    world_grid &new_world = WORLD2;
};

void print_sep(int length)
{
    for (int x=0; x<length; x++) printf("=");
    printf("\n");
}

int main()
{
    const int size     = 100;
    const int num_gens = 100;

    GOL world(size);
    world.print();
    print_sep(size);

    for (int gen_num = 0; gen_num < num_gens; gen_num++)
    {
        world.next_gen();
        world.print();
        print_sep(size);
    }

    return 0;
}
