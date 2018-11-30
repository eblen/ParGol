#include <array>

const int WORLD_SIZE = 100;
const int NUM_GENS   = 100;

using world_grid = std::array<std::array<bool, WORLD_SIZE>, WORLD_SIZE>;

// Assumes input world is empty
void init_world(world_grid &w)
{
    w[49][48] = true;
    w[49][49] = true;
    w[49][50] = true;
    w[49][51] = true;
}

void next_gen(const world_grid &old_world, world_grid &new_world)
{
    for (int x=0; x<WORLD_SIZE; x++)
    {
        for (int y=0; y<WORLD_SIZE; y++)
        {
            // Count neighbors
            int num_nbrs = 0;
            int x_minus = x == 0 ? WORLD_SIZE-1 : x-1;
            int y_minus = y == 0 ? WORLD_SIZE-1 : y-1;
            int x_plus  = (x+1) % WORLD_SIZE;
            int y_plus  = (y+1) % WORLD_SIZE;
            if (old_world[x_minus][y_minus]) num_nbrs++;
            if (old_world[x]      [y_minus]) num_nbrs++;
            if (old_world[x_plus] [y_minus]) num_nbrs++;
            if (old_world[x_minus][y])       num_nbrs++;
            if (old_world[x_plus] [y])       num_nbrs++;
            if (old_world[x_minus][y_plus])  num_nbrs++;
            if (old_world[x]      [y_plus])  num_nbrs++;
            if (old_world[x_plus] [y_plus])  num_nbrs++;

            // Apply rules
            switch(num_nbrs) {
                case 2:
                    new_world[x][y] = old_world[x][y];
                    break;
                case 3:
                    new_world[x][y] = true;
                    break;
                default:
                    new_world[x][y] = false;

            }
        }
    }
}

void print_world(const world_grid &w)
{
    for (int x=0; x<WORLD_SIZE; x++)
    {
        for (int y=0; y<WORLD_SIZE; y++)
        {
            if (w[x][y]) printf("*");
            else printf(" ");
        }
        printf("\n");
    }
}

void print_sep()
{
    for (int x=0; x<WORLD_SIZE; x++) printf("=");
    printf("\n");
}

int main()
{
    static world_grid WORLD1;
    static world_grid WORLD2;
    world_grid &old_world = WORLD1;
    world_grid &new_world = WORLD2;

    init_world(old_world);
    print_world(old_world);
    print_sep();
    for (int gen_num = 0; gen_num < NUM_GENS; gen_num++)
    {
        next_gen(old_world, new_world);
        print_world(new_world);
        print_sep();

        if (gen_num % 2 != 0) {
            old_world = WORLD1;
            new_world = WORLD2;
        }
        else {
            old_world = WORLD2;
            new_world = WORLD1;
        }
    }

    return 0;
}
