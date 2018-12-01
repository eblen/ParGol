#include <array>

template<int WORLD_SIZE>
class GOL
{
    using world_grid = std::array<std::array<bool, WORLD_SIZE>, WORLD_SIZE>;

    public:
    GOL()
    {
        for (int x=0; x<WORLD_SIZE; x++)
        {
            for (int y=0; y<WORLD_SIZE; y++)
            {
                old_world[x][y] = (y==0) || (x % 2 == 0) ? false : true;
            }
        }
    }

    void next_gen()
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

        world_grid &tmp = old_world;
        old_world = new_world;
        new_world = tmp;
    }

    void print()
    {
        for (int x=0; x<WORLD_SIZE; x++)
        {
            for (int y=0; y<WORLD_SIZE; y++)
            {
                if (old_world[x][y]) printf("*");
                else printf(" ");
            }
            printf("\n");
        }
    }

    private:
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
    const int WORLD_SIZE = 100;
    const int NUM_GENS   = 100;

    GOL<WORLD_SIZE> world;
    world.print();
    print_sep(WORLD_SIZE);

    for (int gen_num = 0; gen_num < NUM_GENS; gen_num++)
    {
        world.next_gen();
        world.print();
        print_sep(WORLD_SIZE);
    }

    return 0;
}
