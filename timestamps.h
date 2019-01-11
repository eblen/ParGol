#include <chrono>
#include <vector>

using namespace std::chrono;

class TimeStamps
{
    using time_stamp = std::pair<long, long>;

    static long get_time()
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(steady_clock::now().time_since_epoch()).count();
    }

    public:
    TimeStamps(int num_tasks)
    {
        time_stamps.resize(num_tasks);
        total_times.assign(num_tasks, 0);
    }

    void start(int task_num, long start_time = get_time())
    {
        if (time_stamps[task_num].size() > 0 && time_stamps[task_num].back().second == 0)
        {
            fprintf(stderr, "TimeStamps: Error - task started twice\n");
            return;
        }
        time_stamps[task_num].emplace_back(time_stamp(start_time, 0));
    }

    void stop(int task_num, long stop_time = get_time())
    {
        if (time_stamps[task_num].size() == 0 || time_stamps[task_num].back().second != 0)
        {
            fprintf(stderr, "TimeStamps: Error - task stopped before started\n");
            return;
        }
        time_stamps[task_num].back().second = stop_time;
        total_times[task_num] += (stop_time - time_stamps[task_num].back().first);
    }

    time_stamp get(int task_num, int pos) const {return time_stamps[task_num][pos];}
    long get_start(int task_num, int pos) const {return time_stamps[task_num][pos].first;}
    long get_stop( int task_num, int pos) const {return time_stamps[task_num][pos].second;}

    long get_total_time(int task_num) const {return total_times[task_num];}

    long get_overlap(int task1, int task2) const
    {
        if (time_stamps[task1].size() == 0 || time_stamps[task2].size() == 0) return 0;

        // Small helper struct for iterating through time stamps
        struct time_bound_ptr
        {
            int pos;
            bool isEnd;
            int size;
            time_bound_ptr(int s) :pos(0), isEnd(false), size(s) {};
            bool next()
            {
                if (pos==(size-1) && isEnd == true) return false;
                if (isEnd) pos++;
                isEnd = !isEnd;
                return true;
            }
            bool prev()
            {
                if (pos==0 && isEnd == false) return false;
                if (!isEnd) pos--;
                isEnd = !isEnd;
                return true;
            }
        };

        time_bound_ptr tb_ptr1(time_stamps[task1].size());
        time_bound_ptr tb_ptr2(time_stamps[task2].size());
        long overlap = 0;
        while (true)
        {
            long prev_time_1 = tb_ptr1.isEnd ? time_stamps[task1][tb_ptr1.pos].second : time_stamps[task1][tb_ptr1.pos].first;
            long prev_time_2 = tb_ptr2.isEnd ? time_stamps[task2][tb_ptr2.pos].second : time_stamps[task2][tb_ptr2.pos].first;
            if (!tb_ptr1.next()) break;
            if (!tb_ptr2.next()) break;
            long next_time_1 = tb_ptr1.isEnd ? time_stamps[task1][tb_ptr1.pos].second : time_stamps[task1][tb_ptr1.pos].first;
            long next_time_2 = tb_ptr2.isEnd ? time_stamps[task2][tb_ptr2.pos].second : time_stamps[task2][tb_ptr2.pos].first;
            if (next_time_1 <= next_time_2)
            {
                tb_ptr2.prev();
                if (tb_ptr1.isEnd && !tb_ptr2.isEnd) overlap += (next_time_1 - std::max(prev_time_1, prev_time_2));
            }
            else
            {
                tb_ptr1.prev();
                if (tb_ptr2.isEnd && !tb_ptr1.isEnd) overlap += (next_time_2 - std::max(prev_time_1, prev_time_2));
            }
        }
        return overlap;
    }

    private:
    std::vector<std::vector<time_stamp>> time_stamps;
    std::vector<long> total_times;
};
