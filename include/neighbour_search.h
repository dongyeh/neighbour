#ifndef NEIGHBOUR_SEARCH_H
#define NEIGHBOUR_SEARCH_H

#include <cmath>

class NeighbourSearch
{
public:
    NeighbourSearch(int num_thread);
    ~NeighbourSearch();

    void Init(int maxn);
    void Search(const char *infile, const char *outfile, bool enable_timer);
    long int MemStats() const { return mem_size_ / (1024 * 1024); }

    int npair();
private:
    typedef int IntPair[2];

    enum DisMode {
        HORIZONTAL = 0,
        VERTICAL
    };

    struct Point {
        double x, y, rx, ry;
        int from;
    };

    struct Event {
        int x, type, bottom_point, top_point, which;
    };

    struct Node {
        int x, y, lson, rson, sum, last, sec_last, to_leaf;
    };

    struct ThreadData {
        int *from, *rect_centre_x, *rect_left, *rect_right, *rect_bottom,
        *rect_centre_y, *rect_top, *rect_btm_time, *insert_time, *fir, *next,
        *fir_s, *fir_sp, *next_s, *next_sp;
        double *x, *y, *rx, *ry, *ydiff;
        bool *real;
        int npoint, npair;
        IntPair *pair;
        Event *tmp_event;
        Event **final_event;
        Node *tree;
        int *backtracker, *get_queue;

        double sxtime, gdxtime, bstime, dvtime, dytime, cetime, setime, istime, bttime, sctime;
    };

    const int kNumThread;
    const int kPairFactor;
    const int kTreeIncrement;
    const double kNegativeInfinity;

    void Reset(const int maxn);
    void GetNeighbourPair(const int npnt);
    void Process(ThreadData &td);
    void ParallelDiscretize(int npnt, Point **original, Point **final, int &num_diff);
    void ParallelDivide(int npnt, int num_diff, Point**final);
    void PSRS(int npnt, Point **original, Point **final);
    void MergeSort(int *local_pivot, int segnum, int size, Point **final);
    void MergePass(int *local_pivot, int segnum, Point **a, Point **b);
    void Merge(int low, int mid, int high, Point **a, Point **b);
    int SubFind(Point **point, int min, int max, double target);
    void Discretize(int numPoints, double diffCoord[],
                    int &numDiffCoord, double coord[], double r[],
                    int left[], int point[], int right[], DisMode mod);
    int Find(double diffCoord[], int min, int max, double target);
    void InsertData(ThreadData &data, int from, int x_centre,
                    int x_left, int x_right, double x, double y, double rx, double ry,
                    int x_mod, int x_max, bool real);
    void BucketSort(int num_event, Event tmp_event[], Event **event, int pnum,
                    ThreadData &td);
    void BuildTree(Node tree[], int left, int right);
    void Insert(Node tree[], int root, int fir[], int next[],
                int ins_which, int ins_time, int btm, int backtracker[]);
    void GetPair(Node tree[], int root, ThreadData &td, int ins_which,
                 int insert_time[], int ins_lim, int btm, int top, int fir[],
                 int next[]);
    bool enable_timer() const { return enable_timer_; }
    bool equal(double a, double b);

    struct PointLess {
        bool operator()(Point *a, Point *b) const { return a->x < b->x; }
    };

    double *diff_x_;
    int *rect_left_;
    int *rect_centre_x_;
    int *rect_right_;
    Point **point_;
    Point **final_;
    int **pseudo_prev_point_;
    int **pseudo_next_point_;
    int *num_prev_;
    int *num_next_;
    ThreadData *data_;
    double dxtime_, prtime_, time_;

    long int mem_size_;
    int npair_;
    int npnt_;
    bool enable_timer_;
};

#endif // NEIGHBOUR_SEARCH_H
