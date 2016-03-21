#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include <omp.h>
#include <timsort.hpp>

#include <H5hut.h>

#include "neighbour_search.h"
#include "timer.hpp"

using namespace std;

NeighbourSearch::NeighbourSearch(int num_thread)
  : kNumThread(num_thread),
    kPairFactor(60),
    kTreeIncrement(3),
    kNegativeInfinity(-1e10)
{

}

void NeighbourSearch::Init(int maxn)
{
  mem_size_ = 0;

  diff_x_ = new double[maxn + 1]();
  rect_left_ = new int[maxn + 1]();
  rect_centre_x_ = new int[maxn]();
  rect_right_ = new int[maxn]();
  point_ = new Point*[maxn + 1]();
  final_ = new Point*[maxn + 1]();
  for (int i = 0; i< maxn + 1; i++)
  {
    point_[i] = new Point();
  }
  data_ = new ThreadData[kNumThread]();
  int rsv_size = 3 * maxn / kNumThread;
  int tree_size = 2 * (rsv_size + kTreeIncrement) - 1;
  int tree_depth = ceil(log(tree_size) / log(2));
  for (int i = 0; i < kNumThread; i++)
  {

    data_[i].npoint = 0;
    data_[i].npair = 0;

    data_[i].from = new int[rsv_size]();
    data_[i].rect_left = new int[rsv_size + 1]();
    data_[i].rect_centre_x = new int[rsv_size]();
    data_[i].rect_right = new int[rsv_size]();
    data_[i].rect_bottom = new int[rsv_size]();
    data_[i].rect_centre_y = new int[rsv_size]();
    data_[i].rect_top = new int[rsv_size]();
    data_[i].rect_btm_time = new int[rsv_size]();
    data_[i].insert_time = new int[rsv_size]();
    data_[i].fir = new int[rsv_size]();
    data_[i].next = new int[rsv_size]();
    data_[i].fir_s = new int[4]();
    data_[i].fir_sp = new int[rsv_size + 2]();
    data_[i].next_s = new int[rsv_size * 3]();
    data_[i].next_sp = new int[rsv_size * 3]();

    data_[i].x = new double[rsv_size]();
    data_[i].y = new double[rsv_size]();
    data_[i].rx = new double[rsv_size]();
    data_[i].ry = new double[rsv_size]();
    data_[i].ydiff = new double[rsv_size]();

    data_[i].real = new bool[rsv_size]();

    //    data_[i].pair = new IntPair[kPairFactor * rsv_size]();

    data_[i].tmp_event = new Event[rsv_size * 3]();

    data_[i].final_event = new Event*[rsv_size * 3]();

    data_[i].tree = new Node[tree_size]();
    data_[i].backtracker = new int[tree_depth]();
    data_[i].get_queue = new int[tree_size]();

  }
  pseudo_prev_point_ = new int*[kNumThread]();
  pseudo_next_point_ = new int*[kNumThread]();
  for (int i = 0; i < kNumThread; i++)
  {
    pseudo_prev_point_[i] = new int[maxn / kNumThread]();
    pseudo_next_point_[i] = new int[maxn / kNumThread]();
  }
  num_prev_ = new int[kNumThread]();
  num_next_ = new int[kNumThread]();

  mem_size_ = sizeof(int) * (5 * maxn + 1 + kNumThread * (17 * rsv_size + 2)
                             + 2 * kNumThread + kNumThread * 2 * tree_size
                             + kNumThread * tree_depth)
              + sizeof(double) * (maxn + 1 + kNumThread * (4 * rsv_size))
              + sizeof(bool) * (kNumThread * rsv_size)
              + sizeof(IntPair) * (kNumThread * kPairFactor * rsv_size)
              + sizeof(Point) * (maxn + 1)
              + sizeof(Event) * (kNumThread * rsv_size * 3)
              + sizeof(Node) * kNumThread * (2 * (rsv_size + kTreeIncrement) - 1);

}

void NeighbourSearch::Reset(const int maxn)
{
  npair_ = 0;
#pragma omp parallel
  {
    int rank = omp_get_thread_num();
    int nstart = maxn * rank / kNumThread;
    int nend = maxn * (rank + 1) / kNumThread;
    int nsize = nend - nstart;

    memset(diff_x_ + nstart, 0, nsize * sizeof(double));
    memset(rect_left_ + nstart, 0, nsize * sizeof(int));
    memset(rect_centre_x_ + nstart, 0, nsize * sizeof(int));
    memset(rect_right_ + nstart, 0, nsize * sizeof(int));

    data_[rank].npair = 0;
    data_[rank].npoint = 0;
  }
}

void NeighbourSearch::Search(const char *infile, const char *outfile,
                             bool enable_timer)
{
  omp_set_num_threads(kNumThread);
  enable_timer_ = enable_timer;

  h5_file_t *h5file = H5OpenFile(infile, H5_O_RDONLY, 0 /* MPI_Comm */);

  if (H5CheckFile(h5file) != H5_SUCCESS) // In the implementation, SUCCESS returns 0 !
  {
    cerr << "File: " << infile << " NOT exists!" << endl;
    return;
  }

  H5SetStep(h5file, 0);
  npnt_ = H5PartGetNumParticles(h5file);

  long *id = new long[npnt_]();
  double *x = new double[npnt_]();
  double *y = new double[npnt_]();
  double *rx = new double[npnt_]();
  double *ry = new double[npnt_]();

  //  H5PartReadDataInt64(h5file, "id", id);
  H5PartReadDataFloat64(h5file,"x", x);
  H5PartReadDataFloat64(h5file,"y", y);
  H5PartReadDataFloat64(h5file,"rx", rx);
  H5PartReadDataFloat64(h5file,"ry", ry);

  Init(npnt_);
  Reset(npnt_);

  if (kNumThread == 1)
  {
    data_[0].npoint = npnt_;
  }
  else
  {
    point_[0]->x = kNegativeInfinity;
    point_[0]->y = kNegativeInfinity;
    point_[0]->rx = -1;
    point_[0]->ry = -1;
    point_[0]->from = -1;
  }

  for (int i = 0; i < npnt_; i++)
  {
    if (kNumThread == 1)
    {
      data_[0].x[i] = x[i];
      data_[0].y[i] = y[i];
      data_[0].rx[i] = rx[i];
      data_[0].ry[i] = ry[i];
      data_[0].from[i] = i;
      data_[0].real[i] = true;
    }
    else
    {
      point_[i + 1]->x = x[i];
      point_[i + 1]->y = y[i];
      point_[i + 1]->rx = rx[i];
      point_[i + 1]->ry = ry[i];
      point_[i + 1]->from = i;
    }
  }

  H5CloseFile(h5file);
  delete [] id;
  delete [] x;
  delete [] y;
  delete [] rx;
  delete [] ry;
  cout << "Data Loaded" << endl;

  GetNeighbourPair(npnt_);
  assert(npair() <= kPairFactor * npnt_);

  for (int i = 0; i < kNumThread; i++)
  {
    cout << "********** Thread " << i << " **********" << endl;
    if (kNumThread != 1)
    {
      cout << "Sorting X: " << data_[i].sxtime << " s" << endl;
    }
    cout << "Discretizing X: " << dxtime_ << " s" << endl;
    if (kNumThread != 1)
    {
      cout << "Dividing: " << data_[i].dvtime << " s" << endl;
    }
    cout << "Discretizing Y: " << data_[i].dytime << " s" << endl;
    cout << "Construcing Events: " << data_[i].cetime << " s" << endl;
    cout << "Sorting Events: " << data_[i].setime << " s" << endl;
    cout << "Scanning Initialization: " << data_[i].istime << " s" << endl;
    cout << "Building Tree: " << data_[i].bttime << " s" << endl;
    cout << "Scanning: " << data_[i].sctime << " s" << endl;
    cout << "Point: " << data_[i].npoint << endl;
    cout << "Pair: " << data_[i].npair << endl;
    cout << endl;
  }

  cout << "Pair: " << npair() << endl;
  cout << "Memory: " << MemStats() << " MB" << endl;
  cout << "Total time: " << time_ << " s" << endl;

  //  int *count = new int[npnt_];
  //  memset(count, 0, sizeof(int) * npnt_);
  //  for (int i = 0; i < kNumThread; i++)
  //  {
  //    for (int j = 0; j < data_[i].npair; j++)
  //    {
  //      int a = data_[i].pair[j][0];
  //      int b = data_[i].pair[j][1];
  //      count[a]++;
  //      count[b]++;
  //    }
  //  }
  //  ofstream out("count");
  //  for (int i = 0; i < npnt_; i++)
  //  {
  //    out << i << "\t" << count[i] << endl;
  //  }
  //  out.close();
  //  delete [] count;
}

void NeighbourSearch::GetNeighbourPair(const int npnt)
{
  Timer tclk, clk;
  int num_diff = 0;
  double start = omp_get_wtime();
  if (kNumThread == 1)
  {
    clk.restart();
    Discretize(npnt, diff_x_, num_diff, data_[0].x, data_[0].rx,
               data_[0].rect_left, data_[0].rect_centre_x,
               data_[0].rect_right, HORIZONTAL);
    dxtime_ = clk.elapsed();
  }
  else
  {
    clk.restart();
    ParallelDiscretize(npnt, point_, final_, num_diff);
    dxtime_ = clk.elapsed();
    ParallelDivide(npnt, num_diff, final_);
  }

  clk.restart();
#pragma omp parallel
  {
    int rank = omp_get_thread_num();
    Process(data_[rank]);
  }
  //    for (int i = 0; i < kNumThread; i++)
  //    {
  //        Process(data_[i]);
  //    }
  prtime_ = clk.elapsed();
  time_ = tclk.elapsed();
  double end = omp_get_wtime();
  cout << "OMP time: " << end - start << endl;
}

void NeighbourSearch::ParallelDiscretize(int npnt, Point **original, Point **final,
                                         int &num_diff)
{
  Timer clk;
  int nfinal = npnt + 1;
  PSRS(nfinal, original, final);

  clk.restart();
  diff_x_[num_diff] = kNegativeInfinity;
  for(int i = 1; i < nfinal; i++)
  {
    if(final[i]->x > final[i - 1]->x)
    {
      diff_x_[++num_diff] = final[i]->x;
    }
  }
  //    cout << num_diff << endl;
  //    cout << "Generate diff x: " << clk.elapsed() << " s" << endl;

  clk.restart();
  double x = 0, rx = 0;
#pragma omp parallel for private(x, rx)
  for (int i = 0; i < npnt; i++)
  {
    x = final[i + 1]->x;
    rx = final[i + 1]->rx;


    int tmp = Find(diff_x_, 0, num_diff, x - rx);
    rect_left_[i] = equal(diff_x_[tmp], x - rx) ? tmp - 1 : tmp;
    //rect_left_[i] = Find(diff_x_, 0, num_diff, x - r - FLT_EPSILON * 2);
    rect_centre_x_[i] = Find(diff_x_, 0, num_diff, x);
    tmp = Find(diff_x_, 0, num_diff, x + rx);
    rect_right_[i] = equal(diff_x_[tmp + 1], x + rx) ? tmp +  1 : tmp;
    //        rect_right_[i] = Find(diff_x_, 0, num_diff, x + r);
  }
  //    cout << "Binary Search: " << clk.elapsed() << " s" << endl;
}

void NeighbourSearch::ParallelDivide(int npnt, int num_diff, Point**final)
{
#pragma omp parallel
  {
    Timer clk;
    int rank = omp_get_thread_num();
    int nth = kNumThread;
    int first = rank * npnt / nth;
    int last = (rank + 1) * npnt / nth;
    int div = Find(diff_x_, 0, num_diff, final[first + 1]->x);
    int next_div = rank != nth - 1 ?
                             Find(diff_x_, 0, num_diff, final[last + 1]->x) : num_diff + 1;

    int x_mod = div - 1;
    int x_max = next_div - div + 1;

    int *for_prev = pseudo_prev_point_[rank];
    int *for_next = pseudo_next_point_[rank];
    int nprev = 0;
    int nnext = 0;
    for (int i = first; i < last; i++)
    {
      int from = final[i + 1]->from;
      InsertData(data_[rank], from, rect_centre_x_[i], rect_left_[i],
                 rect_right_[i], final[i + 1]->x, final[i + 1]->y,
                 final[i + 1]->rx, final[i + 1]->ry, x_mod, x_max, true);
      if (rect_left_[i] < x_mod + 1 && rank != 0)
        for_prev[nprev++] = i;
      if (rect_right_[i] - x_mod >= x_max && rank + 1 < nth)
        for_next[nnext++] = i;
    }

    assert(nprev <= 2 * (npnt / kNumThread));
    assert(nnext <= 2 * (npnt / kNumThread));
    num_prev_[rank] = nprev;
    num_next_[rank] = nnext;

#pragma omp barrier

    if (rank != 0)
    {
      int *from_prev = pseudo_next_point_[rank - 1];
      for(int i = 0; i < num_next_[rank - 1]; i++)
      {
        int k = from_prev[i];
        int from = final[k + 1]->from;
        InsertData(data_[rank], from, rect_centre_x_[k], rect_left_[k],
                   rect_right_[k], final[k + 1]->x, final[k + 1]->y,
                   final[k + 1]->rx, final[k + 1]->ry, x_mod, x_max, false);
      }
    }
    if (rank != nth - 1)
    {
      int *from_next = pseudo_prev_point_[rank + 1];
      for(int i = 0; i < num_prev_[rank + 1]; i++)
      {
        int k = from_next[i];
        int from = final[k + 1]->from;
        InsertData(data_[rank], from, rect_centre_x_[k], rect_left_[k],
                   rect_right_[k], final[k + 1]->x, final[k + 1]->y,
                   final[k + 1]->rx, final[k + 1]->ry, x_mod, x_max, false);
      }
    }

    data_[rank].dvtime = clk.elapsed();
  }
}

void NeighbourSearch::Process(ThreadData &td)
{
  Timer clk;
  int pnum = td.npoint;
  int *rect_left = td.rect_left;
  int *rect_centre_x = td.rect_centre_x;
  int *rect_right = td.rect_right;
  int *rect_bottom = td.rect_bottom;
  int *rect_centre_y = td.rect_centre_y;
  int *rect_top = td.rect_top;
  double *y = td.y;
  double *rx = td.rx;
  double *ry = td.ry;
  bool *real = td.real;
  int num_diff = 0;
  double *ydiff = td.ydiff;

  //  cout << pnum << endl;

  clk.restart();
  Discretize(pnum, ydiff, num_diff, y, ry, rect_bottom,
             rect_centre_y, rect_top, VERTICAL);

  td.dytime = clk.elapsed();
  clk.restart();

  Event *tmp_event = td.tmp_event;
  Event **final_event = td.final_event;
  Event *event;

  int num_event = 0;
  for (int i = 0; i < pnum * 3; i++)
  {
    if (i % 3 != 1 || real[i / 3])
    {
      event = &tmp_event[num_event];
      event->x = i % 3 == 0 ? rect_left[i / 3] :
                              (i % 3 == 1 ? rect_centre_x[i / 3] :
                                            rect_right[i / 3]);
      event->type = i % 3;
      event->bottom_point = i % 3 == 1 ? rect_centre_y[i / 3] :
                                         rect_bottom[i / 3];
      event->top_point = i % 3 == 1 ? rect_centre_y[i / 3] :
                                      rect_top[i / 3];
      event->which = i / 3;
      num_event++;
    }
  }
  td.cetime = clk.elapsed();
  clk.restart();

  BucketSort(num_event, tmp_event, final_event, pnum, td);
  td.setime = clk.elapsed();
  clk.restart();

  int *rect_btm_time = td.rect_btm_time;
  int *insert_time = td.insert_time;
  int *fir = td.fir;
  int *next = td.next;

  memset(fir, -1, sizeof(int) * (pnum + 1));

  for (int i = 0; i < num_event; i++)
  {
    if (final_event[i]->type == 1)
    {
      insert_time[final_event[i]->which] = i;
    }
  }
  td.istime = clk.elapsed();
  clk.restart();

  Node *tree = td.tree;
  BuildTree(tree, 1, num_diff + kTreeIncrement);
  td.bttime = clk.elapsed();
  clk.restart();

  int ins_which, ins_time, btm, ins_lim, top;
  int *backtracker = td.backtracker;

  for (int i = 0; i < num_event; i++)
  {
    event = final_event[i];
    if (event->type == 0)
    {
      rect_btm_time[event->which] = i;
    }
    if (event->type == 1)
    {
      ins_which = event->which;
      ins_time = i;
      btm = event->bottom_point;
      Insert(tree, 0, fir, next, ins_which, ins_time, btm, backtracker);
    }
    if (event->type == 2 && tree[0].last >= rect_left[event->which])
    {
      ins_lim = rect_btm_time[event->which];
      btm = event->bottom_point;
      top = event->top_point;
      ins_which = event->which;
      if (tree[0].last >= ins_lim)
        GetPair(tree, 0, td, ins_which, insert_time, ins_lim, btm,
                top, fir, next);
    }
  }
  td.sctime = clk.elapsed();
}

void NeighbourSearch::PSRS(int npnt, Point **original, Point **final)
{
  int nth = kNumThread;
  double samples[nth * nth];
  int pivot[nth][nth + 1];
  int seg_length[nth][nth];
  int new_bound[nth + 1];

#pragma omp parallel
  {
    Timer clk;
    int rank = omp_get_thread_num();
    int start = rank * npnt / nth;
    int end = (rank + 1)  * npnt / nth;
    int size = end - start;
    sort(original + start, original + end, PointLess());

    for (int i = 0; i < nth; i++)
    {
      samples[rank * nth + i] = original[start + size * i / nth]->x;
    }

#pragma omp barrier

#pragma omp single
    {
      sort(samples, samples + nth * nth);
    }

#pragma omp barrier
    pivot[rank][0] = start;
    for (int i = 1; i < nth; i++)
    {
      pivot[rank][i] = SubFind(original, start, start + size - 1,
                               samples[nth * i]);
    }
    pivot[rank][nth] = end;

#pragma omp barrier

#pragma omp single
    {
      new_bound[0] = 0;
      for (int i = 0; i < nth; i++)
      {
        new_bound[i + 1] = new_bound[i];
        for (int j = 0; j < nth; j++)
        {
          seg_length[i][j] = pivot[j][i + 1] - pivot[j][i];
          new_bound[i + 1] += seg_length[i][j];
        }
      }
    }

    int local_start = new_bound[rank];
    int local_pivot[nth];
    local_pivot[0] = 0;
    for (int i = 0; i < nth; i++)
    {
      memcpy(final + local_start + local_pivot[i], original + pivot[i][rank],
             sizeof(Point*) * seg_length[rank][i]);
      local_pivot[i + 1] = local_pivot[i] + seg_length[rank][i];
    }

    //timsort(final + new_bound[rank], final + new_bound[rank+1], PointLess());
    //sort(final + new_bound[rank], final + new_bound[rank+1], PointLess());
    MergeSort(local_pivot, nth, local_pivot[nth] - local_pivot[0],
              final + new_bound[rank]);
    data_[rank].sxtime = clk.elapsed();
  }

}

void NeighbourSearch::MergeSort(int *local_pivot, int segnum, int size,
                                Point **final)
{
  Point **buffer = new Point*[size]();
  while (segnum > 1)
  {
    MergePass(local_pivot, segnum, final, buffer);
    segnum = segnum % 2 == 0 ? segnum / 2 : segnum / 2 + 1;

    MergePass(local_pivot, segnum, buffer, final);
    segnum = segnum % 2 == 0 ? segnum / 2 : segnum / 2 + 1;
  }
  delete [] buffer;
}

void NeighbourSearch::MergePass(int *local_pivot, int segnum, Point **a, Point **b)
{
  int i = 0;
  while (i < segnum - 1)
  {
    Merge(local_pivot[i], local_pivot[i + 1], local_pivot[i + 2], a, b);
    local_pivot[i / 2] = local_pivot[i];
    local_pivot[i / 2 + 1] = local_pivot[i + 2];
    i += 2;
  }

  if (segnum % 2)
  {
    memcpy(b + local_pivot[i], a + local_pivot[i],
           (local_pivot[i + 1] - local_pivot[i]) * sizeof(Point*));
    local_pivot[i / 2 + 1] = local_pivot[i + 1];
  }
}

void NeighbourSearch::Merge(int low, int mid, int high, Point **a, Point **b)
{
  int i = low;
  int left_min = low;
  int left_max = mid;
  int right_min = mid;
  int right_max = high;

  while (i < high)
  {
    if (left_min >= left_max)
    {
      b[i++] = a[right_min++];
    }
    else if (right_min >= right_max)
    {
      b[i++] = a[left_min++];
    }
    else
    {
      b[i++] = a[left_min]->x <= a[right_min]->x ? a[left_min++] : a[right_min++];
    }
  }
}

int NeighbourSearch::SubFind(Point **point, int min, int max, double target)
{
  while (min < max)
  {
    if (point[(min + max) / 2 + 1]->x <= target)
    {
      min = (min + max) / 2 + 1;
    }
    else
    {
      max = (min + max) / 2;
    }
  }
  return point[min]->x <= target ? min + 1 : min;
}

int NeighbourSearch::Find(double diffCoord[], int min, int max, double target)
{
  int lowerBnd = min;
  int upperBnd = max;

  //  if (fabs(diffCoord[(min + max) / 2] - target) <= DBL_EPSILON)
  //  {
  //    lowerBnd = (min + max) / 2;
  //  }
  //  else
  //  {
  while (lowerBnd < upperBnd)
  {
    if (diffCoord[(lowerBnd + upperBnd) / 2 + 1] > target)
    {
      upperBnd = (lowerBnd + upperBnd) / 2;
    }
    else
    {
      lowerBnd = (lowerBnd + upperBnd) / 2 + 1;
    }
  }
  //  }
  return lowerBnd;
}

void NeighbourSearch::InsertData(ThreadData &data, int from, int x_centre,
                                 int x_left, int x_right, double x, double y,
                                 double rx, double ry, int x_mod, int x_max, bool real)
{
  int id = data.npoint;

  data.rect_centre_x[id] = x_centre - x_mod;
  data.rect_left[id] = max(x_left - x_mod, 0);
  data.rect_right[id] = min(x_right - x_mod, x_max);
  data.rx[id] = rx;
  data.ry[id] = ry;
  data.y[id] = y;
  data.x[id] = x;
  data.from[id] = from;
  data.real[id] = real;
  data.npoint++;
}

void NeighbourSearch::Discretize(int npnt, double diff[], int &ndiff,
                                 double coord[], double r[], int left[],
                                 int centre[], int right[], DisMode mod)
{
  Timer clk;
  memcpy(diff + 1, coord, sizeof(double) * npnt);

  diff[0] = kNegativeInfinity;
  sort(diff, diff + npnt + 1);


  for (int i = 1; i <= npnt; i++)
  {
    if (diff[i] > diff[i - 1])
    {
      diff[++ndiff] = diff[i];
    }
  }
  //    cout << ndiff << endl;

  clk.restart();
  double pos, radius;
  for (int i = 0; i < npnt; i++)
  {
    pos = coord[i];
    radius = r[i];

    int tmp = Find(diff, 0, ndiff, pos - radius);
    left[i] = (equal(diff[tmp], pos - radius) ? tmp - 1 : tmp) + mod;
    //left[i] = Find(diff, 0, ndiff, pos - radius - FLT_EPSILON * 2) + mod;
    centre[i] = Find(diff, 0, ndiff, pos);
    //        right[i] = Find(diff, 0, ndiff, pos + radius);
    tmp = Find(diff, 0, ndiff, pos + radius);
    right[i] = equal(diff[tmp + 1], pos + radius) ? tmp + 1 : tmp;
  }
  //  cout << clk.elapsed() << endl;
}

void NeighbourSearch::BucketSort(int num_event, Event tmp_event[], Event **event,
                                 int pnum, ThreadData &td)
{
  int *fir_s = td.fir_s;
  int *fir_sp = td.fir_sp;
  int *next_s = td.next_s;
  int *next_sp = td.next_sp;
  fir_s[1] = -1;
  fir_s[2] = -1;

  for (int i = 0; i < num_event; i++)
  {
    int tmp = tmp_event[i].type == 1 ? 1 : 2;
    next_s[i] = fir_s[tmp];
    fir_s[tmp] = i;
  }

  memset(fir_sp, -1, sizeof(int) * (pnum + 1));

  for (int i = 2; i > 0; i--)
  {
    for (int j = fir_s[i]; j != -1; j = next_s[j])
    {
      int tmp = tmp_event[j].x;
      //            if(tmp<0)
      //                cout << tData.from[tmpEvent[j].which] << endl;
      next_sp[j] = fir_sp[tmp];
      fir_sp[tmp] = j;
    }
  }

  int sum = 0;
  for(int i = 0; i <= pnum; i++)
  {
    for(int j = fir_sp[i]; j != -1; j = next_sp[j])
    {
      event[sum++] = &tmp_event[j];
    }
  }
}

void NeighbourSearch::BuildTree(Node tree[], int left, int right)
{
  int now = 0;
  int sum = 1;

  tree[now].x = left;
  tree[now].y = right;

  while (now < sum)
  {
    tree[now].sum = 0;
    tree[now].last = -1;
    tree[now].sec_last = -1;
    tree[now].to_leaf = -1;

    if (tree[now].x == tree[now].y)
    {
      now++;
      continue;
    }
    tree[sum].x = tree[now].x;
    tree[sum].y = (tree[now].x + tree[now].y) / 2;
    tree[now].lson = sum++;
    tree[sum].x = (tree[now].x + tree[now].y) / 2 + 1;
    tree[sum].y = tree[now].y;
    tree[now].rson = sum++;
    now++;
  }
  //  int now = ++tsum;
  //  tree[now].x = left;
  //  tree[now].y = right;
  //  tree[now].sum = 0;
  //  tree[now].last = -1;
  //  tree[now].sec_last = -1;
  //  tree[now].to_leaf = -1;

  //  if (left == right)
  //  {
  //      return;
  //  }

  //  tree[now].lson = tsum + 1;
  //  BuildTree(tree, left, (left + right) / 2, queue, tsum);

  //  tree[now].rson = tsum + 1;
  //  BuildTree(tree, (left + right) / 2 + 1, right, queue, tsum);

}

void NeighbourSearch::Insert(Node tree[], int root, int fir[], int next[],
                             int ins_which, int ins_time, int btm,
                             int backtracker[])
{
  int depth = 0;
  int now = root;
  while (tree[now].x != tree[now].y)
  {
    tree[now].sum++;
    tree[now].sec_last = tree[now].last;
    tree[now].last = ins_time;
    backtracker[depth++] = now;
    if (btm <= tree[tree[now].lson].y)
    {
      now = tree[now].lson;
    }
    else
    {
      now = tree[now].rson;
    }
  }
  tree[now].sum++;
  tree[now].sec_last = tree[now].last;
  tree[now].last = ins_time;

  next[ins_which] = fir[tree[now].x];
  fir[tree[now].x] = ins_which;
  tree[now].to_leaf = now;

  for (int i = 0; i < depth; i++)
  {
    tree[backtracker[i]].to_leaf = now;
  }

  //  tree[now].sum++;
  //  tree[now].sec_last = tree[now].last;
  //  tree[now].last = ins_time;

  //  if (tree[now].x == tree[now].y)
  //  {
  //    next[ins_which] = fir[tree[now].x];
  //    fir[tree[now].x] = ins_which;
  //    tree[now].to_leaf = now;
  //    target = now;
  //    return;
  //  }

  //  if (btm <= tree[tree[now].lson].y)
  //  {
  //    Insert(tree, tree[now].lson, target, fir, next, ins_which, ins_time, btm);
  //  }
  //  else
  //  {
  //    Insert(tree, tree[now].rson, target, fir, next, ins_which, ins_time, btm);
  //  }

  //  tree[now].to_leaf = target;
}

void NeighbourSearch::GetPair(Node tree[], int root, ThreadData &td,
                              int ins_which, int insert_time[], int ins_lim,
                              int btm, int top, int fir[], int next[])
{
  int now = -1;
  int sum = 0;

//  static int get_queue[100];
  int *get_queue = td.get_queue;
  get_queue[sum] = root;

  while (now++ < sum)
  {
    if (tree[get_queue[now]].x == tree[get_queue[now]].y)
    {
      for (int tmp = fir[tree[get_queue[now]].x]; (tmp != -1)
           && (insert_time[tmp] >= ins_lim); tmp = next[tmp])
      {
        int a = td.from[ins_which];
        int b = td.from[tmp];
        if (a != b)
        {
          int id = td.npair;
          //          td.pair[id][0] = a;
          //          td.pair[id][1] = b;
          td.npair++;
        }
      }
    }
    else
    {
      if (tree[get_queue[now]].sec_last < ins_lim
          && tree[tree[get_queue[now]].to_leaf].x >= btm
          && tree[tree[get_queue[now]].to_leaf].x <= top)
      {
        get_queue[++sum] = tree[get_queue[now]].to_leaf;
        continue;
      }

      if (tree[tree[get_queue[now]].lson].last >= ins_lim
          && tree[tree[get_queue[now]].lson].y >= btm)
      {
        get_queue[++sum] = tree[get_queue[now]].lson;
      }

      if (tree[tree[get_queue[now]].rson].last >= ins_lim
          && tree[tree[get_queue[now]].rson].x <= top)
      {
        get_queue[++sum] = tree[get_queue[now]].rson;
      }
    }

  }

//  if (ins_which % 1000 == 0)
//  cout << sum << endl;

  //    if (tree[now].x == tree[now].y)
  //    {
  //      for (int tmp = fir[tree[now].x]; (tmp != -1)
  //           && (insert_time[tmp] >= ins_lim); tmp = next[tmp])
  //      {
  //        int a = td.from[ins_which];
  //        int b = td.from[tmp];
  //        if (a != b)
  //        {
  //          int id = td.npair;
  //          td.pair[id][0] = a;
  //          td.pair[id][1] = b;
  //          td.npair++;
  //        }
  //      }
  //      return;
  //    }

  //    if (tree[now].sec_last < ins_lim && tree[tree[now].to_leaf].x >= btm
  //        && tree[tree[now].to_leaf].x <= top)
  //    {
  //      GetPair(tree, tree[now].to_leaf, td, ins_which, insert_time, ins_lim,
  //              btm, top, fir, next);
  //      return;
  //    }

  //    if (tree[tree[now].lson].last >= ins_lim && tree[tree[now].lson].y >= btm)
  //    {
  //      GetPair(tree, tree[now].lson, td, ins_which, insert_time, ins_lim,
  //              btm, top, fir, next);
  //    }

  //    if(tree[tree[now].rson].last>=ins_lim && tree[tree[now].rson].x<=top)
  //    {
  //      GetPair(tree, tree[now].rson, td, ins_which, insert_time, ins_lim,
  //              btm, top, fir, next);
  //    }
}

int NeighbourSearch::npair()
{
  npair_ = 0;
  for (int i = 0; i < kNumThread; i++)
  {
    npair_ += data_[i].npair;
  }
  return npair_;
}

bool NeighbourSearch::equal(double a, double b)
{
  return (fabs(a - b) <= FLT_EPSILON * max(fabs(a), fabs(b)));
}

NeighbourSearch::~NeighbourSearch()
{
  delete [] diff_x_;
  delete [] rect_left_;
  delete [] rect_centre_x_;
  delete [] rect_right_;
  for (int i = 0; i < kNumThread; i++)
  {
    delete [] data_[i].from;
    delete [] data_[i].rect_left;
    delete [] data_[i].rect_centre_x;
    delete [] data_[i].rect_right;
    delete [] data_[i].rect_bottom;
    delete [] data_[i].rect_centre_y;
    delete [] data_[i].rect_top;
    delete [] data_[i].rect_btm_time;
    delete [] data_[i].insert_time;
    delete [] data_[i].fir;
    delete [] data_[i].next;
    delete [] data_[i].fir_s;
    delete [] data_[i].fir_sp;
    delete [] data_[i].next_s;
    delete [] data_[i].next_sp;

    delete [] data_[i].x;
    delete [] data_[i].y;
    delete [] data_[i].rx;
    delete [] data_[i].ry;
    delete [] data_[i].ydiff;

    delete [] data_[i].real;

    //    delete [] data_[i].pair;

    delete [] data_[i].tmp_event;
    delete [] data_[i].final_event;

    delete [] data_[i].tree;
    delete [] data_[i].backtracker;
    delete [] data_[i].get_queue;
  }
  delete [] data_;
  for(int i = 0; i < npnt_ + 1; i++)
  {
    delete point_[i];
  }
  delete [] point_;
  delete [] final_;
  for (int i = 0; i < kNumThread; i++)
  {
    delete [] pseudo_prev_point_[i];
    delete [] pseudo_next_point_[i];
  }
  delete [] pseudo_prev_point_;
  delete [] pseudo_next_point_;
  delete [] num_prev_;
  delete [] num_next_;
}
