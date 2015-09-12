#ifndef SEARCH_H
#define SEARCH_H

#include "global.h"

class Search
{
public:

  Search();
  ~Search();

  void start(const char *inputFile, const char *outputFile, bool enableTime);
  int pairNum();

private:
  enum DiscreteMode
  {
    Horizontal,
    Vertical
  };

  struct Event {
    int X, Type, bottomPoint, topPoint, Which;
  };
  struct Node {
    int x, y, lson, rson, sum, Last, SecL, toLeaf;
  };

  void initialize(int npnt);
  void GetNeighbourPair(int n, int corr[], double Datx[], double Daty[], double DatR[], int Pair[][2]);
  void discretize(double coord[], double r[], int left[], int point[], int right[], DiscreteMode Mod);
  int Find(int min,int max,double target);
  bool equal(double a, double b);
  bool lessequal(double a, double b);
  void bucketSort();
  void BuildTree(int x, int y);
  void Insert(int now, int &target);
  void GetPair(int now);
  double Dist(double x,double y);

  Node *tree;
  Event *tmpEvent;
  Event **event;
  double *diff;
  double *Datx, *Daty, *DatR;
  int *corr, *rectLeft, *rectBottom, *rectRight, *rectTop, *pointX, *pointY, *rectBtmTime, *Fir, *Next, *ITime, *FirSp, *NextSp, *NextS;
  int FirS[4], n, sum, numDiffCoord, Tsum, Psum;
  NeighbourPair *Pair;
  int IWhi,ITT,Ix,Iy,ILim,ITar;
  int _maxn, _maxpair;
  double dxTime, dyTime, ceTime, seTime, iaTime, btTime, scTime, nxTime, ttTime;
  bool enableProcessTime;
};

#endif // SEARCH_H
