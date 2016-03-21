#include "direct.h"

#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>

#include <H5hut.h>

#include "timer.hpp"

using namespace std;

Direct::Direct()
{

}

void Direct::search(const char *inputFile, const char *outputFile, bool enableTime)
{
  Timer clk;

  h5_file_t *h5file = H5OpenFile(inputFile, H5_O_RDONLY, 0 /* MPI_Comm */);

  if (H5CheckFile(h5file) != H5_SUCCESS) // In the implementation, SUCCESS returns 0 !
  {
    cerr << "File: " << inputFile << " NOT exists!" << endl;
    return;
  }

  H5SetStep(h5file, 0);
  int n = H5PartGetNumParticles(h5file);
  initialize(n);

  H5PartReadDataFloat64(h5file, "x", x_vec_);
  H5PartReadDataFloat64(h5file, "y", y_vec_);
  H5PartReadDataFloat64(h5file, "rx", rx_vec_);
  H5PartReadDataFloat64(h5file, "ry", ry_vec_);

  H5CloseFile(h5file);

  clk.restart();
  int npair = 0;
  for (int i = 0; i < n; i++)
  {
    double xa = x_vec_[i];
    double ya = y_vec_[i];
    double rx_a = rx_vec_[i];
    double ry_a = ry_vec_[i];

    for (int j = 0; j < n; j++)
    {
      if (i != j)
      {
        double distx = fabs(xa - x_vec_[j]);
        double disty = fabs(ya - y_vec_[j]);
        if (distx <= rx_a + numeric_limits<float>::epsilon() * max(fabs(distx), fabs(rx_a))
            && disty <= ry_a + numeric_limits<float>::epsilon() * max(fabs(disty), fabs(ry_a)))
        {
          npair++;
        }
      }
    }
  }

  double time = clk.elapsed();
  cout << "Time: " << time << endl;
  cout << "Pair: " << npair << endl;
}

void Direct::initialize(int npnt)
{
  x_vec_ = new double[npnt]();
  y_vec_ = new double[npnt]();
  rx_vec_ = new double[npnt]();
  ry_vec_ = new double[npnt]();

  double size = sizeof(double) * (3 * npnt);
  cout << "**********Memory**********" << endl;
  cout << size/(1024*1024) << " MB" << endl;
  cout << endl;
}

Direct::~Direct()
{
  delete [] x_vec_;
  delete [] y_vec_;
  delete [] rx_vec_;
  delete [] ry_vec_;
}
