#include "link_list_algorithm.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <limits>

#include <H5hut.h>

#include "timer.hpp"

using namespace std;

LinkListAlgorithm::LinkListAlgorithm()
  : num_of_pair_(0),
    num_of_particle_(0)
{
}


void LinkListAlgorithm::Search(const char *infile)
{
  LoadParticleFile(infile);

  Timer clk;
  GenerateCell();
  LinkCellParticle();
  FindNeighbours();
  cout << "Time: " << clk.elapsed() << endl;

  cout << "Num of pairs: " << num_of_pair() << endl;
}

void LinkListAlgorithm::LoadParticleFile(const char *infile)
{

  h5_file_t *h5file = H5OpenFile(infile, H5_O_RDONLY, 0 /* MPI_Comm */);

  if (H5CheckFile(h5file) != H5_SUCCESS) // In the implementation, SUCCESS returns 0 !
  {
    cerr << "File: " << infile << " NOT exists!" << endl;
    return;
  }

  H5SetStep(h5file, 0);
  num_of_particle_ = H5PartGetNumParticles(h5file);

  x_vec_ = new double[num_of_particle_]();
  y_vec_ = new double[num_of_particle_]();
  rx_vec_ = new double[num_of_particle_]();
  ry_vec_ = new double[num_of_particle_]();

  H5PartReadDataFloat64(h5file,"x", x_vec_);
  H5PartReadDataFloat64(h5file,"y", y_vec_);
  H5PartReadDataFloat64(h5file,"rx", rx_vec_);
  H5PartReadDataFloat64(h5file,"ry", ry_vec_);

  H5CloseFile(h5file);
}

void LinkListAlgorithm::GenerateCell()
{
  min_x_ = *min_element(x_vec_, x_vec_ + num_of_particle_);
  max_x_ = *max_element(x_vec_, x_vec_ + num_of_particle_);
  min_y_ = *min_element(y_vec_, y_vec_ + num_of_particle_);
  max_y_ = *max_element(y_vec_, y_vec_ + num_of_particle_);

  double sum_r = 0;
  for (int i = 0 ; i < num_of_particle_; i++)
  {
    sum_r += rx_vec_[i];
    sum_r += ry_vec_[i];
  }
  ref_r_ = sum_r / (num_of_particle_ * 2);

  double x_range = max_x_ - min_x_;
  double y_range = max_y_ - min_y_;

  num_x_cell_ = ceil(x_range / ref_r_);
  num_y_cell_ = ceil(y_range / ref_r_);
  cout << "Cells: " << num_x_cell_ << " x " << num_y_cell_ << endl;

  cells_ = new Cell*[num_x_cell_]();
  for (int i = 0; i < num_x_cell_; i++)
  {
    cells_[i] = new Cell[num_y_cell_]();
  }
}

void LinkListAlgorithm::LinkCellParticle()
{
  for (int i = 0; i < num_of_particle_; i++)
  {
    int x_id = static_cast<int>((x_vec_[i] - min_x_) / ref_r_);
    int y_id = static_cast<int>((y_vec_[i] - min_y_) / ref_r_);
    cells_[x_id][y_id].particle_list.push_back(i);
  }
}

void LinkListAlgorithm::FindNeighbours()
{
  AdjCell adj_cell[5] = {{-1, 1}, {0, 1}, {1, 1}, {0, 0}, {1, 0}};
  AdjCell complete_adj_cell[9] = {{-1, 1}, {0, 1}, {1, 1},
                                  {-1, 0}, {0, 0}, {1, 0},
                                  {-1, -1}, {0, -1}, {1, -1}};
  for (int i = 0; i < num_y_cell_; i++)
  {
    for (int j = 0; j < num_x_cell_; j++)
    {
      //      for (int k = 0; k < 5; k++)
      //      {
      //        int cur_x = j + adj_cell[k].x;
      //        int cur_y = i + adj_cell[k].y;
      for (int k = 0; k < 9; k++)
      {
        int cur_x = j + complete_adj_cell[k].x;
        int cur_y = i + complete_adj_cell[k].y;
        if (cur_x >= 0 && cur_x < num_x_cell_ && cur_y >= 0 && cur_y < num_y_cell_)
        {
          for (list<int>::iterator iter = cells_[j][i].particle_list.begin();
               iter != cells_[j][i].particle_list.end(); iter++)
          {
            for (list<int>::iterator adj_iter = cells_[cur_x][cur_y].particle_list.begin();
                 adj_iter != cells_[cur_x][cur_y].particle_list.end(); adj_iter++)
            {
              int a = *iter;
              int b = *adj_iter;

              if (a != b)
              {
                double distx = fabs(x_vec_[a] - x_vec_[b]);
                double disty = fabs(y_vec_[a] - y_vec_[b]);
                if (distx <= rx_vec_[a] + numeric_limits<float>::epsilon() * max(fabs(distx), fabs(rx_vec_[a]))
                    && disty <= ry_vec_[a] + numeric_limits<float>::epsilon() * max(fabs(disty), fabs(ry_vec_[a])))
                  //                double distx = xa - Datx[j];
                  //                double disty = ya - Daty[j];
                  //                if (fabs(distx) <= ra + FLT_EPSILON * max(fabs(distx), fabs(ra))
                  //                    && fabs(disty) <= ra + FLT_EPSILON * max(fabs(disty), fabs(ra)))
                {
                  num_of_pair_++;
                }
              }
            }
          }
        }
      }
    }
  }
}

LinkListAlgorithm::~LinkListAlgorithm()
{ 
  for (int i = 0; i < num_x_cell_; i++)
  {
    delete [] cells_[i];
  }
  delete [] cells_;

  delete [] x_vec_;
  delete [] y_vec_;
  delete [] rx_vec_;
  delete [] ry_vec_;
}
