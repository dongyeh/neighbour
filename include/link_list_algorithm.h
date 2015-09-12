#ifndef LINK_LIST_ALGORITHM_H
#define LINK_LIST_ALGORITHM_H

#include <vector>
#include <list>

class LinkListAlgorithm
{
public:
  LinkListAlgorithm();
  ~LinkListAlgorithm();

  void Search(const char *infile);

  int num_of_pair() const {return num_of_pair_;}

private:
  struct Cell {
    std::list<int> particle_list;
  };

  struct AdjCell {
    int x, y;
  };

  void LoadParticleFile(const char *infile);
  void GenerateCell();
  void LinkCellParticle();
  void FindNeighbours();

  int num_of_pair_;
  int num_of_particle_;
  int num_x_cell_;
  int num_y_cell_;
  double min_x_;
  double max_x_;
  double min_y_;
  double max_y_;
  double ref_r_;
  double *x_vec_;
  double *y_vec_;
  double *rx_vec_;
  double *ry_vec_;
  Cell **cells_;
};

#endif // LINK_LIST_ALGORITHM_H
