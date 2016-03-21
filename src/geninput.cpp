#include "geninput.h"

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#include <string>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <H5hut.h>

using namespace std;

GenInput::GenInput()
{
}

void GenInput::generate(int L, int W, string modeStr)
{
  long numParticle = L * W;

  id_ = new long[numParticle]();
  x_ = new double[numParticle]();
  y_ = new double[numParticle]();
  rx_ = new double[numParticle]();
  ry_ = new double[numParticle]();

  char filename[100];

  bool gen_success = false;

  if (modeStr == "randrand")
  {
    sprintf(filename, "input%ldrandrand.h5part", numParticle);
    gen_success = GenerateRandRand(filename, numParticle, L, W);
  }
  else if (modeStr == "randunif")
  {
    sprintf(filename, "input%ldrandunif.h5part", numParticle);
    gen_success = GenerateRandUnif(filename, numParticle, L, W);
  }
  else if (modeStr == "unifrand")
  {
    sprintf(filename, "input%ldunifrand.h5part", numParticle);
    gen_success = GenerateUnifRand(filename, numParticle, L, W);
  }
  else if (modeStr == "unifunif")
  {
    sprintf(filename, "input%ldunifunif.h5part", numParticle);
    gen_success = GenerateUnifUnif(filename, numParticle, L, W);
  }
  else
  {
    cerr << "Illegal mode!" << endl;
    return;
  }
}

bool GenInput::GenerateRandRand(const char *filename, int numParticle,
                                int L, int W)
{
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<> l_dist(0, L);
  boost::random::uniform_real_distribution<> w_dist(0, W);
  boost::random::uniform_real_distribution<> r_dist(0, 2);

  h5_file_p outFile = H5OpenFile(filename, H5_O_WRONLY, 0 /* MPI_Comm */);

  if (H5CheckFile(outFile) == H5_SUCCESS) // In this implementation, SUCCESS returns 0 !
  {
    cerr << "File: " << filename << " NOT exists!" << endl;
    return false;
  }

  H5SetStep(outFile, 0);
  H5PartSetNumParticles(outFile, numParticle);

  double scale = 0.01;
  double r_mean = 3 * scale;

  for(long int i = 0; i < numParticle; i++)
  {
    id_[i] = i;
    x_[i] = l_dist(gen) * scale;
    y_[i] = w_dist(gen) * scale;
    rx_[i] = (r_dist(gen) - 1) * scale * 0.2 + r_mean;
    ry_[i] = (r_dist(gen) - 1) * scale * 0.2 + r_mean;
  }

  WriteFile(outFile);
  H5CloseFile(outFile);

  return true;
}

bool GenInput::GenerateRandUnif(const char *filename,
                                int numParticle, int L, int W)
{
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<> l_dist(0, L);
  boost::random::uniform_real_distribution<> w_dist(0, W);

  h5_file_p outFile = H5OpenFile(filename, H5_O_WRONLY, 0 /* MPI_Comm */);

  if (H5CheckFile(outFile) != H5_SUCCESS) // In the implementation, SUCCESS returns 0 !
  {
    cerr << "File: " << filename << " NOT exists!" << endl;
    return false;
  }

  H5SetStep(outFile, 0);
  H5PartSetNumParticles(outFile, numParticle);

  double scale = 0.01;
  double r_mean = 3 * scale;

  for(long int i = 0; i < numParticle; i++)
  {
    id_[i] = i;
    x_[i] = l_dist(gen) * scale;
    y_[i] = w_dist(gen) * scale;
    rx_[i] = r_mean;
    ry_[i] = r_mean;
  }

  WriteFile(outFile);
  H5CloseFile(outFile);
  return true;
}

bool GenInput::GenerateUnifRand(const char *filename,
                                int numParticle, int L, int W)
{
  boost::random::mt19937 gen;
  boost::random::uniform_real_distribution<> r_dist(0, 2);

  h5_file_p outFile = H5OpenFile(filename, H5_O_WRONLY, 0 /* MPI_Comm */);

  if (H5CheckFile(outFile) != H5_SUCCESS) // In the implementation, SUCCESS returns 0 !
  {
    cerr << "File: " << filename << " NOT exists!" << endl;
    return false;
  }

  H5SetStep(outFile, 0);
  H5PartSetNumParticles(outFile, numParticle);

  double scale = 0.01;
  double r_mean = 3 * scale;

  for(int i = 0; i < L; i++)
  {
    for(int j = 0; j < W; j++)
    {
      long int index = i * L + j;
      id_[index] = index;
      x_[index] = i * scale;
      y_[index] = j * scale;
      rx_[index] = (r_dist(gen) - 1) * scale * 0.2 + r_mean;
      ry_[index] = (r_dist(gen) - 1) * scale * 0.2 + r_mean;
    }
  }
  WriteFile(outFile);
  H5CloseFile(outFile);
  return true;
}

bool GenInput::GenerateUnifUnif(const char *filename,
                                int numParticle, int L, int W)
{
  h5_file_p outFile = H5OpenFile(filename, H5_O_WRONLY, 0 /* MPI_Comm */);

  if (H5CheckFile(outFile) != H5_SUCCESS) // In the implementation, SUCCESS returns 0 !
  {
    cerr << "File: " << filename << " NOT exists!" << endl;
    return false;
  }

  H5SetStep(outFile, 0);
  H5PartSetNumParticles(outFile, numParticle);

  double scale = 0.01;
  double r_mean = 3 * scale;
  long int index = 0;

  for(int i = 0; i < L; i++)
  {
    for(int j = 0; j < W; j++)
    {
      id_[index] = index;
      x_[index] = i * scale;
      y_[index] = j * scale;
      rx_[index++] = r_mean;
      ry_[index++] = r_mean;
    }
  }
  WriteFile(outFile);
  H5CloseFile(outFile);
  return true;
}

void GenInput::WriteFile(h5_file_t *file)
{
  H5PartWriteDataInt64(file, "id", id_);
  H5PartWriteDataFloat64(file, "x", x_);
  H5PartWriteDataFloat64(file, "y", y_);
  H5PartWriteDataFloat64(file, "rx", rx_);
  H5PartWriteDataFloat64(file, "ry", ry_);
}

GenInput::~GenInput()
{
  delete [] id_;
  delete [] x_;
  delete [] y_;
  delete [] rx_;
  delete [] ry_;
}
