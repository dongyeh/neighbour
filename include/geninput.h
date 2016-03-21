#ifndef GENINPUT_H
#define GENINPUT_H

#include <string>

#include <H5hut.h>

class GenInput
{
public:
  GenInput();
  ~GenInput();

  void generate(int L, int W, std::string modeStr);

private:
  bool GenerateRandRand(const char *filename,
                        int numParticle, int L, int W);
  bool GenerateRandUnif(const char *filename,
                        int numParticle, int L, int W);
  bool GenerateUnifRand(const char *filename,
                        int numParticle, int L, int W);
  bool GenerateUnifUnif(const char *filename,
                        int numParticle, int L, int W);
  void WriteFile(h5_file_t *file);

  long *id_;
  double *x_;
  double *y_;
  double *rx_;
  double *ry_;
};

#endif // GENINPUT_H
