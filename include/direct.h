#ifndef DIRECT_H
#define DIRECT_H

class Direct
{
public:
  Direct();
  ~Direct();

  void search(const char *inputFile, const char *outputFile, bool enableTime);
  void initialize(int npnt);

private:
  double *x_vec_;
  double *y_vec_;
  double *rx_vec_;
  double *ry_vec_;
};

#endif // DIRECT_H
