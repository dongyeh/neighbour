#ifndef DIRECT_H
#define DIRECT_H

class Direct
{
public:
  Direct();
  ~Direct();

  void search(const char *inputFile, const char *outputFile, bool enableTime);
  void initialize(int npnt);
};

#endif // DIRECT_H
