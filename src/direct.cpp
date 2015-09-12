#include "direct.h"
#include "timer.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

Direct::Direct()
{

}

void Direct::search(const char *inputFile, const char *outputFile, bool enableTime)
{
  Timer clk;

  ifstream inFile(inputFile);

  if(!inFile)
  {
      cerr << "Error while loading input file" << endl;
      return;
  }

  string line;
  int id;
  double x, y, r;

  getline(inFile, line);
  istringstream stream(line);
  stream >> n;

  initialize(n);

  for(int i=0; i<n; i++)
  {
      getline(inFile, line);
      istringstream dataStream(line);
      dataStream >> id >> x >> y >> r;
      Datx[i]=x;
      Daty[i]=y;
      DatR[i]=r;
  }

  inFile.close();

  for (int i = 0; i < n; i++)
  {
    double xa = Datx[i];
    double ya = Daty[i];
    double ra = DatR[i];
    for (int j = 0; j < n; j++)
    {
      if (i != j)
      {
      double distx = fabs(xa - Datx[j]);
      double disty = fabs(ya - Daty[j]);
      if (distx )
      }
    }
  }
}

void Direct::initialize(int npnt)
{
    Datx=new double[npnt]();
    Daty=new double[npnt]();
    DatR=new double[npnt]();

    double size = sizeof(double) * (3 * npnt);
    cout << "**********Memory**********" << endl;
    cout << size/(1024*1024) << " MB" << endl;
    cout << endl;
}

Direct::~Direct()
{
}
