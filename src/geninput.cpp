#include "geninput.h"

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

GenInput::GenInput()
{
}

void GenInput::generate(int numParticle, bool genRandom)
{

  fstream outFile;

  char filename[100];
  if (genRandom)
  {
    sprintf(filename, "input%drand.txt", numParticle);

    outFile.open(filename, fstream::out);

    outFile << numParticle << endl;

    srand(time(NULL));
    double max=static_cast<double>(RAND_MAX);
    double range=sqrt(numParticle);
    double radius_range=10;

    for(int i=0; i<numParticle; i++)
    {
      outFile << i << "\t" << rand()*range*2/max << "\t" << rand()*range*2/max << "\t" << rand()*radius_range/max << endl;
    }
  }
  else
  {
    int size = static_cast<int>(sqrt(numParticle));
    sprintf(filename, "input%d.txt", numParticle);

    outFile.open(filename, fstream::out);

    outFile << numParticle << endl;

    for(int i=0; i<size; i++)
    {
      for(int j=0; j<size; j++)
      {
        outFile << i*size+j << "\t" << i << "\t" << j << "\t" << 3 << endl;
      }
    }
  }

  outFile.close();
}

GenInput::~GenInput()
{

}
