/*

File: parse_particles.cpp
Brief: Script to extract particles from binary simulation output file
Author: Andrea Klein       <alklein@alumni.stanford.edu>

Usage: g++ parse_particles.cpp -o parser
       ./parser > particles.txt
Note: the infile halo_part.z=00.0000 must be in the same directory.

*/

#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

int Np = 128*128*128;

struct _particle
{
  int ip;
  int ih;
  float x, y, z;
  float vx, vy, vz;
};

int readsim(string fname, _particle *pp)
{

  ifstream infile;
  infile.open(fname.c_str(), ios::binary|ios::in);
  double * dummy = new double[1];
  infile.read( (char *)&dummy[0], 8 ); // read first 8 bytes into dummy

  for(int i=0; i<Np; i++)
  {
    infile.read( (char *)&pp[i].ip, 8 );
    infile.read( (char *)&pp[i].ih, 8 );
    infile.read( (char *)&pp[i].x, sizeof(float) );
    infile.read( (char *)&pp[i].y, sizeof(float) );
    infile.read( (char *)&pp[i].z, sizeof(float) );
    infile.read( (char *)&pp[i].vx, sizeof(float) );
    infile.read( (char *)&pp[i].vy, sizeof(float) );
    infile.read( (char *)&pp[i].vz, sizeof(float) );
  }

  return 0;
}

int main()
{
  _particle *pp = new _particle[Np];
  readsim("halo_part.z=00.0000", pp);

  for (int i=0; i<Np; i++) {
    cout << pp[i].ip << ' ' 
         << pp[i].ih << ' '
	 << pp[i].x  << ' ' << pp[i].y << ' ' << pp[i].z << ' '
	 << pp[i].vx << ' ' << pp[i].vy << ' ' << pp[i].vz << endl;
  }

  return 0;
}
