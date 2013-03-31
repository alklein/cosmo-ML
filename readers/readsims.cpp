#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

int Np = 128*128*128;

struct _particle
{
  float x, y, z;
  float vx, vy, vz;
  int id;
};

int readsim(string fname, _particle *pp)
{

  ifstream infile;
  infile.open( fname.c_str(), ios::binary|ios::in);

  for(int i=0;i<Np;i++)
  {
    pp[i].id = i;
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
  _particle *pp;
  pp = new _particle[Np];

  //  readsim("sim0000/xv_dm.init",pp);
  readsim("xv_dm.z=00.0000",pp);


  for (int i=0; i<Np; i++) {
    cout << pp[i].id << ' ' 
	 << pp[i].x << ' ' << pp[i].y << ' ' << pp[i].z << ' '
	 << pp[i].vx << ' ' << pp[i].vy << ' ' << pp[i].vz << endl;
  }
  /*
  cout << pp[0].id << ' ' 
       << pp[0].x << ' ' << pp[0].y << ' ' << pp[0].z << ' '
       << pp[0].vx << ' ' << pp[0].vy << ' ' << pp[0].vz << endl;
  cout << pp[Np-1].id << ' ' 
       << pp[Np-1].x << ' ' << pp[Np-1].y << ' ' << pp[Np-1].z << ' ' 
       << pp[Np-1].vx << ' '  << pp[Np-1].vy << ' '  << pp[Np-1].vz << endl;
  */

  return 0;
}
