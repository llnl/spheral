#include "libs.hh"
#include "classes.hh"
#include "headers.hh"

#include <random>

namespace FractalSpace
{
  typedef deque<double>::iterator _ITD__;
  template <class ForwardIterator>
  void make_me_a_galaxy(Fractal_Memory* PFM,int numbers,double total_mass,double G,
			ForwardIterator posxb,ForwardIterator posyb,ForwardIterator poszb,
			ForwardIterator velxb,ForwardIterator velyb,ForwardIterator velzb)
  {
    // ofstream& FHT{PFM->p_file->DUMPS};
    const int FractalRank=PFM->p_mess->FractalRank;
    int seed=9973+256*FractalRank;
    std::mt19937 genRAN(seed);
    std::uniform_real_distribution<double> distRAD(1.0e-10,1.0);
    std::uniform_real_distribution<double> distTHETA(-0.99999999,0.99999999);
    std::uniform_real_distribution<double> distPHI(0.0,8.0*atan(1.0));
    std::normal_distribution<double> distNORM(0.0,1.0);
    const double rmax=30.0;
    const double x_off=-10.0;
    const double y_off=0.0;
    const double z_off=10.0;
    const double slope=-0.9;
    double velratio=0.7;
    double sigratio=0.002;
    velratio*=sqrt(G);
    sigratio*=sqrt(G);
    const double slope3=slope+3.0;
    const double expo=1.0/(3.0+slope);
    // const double mx=1.0e30;
    // const double my=1.0e30;
    // const double mz=1.0e30;
    // const double nx=-1.0e30;
    // const double ny=-1.0e30;
    // const double nz=-1.0e30;
    // const double st1=0.0;
    // const double sp1=0.0;
    for(int ni=0;ni<numbers;ni++)
      {
	double r=rmax*pow(distRAD(genRAN),expo);
	double ctheta=distTHETA(genRAN);
	// st1+=ctheta;
	double R=r*sqrt(1.0-ctheta*ctheta);
	*poszb=r*ctheta+z_off;
	double phi=distPHI(genRAN);
	// sp1+=phi;
	double cphi=cos(phi);
	double sphi=sin(phi);
	*posxb=R*cphi+x_off;
	*posyb=R*sphi+y_off;
	double sqmr=sqrt(total_mass*pow((r/rmax),slope3)/r);
	double vt=velratio*sqmr;
	double vs=sigratio*sqmr;
	*velxb=-sphi*vt+vs*distNORM(genRAN);
	*velyb=cphi*vt+vs*distNORM(genRAN);
	*velzb=vs*distNORM(genRAN);
	// mx=min(mx,*posxb);
	// my=min(my,*posyb);
	// mz=min(mz,*poszb);
	// nx=max(nx,*posxb);
	// ny=max(ny,*posyb);
	// nz=max(nz,*poszb);
	// FHT << "A GALAXY " << FractalRank << " " << ni << " " << *posxb << " " << *posyb << " " << *poszb;
	// FHT << " " << *velxb << " " << *velyb << " " << *velzb << "\n";
	posxb++;
	posyb++;
	poszb++;
	velxb++;
	velyb++;
	velzb++;
      }
    // FHT << " MINMAX " << FractalRank << " " << mx  << " " << my  << " " << mz  << " " << nx  << " " << ny  << " " << nz << endl;
    // FHT << " AVER " << FractalRank << " " << st1/numbers << " " << sp1/numbers << endl;
  }
}
namespace FractalSpace
{
  template 
  void make_me_a_galaxy(Fractal_Memory* PFM,int numbers,double total_mass,double G,
			_ITD__ posxb,_ITD__ posyb,_ITD__ poszb,
			_ITD__ velxb,_ITD__ velyb,_ITD__ velzb);
}
namespace FractalSpace
{
  void make_me_a_galaxy(int FractalRank,int numbers,double total_mass,double G,
			vector <double>& posx,vector <double>& posy,vector <double>& posz,
			vector <double>& velx,vector <double>& vely,vector <double>& velz)
  {
    int seed=9973+256*FractalRank;
    srand(seed);
    double rand_max=(double)RAND_MAX;
    double rmax=30.0;
    double x_off=-1.0;
    double y_off=1.5;
    double z_off=0.5;
    double slope=-0.9;
    double velratio=0.7;
    double sigratio=0.002;
    velratio*=sqrt(G);
    sigratio*=sqrt(G);
    double slope3=slope+3.0;
    double expo=1.0/(3.0+slope);
    double twopi=8.0*atan(1.0);
    bool allok=true;
    for(int ni=0;ni<numbers;ni++)
      {
	//	double r1=std::distribution(generator);
	double r1=Fractal::my_rand_not_zero(rand_max);
	r1=pow(r1,expo);
	//	double ctheta=2.0*std::distribution(generator)-1.0;
	double ctheta=0.99999*(2.0*Fractal::my_rand(rand_max)-1.0);
	//	double phi=twopi*std::distribution(generator);
	double phi=twopi*Fractal::my_rand(rand_max);
	double r=rmax*r1;
	posz[ni]=r*ctheta+z_off;
	double stheta=sqrt(1.0-ctheta*ctheta);
	double cphi=cos(phi);
	double sphi=sin(phi);
	posx[ni]=r*stheta*cphi+x_off;
	posy[ni]=r*stheta*sphi+y_off;
	double massr=total_mass*pow((r/rmax),slope3);
	double vt=velratio*sqrt(massr/r);
	double vs=sigratio*sqrt(massr/r);
	double v2=vs*sqrt(-2.0*log(Fractal::my_rand_not_zero(rand_max)));
	double ang=twopi*Fractal::my_rand(rand_max);
	velx[ni]=-sphi*vt+v2*cos(ang);
	vely[ni]=cphi*vt+v2*sin(ang);
	velz[ni]=vs*sqrt(-2.0*log(Fractal::my_rand_not_zero(rand_max)))*cos(twopi*Fractal::my_rand(rand_max));
	if(isfinite(posx[ni]) && isfinite(posy[ni]) && isfinite(posz[ni]) && isfinite(velx[ni]) && isfinite(vely[ni]) && isfinite(velz[ni]))
	  continue;
	cerr << "BAD GALAXY " << FractalRank << " " << ni << " " << posx[ni] << " " << posy[ni] << " " << posz[ni];
	cerr << " " << velx[ni] << " " << vely[ni] << " " << velz[ni];
	cerr << " " << r1 << " " << phi << " " << ctheta << " " << v2 << " " << ang << endl;
	allok=false;
      }
    assert(allok);
  }
}
