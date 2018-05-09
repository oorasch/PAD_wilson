#include <gsl/gsl_sf_bessel.h>
#include "initialize.h"

const unsigned seed1(std::chrono::system_clock::now().time_since_epoch().count());
//std::default_random_engine generator(seed);
std::ranlux48 generator1(seed1);

using namespace std;


void neib_init (const int &N1, const int &N2, int neib[][4])
{
    int k,xm,xp,yp,ym;

    for(int x=0; x < N1; x++)
    {
        xp= x+1;
        xm= x-1;

        if (xp == N1) xp = 0;
        if (xm == -1) xm = N1 - 1;

        for(int y=0; y < N2; y++)
        {
            yp = y + 1;
            ym = y - 1;

            if (yp == N2) yp = 0;
            if (ym == -1) ym = N2 - 1;

            k = x + y*N1;

            neib[k][0]= x + yp*N1;
            neib[k][1]= xp + y*N1;
            neib[k][2]= x + ym*N1;
            neib[k][3]= xm + y*N1;
        }
    }

    return;
}


void sig_init (vector<vector<int>>& sig_link)
{
	std::uniform_int_distribution<int> distribution(0,1);

  for(int i = 0; i < constants::V; i++)
  {
    sig_link[i][0] = 2*distribution(generator1)-1;// -1 or 1
    sig_link[i][1] = 2*distribution(generator1)-1;
  }
	return;
}


void calc_M_array(const int neib [][4], const vector<vector<int>>& sig_link, const vector<int>& s_site, vector<double>& M_prime)
{
	//vector<double> M_prime (V*V, 0.0);
	
	for(int x = 0; x < constants::V; x++)
	{
		for(int y = 0; y < constants::V; y++)
		{
			int k = y + x*constants::V;
			
			M_prime[k] = 0.0;
			
			if(x == y)
			{
				M_prime[k] += 1.0 + (constants::m - 1.0)*s_site.at(x);
				
			}
			else
			{
				if(y == neib[x][0])
				{
					M_prime[k] += s_site.at(x)*s_site.at(y)*sig_link.at(x).at(0);
				}
			
				if(y == neib[x][1])
				{
					M_prime[k] += s_site.at(x)*s_site.at(y)*sig_link.at(x).at(1);
				}
			
				if(x == neib[y][0])
				{
					M_prime[k] -= s_site.at(x)*s_site.at(y)*sig_link.at(y).at(0);
				}
			
				if(x == neib[y][1])
				{
					M_prime[k] -= s_site.at(x)*s_site.at(y)*sig_link.at(y).at(1);
				}
				
				if(M_prime[k] != 0) M_prime[k] *= 0.5;
				
				//cout << M_prime[k] << endl;
			}
		}
		//cout << endl;
	}
	
	//cout << endl;
	
	return;
}

void bessel_init(vector<double>& vec)
{
	vec.clear();
	double aux = 2.0*sqrt(constants::eta*constants::eta_bar);
	for(int i = 0; i < constants::max_bessel; i++)
	{
		vec.push_back(gsl_sf_bessel_In(i, aux));
	}
}


