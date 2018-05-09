#include <iostream>
#include <vector>
#include <cmath>
#include "constants.h"
#include "debug.h"

using namespace std;

void pauli_check(const int neib[][4], const vector<vector<int>>& k_link)
{
	for(int i = 0; i < constants::V; i++)
	{
		//check if loops touch
		if(abs(k_link[i][0]) > 1 || abs(k_link[i][1]) > 1) cout<< i << ": Uuuuuups ... there goes the Pauli principle :'(" << endl;
		
		//check if loops cross
		int count = 0;
		
		if(abs(k_link[i][0]) == 1)	count++;
		if(abs(k_link[i][1]) == 1)	count++;
		if(abs(k_link[neib[i][2]][0]) == 1)	count++;
		if(abs(k_link[neib[i][3]][1]) == 1)	count++;
		
		if(count > 2) cout << i << ": Uuuuuups ... there goes the Pauli principle :|" << endl;
	}
}

//check the fermion flux conservation
void flux_check(const int neib[][4], const vector<vector<int>>& k_link)
{
	for(int i = 0; i < constants::V; i++)
	{		
		if((k_link[i][0] + k_link[i][1] - k_link[neib[i][2]][0] - k_link[neib[i][3]][1]) != 0) cout << i << ": Uuuuuups ... there goes the flux conservation :0" << endl;
	}
}

//check if the plaquette occupation numbers add up correctly
void plaquette_check(const int neib[][4], const vector<vector<int>>& k_link, const vector<int>& plaq_occ)
{
	for(int i = 0; i < constants::V; i++)
	{
		if(plaq_occ.at(neib[i][3]) != plaq_occ.at(i) - k_link[i][0]) cout<< i << ": Uuuuuups ... plaquette occupation wrong!" << endl;
		if(plaq_occ.at(neib[i][2]) != plaq_occ.at(i) + k_link[i][1]) cout<< i << ": Uuuuuups ... plaquette occupation wrong!" << endl;
	}
}
