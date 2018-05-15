#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <vector>


using namespace std;

int main (int argc, char* argv[])
{
    cout << "Hello world - this is a statistics program" << endl << endl;
    cout << "Currently the method is: Blocking" << endl;

    const int Nt = 4, Ns = 4;
    const int nmeas = 1e4;
    const double pi(3.14159265359);


		string VERSION("_full1_newexp");
//		string VERSION("_gpdecs1");
//		string PARAM("_theta");austria
		string PARAM("_beta");
		string OBS(argv[1]);
		
		string FILE_NAME(OBS+PARAM+VERSION+".dat");

    vector<double> data1;

//    stringstream s1, s2, s3, s4;
//    s1 << Nt;
//    s2 << Ns;
//    s3 << eta;
//    s4 << lambda;
//    string s_Nt = s1.str(), s_Ns = s2.str(), s_eta = s3.str(), s_lambda = s4.str();
    //ifstream infile ("phi_sq("+s_Nt+","+s_Ns+","+s_eta+","+s_lambda+").dat");
    ifstream infile1;
    infile1.open(FILE_NAME);

    if(infile1.is_open())
    {
        double tmp;
        while(!infile1.eof())
        {
            infile1 >> tmp;
            data1.push_back(tmp);
        }

        data1.pop_back();

        if((data1.size()%nmeas)) cout << "Something went wrong!!!" << endl;

        cout << data1.size() << " " << nmeas << endl;

    }
    else
    {
        cout << "Unable to open file";
        abort();
    }
    infile1.close();
    
    ofstream outfile1;
    outfile1.open(OBS+PARAM+VERSION+"_stats.dat");

    
    double th = 0.1;
    for(int i = 0; i < data1.size(); i+= nmeas)
		{
			//calc mean
    	double avg1 = 0.0;
			for(int j = i; j < i+nmeas; j++)
			{
				avg1 += data1.at(j);
			}
			
			//cout << avg1 << " " << avg2 << endl;
			
			avg1 /= nmeas;
			
			//cout << avg1 << " " << avg2 << endl;
			
			double sigma1 = 0.0;
			for(int j = i; j < i+nmeas; j++)
			{
				sigma1 += (avg1 - data1.at(j))*(avg1 - data1.at(j));
			}
			//cout << sigma1 << " " << sigma2 << endl;
			sigma1 /= nmeas - 1.0;
			sigma1 = sqrt(sigma1);
			//cout << sigma1 << " " << sigma2 << endl;
		
			outfile1 << th << " " << avg1 << " " << sigma1 << endl;
			th += 0.5;
    }
    outfile1.close();
    
    

//    for(auto it: data) avg += it;

//    avg /= data.size();

//    cout << avg << endl;

//    //calc statistical error with blocking

//    double sigma = 0.0;
//    int binNum = 2;
//    int binSize = data.size()/binNum;
//    vector<double> block_avg;

//    cout << binSize << endl;

//    int n, b;
//    for(n = 0; n < binNum; n++)
//    {
//        double block_avg_tmp = 0.0;

//        for(b = binSize*n; b < binSize*(n+1); b++)
//        {
//            block_avg_tmp += data.at(b);
//        }

//        block_avg.push_back(block_avg_tmp/(1.0*binSize));
//    }

//    for(auto it: block_avg) sigma += (avg - it)*(avg - it);

//    sigma /= 1.0*(binNum - 1.0);
//    sigma = sqrt(sigma/(1.0*binNum));

//    cout << sigma << endl;



    return 0;
}
