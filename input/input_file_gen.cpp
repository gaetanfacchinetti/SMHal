#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>

int main()
{
    std::ofstream outfile;
    outfile.open("./my_input.dat");
    
    int N = 10000;
    
    double mA1_min = 0.01;
    double mA1_max = 1.5;
    
    double d_mA1 = (mA1_max - mA1_min)/N;
    double mA1;
    
    for(int i = 0; i<N+1; i++)
    {
        mA1 = mA1_min + i*d_mA1;
        outfile << "0.3 0.3 10 125 " << mA1 << " 1 10 10 0.3 1 0 0 1 1 0 0 1 0.03 0" << std::endl;
    }
    
    outfile.close();
    
    return 0;
}
