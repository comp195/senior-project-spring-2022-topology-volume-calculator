#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <unordered_set>
#include <random>
#include <chrono>
using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 8)
    {
        cout<<"Usage: <Points File> <Triangles File> <X size> <Y size> <Step size X> <Step size Y> <Noise>"<<endl;
        return -1;
    }
    float pointHeight = 1;
    long long int x = atoi(argv[3]);
    long long int y = atoi(argv[4]);
    float stepX = atof(argv[5]);
    float stepY = atof(argv[6]);
    long long int numPointsX = x/stepX;
    if(numPointsX*stepX < x)
        numPointsX += 1;
    long long int numPointsY = y/stepY;
    if(numPointsY*stepY < y)
        numPointsY += 1;

    ofstream pointsFile;
    pointsFile.open(argv[1]);
    float noisePower = atof(argv[7]);
    default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    uniform_real_distribution<float> distribution(-noisePower,noisePower);

    for(long double j=0; j< numPointsY; j++)
    {
        for(long double i=0; i< numPointsX; i++)
        {
            long double x = i*stepX+distribution(generator);
            long double y = j*stepY+distribution(generator);
            pointsFile<<x<<" "<<y<<" "<<pointHeight<<endl;
        }
    }
    pointsFile.close();

    
    ofstream trianglesFile;
    trianglesFile.open(argv[2]);
    cout<<"x/step = "<<(x/stepX)<<",  y/step = "<<(y/stepY)<<endl;

    cout<<"x: "<<x<<",  stepX: "<<stepX<<",  numX: "<<numPointsX<<"  -  ";
    cout<<"y: "<<y<<",  stepy: "<<stepY<<",  numY: "<<(y/stepY)<<" ";
    for(long long j=0; j< numPointsY-1; j++)
    {
        for(long long i=0; i< numPointsX-1; i++)
        {
            long long a = j*numPointsX+i;
            //trianglesFile<<a<<" "<<a+1<<" "<<a+numPointsX<<endl;
            //trianglesFile<<a+1<<" "<<a+1+numPointsX<<" "<<a+numPointsX<<endl;

            trianglesFile<<a<<" "<<a+1+numPointsX<<" "<<a+numPointsX<<endl;
            trianglesFile<<a<<" "<<a+1<<" "<<a+numPointsX+1<<endl;
        }
    }
    trianglesFile.close();
}