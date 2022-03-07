#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <unordered_set>
using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 7)
    {
        cout<<"Usage: <Points File> <Triangles File> <X size> <Y size> <Step size X> <Step size Y>"<<endl;
        return -1;
    }
    float pointHeight = 0;
    int x = atoi(argv[3]);
    int y = atoi(argv[4]);
    float stepX = atof(argv[5]);
    float stepY = atof(argv[6]);

    ofstream pointsFile;
    pointsFile.open(argv[1]);

    for(float j=0; j< y; j+=stepY)
    {
        for(float i=0; i< x; i+=stepX)
        {
            pointsFile<<i<<" "<<j<<" "<<pointHeight<<endl;
        }
    }
    pointsFile.close();

    
    ofstream trianglesFile;
    trianglesFile.open(argv[2]);
    cout<<"x/step = "<<(x/stepX)<<",  y/step = "<<(y/stepY)<<endl;

    int numPointsX = x/stepX;
    if(numPointsX*stepX < x)
        numPointsX += 1;
    cout<<"x: "<<x<<",  stepX: "<<stepX<<",  numX: "<<numPointsX<<"  -  ";
    cout<<"y: "<<y<<",  stepy: "<<stepY<<",  numY: "<<(y/stepY)<<" ";
    for(int j=0; j< (y/stepY)-1; j++)
    {
        for(int i=0; i< numPointsX-1; i++)
        {
            int a = j*numPointsX+i;
            //trianglesFile<<a<<" "<<a+1<<" "<<a+numPointsX<<endl;
            //trianglesFile<<a+1<<" "<<a+1+numPointsX<<" "<<a+numPointsX<<endl;

            trianglesFile<<a<<" "<<a+1+numPointsX<<" "<<a+numPointsX<<endl;
            trianglesFile<<a<<" "<<a+1<<" "<<a+numPointsX+1<<endl;
        }
    }
    trianglesFile.close();
}