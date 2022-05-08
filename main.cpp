// C++ program to demonstrate the
// boilerplate code
#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include <unordered_map>
#include <math.h>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
#include <omp.h>

using namespace std;

//1 = No triangle march
//2 = Triangle march
#define VERSION 2

clock_t overall_time;
clock_t time_req;
double pi = 3.141592653589793238462643383279502884197169399375105820974944;
double root2 = sqrt(2);
int SYSTEM_NUM_BITS = 32;

struct Point;
struct LineSegmentKey
{
    long points[2];

    LineSegmentKey(long a, long b)
    {
        if(a < b)
        {
            points[0] = a;
            points[1] = b;
        }
        else
        {
            points[0] = b;
            points[1] = a;
        }
    }

    bool operator< (const LineSegmentKey & otherSegment) const
    {
        bool t = points[0]<otherSegment.points[0];
        if(t)
            return t;
        return points[1] < otherSegment.points[1];
    }

    bool operator== (const LineSegmentKey & otherSegment) const
    {
        //cout<<"\thacking:   "<<points[0]<<" - "<<points[1]<<"      vs       "<<otherSegment.points[0]<<" - "<<otherSegment.points[1]<<"\t";

        return (points[0] == otherSegment.points[0] && points[1] == otherSegment.points[1]);
    }

    bool operator!=(const LineSegmentKey & other) const
    {
        return (points[0] != other.points[0] || points[1] != other.points[1]);
    }

    struct HashFunction
    {
        size_t operator()(const LineSegmentKey& segment) const
        {
            return (segment.points[0]<<(SYSTEM_NUM_BITS/2))+segment.points[1];
        }
    };
};
struct LineSegmentData
{
    vector<long> triangles;
};

struct Triangle;
struct Square;
struct Circle;
struct Intersection;

//Data Variables
Point *points;
unordered_map<LineSegmentKey, LineSegmentData, LineSegmentKey::HashFunction> segments;
Triangle *triangles;
Circle *circles;
vector<long double> radii;
long numPoints;
long numSegments;
long numTriangles;
long numCircles;

//Square Algorithm Variables
Square **squares;
//double largestX;
//double largestY;
//double squareSize;
//long maxSquareX;
//long maxSquareY;
int numThreads = 1;

bool debug = true;
string debugSquaresFile = "DebugSquares.txt";
string debugIntersectionsFile = "DebugIntersections.txt";


struct Point
{
    long double x;
    long double y;
    long double h;

    bool operator==(const Point& otherPoint) const
    {
        if (this->x == otherPoint.x && this->y == otherPoint.y) return true;
        else return false;
    }
    
    bool operator<(const Point& otherPoint) const
    {
        if(x == otherPoint.x)
            return y < otherPoint.y;
        return x < otherPoint.x;
    }

    struct HashFunction
    {
        size_t operator()(const Point& point) const
        {
            size_t xHash = std::hash<int>()(point.x);
            size_t yHash = std::hash<int>()(point.y) << 1;
            return xHash ^ yHash;
        }
    };
};

struct Triangle
{
    long int points[3];
    long double center[2];
    //z = -mx - ny + h
    double m;
    double n;
    double h;

    vector<LineSegmentKey> GetSegments()
    {
        vector<LineSegmentKey> trigSegments;
        trigSegments.push_back(LineSegmentKey(points[0],points[1]));
        trigSegments.push_back(LineSegmentKey(points[0],points[2]));
        trigSegments.push_back(LineSegmentKey(points[1],points[2]));
        return trigSegments;
    }
};
struct Square
{
    vector<LineSegmentKey> intersectingSegments;
};
struct Circle
{
    //x, y, h = radius
    long double x;
    long double y;
    long double lineIntegral;
};
struct Intersection
{
    double theta;
    LineSegmentKey ls;
    Intersection(double t, LineSegmentKey segment): theta(t), ls(segment)
    {
    }

    bool operator< (const Intersection & other) const
    {
        return theta < other.theta;
    }

    bool operator==(const Intersection & other) const
    {
        return other.theta == this->theta;
    }
};

//Debugging Methods
void PrintPoint(long a)
{
    cout<<"("<<points[a].x<<", "<<points[a].y<<")";
}
void PrintSegment(LineSegmentKey s)
{
    //cout<<"Hash: "<<LineSegmentKey::HashFunction().operator(s);
    PrintPoint(s.points[0]);
    cout<<" - ";
    PrintPoint(s.points[1]);
}
void PrintSegmentSimple(LineSegmentKey s)
{
    //cout<<"Hash: "<<LineSegmentKey::HashFunction().operator(s);
    cout<<"("<<s.points[0]<<" - "<<s.points[1]<<")";
}
/*void PrintSquares()
{
    ofstream debugFile;
    if(debug)
    {
        debugFile.open(debugSquaresFile);
        debugFile<<maxSquareX<<" "<<maxSquareY<<endl;
    }
    std::cout<<endl<<"Printing Squares: ["<<maxSquareX<<", "<<maxSquareY<<"]"<<endl;
    for(int j=0; j<maxSquareY; j++)
    {
        std::cout<<"      - y: "<<j<<endl;
        for(int i=0; i<maxSquareX; i++)
        {
            std::cout<<"   x: "<<i<<endl;
            if(debug)
            {
                debugFile<<squares[i][j].intersectingSegments.size()<<endl;
            }
            for(auto k = squares[i][j].intersectingSegments.begin(); k!= squares[i][j].intersectingSegments.end(); k++)
            {
                std::cout<<"            - seg: ("<<points[k->points[0]].x<<", "<<points[k->points[0]].y<<") - ("<<points[k->points[1]].x<<", "<<points[k->points[1]].y<<")"<<endl;
                if(debug)
                {
                    debugFile<<points[k->points[0]].x<<" "<<points[k->points[0]].y<<" "<<points[k->points[1]].x<<" "<<points[k->points[1]].y<<endl;
                }
            }
        }
    }
}*/
void PrintTriangles()
{
    cout<<"\nPrinting "<<numTriangles<<" Trianges:"<<endl;
    for(long i=0; i<numTriangles; i++)
    {
        cout<<"\tTriangle "<<&triangles[i]<<endl<<"\t\t";
        PrintSegment(LineSegmentKey(triangles[i].points[0], triangles[i].points[1]));
        cout<<"  :  ";
        PrintSegment(LineSegmentKey(triangles[i].points[0], triangles[i].points[2]));
        cout<<"  :  ";
        PrintSegment(LineSegmentKey(triangles[i].points[1], triangles[i].points[2]));
        cout<<endl;
    }
}

long FileNumLines(string fileName)
{
    time_req = clock();
    string line;
    long count = 0;
    ifstream readFile(fileName);
    if(readFile.is_open())
    {
        while(getline(readFile,line)){
            count++;
        }
    }
    else
    {
        std::cout<<"File Not Open"<<endl;
        exit(EXIT_FAILURE);
    }
	time_req = clock() - time_req;
    cout<<"Checking File Size took "<<(float)time_req/CLOCKS_PER_SEC<<" seconds"<<endl;
    return count;
}
void ReadPoints(char *pointsFile)
{
    ifstream readFile(pointsFile);
    if(readFile)
    {
        numPoints = FileNumLines(pointsFile);
        //std::cout<<"Lines in points file: "<<numPoints<<endl;
        points = new Point[numPoints];
        double input;

        for(long i = 0; i < numPoints; i++)
        {
            readFile >> input;
            points[i].x = input;

            readFile >> input;
            points[i].y = input;

            readFile >> input;
            points[i].h = input;

            //std::cout<<"("<<points[i].x<<", "<<points[i].y<<", "<<points[i].h<<")"<<endl;
        }
        std::cout<<endl;
    }
    else
    {
        std::cout<<"Could Not Find Points File: "<<pointsFile<<endl;
        exit(EXIT_FAILURE);
    }
}
void ReadTriangles(char *triangleFile)
{
    cout<<"READING TRIANGLES: "<<endl;
    ifstream readFile(triangleFile);
    if(readFile)
    {
        numTriangles = FileNumLines(triangleFile);
        cout<<"NumTriangles: "<<numTriangles<<endl;
        triangles = new Triangle[numTriangles];
        double input;
        numSegments = 0;
        //squareSize = 0;

        time_req = clock();
        for(long long i = 0; i < numTriangles; i++)
        {
            long ind1, ind2, ind3;
            readFile >> ind1 >> ind2 >> ind3;
            Point a = points[ind1];
            Point b = points[ind2];
            Point c = points[ind3];
            LineSegmentKey s1 = LineSegmentKey(ind1,ind2);
            LineSegmentKey s2 = LineSegmentKey(ind1,ind3);
            LineSegmentKey s3 = LineSegmentKey(ind2,ind3);
            
            triangles[i].center[0] = (a.x+b.x+c.x)/3;
            triangles[i].center[1] = (a.y+b.y+c.y)/3;

            //Calculate equation for plane
            double AB[3];
            double AC[3];

            AB[0] = b.x-a.x;
            AB[1] = b.y-a.y;
            AB[2] = b.h-a.h;

            AC[0] = c.x-a.x;
            AC[1] = c.y-a.y;
            AC[2] = c.h-a.h;

            double normalVector[3];
            normalVector[0] = AB[1]*AC[2] - AB[2]*AC[1];
            normalVector[1] = AB[0]*AC[2] - AB[2]*AC[0];
            normalVector[2] = AB[0]*AC[1] - AB[1]*AC[0];

            double m = normalVector[0]/normalVector[2];
            double n = normalVector[1]/normalVector[2];
            double h = a.x*m + a.y*n + a.h;

            triangles[i].m = m;
            triangles[i].n = n;
            triangles[i].h = h;

            segments[s1].triangles.push_back(i);
            segments[s2].triangles.push_back(i);
            segments[s3].triangles.push_back(i);

            triangles[i].points[0] = ind1;
            triangles[i].points[1] = ind2;
            triangles[i].points[2] = ind3;
        }
        //SEGMENTS ARE GOING TO BE A DOUBLE KEY MAP

	    time_req = clock() - time_req;
        cout<<"Reading Triangles took "<<(float)time_req/CLOCKS_PER_SEC<<" seconds"<<endl;
    }
}
void ReadCircles(char *circleFile)
{
    ifstream readFile(circleFile);
    if(readFile)
    {
        numCircles = FileNumLines(circleFile) - 1;
        circles = new Circle[numCircles];
        long double radiStart;
        long double radiEnd;
        double radiStep;
        readFile >> radiStart >> radiEnd >> radiStep;

        for(long double i = radiStart; i<=radiEnd; i += radiStep)
        {
            //cout<<"Reading Circle Radii: "<<i<<endl;
            radii.push_back(i);
        }

        double input;

        for(long i = 0; i < numCircles; i++)
        {
            readFile >> input;
            circles[i].x = input;

            readFile >> input;
            circles[i].y = input;
        }
        std::cout<<endl;
    }
    else
    {
        std::cout<<"Could Not Find Circles File: "<<circleFile<<endl;
        exit(EXIT_FAILURE);
    }
}

vector<Intersection> CircleLineIntersection(LineSegmentKey ls, long int circleIndex, long int radiiIndex)
{
    //PrintSegment(ls);
    vector<Intersection> intersections;
    double AB[2];
    double AC[2];
    double IC[2];

    AB[0] = points[ls.points[1]].x - points[ls.points[0]].x;
    AB[1] = points[ls.points[1]].y - points[ls.points[0]].y;

    AC[0] = points[ls.points[0]].x-circles[circleIndex].x;
    AC[1] = points[ls.points[0]].y-circles[circleIndex].y;

    double a = AB[0]*AB[0] + AB[1]*AB[1];
    double bHalf = (AC[0]*AB[0] + AC[1]*AB[1]);
    double c = AC[0]*AC[0] + AC[1]*AC[1] - radii[radiiIndex]*radii[radiiIndex];

    double discriminant = bHalf*bHalf - a*c;
    if(discriminant < 0)
    {
        //cout<<"\t == No Hit"<<endl;
        return intersections;
    }
    else if(discriminant == 0)
    {
        //cout<<"\t== Single hit"<<endl;
        double n = -bHalf/a;
        //cout<<"\t\t";
        //PrintSegment(ls);
        //cout<<" - n: "<<n;
        if(n>=0&&n<=1)
        {
            IC[0] = AC[0] + n*AB[0];
            IC[1] = AC[1] + n*AB[1];

            double theta = asin(IC[1]/radii[radiiIndex]);
            if(IC[0]<0)
                theta = pi-theta;
            if(theta < 0)
                theta += 2*pi;

            intersections.push_back(Intersection(theta, ls));
        }
        return intersections;
    }
    else
    {
        //cout<<"\t== Double hit"<<endl;
        double n = (-bHalf + sqrt(discriminant))/a;
        //cout<<"\t\t";
        //PrintSegment(ls);
        //cout<<" - n: "<<n;
        if(n>=0&&n<=1)
        {
            IC[0] = AC[0] + n*AB[0];
            IC[1] = AC[1] + n*AB[1];
            //double theta = atan(IC[0]/IC[1]);
            double numBeforeAsin = IC[1]/radii[radiiIndex];
            if(numBeforeAsin > 1)
                numBeforeAsin = 1;
            if(numBeforeAsin < -1)
                numBeforeAsin = -1;
            double theta = asin(numBeforeAsin);
            //cout<<"NumBeforeAsin: "<<IC[1]/radii[radiiIndex]<<",   Theta: "<<theta<<endl;
            if(IC[0]<0)
                theta = pi-theta;
            if(theta < 0)
                theta += 2*pi;

            intersections.push_back(Intersection(theta, ls));
        }
        n = (-bHalf - sqrt(discriminant))/a;
        //cout<<"\t\t";
        //PrintSegment(ls);
        //cout<<" - n: "<<n;
        if(n>=0&&n<=1)
        {
            IC[0] = AC[0] + n*AB[0];
            IC[1] = AC[1] + n*AB[1];
            //theta = atan(IC[0]/IC[1]);
            double numBeforeAsin = IC[1]/radii[radiiIndex];
            if(numBeforeAsin > 1)
                numBeforeAsin = 1;
            if(numBeforeAsin < -1)
                numBeforeAsin = -1;
            double theta = asin(numBeforeAsin);
            //cout<<"NumBeforeAsin: "<<IC[1]/radii[radiiIndex]<<",   Theta: "<<theta<<endl;
            if(IC[0]<0)
                theta = pi-theta;
            if(theta < 0)
                theta += 2*pi;

            intersections.push_back(Intersection(theta, ls));
        }
        //cout<<endl;
        return intersections;
    }
}

ofstream debugIntersectionsStream;
void SingleCircleIntegral(long int circleIndex, long int radiiIndex)
{
        /*
         * Optimizations:
         *  - Save intersection values of lines already calculated = x2 speed
        */
    vector<Intersection> intersections;
    unordered_set<LineSegmentKey, LineSegmentKey::HashFunction> checkedSegments;
    vector<LineSegmentKey> queue;
    vector<LineSegmentKey> nextQueue;
#if VERSION == 1
    //std::cout<<"SingleCircleIntegral: <"<<circles[circleIndex].x<<", "<<circles[circleIndex].y<<"> : "<<radii[radiiIndex]<<endl;
    for(auto [key, value] : segments)
    {
        //cout<<"\t";
        //PrintSegment(key);
        //cout<<endl;
        vector<Intersection> temp = CircleLineIntersection(key, circleIndex, radiiIndex);
        if(!temp.empty())
            //cout<<"\t - Hit!"<<endl;
            intersections.insert(intersections.end(), temp.begin(), temp.end());
    }
#else
    //std::cout<<"SingleCircleIntegral: <"<<circles[circleIndex].x<<", "<<circles[circleIndex].y<<"> : "<<radii[radiiIndex]<<endl;
    for(auto [key, value] : segments)
    {
        //cout<<"\t";
        //PrintSegment(key);
        //cout<<endl;
        vector<Intersection> temp = CircleLineIntersection(key, circleIndex, radiiIndex);
        if(!temp.empty())
        {
            //cout<<"\t - Hit!"<<endl;
            checkedSegments.insert(key);
            for(long int trig : value.triangles)
            {
                Triangle t = triangles[trig];
                LineSegmentKey ls1 = LineSegmentKey(t.points[0], t.points[1]);
                LineSegmentKey ls2 = LineSegmentKey(t.points[0], t.points[2]);
                LineSegmentKey ls3 = LineSegmentKey(t.points[1], t.points[2]);
                if(checkedSegments.find(ls1) != checkedSegments.end())
                {
                    checkedSegments.insert(ls1);
                    queue.push_back(ls1);
                }
                if(checkedSegments.find(ls2) != checkedSegments.end())
                {
                    checkedSegments.insert(ls2);
                    queue.push_back(ls2);
                }
                if(checkedSegments.find(ls3) != checkedSegments.end())
                {
                    checkedSegments.insert(ls3);
                    queue.push_back(ls3);
                }
            }
            intersections.insert(intersections.end(), temp.begin(), temp.end());
            break;
        }
    }
    long int size = queue.size();
    for(long int i=0; i<size; i++)
    {
        vector<Intersection> temp = CircleLineIntersection(queue[i], circleIndex, radiiIndex);
        if(!temp.empty())
        {
            //cout<<"\t - Hit!"<<endl;
            for(long int trig : segments[queue[i]].triangles)
            {
                Triangle t = triangles[trig];
                LineSegmentKey ls1 = LineSegmentKey(t.points[0], t.points[1]);
                LineSegmentKey ls2 = LineSegmentKey(t.points[0], t.points[2]);
                LineSegmentKey ls3 = LineSegmentKey(t.points[1], t.points[2]);
                if(checkedSegments.find(ls1) == checkedSegments.end())
                {
                    checkedSegments.insert(ls1);
                    queue.push_back(ls1);
                    size++;
                }
                if(checkedSegments.find(ls2) == checkedSegments.end())
                {
                    checkedSegments.insert(ls2);
                    queue.push_back(ls2);
                    size++;
                }
                if(checkedSegments.find(ls3) == checkedSegments.end())
                {
                    checkedSegments.insert(ls3);
                    queue.push_back(ls3);
                    size++;
                }
            }
            intersections.insert(intersections.end(), temp.begin(), temp.end());
        }
    } 
#endif
    sort(intersections.begin(), intersections.end());
    //Intersections should be finished Calculating, print them here:
    if(debug)
    {
    #pragma omp critical
        debugIntersectionsStream<<circles[circleIndex].x<<" "<<circles[circleIndex].y<<" "<<radii[radiiIndex]<<" ";
    }
    //cout<<"Intersections for circle at <"<<circles[circleIndex].x<<", "<<circles[circleIndex].y<<"> with r: "<<radii[radiiIndex]<<"\n";
    
    //Calculate Integral
    long double integral = 0;
    long int numIntersections = intersections.size();
    for(long int i=0; i<numIntersections; i++)
    {
        long int j = i+1;
        double theta1 = intersections[i].theta;
        double theta2;
        if(i == numIntersections-1)
        {
            j = 0;
            theta2 = intersections[j].theta + 2*pi;
        }
        else
        {
            theta2 = intersections[j].theta;
        }
        //cout<<"i: "<<i<<"   j: "<<j<<"\t\t"<<theta1<<" "<<theta2<<endl;
        long int integralTriangle = -1;
        //Normal Case
        if(intersections[i].ls != intersections[j].ls)
        {
            LineSegmentData ls1 = segments[intersections[i].ls];
            LineSegmentData ls2 = segments[intersections[j].ls];
            //cout<<"before set intersection"<<endl;

            for(long int a : ls1.triangles)
            {
                for(long int b : ls2.triangles)
                {
                    //cout<<"a: "<<a<<", b: "<<b<<endl;
                    if(a == b)
                    {
                        integralTriangle = a;
                        goto foundTriangle;
                    }
                }
            }
        }
        //Case where both intersections are on the same line segment, need extra logic to find which triangle to use
        else
        {
            integralTriangle = segments[intersections[i].ls].triangles[0];
            //Only one available triangle
            if(segments[intersections[i].ls].triangles.size() == 1)
            {
                goto foundTriangle;
            }

            //Two available triangles
            double A[2];
            double B[2];
            A[0] = cos(intersections[i].theta);
            A[1] = sin(intersections[i].theta);
            B[0] = cos(intersections[j].theta);
            B[1] = sin(intersections[j].theta);

            //Left > 0, on line = 0, Right < 0
            double side = (B[0] - A[0]) * (triangles[integralTriangle].center[1] - A[1])
                            - (B[1] - A[1]) * (triangles[integralTriangle].center[0] - A[0]);
            if(side > 0)
                integralTriangle = segments[intersections[i].ls].triangles[1];
        }
foundTriangle:
        if(integralTriangle != -1)
        {
            //cout<<"Common Triangle: "<<commonTriangle<<endl;
            double arcIntegral = - radii[radiiIndex]*triangles[integralTriangle].m*(sin(theta2) - sin(theta1)) 
                                    + radii[radiiIndex]*triangles[integralTriangle].n * (cos(theta2) - cos(theta1))
                                    + triangles[integralTriangle].h * (theta2 - theta1);
            integral+=arcIntegral;
        }
        //cout<<"\tt1: "<<theta1<<", t2:"<<theta2<<"\t    arcInt: "<<arcIntegral<<",  integral:"<<integral<<endl;

        #pragma omp critical
            debugIntersectionsStream<<intersections[i].theta<<" ";
    }
    integral *= radii[radiiIndex];
    #pragma omp critical
        cout<<"RADII: "<<radii[radiiIndex]<<"\t\tINTEGRAL VALUE: "<<integral<<endl;
    debugIntersectionsStream<<endl;
    intersections.clear();
}

void CalculateAllCircleIntegrals()
{
    if(debug)
        debugIntersectionsStream.open(debugIntersectionsFile);

    #pragma omp parallel for num_threads(numThreads)
        for(long int i=0; i<numCircles; i++)
            for(long int radiiIndex = 0; radiiIndex < radii.size(); radiiIndex++)
                SingleCircleIntegral(i, radiiIndex);

    if(debug)
        debugIntersectionsStream.close();
}

// Driver Code
int main(int argc, char *argv[])
{
    if(argc != 5)
    {
        std::cout<<"Usage: <Points.txt> <Trianges.txt> <Circles.txt> <NumThreads>"<<endl;
        return 0;
    }
    overall_time = clock();
    cout<<"Reading In Points\n";
    ReadPoints(argv[1]);

    cout<<"Reading in Triangles\n";
    ReadTriangles(argv[2]);
    cout<<"NumSegments = "<<segments.size()<<endl;

    //for(int i=0; i<numPoints; i++)
    //{
    //    cout<<"\nPrinting Point: "<<&points[i]<<",    ";
    //    PrintPoint(i);
    //}
    cout<<"\n"<<endl;
    /*for(auto [key, value] : segments)
    {
        cout<<"\nPrinting Segment: ";
        PrintSegmentSimple(key);
        cout<<" - Trigs: ";
        
        cout<<" - count: "<<segments.count(key)<<"  \t||   ";
        for(int t:segments[key].triangles)
            cout<<t<<", ";
    }*/
    //PrintTriangles();
    cout<<"Reading in Circles\n";
    ReadCircles(argv[3]);
    //cout<<"Starting Squarify"<<endl;

    //DivideSquares();
    //PrintSquares();
    numThreads = atoi(argv[4]);

    cout<<"Starting Circle Calculations"<<endl;
    CalculateAllCircleIntegrals();
	
    overall_time = clock() - overall_time;
    cout<<"Program execution took "<<(float)overall_time/CLOCKS_PER_SEC<<" seconds.\t\tPoints: "<<numPoints<<",  Triangles: "<<numTriangles<<",   Circles:  "<<numCircles*radii.size()<<endl;

    return 0;
}

/*
    ToDo:   Triangle March Optimization
            Edge case for small triangles: 
                - Find a way to deal with triangles that are too small to intersect a triangle
                 

                - DONE: NEEDS TO BE TESTED!!!!!!
                        -When two intersection on the same line right next to each other, the program will need to check if 
                        line segment is the same for both intersections, then it will find the triangle whose center is 
                        on the same side of the segment as the arc

    Bugs:   Certain circles not having intersections calculated
Circle File = 
0.2 0.7 0.2
0 1.5
0 2

*/
