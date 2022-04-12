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

clock_t time_req;
double pi = 3.141592653589793238462643383279502884197169399375105820974944;
double root2 = sqrt(2);
int SYSTEM_NUM_BITS = 32;

struct Point;
struct LineSegmentKey
{
    long points[2];

    LineSegmentKey(int a, int b)
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
    vector<int> triangles;
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
vector<double> radii;
long numPoints;
long numSegments;
long numTriangles;
long numCircles;

//Square Algorithm Variables
Square **squares;
double largestX;
double largestY;
//double squareSize;
long maxSquareX;
long maxSquareY;
int numThreads = 1;

bool debug = true;
string debugSquaresFile = "DebugSquares.txt";
string debugIntersectionsFile = "DebugIntersections.txt";


struct Point
{
    double x;
    double y;
    double h;

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
    int points[3];
    //z = -mx - ny + h
    double m;
    double n;
    double h;
};
struct Square
{
    vector<LineSegmentKey> intersectingSegments;
};
struct Circle
{
    //x, y, h = radius
    double x;
    double y;
    double lineIntegral;
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

    bool operator==(const Intersection& other) const
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
void PrintSquares()
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
}
void PrintTriangles()
{
    cout<<"\nPrinting "<<numTriangles<<" Trianges:"<<endl;
    for(int i=0; i<numTriangles; i++)
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

int FileNumLines(string fileName)
{
    time_req = clock();
    string line;
    int count = 0;
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
    largestX = 0;
    largestY = 0;
    ifstream readFile(pointsFile);
    if(readFile)
    {
        numPoints = FileNumLines(pointsFile);
        //std::cout<<"Lines in points file: "<<numPoints<<endl;
        points = new Point[numPoints];
        double input;

        for(int i = 0; i < numPoints; i++)
        {
            readFile >> input;
            points[i].x = input;
            if(input > largestX)
                largestX = input;

            readFile >> input;
            points[i].y = input;
            if(input > largestY)
                largestY = input;

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
        triangles = new Triangle[numTriangles];
        double input;
        numSegments = 0;
        //squareSize = 0;

        time_req = clock();
        for(int i = 0; i < numTriangles; i++)
        {
            int ind1, ind2, ind3;
            readFile >> ind1 >> ind2 >> ind3;
            Point a = points[ind1];
            Point b = points[ind2];
            Point c = points[ind3];
            LineSegmentKey s1 = LineSegmentKey(ind1,ind2);
            LineSegmentKey s2 = LineSegmentKey(ind1,ind3);
            LineSegmentKey s3 = LineSegmentKey(ind2,ind3);

            //cout<<"\tTriangle "<<i<<endl<<"\t\t";
            //PrintSegment(s1);
           // cout<<"  :  ";
            //PrintSegment(s2);
            //cout<<"  :  ";
            //PrintSegment(s3);
            //cout<<endl;

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

            //cout<<"<"<<AB[0]<<", "<<AB[1]<<", "<<AB[2]<<">     ";
            //cout<<"<"<<AC[0]<<", "<<AC[1]<<", "<<AC[2]<<">     ";
            //cout<<"<"<<normalVector[0]<<", "<<normalVector[1]<<", "<<normalVector[2]<<">"<<endl;

            double m = normalVector[0]/normalVector[2];
            double n = normalVector[1]/normalVector[2];
            double h = a.x*m + a.y*n + a.h;
            cout<<"m: "<<m<<"   n: "<<n<<"   h: "<<h<<endl;
            triangles[i].m = m;
            triangles[i].n = n;
            triangles[i].h = h;


            //Determine line segment x and y length, for use in finding squareSize - square size will be largest segment length size
            /*double xLen1 = abs(points[ind1].x - points[ind2].x);            
            double xLen2 = abs(points[ind1].x - points[ind3].x); 
            double xLen3 = abs(points[ind2].x - points[ind3].x); 
            double xLen = max(max(xLen1, xLen2), xLen3);

            double yLen1 = abs(points[ind1].y - points[ind2].y);
            double yLen2 = abs(points[ind1].y - points[ind3].y);
            double yLen3 = abs(points[ind2].y - points[ind3].y);
            double yLen = max(max(yLen1, yLen2), yLen3);

            squareSize = max(max(xLen, yLen), squareSize);*/

            /*cout<<"\t"<<i<<":\n";
            PrintSegmentSimple(s1);
            cout<<" - count: "<<segments.count(s1)<<endl;
            PrintSegmentSimple(s2);
            cout<<" - count: "<<segments.count(s2)<<endl;
            PrintSegmentSimple(s3);
            cout<<" - count: "<<segments.count(s3)<<endl;

            for(auto [key, value] : segments)
            {
                cout<<"\n\tPrinting Segment: ";
                PrintSegmentSimple(key);
                cout<<" - Trigs: ";
                
                cout<<" - count: "<<segments.count(key)<<"  \t||   ";
                for(int t:segments[key].triangles)
                    cout<<t<<", ";
            }
            cout<<"\ns1: --------------------------------------";*/
            
            segments[s1].triangles.push_back(i);

            
            /*for(auto [key, value] : segments)
            {
                cout<<"\n\tPrinting Segment: ";
                PrintSegmentSimple(key);
                cout<<" - Trigs: ";
                
                cout<<" - count: "<<segments.count(key)<<"  \t||   ";
                for(int t:segments[key].triangles)
                    cout<<t<<", ";
            }
            cout<<"\ns2: --------------------------------------";*/

            segments[s2].triangles.push_back(i);

            
            /*for(auto [key, value] : segments)
            {
                cout<<"\n\tPrinting Segment: ";
                PrintSegmentSimple(key);
                cout<<" - Trigs: ";
                
                cout<<" - count: "<<segments.count(key)<<"  \t||   ";
                for(int t:segments[key].triangles)
                    cout<<t<<", ";
            }
            cout<<"\ns3: --------------------------------------";*/

            segments[s3].triangles.push_back(i);

            
            /*for(auto [key, value] : segments)
            {
                cout<<"\n\tPrinting Segment: ";
                PrintSegmentSimple(key);
                cout<<" - Trigs: ";
                
                cout<<" - count: "<<segments.count(key)<<"  \t||   ";
                for(int t:segments[key].triangles)
                    cout<<t<<", ";
            }
            cout<<"\n----------------------------------------";*/


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
        double radiStart;
        double radiEnd;
        double radiStep;
        readFile >> radiStart >> radiEnd >> radiStep;

        for(double i = radiStart; i<=radiEnd; i += radiStep)
        {
            cout<<"Reading Circle Radii: "<<i<<endl;
            radii.push_back(i);
        }

        double input;

        for(int i = 0; i < numCircles; i++)
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

/*void AddSegmentToSquares(LineSegmentKey ls)
{
    //Line segment in square coordinates
    double x = points[ls.points[0]].x / squareSize;
    double y = points[ls.points[0]].y / squareSize;
    double targetX = points[ls.points[1]].x / squareSize;
    double targetY = points[ls.points[1]].y / squareSize;

    //Delta x,y
    double dX = targetX - x;
    double dY = targetY - y;
    //Is x increasing or decreasing
    bool dxPositive = dX>0? true : false;
    bool dyPositive = dY>0? true : false;
    //std::cout<<"\tdX: "<<dX<<",  dY: "<<dY<<endl;

    if(dY == 0)
    {
        for(int i = x; dxPositive? i<targetX : i>targetX; dxPositive? i++ : i--)
        {
            squares[i][(int)y].intersectingSegments.push_back(ls);
            //std::cout<<"\t\tADDED TO SQUARE: ["<<i<<", "<<(int)y<<"] - <"<<(int)i*squareSize<<", "<<(int)y*squareSize<<">"<<endl;
        }
        //squares[(int)targetX][(int)targetY].intersectingSegments.insert(ls);
    }
    else if(dX == 0)
    {
        for(int j = y; dyPositive? j<targetY : j>targetY; dyPositive? j++ : j--)
        {
            squares[(int)x][j].intersectingSegments.push_back(ls);
            //std::cout<<"\t\tADDED TO SQUARE: ["<<(int)x<<", "<<j<<"] - <"<<(int)x*squareSize<<", "<<(int)j*squareSize<<">"<<endl;
        }
        //squares[(int)targetX][(int)targetY].intersectingSegments.insert(ls);
    }
    else
    {
        double slope = dY/dX;
        double b = y - slope*x;
        //cout<<"SLOPE = "<<slope<<endl;

        //cout<<"\tslope: "<<slope<<",  b: "<<b<<endl;
        if(slope*slope < 1)
        {
            int i = x;
            double j = slope*i+b;
            //cout<<"\t\ti: "<<i<<",  j: "<<j<<endl;
            
            int incrementX = 1;
            if(!dxPositive)
            {
                incrementX = -1;
                slope *=-1;
            }
            //cout<<"  -  incX: "<<incrementX<<", newSlope: "<<slope<<endl;

            while(i != targetX || (int)j != targetY)
            {
                //                            cout<<"   --- i: "<<i<<",  j: "<<j<<endl;
                squares[i][(int)j].intersectingSegments.push_back(ls);
                //                            std::cout<<"\t\tADDED TO SQUARE: ["<<(int)i<<", "<<(int)j<<"] - <"<<(int)i*squareSize<<", "<<(int)j*squareSize<<">"<<endl;
                if(slope > 0)
                    i += incrementX;
                else 
                    j += slope;
                //                            cout<<"\ti: "<<i<<",  j: "<<j<<",    target: ("<<targetX<<", "<<targetY<<")"<<endl;
                if(i == targetX && (int)j == targetY)
                    break;
                
                squares[i][(int)j].intersectingSegments.push_back(ls);
                //                            std::cout<<"\t\tADDED TO SQUARE: ["<<(int)i<<", "<<(int)j<<"] - <"<<(int)i*squareSize<<", "<<(int)j*squareSize<<">"<<endl;
                if(slope <= 0)
                    i += incrementX;
                else 
                    j += slope;
                //                            cout<<"\ti: "<<i<<",  j: "<<j<<",    target: ("<<targetX<<", "<<targetY<<")"<<endl;
            }
        }
        else
        {
            int j = y;
            double i = (j-b)/slope;
            //cout<<"\t\ti: "<<i<<",  j: "<<j;

            slope = 1/slope;
            int incrementY= 1;
            if(!dyPositive)
            {
                incrementY = -1;
                slope *=-1;
            }
            //cout<<"   -   incY: "<<incrementY<<", newSlope: "<<slope<<endl;

            while(i != targetX || (int)j != targetY)
            {
                                                    //cout<<"   --- i: "<<i<<",  j: "<<j;
                squares[(int)i][j].intersectingSegments.push_back(ls);
                //                                    std::cout<<"\t\tADDED TO SQUARE: ["<<(int)i<<", "<<(int)j<<"] - <"<<(int)i*squareSize<<", "<<(int)j*squareSize<<">"<<endl;
                if(slope > 0)
                    j += incrementY;
                else 
                    i += slope;
                                                    //cout<<"\ti: "<<i<<",  j: "<<j<<",    target: ("<<targetX<<", "<<targetY<<")"<<endl;
                if(i == targetX && (int)j == targetY)
                    break;
                squares[(int)i][j].intersectingSegments.push_back(ls);
                //                                    std::cout<<"\t\tADDED TO SQUARE: ["<<(int)i<<", "<<(int)j<<"] - <"<<(int)i*squareSize<<", "<<(int)j*squareSize<<">"<<endl;
                
                squares[(int)i][j].intersectingSegments.push_back(ls);
                if(slope <= 0)
                    j += incrementY;
                else 
                    i += slope;
                                                    //cout<<"\ti: "<<i<<",  j: "<<j<<",    target: ("<<targetX<<", "<<targetY<<")"<<endl;
            }
        }
    }
    squares[(int)targetX][(int)targetY].intersectingSegments.push_back(ls);
    //std::cout<<"\t\tTarget ADDED TO SQUARE: ["<<(int)targetX<<", "<<(int)targetY<<"] - <"<<(int)targetX*squareSize<<", "<<(int)targetY*squareSize<<">"<<endl;
}
void DivideSquares()
{
    //Determine Square Size and number of squares
    //squareSize = sqrt(min((double)(largestX)/2.0,(double)(largestY)/2.0));
    maxSquareX = (int)(largestX/squareSize)+1;
    maxSquareY = (int)(largestY/squareSize)+1;
    std::cout<<endl<<"Squarify: size = "<<squareSize<<",  ["<<maxSquareX<<" x "<<maxSquareY<<"] - ["<<maxSquareX*squareSize<<" x "<<maxSquareY*squareSize<<"] ----"<<endl;

    //Populate Squares
    squares = new Square*[maxSquareX];
    for(int i=0; i<maxSquareX; i++)
        squares[i] = new Square[maxSquareY];

    //std::cout<<"Segments"<<endl;

    //Add Segments to squares
    for(auto i = segments.begin(); i != segments.end(); i++)
    {
        cout<<" - Divide Squares on segment: <"<<points[i->first.points[0]].x<<", "<<points[i->first.points[0]].y<<"> - <"<<points[i->first.points[1]].x<<", "<<points[i->first.points[1]].y<<">"<<endl;
        //std::cout<<"("<<i->points[0]->x<<", "<<i->points[0]->y<<") - ("<<i->points[1]->x<<", "<<i->points[1]->y<<")"<<endl;
        AddSegmentToSquares(i->first);
    }
    cout<<"Divide Squares Finishing"<<endl;
}*/

vector<Intersection> CircleLineIntersection(LineSegmentKey ls, int circleIndex, int radiiIndex)
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
void SingleCircleIntegral(int circleIndex, int radiiIndex)
{
        /*
         * Optimizations:
         *  - Save intersection values of lines already calculated = x2 speed
        */
    vector<Intersection> intersections;
    unordered_set<LineSegmentKey, LineSegmentKey::HashFunction> checkedSegments;
    vector<LineSegmentKey> queue;
    vector<LineSegmentKey> nextQueue;

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
    sort(intersections.begin(), intersections.end());
    //Intersections should be finished Calculating, print them here:
    if(debug)
    {
        debugIntersectionsStream<<circles[circleIndex].x<<" "<<circles[circleIndex].y<<" "<<radii[radiiIndex]<<" ";
    }
    //cout<<"Intersections for circle at <"<<circles[circleIndex].x<<", "<<circles[circleIndex].y<<"> with r: "<<radii[radiiIndex]<<"\n";
    double integral = 0;
    int numIntersections = intersections.size();
    for(int i=0; i<numIntersections; i++)
    {
        int j = i+1;
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

        LineSegmentData ls1 = segments[intersections[i].ls];
        LineSegmentData ls2 = segments[intersections[j].ls];
        //cout<<"before set intersection"<<endl;

        int commonTriangle = -1;
        for(int a : ls1.triangles)
        {
            for(int b : ls2.triangles)
            {
                //cout<<"a: "<<a<<", b: "<<b<<endl;
                if(a == b)
                {
                    commonTriangle = a;
                    goto foundTriangle;
                }
            }
        }
foundTriangle:
        //cout<<"Common Triangle: "<<commonTriangle<<endl;
        double arcIntegral = - radii[radiiIndex]*triangles[commonTriangle].m*(sin(theta2) - sin(theta1)) 
                                + radii[radiiIndex]*triangles[commonTriangle].n * (cos(theta2) - cos(theta1))
                                + triangles[commonTriangle].h * (theta2 - theta1);
        integral+=arcIntegral;
        //cout<<"\tt1: "<<theta1<<", t2:"<<theta2<<"\t    arcInt: "<<arcIntegral<<",  integral:"<<integral<<endl;

        debugIntersectionsStream<<intersections[i].theta<<" ";
    }
    integral *= radii[radiiIndex];
    cout<<"RADII: "<<radii[radiiIndex]<<"\t\tINTEGRAL VALUE: "<<integral<<endl;
    debugIntersectionsStream<<endl;
    intersections.clear();
}

void CalculateAllCircleIntegrals()
{
    if(debug)
        debugIntersectionsStream.open(debugIntersectionsFile);

    #pragma omp parallel for num_threads(numThreads)
        for(int i=0; i<numCircles; i++)
            for(int radiiIndex = 0; radiiIndex < radii.size(); radiiIndex++)
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

    return 0;
}

/*
    ToDo:   Triangle March Optimization


    Bugs:   Certain circles not having intersections calculated
Circle File = 
0.2 0.7 0.2
0 1.5
0 2

*/
