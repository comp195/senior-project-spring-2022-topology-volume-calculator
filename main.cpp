// C++ program to demonstrate the
// boilerplate code
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <math.h>
#include <ctime>
using namespace std;

struct LineSegment;
clock_t time_req;
double pi = 3.141592653589793238462643383279502884197169399375105820974944;
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
    const LineSegment* segments[3];
    //z = -mx - ny + h
    double m;
    double n;
    double h;
};

struct LineSegment
{
    Point* points[2];
    size_t hashNum;
    set<int> triangles;

    LineSegment(Point *a, Point *b)
    {
        points[0] = a;
        points[1] = b;

        size_t x1Hash = std::hash<int>()(a->x);
        size_t x2Hash = std::hash<int>()(b->x);
        size_t y1Hash = std::hash<int>()(a->y)/2;
        size_t y2Hash = std::hash<int>()(b->y)/2;
        hashNum = x1Hash ^ y1Hash + x2Hash ^ y2Hash;
    }

    bool operator==(const LineSegment& otherSegment) const
    {
        if ((points[0] == otherSegment.points[0] && points[1] == otherSegment.points[1]) || 
                points[0] == otherSegment.points[1] && points[1] == otherSegment.points[0])
            return true;
        else return false;
    }

    struct HashFunction
    {
        size_t operator()(const LineSegment& segment) const
        {
            return segment.hashNum;
        }
    };
};

struct Intersection
{
    double theta;
    const LineSegment* ls;
    Intersection(double t, const LineSegment *segment)
    {
        theta = t;
        ls = segment;
    }

    bool operator==(const Intersection& other) const
    {
        return other.theta == this->theta;
    }
};
struct compareIntersection {
    bool operator() (Intersection a, Intersection b) const {
        return a.theta<b.theta;
    }
};

struct Square
{
    unordered_set<LineSegment, LineSegment::HashFunction> intersectingSegments;
};

//x, y, h = radius
struct Circle
{
    double x;
    double y;
    double lineIntegral;
};

void PrintPoint(Point a)
{
    cout<<"("<<a.x<<", "<<a.y<<")";
}
void PrintSegment(LineSegment s)
{
    PrintPoint(*s.points[0]);
    cout<<" - ";
    PrintPoint(*s.points[1]);
}

Point *points;
int numPoints;
unordered_set<LineSegment, LineSegment::HashFunction> segments;
int numSegments;
Triangle *triangles;
int numTriangles;
Circle *circles;
int numCircles;
vector<double> radii;

double largestX;
double largestY;
double squareSize;
Square **squares;
int maxSquareX;
int maxSquareY;
string debugSquaresFile = "DebugSquares.txt";
bool debug = true;


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
                std::cout<<"            - seg: ("<<k->points[0]->x<<", "<<k->points[0]->y<<") - ("<<k->points[1]->x<<", "<<k->points[1]->y<<")"<<endl;
                if(debug)
                {
                    debugFile<<k->points[0]->x<<" "<<k->points[0]->y<<" "<<k->points[1]->x<<" "<<k->points[1]->y<<endl;
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
        cout<<"\tTriangle "<<&triangles[i]<<"   "<<triangles[i].segments[0]<<endl<<"\t\t";
        PrintSegment(*triangles[i].segments[0]);
        cout<<"  :  ";
        PrintSegment(*triangles[i].segments[1]);
        cout<<"  :  ";
        PrintSegment(*triangles[i].segments[2]);
        cout<<endl;
    }
}

int FileNumLines(string fileName)
{
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
        squareSize = 0;

        time_req = clock();
        for(int i = 0; i < numTriangles; i++)
        {
            int ind1, ind2, ind3;
            readFile >> ind1 >> ind2 >> ind3;
            Point* a = &points[ind1];
            Point* b = &points[ind2];
            Point* c = &points[ind3];
            LineSegment s1 = LineSegment(a,b);
            LineSegment s2 = LineSegment(a,c);
            LineSegment s3 = LineSegment(b,c);

            cout<<"\tTriangle "<<i<<endl<<"\t\t";
            PrintSegment(s1);
            cout<<"  :  ";
            PrintSegment(s2);
            cout<<"  :  ";
            PrintSegment(s3);
            cout<<endl;

            //Calculate equation for plane
            double AB[3];
            double AC[3];

            AB[0] = b->x-a->x;
            AB[1] = b->y-a->y;
            AB[2] = b->h-a->h;

            AC[0] = c->x-a->x;
            AC[1] = c->y-a->y;
            AC[2] = c->h-a->h;

            double normalVector[3];
            normalVector[0] = AB[1]*AC[2] - AB[1]*AC[2];
            normalVector[1] = AB[0]*AC[2] - AB[0]*AC[2];
            normalVector[2] = AB[0]*AC[1] - AB[0]*AC[1];

            double m = normalVector[0]/normalVector[2];
            double n = normalVector[1]/normalVector[2];
            double h = a->x*m + a->y*n + a->h;
            triangles[i].m = m;
            triangles[i].n = n;
            triangles[i].h = h;


            //Determine line segment x and y length, for use in finding squareSize - square size will be largest segment length size
            double xLen1 = abs(points[ind1].x - points[ind2].x);            
            double xLen2 = abs(points[ind1].x - points[ind3].x); 
            double xLen3 = abs(points[ind2].x - points[ind3].x); 
            double xLen = max(max(xLen1, xLen2), xLen3);

            double yLen1 = abs(points[ind1].y - points[ind2].y);
            double yLen2 = abs(points[ind1].y - points[ind3].y);
            double yLen3 = abs(points[ind2].y - points[ind3].y);
            double yLen = max(max(yLen1, yLen2), yLen3);

            squareSize = max(max(xLen, yLen), squareSize);


            unordered_set<LineSegment, LineSegment::HashFunction>::iterator got = segments.find(s1);
            //Line Segment not created yet
            if(got == segments.end())
            {
                segments.insert(s1);
                numSegments++;
            }
            got = segments.find(s1);
            s1 = *got;

            got = segments.find(s2);
            if(got == segments.end())
            {
                segments.insert(s2);
                numSegments++;
            }
            got = segments.find(s2);
            s2 = *got;

            got = segments.find(s3);
            if(got == segments.end())
            {
                segments.insert(s3);
                numSegments++;
            }
            got = segments.find(s3);
            s3 = *got;

            triangles[i].segments[0] = &*segments.find(s1);
            triangles[i].segments[1] = &*segments.find(s2);
            triangles[i].segments[2] = &*segments.find(s3);
            s1.triangles.insert(i);
            s2.triangles.insert(i); 
            s3.triangles.insert(i);
            cout<<"\t\t";
            PrintSegment(s1);
            cout<<"  :  ";
            PrintSegment(s2);
            cout<<"  :  ";
            PrintSegment(s3);
            cout<<endl;
        }
	    time_req = clock() - time_req;
        cout<<"Reading Triangles took "<<(float)time_req/CLOCKS_PER_SEC<<" seconds"<<endl;
    }
}
void ReadCircles(char *circleFile)
{
    ifstream readFile(circleFile);
    if(readFile)
    {
        numCircles = FileNumLines(circleFile);
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

void AddSegmentToSquares(const LineSegment &ls)
{
    //Line segment in square coordinates
    double x = ls.points[0]->x / squareSize;
    double y = ls.points[0]->y / squareSize;
    double targetX = ls.points[1]->x / squareSize;
    double targetY = ls.points[1]->y / squareSize;
    //std::cout<<"Adding Segment: ("<<x<<","<<y<<") - ("<<targetX<<","<<targetY<<")"<<endl;

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
            squares[i][(int)y].intersectingSegments.insert(ls);
            //std::cout<<"\t\tADDED TO SQUARE: ["<<i<<", "<<(int)y<<"] - <"<<(int)i*squareSize<<", "<<(int)y*squareSize<<">"<<endl;
        }
        //squares[(int)targetX][(int)targetY].intersectingSegments.insert(ls);
    }
    else if(dX == 0)
    {
        for(int j = y; dyPositive? j<targetY : j>targetY; dyPositive? j++ : j--)
        {
            squares[(int)x][j].intersectingSegments.insert(ls);
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
                squares[i][(int)j].intersectingSegments.insert(ls);
                //                            std::cout<<"\t\tADDED TO SQUARE: ["<<(int)i<<", "<<(int)j<<"] - <"<<(int)i*squareSize<<", "<<(int)j*squareSize<<">"<<endl;
                if(slope > 0)
                    i += incrementX;
                else 
                    j += slope;
                //                            cout<<"\ti: "<<i<<",  j: "<<j<<",    target: ("<<targetX<<", "<<targetY<<")"<<endl;
                if(i == targetX && (int)j == targetY)
                    break;
                
                squares[i][(int)j].intersectingSegments.insert(ls);
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
                squares[(int)i][j].intersectingSegments.insert(ls);
                //                                    std::cout<<"\t\tADDED TO SQUARE: ["<<(int)i<<", "<<(int)j<<"] - <"<<(int)i*squareSize<<", "<<(int)j*squareSize<<">"<<endl;
                if(slope > 0)
                    j += incrementY;
                else 
                    i += slope;
                                                    //cout<<"\ti: "<<i<<",  j: "<<j<<",    target: ("<<targetX<<", "<<targetY<<")"<<endl;
                if(i == targetX && (int)j == targetY)
                    break;
                squares[(int)i][j].intersectingSegments.insert(ls);
                //                                    std::cout<<"\t\tADDED TO SQUARE: ["<<(int)i<<", "<<(int)j<<"] - <"<<(int)i*squareSize<<", "<<(int)j*squareSize<<">"<<endl;
                
                squares[(int)i][j].intersectingSegments.insert(ls);
                if(slope <= 0)
                    j += incrementY;
                else 
                    i += slope;
                                                    //cout<<"\ti: "<<i<<",  j: "<<j<<",    target: ("<<targetX<<", "<<targetY<<")"<<endl;
            }
        }
    }
    squares[(int)targetX][(int)targetY].intersectingSegments.insert(ls);
    //std::cout<<"\t\tTarget ADDED TO SQUARE: ["<<(int)targetX<<", "<<(int)targetY<<"] - <"<<(int)targetX*squareSize<<", "<<(int)targetY*squareSize<<">"<<endl;
}
void DivideSquares()
{
    //Determine Square Size and number of squares
    //squareSize = sqrt(min((double)(largestX)/2.0,(double)(largestY)/2.0));
    maxSquareX = (int)(largestX/squareSize)+1;
    maxSquareY = (int)(largestY/squareSize)+1;
    std::cout<<endl<<"Squarify: size = "<<squareSize<<",  ["<<maxSquareX<<" x "<<maxSquareY<<"] - ["<<maxSquareX*squareSize<<" x "<<maxSquareY*squareSize<<"] ----"<<endl;

    squares = new Square*[maxSquareX];
    for(int i=0; i<maxSquareX; i++)
        squares[i] = new Square[maxSquareY];

    //std::cout<<"Segments"<<endl;

    for(auto i = segments.begin(); i != segments.end(); i++)
    {
        cout<<" - Divide Squares on segment: <"<<(*i).points[0]->x<<", "<<(*i).points[0]->y<<"> - <"<<(*i).points[1]->x<<", "<<(*i).points[1]->y<<">"<<endl;
        //std::cout<<"("<<i->points[0]->x<<", "<<i->points[0]->y<<") - ("<<i->points[1]->x<<", "<<i->points[1]->y<<")"<<endl;
        AddSegmentToSquares(*i);
    }
    cout<<"Divide Squares Finishing"<<endl;
}

set<Intersection,compareIntersection> CircleLineIntersection(const LineSegment &ls, int circleIndex, int radiiIndex)
{
    set<Intersection,compareIntersection> intersections;
    double AB[2];
    double AC[2];
    double IC[2];

    AB[0] = ls.points[1]->x-ls.points[0]->x;
    AB[1] = ls.points[1]->y-ls.points[0]->y;

    AC[0] = ls.points[0]->x-circles[circleIndex].x;
    AC[1] = ls.points[0]->y-circles[circleIndex].y;

    double a = AB[0]*AB[0] + AB[1]*AB[1];
    double bHalf = (AC[0]*AB[0] + AC[1]*AB[1]);
    double c = AC[0]*AC[0] + AC[1]*AC[1] - radii[radiiIndex]*radii[radiiIndex];

    double discriminant = bHalf*bHalf - a*c;
    if(discriminant < 0)
    {
        return intersections;
    }
    else if(discriminant == 0)
    {
        double n = -bHalf/a;
        IC[0] = AC[0] + n*AB[0];
        IC[1] = AC[1] + n*AB[1];

        double theta = atan(IC[0]/IC[1]);
        if(IC[0]<0)
            theta = pi-theta;

        intersections.insert(Intersection(theta, &ls));
        /* NEED NEW ALGORITHM: OPERATE ON A FULL TRIANGLE, FIND ALL INTERSECTIONS, CALCULATE LINE INTEGRALS, QUEUE TWO MORE TRIANGLES FROM THE UNUSED SIDES
         *    Find Intersecting Segment using squares -> queue triangles -> calculate integrals -> queue triangle on other side of segment if has intersection
                      -> stop when queue'd triangles don't queue any more
         *  
         * Optimizations:
         *  - Save intersection values of lines already calculated = x2 speed
        */
        return intersections;
    }
    else
    {
        double n = (-bHalf + sqrt(discriminant))/a;
        IC[0] = AC[0] + n*AB[0];
        IC[1] = AC[1] + n*AB[1];
        double theta = atan(IC[0]/IC[1]);
        if(IC[0]<0)
            theta = pi-theta;

        intersections.insert(Intersection(theta, &ls));

        n = (-bHalf - sqrt(discriminant))/a;
        IC[0] = AC[0] + n*AB[0];
        IC[1] = AC[1] + n*AB[1];
        theta = atan(IC[0]/IC[1]);
        if(IC[0]<0)
            theta = pi-theta;

        intersections.insert(Intersection(theta, &ls));
        
        return intersections;
    }
}

void SingleCircleIntegral(int circleIndex, int radiiIndex)
{
    set<Intersection,compareIntersection> intersections;
    set<int> triangleQueue;
    set<int> nextTriangleQueue;
    set<int> coveredTriangles;

    //Find square somehow ----------------------------------------------------NEED TO DO THIS ASAP
    Square s = squares[0][0];
    auto start = s.intersectingSegments.begin();
    auto end = s.intersectingSegments.end();

    //Find First line segment that intersects circle,. then transition to marching algorithm
    for(auto i = start; i != end; i++)
    {
        LineSegment ls = *i;
        intersections = CircleLineIntersection(ls, circleIndex, radiiIndex);
        if(!intersections.empty())
        {
            triangleQueue = ls.triangles;
            break;
        }
    }

    //While there are triangles left to march around the circle
    while(!triangleQueue.empty())
    {
        for(int currTriangle:triangleQueue)
        {
            coveredTriangles.insert(currTriangle);
            set<Intersection,compareIntersection> temp;
            for(int segmentIndex = 0; segmentIndex < 3; segmentIndex++)
            {
                temp = CircleLineIntersection(*(triangles[currTriangle].segments[segmentIndex]), currTriangle, radiiIndex);
                if(temp.empty())
                {
                    //Save the calculated intersection points
                    intersections.insert(temp.begin(), temp.end());
                    //Queue triangles on both sides of this line segment for next round (duplicate triangle removed through subtract covered triangles set)
                    //This makes line segment intersections get calculated twice!! try memoization optimization later
                    nextTriangleQueue.insert(triangles[currTriangle].segments[segmentIndex]->triangles.begin(), triangles[currTriangle].segments[segmentIndex]->triangles.end());
                }
            }
        }

        //Make sure we don't re-do triangles that we have already covered
        triangleQueue.clear();
        set_difference(nextTriangleQueue.begin(), nextTriangleQueue.end(), coveredTriangles.begin(), coveredTriangles.end(), inserter(triangleQueue,triangleQueue.end()));
    }

    //Intersections should be finished Calculating, print them here:
    cout<<"Intersections for circle at <"<<circles[circleIndex].x<<", "<<circles[circleIndex].y<<"> with r: "<<radii[radiiIndex]<<"\n";
}

void CalculateAllCircleIntegrals()
{
    for(int i=0; i<numCircles; i++)
        for(int radiiIndex: radii)
            SingleCircleIntegral(i, radiiIndex);
}

// Driver Code
int main(int argc, char *argv[])
{
    if(argc != 4)
    {
        std::cout<<"Usage: <Points.txt> <Trianges.txt> <Circles.txt>"<<endl;
        return 0;
    }

    cout<<"Reading In Points\n";
    ReadPoints(argv[1]);

    cout<<"Reading in Triangles\n";
    ReadTriangles(argv[2]);
    cout<<"NumSegments = "<<numSegments<<endl;

    for(int i=0; i<numPoints; i++)
    {
        cout<<"\nPrinting Point: "<<&points[i]<<",    ";
        PrintPoint(points[i]);
    }
    cout<<"\n"<<endl;
    for(auto i = segments.begin(); i!=segments.end(); i++)
    {
        cout<<"\nPrinting Segment: "<<&*i<<",    ";
        PrintSegment(*i);
    }
    PrintTriangles();
    ReadCircles(argv[3]);
    cout<<"Starting Squarify"<<endl;

    DivideSquares();
    //PrintSquares();

    cout<<"Starting Circle Calculations"<<endl;
    CalculateAllCircleIntegrals();

    return 0;
}


