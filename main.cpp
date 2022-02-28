// C++ program to demonstrate the
// boilerplate code
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <unordered_set>
using namespace std;

struct LineSegment;

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
    LineSegment* segments[3];
};

struct LineSegment
{
    Point* points[2];
    size_t hashNum;
    vector<Triangle*> triangles;

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

struct Square
{
    unordered_set<LineSegment, LineSegment::HashFunction> intersectingSegments;
};

struct Circle
{
    double radius;
    Point center;
};

Point *points;
int numPoints;
unordered_set<LineSegment, LineSegment::HashFunction> segments;
int numSegments;
Triangle *triangles;

double largestX;
double largestY;
double squareSize;
Square **squares;
int maxSquareX;
int maxSquareY;


void PrintSquares()
{
    cout<<endl<<"Printing Squares: ["<<maxSquareX<<", "<<maxSquareY<<"]"<<endl;
    for(int i=0; i<maxSquareX; i++)
    {
        cout<<"   x: "<<i<<endl;
        for(int j=0; j<maxSquareY; j++)
        {
            cout<<"      - y: "<<j<<endl;
            for(auto k = squares[i][j].intersectingSegments.begin(); k!= squares[i][j].intersectingSegments.end(); k++)
            {
                cout<<"            - seg: ("<<k->points[0]->x<<", "<<k->points[0]->y<<") - ("<<k->points[1]->x<<", "<<k->points[1]->y<<")"<<endl;
            }
        }
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
        cout<<"File Not Open"<<endl;
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
        cout<<"Lines in points file: "<<numPoints<<endl;
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

            cout<<"("<<points[i].x<<", "<<points[i].y<<", "<<points[i].h<<")"<<endl;
        }
        cout<<endl;
    }
    else
    {
        cout<<"Could Not Find Points File: "<<pointsFile<<endl;
        exit(EXIT_FAILURE);
    }
}
void ReadTriangles(char *triangleFile)
{
    ifstream readFile(triangleFile);
    if(readFile)
    {
        int fileLength = FileNumLines(triangleFile);
        triangles = new Triangle[fileLength];
        double input;
        numSegments = 0;
        squareSize = 0;

        for(int i = 0; i < fileLength; i++)
        {
            int ind1, ind2, ind3;
            readFile >> ind1 >> ind2 >> ind3;
            LineSegment s1 = LineSegment(&points[ind1], &points[ind2]);
            LineSegment s2 = LineSegment(&points[ind1], &points[ind3]);
            LineSegment s3 = LineSegment(&points[ind2], &points[ind3]);

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

            unordered_set<LineSegment, LineSegment::HashFunction>::const_iterator got = segments.find(s1);
            //Line Segment not created yet
            if(got == segments.end())
            {
                segments.insert(s1);
                numSegments++;
            }
            else
            {
                s1 = *got;
            }

            got = segments.find(s2);
            if(got == segments.end())
            {
                segments.insert(s2);
                numSegments++;
            }
            else
            {
                s2 = *got;
            }

            got = segments.find(s3);
            if(got == segments.end())
            {
                segments.insert(s3);
                numSegments++;
            }
            else
            {
                s3 = *got;
            }

            triangles[i].segments[0] = &s1;
            triangles[i].segments[1] = &s2;
            triangles[i].segments[2] = &s3;
            s1.triangles.push_back(&triangles[i]); 
            s2.triangles.push_back(&triangles[i]); 
            s3.triangles.push_back(&triangles[i]);
        }
    }
}
void ReadCircles(char *circleFile)
{
    ifstream readFile(circleFile);
    if(readFile)
    {

    }
    else
    {
        cout<<"Could Not Find Circle File: "<<circleFile<<endl;
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

    //Set highest point as x/y, lower point as targetX/targetY
    if(y < targetY)
    {
        double temp = y;
        y = targetY;
        targetY = temp;

        temp = x;
        x = targetX;
        targetX = temp;
    }
    cout<<"     : <"<<x<<", "<<y<<"> - <"<<targetX<<", "<<targetY<<">"<<endl;

    //Delta x,y
    double dX = targetX - x;
    double dY = targetY - y;
    //Is x increasing or decreasing
    int xinc = dX > 0? 1 : -1;
    //y is always decreasing

    //If ls is a horizontal line, can't calculate slope so special case
    if(dY == 0)
    {
        while((int)x != (int)targetX)
        {
            //Add square[(int)x][(int)y];
            squares[(int)x][(int)y].intersectingSegments.insert(ls);
            cout<<"\t\tADDED TO SQUARE: ["<<(int)x<<", "<<(int)y<<"]"<<endl;
            x += xinc;
        }
    }
    else
    {
        //Change in Y per x, negative/positive adjusted for x increase/decrease
        double yPerX = dY/dX * xinc;
        //Number of X squares per one Y square;
        double xPerY = dX/dY;
        //Percent square that y is at (y is descending so counting from bottom) (it is already out of 1)
        double initialYPercent = (y - (int)y);
        cout<<"yPx: "<<yPerX<<",  xPy: "<<xPerY<<",  initYP: "<<initialYPercent<<endl;

        double nextXInt = x - xPerY*initialYPercent;
        while((int)x != (int)targetX || (int)y != (int)targetY)
        {
            cout<<"pos: ("<<x<<", "<<y<<")  ->  ("<<targetX<<", "<<targetY<<")"<<endl;
            cout<<"   -- "<<nextXInt<<endl;
            squares[(int)x][(int)y].intersectingSegments.insert(ls);
            if((int)x == (int)nextXInt)
            {
                nextXInt += xPerY;
                y--;
            }
            else
                x++;
        }
    }
    squares[(int)targetX][(int)targetY].intersectingSegments.insert(ls);
    cout<<"\t\tADDED TO SQUARE: ["<<(int)targetX<<", "<<(int)targetY<<"]"<<endl;
}

void DivideSquares()
{
    //Determine Square Size and number of squares
    //squareSize = sqrt(min((double)(largestX)/2.0,(double)(largestY)/2.0));
    maxSquareX = (int)(largestX/squareSize)+1;
    maxSquareY = (int)(largestY/squareSize)+1;
    cout<<endl<<"Squarify: size = "<<squareSize<<",  ["<<maxSquareX<<" x "<<maxSquareY<<"] - ["<<maxSquareX*squareSize<<" x "<<maxSquareY*squareSize<<"] ----"<<endl;

    squares = new Square*[maxSquareX];
    for(int i=0; i<maxSquareX; i++)
        squares[i] = new Square[maxSquareY];

    cout<<"Segments"<<endl;

    for(auto i = segments.begin(); i != segments.end(); i++)
    {
        cout<<"("<<i->points[0]->x<<", "<<i->points[0]->y<<") - ("<<i->points[1]->x<<", "<<i->points[1]->y<<")"<<endl;
        AddSegmentToSquares(*i);
    }
}

void CircleLineIntersection(LineSegment &ls, Circle &c)
{
}

// Driver Code
int main(int argc, char *argv[])
{
    if(argc != 4)
    {
        cout<<"Usage: <Points.txt> <Trianges.txt> <Circles.txt>"<<endl;
        return 0;
    }

    ReadPoints(argv[1]);
    ReadTriangles(argv[2]);
    //ReadCircles(argv[3]);

    DivideSquares();
    PrintSquares();

     
    /*
    

    maxSquareX = (int)(largestX/squareSize)+1;
    maxSquareY = (int)(largestY/squareSize)+1;
    cout<<endl<<"Squarify: size = "<<squareSize<<",  ["<<maxSquareX<<" x "<<maxSquareY<<"] - ["<<maxSquareX*squareSize<<" x "<<maxSquareY*squareSize<<"] ----"<<endl;

    squares = new Square*[maxSquareX];
    for(int i=0; i<maxSquareX; i++)
        squares[i] = new Square[maxSquareY];

    LineSegment s1 = LineSegment(&points[0], &points[1]);
    AddSegmentToSquares(s1);*/


    return 0;
}


/*
points
0 0 1
0 1 2
1 0 1.5
1 1 2
1.5 1.5 3
2 1 2
1 2 2
*/

