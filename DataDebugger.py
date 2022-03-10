#!/usr/bin/env python3

from tkinter import *
import sys
import math

class Point:
    def __init__(self, x, y) -> None:
        self.x = x
        self.y = y
    
    def __repr__(self):
        return str(self)

    def __str__(self) -> str:
        return "("+str(self.x)+", "+str(self.y)+")"

def main():
    n = len(sys.argv)
    if(n != 5):
        print("Usage: <Points.txt> <Triangles.txt> <DebugSquares.txt> <DebugIntersections.txt>")
        return

    root = Tk()
    root.title("Data Visualization")
    root.geometry("1000x700")
    canvas = Canvas(root, width = 900, height = 600)
    canvas.pack()

    spacing = 50
    offset = 10


    pointsList = []
    with open(sys.argv[1]) as f:
        line = f.readline().split()
        while line:
            #print(line)
            lineAsInt = list(map(float, line))
            #print(lineAsInt)
            pointsList.append( Point(lineAsInt[0]*spacing + offset, lineAsInt[1]*spacing+ offset ) )
            line = f.readline().split()
    
    #print(pointsList)
    with open(sys.argv[2]) as f:
        line = f.readline().split()
        while line:
            #print(line)
            lineAsInt = list(map(int, line))
            #print(lineAsInt)
            p1 = pointsList[lineAsInt[0]]
            p2 = pointsList[lineAsInt[1]]
            p3 = pointsList[lineAsInt[2]]
            #print(p1)
            #print(p2)
            #print(p3)

            canvas.create_line(p1.x,p1.y,p2.x,p2.y, p3.x, p3.y, p1.x, p1.y)
            line = f.readline().split()
    
    #DebugCircles/Intersections
    with open(sys.argv[4]) as f:
        line = f.readline().split()
        while line:
            #print(line)
            lineAsFloat = list(map(float, line))
            #print(lineAsInt)
            centerX = lineAsFloat[0]*spacing+offset
            centerY = lineAsFloat[1]*spacing+offset
            radius = lineAsFloat[2]*spacing

            canvas.create_oval(centerX - radius, centerY - radius, centerX+radius, centerY+radius)

            for theta in lineAsFloat[3:]:
                xOffset = radius*math.cos(theta)
                yOffset = radius*math.sin(theta)
                canvas.create_oval(centerX +xOffset - 2, centerY + yOffset-2, centerX+xOffset+2, centerY+yOffset+2, fill="RED", outline="RED")
            line = f.readline().split()

    root.mainloop()


if __name__ == '__main__':
    main()