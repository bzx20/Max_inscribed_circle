# Calculating the maximum inscribed circle of an arbitrary convex polygon
# in 2D space
# input: a list of points coordinates which are the vertices of the polygon
# output: the center and radius of the maximum inscribed circle

import math
import numpy as np
from queue import LifoQueue
 

eps = 1e-4

def max_inscribed_circle(points):
    """
    Calculate the center and radius of the maximum inscribed circle of an
    arbitrary convex polygon in 2D space
    :param points: a list of points coordinates which are the vertices of the polygon
    :return: the center and radius of the maximum inscribed circle
    """
    init_points = np.array(points)
    # if len(init_points) == 3:
    #     p1 = init_points[0]
    #     p2 = init_points[1]
    #     p3 = init_points[2]
    #     center = cal_inscribed_center_triangle(p1, p2, p3)
    #     d = distance_point_polygon(center, init_points)
    #     return center, d
    while True:
        # select a point inside the polygon
        center = center_point_inside_polygon(points)
        print("center:", center)
        # calculate the min distance between the center and the polygon
        d = distance_point_polygon(center, points)
        # Calculate the line parallel to each edge in the polygon with distance d
        # from the center
        # see if the lines construct a triangle
        # if yes, return the center and radius
        # if not, repeat calculating the line parallel to each edge in the new polygon
        
        lines = []
        for i in range(len(points)):
            p1 = points[i]
            p2 = points[(i + 1) % len(points)]
            a, b, c = cal_abc_from_line(p1[0],p1[1],p2[0],p2[1])
            c1= c-d * np.sqrt(a*a+b*b)
            c2= c+d * np.sqrt(a*a+b*b)
            d1 = abs(a*center[0]+b*center[1]+c1)
            d2 = abs(a*center[0]+b*center[1]+c2)
            if(d1<d2):
                c=c1
            else:
                c=c2
            lines.append([a, b, c])

        new_points = cal_new_polygon(lines)
        print("new_points:",new_points)
        if len(new_points) < 3:
            d = distance_point_polygon(center, init_points)
            return center, d
        # TODO...CHECK
        elif len(new_points) == 3:
            p1 = new_points[0]
            p2 = new_points[1]
            p3 = new_points[2]
            center = cal_inscribed_center_triangle(p1, p2, p3)
            d = distance_point_polygon(center, init_points)
            return center, d
        else:
            new_points = np.array(new_points)
            
            points = new_points
    # repeat calculating the center and radius of the maximum inscribed circle



def cal_abc_from_line(x1, y1, x2, y2):
    """
    Calculate the parameters a, b, c of the line ax+by+c=0
    :param x1, y1: coordinate of the first point
    :param x2, y2: coordinate of the second point
    :return: a, b, c
    """
    a = y1 - y2
    b = x2 - x1
    c = x1 * y2 - x2 * y1
    return a, b, c
    
    
    
# calculate the distance between point p and line segment ab
def distance_point_line(p, a, b):
    """
    Calculate the distance between point p and line segment ab
    :param p: vector point p
    :param a: vector point a
    :param b: vector point b
    :return: the distance between point p and line segment ab
    """
    pa = p - a
    ba = b - a
    h = np.dot(pa, ba) / np.dot(ba, ba)
    h = np.clip(h, 0, 1)
    return np.linalg.norm(pa - h * ba)


# calculate the min distance between point p and polygon
def distance_point_polygon(p, polygon):
    """
    Calculate the min distance between point p and polygon
    :param p: vector point p
    :param polygon: a list of points coordinates which are the vertices of the
    polygon
    :return: the min distance between point p and polygon
    """
    min_distance = math.inf
    for i in range(len(polygon)):
        a = polygon[i]
        b = polygon[(i + 1) % len(polygon)]
        distance = distance_point_line(p, a, b)
        if distance < min_distance:
            min_distance = distance
    return min_distance


# select a point inside the polygon
def center_point_inside_polygon(polygon):
    """
    Select the center point inside the polygon
    :param polygon: a list of points coordinates which are the vertices of the
    polygon
    :return: the center point inside the polygon
    """
    # calculate the center of the polygon
    center = np.mean(polygon, axis=0)
    # select a point inside the polygon
    for i in range(len(polygon)):
        a = np.array(polygon[i])
        b = np.array(polygon[(i + 1) % len(polygon)])
        c = np.array(polygon[(i + 2) % len(polygon)])
        # calculate the angle between line segment ab and bc
        ba = a - b
        bc = c - b
        angle = np.arccos(np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc)))
        # if the angle is larger than 180 degree, the point b is inside the
        # polygon
        if angle > math.pi:
            return b
    return center


# orientation between points p,q,r
def orientation(p,q,r):
    val = (q[1]-p[1])*(r[0]-q[0])-(q[0]-p[0])*(r[1]-q[1])
    # colinear
    if (abs(val)<=eps): return 0
    # counterclockwise
    elif(val>eps): return 1
    # clockwise
    else: return 2

def cal_new_polygon(lines):
    new_points=LifoQueue(maxsize=len(lines))
    line_stack=LifoQueue(maxsize=len(lines))
    # line_stack.put(lines[0])
    line_stack.put(lines[1])
    line_stack.put(lines[2])
    # line_stack.append(lines[2])
    first_point=cal_intersection_point(lines[0],lines[1])
    new_points.put(first_point)
    new_points.put(cal_intersection_point(lines[1],lines[2]))
    k = 3
    n = len(lines)
    while(k <= n):
        top = line_stack.get()
        [nx,ny] = cal_intersection_point(top,lines[k%n])
        [cx,cy] = new_points.get()
        [px,py] = new_points.get()
        new_points.put([px,py])
        if(orientation([px,py],[cx,cy],[nx,ny])==2):
            new_points.put([cx,cy])
            new_points.put([nx,ny])
            line_stack.put(top)
        elif(orientation([px,py],[cx,cy],[nx,ny])==1):
            # select the former line to calculate the intersection point
            top = line_stack.get()
            [nx,ny]=cal_intersection_point(top,lines[k%n])
            new_points.put([nx,ny])
            line_stack.put(top)
        else:
            line_stack.put(top)
            # new_points.put([cx,cy])
            new_points.put([nx,ny])
        line_stack.put(lines[k%n])
        k+=1
    points = []
    while not new_points.empty():
        points.append(new_points.get())
    points.reverse()
    return points


def cal_intersection_point(l1,l2):
    a1,b1,c1 = l1
    a2,b2,c2 = l2
    D = a1 * b2 - a2 * b1
    # if abs(D) < eps:
    # else:
    x = (b1 * c2 - b2 * c1) / D
    y = (a2 * c1 - a1 * c2) / D
    return [x,y]

# calculate the inscribed center of triangle
def cal_inscribed_center_triangle(p1, p2, p3):
    """
    Calculate the inscribed center of triangle
    :param p1, p2, p3: the vertices of the triangle
    :return: the inscribed center of triangle
    """
    v1 = np.array([p2[0] - p3[0], p2[1] - p3[1]])
    v2 = np.array([p3[0] - p1[0], p3[1] - p1[1]])
    v3 = np.array([p1[0] - p2[0], p1[1] - p2[1]])
    d1 = np.linalg.norm(v1)
    d2 = np.linalg.norm(v2)
    d3 = np.linalg.norm(v3)
    d = d1 + d2 + d3
    d1 = d1/d
    d2 = d2/d
    d3 = d3/d
    # calculate the inscribed center of triangle
    x = d1 * p1[0] + d2 * p2[0] + d3 * p3[0]
    y = d1 * p1[1] + d2 * p2[1] + d3 * p3[1]
    return [x, y]




if __name__ == '__main__':
    # test
    # points = [[0, 0], [2, 0], [2, 1], [0, 3]]
    # points = [[2.0, 0.0], [1.0, 1.732],[-1.0,1.732],[-2.0,0],[-1.0,-1.732],[1.0,-1.732]]
    # points = [[245,-241],[212,-495],[492,-461],[560.36,-184.65]]
    # points = [[417.39,-63.4],[384,-200],[483,-338],[662,-292],[738,-144],[626,-20]]
    # points = [[72,-428],[92,-539],[261.5,-539],[178,-464]]
    # points = [[135,-61],[93,-249],[135,-387],[304.6,-495.5],[364.65,-400.22],[364.65,-228],[345,-140],[233,-34]]
    # points = [[52.17,-420.5],[70.11,-564.5],[222.17,-564.5],[83.67,-380]]
    # points = [[0.0, 0.0], [1.0, 0.0], [2.0, 1.0], [0.0, 3.0], [-1.0, 1.0]]
    # points = [[2.0, 0.0], [1.0, 1.732],[-1.0, 1.732],[-2.0, 0],[-1.0, -1.732]]
    # points = [[0,0],[1,1.732],[2,0]]
    points=[[472,-370],[455.81,-501.5],[668.19,-501.5]]
    # points=[[0,0],[1,1.732],[2,0]]
    points=np.array(points)
    center, radius = max_inscribed_circle(points)
    # print the results, 3 decimal places
    print("center:", "(", '{:.3f}'.format(center[0]),",", '{:.3f}'.format(center[1]),")")
    print("radius:", '{:.3f}'.format(radius))
    print(center[0]-radius,-center[1]-radius)


    
    # [x,y]=cal_inscribed_center_triangle([472,-370],[455.81,-501.5],[668.19,-501.5])
    # print(x,y)
    # s = LifoQueue(maxsize=5)
    # s.put(1)
    # s.put(2)
    # s.put(3)
    # print(s.get())
    # print("size:",s.qsize())
    # print(s.get())
    # print(s.get())
    # a,b,c = cal_abc_from_line(2,0,1,1.732)
    # print(a,b,c)
    # [x,y]=cal_intersection_point([1,2,3],[2,4,6])
    # print([x,y])

    # center = center_point_inside_polygon(points)
    # d = distance_point_polygon(center, points)
    # lines = []
    # for i in range(len(points)):
    #     p1 = points[i]
    #     p2 = points[(i + 1) % len(points)]
    #     a, b, c = cal_abc_from_line(p1[0],p1[1],p2[0],p2[1])
    #     c1= c-d * np.sqrt(a*a+b*b)
    #     c2= c+d * np.sqrt(a*a+b*b)
    #     d1 = abs(a*center[0]+b*center[1]+c1)
    #     d2 = abs(a*center[0]+b*center[1]+c2)
    #     if(d1<d2):
    #         c=c1
    #     else:
    #         c=c2
    #     lines.append([a, b, c])
    # print(lines)
    # new_points = cal_new_polygon(lines)
    # print(new_points)
