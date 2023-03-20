# Calculating the maximum inscribed circle of an arbitrary convex polygon
# in 2D space
# input: a list of points coordinates which are the vertices of the polygon
# output: the center and radius of the maximum inscribed circle

import math
import numpy as np
eps = 1e-4

def max_inscribed_circle(points):
    """
    Calculate the center and radius of the maximum inscribed circle of an
    arbitrary convex polygon in 2D space
    :param points: a list of points coordinates which are the vertices of the polygon
    :return: the center and radius of the maximum inscribed circle
    """
    init_points = np.array(points)
    while True:
        # select a point inside the polygon
        center = center_point_inside_polygon(points)
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
            d0 = distance_point_line(center, p1, p2)
            cp1 = p1-center
            cp2 = p2-center
            new_p1 = center + (1 - d/d0)*cp1
            new_p2 = center + (1 - d/d0)*cp2
            # save lines constructed by new_p1 and new_p2
            a, b, c = cal_abc_from_line(new_p1[0], new_p1[1], new_p2[0], new_p2[1])
            lines.append([a, b, c])
        # find the intersection of the lines
        # print("lines:",lines)

        new_points = []
        for i in range(len(lines)):
            a1,b1,c1 = lines[i]
            a2,b2,c2 = lines[(i+1)%len(lines)]
            # calculate the intersection of the lines
            # if the lines are parallel, the intersection is the center
            D = a1 * b2 - a2 * b1
            if abs(D) < eps:
                d = distance_point_polygon(center, init_points)
                return center, d
            x = (b2 * c1 - b1 * c2) / D
            y = (a1 * c2 - a2 * c1) / D
            # calculate the new polygon, remove the same points
            flag_same = False
            for p in new_points:
                if abs(p[0]-x)<eps and abs(p[1] - y)<eps:
                    flag_same = True
                    break
            if not flag_same:
                new_points.append([x, y])
        if len(new_points) <= 3:
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


if __name__ == '__main__':
    # test
    # points = [[0, 0], [1, 0], [1, 1], [0, 1]]
    # points = [[2.0, 0.0], [1.0, 1.732],[-1.0,1.732],[-2.0,0],[-1.0,-1.732],[1.0,-1.732]]
    # points = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 2.0], [-1.0, 1.0]]
    points = [[2.0, 0.0], [1.0, 1.732],[-1.0,1.732],[-2.0,0],[-1.0,-1.732],[1.0,-1.732]]
    points=np.array(points)
    center, radius = max_inscribed_circle(points)
    # print the results, 3 decimal places
    print("center:", "(", '{:.3f}'.format(center[0]),",", '{:.3f}'.format(center[1]),")")
    print("radius:", '{:.3f}'.format(radius))