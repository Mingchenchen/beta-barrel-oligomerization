"""Monte carlo methods and supporting functions.
point_within_circle(center,radius): gives a random point within a given
    circle.
circle_search(starting_point, starting_range, function, iterations,
              return_type='vector'): minimizes a function of a vector.
find_centers(groupdict, time, selectionname = 'relevant'): minimizes basis
    moment for all groups in a dictionary."""

import random
import math
import moments
from moments import Vector
import time
import csv
import functools

def point_within_circle(center,radius):
    """Returns the position vector of a point within a given radius of a given
    center, in the plane z = k where k is the z component of the center.
    ARGUMENTS:
    an iterable representing the x, y and z components of the center
    a number representing the radius
    RETURNS:
    a vector"""
    center = Vector(center)
    angle = random.random() * 2*math.pi
    reach = random.random() * radius
    return center + reach * Vector((math.cos(angle),math.sin(angle),0))
    
def circle_search(starting_point, starting_range, function, iterations,
                  return_type='vector', target = None, max_time = None):
    """Minimizing a function of a vector by trying random points in the
    vicinity of a given starting vector, then in a smaller vicinity around
    whatever closer point it finds in that vicintiy, ad nauseum.
    ARGUMENTS:
    starting_point: the vector to start from.
    starting_range: the radius around the starting vector to search for
        the second vector.
    function: the function to minimize.
    iterations: number of iterations to run.
    return_type="vector": return type. Alternative is "list", which will give
        a list of every intermediate point on the way to the final minimum.
        As always with my option arguments, you can just use "l" or "v"
    target=None: if set, the search will break if it finds a value below
        target.
    max_time=None: if set, the seach will break if it exceeds the maxtime
        given in seconds.
    RETURNS:
    A list or vector, as specified by "return_type". Vector would be final
    vector that gave the minimum value of the given function, list would be
    every vector that was found along the way that the function mapped to a
    lower value than had ever been found before at that time."""
    starting_point = Vector(starting_point)
    current_min = function(starting_point)
    
    # NOTETOSELF: before you use this, make time the default option, and make
    # it givable in minutes
    
    # Multiplying by standard_radius is like dividing by the first value of
    # the function and then multiplying by the original range. If the next
    # value of the function is half the first, then the next range will be
    # half the starting range.
    standard_radius = starting_range / current_min
    current_point = starting_point
    if return_type[0] == 'l': # for list
        output_list = []
    if max_time is not None:
        start_time = time.clock()
    for iteration in xrange(iterations):
        if return_type == 'l':
            output_list.append(current_point)
        next_point = point_within_circle(current_point,
                        function(current_point * standard_radius))
        next_min = function(next_point)
        if next_min < current_min:
            current_min = next_min
            current_point = next_point
        if target is not None and current_min < target:
            break
        if max_time is not None and time.clock() - start_time > max_time:
            break
    if return_type[0] == 'l':
        return output_list
    elif return_type[0] == 'v':
        return current_point
        

def basis_moment_norm(sele, vector):
    return moments.basis_moment(sele,vector).norm()
    
def find_centers(groupdict, time, path, selectionname = 'relevant'):
    """Minimize basis moment.
    ARGUMENTS:
    Dictionary of Group objects
    Number of seconds to spend on minimization IN TOTAL (time per group is
        this arg divided by number of groups)
    Path to place csv files with centers
    selectionname = 'relevant'. So by default, it'll minimize the basis moment
        for the selection groupname.relevant for each group"""
    
    count = 0
    allrows = []
    for group in groupdict.itervalues():
        sele = group.name + '.' + selectionname
        #print 'Function is functools.partial(basis_moment_norm, ' + sele + ')'
        #print 'Value at (0,0) is ' + str(functools.partial(basis_moment_norm, sele)((0,0)))
        group.center = circle_search((0,0), 5,
            functools.partial(basis_moment_norm, sele), 10000,
            max_time = float(time) / len(groupdict))
        if count < 10:
            countstring = '0' + str(count)
        else:
            countstring = str(count)
        allrows.append([group.name, group.center])
        csvfile = open(path + '/' + countstring + '.csv', 'w')
        writer = csv.writer(csvfile)
        for row in allrows:
            writer.writerow(row)
        csvfile.close()
        count += 1