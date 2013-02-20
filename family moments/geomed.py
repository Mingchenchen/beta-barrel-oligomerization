'''Calculate and test the geometric median'''
import numpy as np

def geomed(starlist, quiet=False):
    '''Calculate the geometric median aka Fermat point using Weiszfeld's
    algorithm. You can learn more about this from the Wikipedia article
    "Geometric Median"'''
    # Throughout this function's comments, the metaphor is maintained
    # that this function moves a spaceship from the origin to the
    # geometric median of the surrounding stars by a series of smaller steps

    # Start with the spaceship at the origin
    ship = np.zeros(len(starlist[0]))

    for i in range(1000):
        # Create a matrix with the star positions as columns
        pos_mat = np.array(starlist).transpose()
    
        # Create a vector where the nth entry is the reciprocal of the
        # distance between the ship and the nth star
        inverse_distances = np.array([np.linalg.norm(star - ship)**-1 \
                                      for star in starlist])

        # Find the next spot dictated by Weiszfeld's algorithm
        destination = inverse_distances.sum()**-1 \
                      * pos_mat.dot(inverse_distances)
    
        # Break the loop if we're moving a billionth of an angstrom at a
        # time
        if np.linalg.norm(ship - destination) < 1e-12:
            break

        # Move the ship to the next spot
        ship = destination

    if not quiet:
        print(str(i) + ' iterations')

    return ship

