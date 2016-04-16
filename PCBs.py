"""
This script provides the IUPAC names of congeners of PCBs.
"""

import itertools

"""
A set of four functions that can yield all possible rotations of the diphenyl.
The molecule is assumed to be oriented as: Ph-Ph (i.e. horizontally)
For all four functions, the input is a tuple and the output is a tuple.
"""
def x_flip(c):
    x_axis_rotation = (c[4], c[3], c[2], c[1], c[0], c[9], c[8], c[7], c[6], c[5])
    return x_axis_rotation

def y_flip(c):
    y_axis_rotation = (c[5], c[6], c[7], c[8], c[9], c[0], c[1], c[2], c[3], c[4])
    return y_axis_rotation

def xy_flip(c):
    xy_axes_rotation = (c[9], c[8], c[7], c[6], c[5], c[4], c[3], c[2], c[1], c[0])
    return xy_axes_rotation

def half_flip(c):
    half_flip_rotation = (c[4], c[3], c[2], c[1], c[0], c[5], c[6], c[7], c[8], c[9])
    return half_flip_rotation


def iupac_name(c):
    """
    Produces the IUPAC name of a given congener of PCB.
    :param c: A tuple, the congener
    :return: A string, the IUPAC name
    """
    iupac_numbers = {
        0: '2',
        1: '3',
        2: '4',
        3: '5',
        4: '6',
        5: '2\'',
        6: '3\'',
        7: '4\'',
        8: '5\'',
        9: '6\''
    }
    iupac_prefix = {
        1: '',
        2: 'di',
        3: 'tri',
        4: 'tetra',
        5: 'penta',
        6: 'hexa',
        7: 'hepta',
        8: 'octa',
        9: 'nona'
    }

    name = ''
    prefix = 0

    for i, bonded in enumerate(c):
        if bonded:
            prefix += 1
            name += iupac_numbers[i] + ','
    name = name[:-1]
    name += '-' + iupac_prefix[prefix] + 'chlorobiphenyl'

    return name


num_binding_carbons = 10
all_sets = itertools.product(range(2), repeat=num_binding_carbons)
unfiltered_congeners = [i for i in all_sets]
all_congeners = []

for congener in unfiltered_congeners:
    if x_flip(congener) not in all_congeners \
            and y_flip(congener) not in all_congeners \
            and xy_flip(congener) not in all_congeners \
            and x_flip(half_flip(congener)) not in all_congeners \
            and y_flip(half_flip(congener)) not in all_congeners \
            and xy_flip(half_flip(congener)) not in all_congeners \
            and half_flip(congener) not in all_congeners \
            and 1 in congener:
                all_congeners.append(congener)

for i, c in enumerate(all_congeners):
    congener = max(c, x_flip(c), y_flip(c), xy_flip(c), half_flip(c),
               x_flip(half_flip(c)), y_flip(half_flip(c)), xy_flip(half_flip(c)))
    all_congeners[i] = congener

num_cl_of_interest = [2, 8]  # Edit this list to get congeners with different numbers of Cl's
congeners_output = []

for c in all_congeners:
    count = 0
    for cl in c:
        count += cl
    if count in num_cl_of_interest:
        congeners_output.append(iupac_name(c))

congeners_output = sorted(congeners_output)
for c in congeners_output:
    print(c)
