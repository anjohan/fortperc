# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:30:51 2013

@author: svenni

Modified by Anders Johansson
"""

import numpy as np


def walk(z):
    #
    # Left turning walker
    #
    # Returns left: nr of times walker passes a site
    #
    # First, ensure that array only has one contact point at left and
    # right : topmost points chosen
    #
    nx = z.shape[0]
    ny = z.shape[1]
    i = np.where(z > 0)
    ix0 = 0  # starting row for walker is always 0
    # starting row
    # (first element where there is a matching column which is zero)
    iy0 = i[1][np.where(i[0] == 0)][0]
    # First do left-turning walker
    directions = np.zeros((4, 2), int)
    directions[0, 0] = -1  # west
    directions[0, 1] = 0
    directions[1, 0] = 0  # south
    directions[1, 1] = -1
    directions[2, 0] = 1  # east
    directions[2, 1] = 0
    directions[3, 0] = 0  # north
    directions[3, 1] = 1
    nwalk = 1
    ix = ix0
    iy = iy0
    direction = 0  # 0 = west, 1 = south, 2 = east, 3 = north
    left = np.zeros((nx, ny), int)
    right = np.zeros((nx, ny), int)
    while (nwalk > 0):
        left[ix, iy] = left[ix, iy] + 1
        # Turn left until you find an occupied site
        nfound = 0
        while (nfound == 0):
            direction = direction - 1
            if (direction < 0):
                direction = direction + 4

            # Check this direction
            iix = ix + directions[direction, 0]
            iiy = iy + directions[direction, 1]
            if (iix >= nx):
                nwalk = 0  # Walker escaped
                nfound = 1
                iix = nx
                ix1 = ix
                iy1 = iy

            # Is there a site here?
            elif (iix >= 0):
                if (iiy >= 0):
                    if (iiy < ny):
                        if (z[iix, iiy] >
                                0):  # there is a site here, move here
                            ix = iix
                            iy = iiy
                            nfound = 1
                            direction = direction + 2
                            if (direction > 3):
                                direction = direction - 4

    # left
    nwalk = 1
    ix = ix0
    iy = iy0
    direction = 1  # 1=left, 2 = down, 3 = right, 4 = up
    while (nwalk > 0):
        right[ix, iy] = right[ix, iy] + 1
        # ix,iy
        # Turn right until you find an occupied site
        nfound = 0
        while (nfound == 0):
            direction = direction + 1
            if (direction > 3):
                direction = direction - 4

            # Check this directionection
            iix = ix + directions[direction, 0]
            iiy = iy + directions[direction, 1]
            if (iix >= nx):
                if (iy >= iy1):
                    nwalk = 0  # Walker escaped
                    nfound = 1
                    iix = nx

            # Is there a site here?
            elif (iix >= 0):
                if (iiy >= 0):
                    if (iiy < ny):
                        if (iix < nx):
                            if (z[iix, iiy] >
                                    0):  # there is a site here, move here
                                ix = iix
                                iy = iiy
                                nfound = 1
                                direction = direction - 2
                                if (direction < 0):
                                    direction = direction + 4

    return left, right
