"""
EDENetworks, a genetic network analyzer
Copyright (C) 2011  Aalto University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
# Files for calculating information theoretic values:
#   entropy_X_Y(P)
#   entropy_XY_Z(P)
#   mutual_info_XY_Z(P)

import numpy as np

def entropy_X(P):
    """Return the entropy of P. P must sum to 1."""
    P = np.array(P)
    P = np.ma.masked_array(P, P==0)
    return -np.inner(np.log2(P), P)

def entropy_X_Y(P):
    """
    Calculate entropy H(X|Y) with probabilites p(x,y) = P[x,y].
    """
    P = np.ma.masked_array(P, P==0)
    dim = np.shape(P)
    P_Y = np.sum(P, 0)

    return -np.sum( P * np.log2( P/np.tile(P_Y, (dim[0], 1)) ) )

def entropy_XY_Z(P):
    """
    Calculate entropy H(X,Y|Z) with probabilites p(x,y,z) = P[x,y,z].
    """
    P = np.ma.masked_array(P, P==0)
    dim = np.shape(P)
    P_Z = np.sum(np.sum(P, 0), 0)

    return -np.sum( P * np.log2( P/np.tile(P_Z, (dim[0], dim[1], 1)) ) )


def mutual_info_XY_Z(P):
    """
    Calculate mutual information I(x,y|z) with probabilities
    p(x,y,z) = P[x,y,z].
    """
    # Calculate necessary entropies
    H_X_Z = entropy_X_Y(np.sum(P, 1))
    H_Y_Z = entropy_X_Y(np.sum(P, 0))
    H_XY_Z = entropy_XY_Z(P)

    #print H_X_Z, H_Y_Z, H_XY_Z
    return H_X_Z + H_Y_Z - H_XY_Z


if __name__ == '__main__':
    """Run unit tests if called."""
    from tests.test_info_theory import *
    unittest.main()
