#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 14:05:59 2023

@author: reidmarkland
"""

import temul.api as tml
import os
import numpy as np
import hyperspy.api as hs
import atomap.api as am
import atomap.tools 
from PolarizationFunctionality import plot_polarisation_vectors
import matplotlib.pyplot as plt

path = r'/Users/reidmarkland/Documents/Research/SULI23/Work For Nick/Electron Microscopy Analysis'
file = r'Slide3_1_bw.png'

image = hs.load(os.path.join(path,file))
fig = plt.figure(figsize=(5,5))

def GetLatticeA(image, separation):
    atom_positions_A = am.get_atom_positions(image, separation)
    sublattice_A = am.Sublattice(atom_positions_A, image = image.data)
    
    sublattice_A.find_nearest_neighbors()
    sublattice_A.refine_atom_positions_using_center_of_mass()
    sublattice_A.refine_atom_positions_using_2d_gaussian()
    
    sublattice_A.construct_zone_axes()
    
    mat_noA = atomap.tools.remove_atoms_from_image_using_2d_gaussian(sublattice_A.image, sublattice_A)
    image_noA = hs.signals.BaseSignal(mat_noA)
    
    return sublattice_A, image_noA



def GetLatticeB(image_noA, separation):

    atom_positions_B = am.get_atom_positions(image_noA, separation = separation)
    sublattice_B = am.Sublattice(atom_positions_B, image = image_noA.data, color='green')
    
    sublattice_B.construct_zone_axes()
    sublattice_B.refine_atom_positions_using_center_of_mass()
    sublattice_B.refine_atom_positions_using_2d_gaussian()
    
    # mat_bkgnd = atomap.tools.remove_atoms_from_image_using_2d_gaussian(sublattice_B.image, sublattice_B)
    # image_bkgnd = hs.signals.BaseSignal(mat_bkgnd)

    return sublattice_B#, image_bkgnd



def PolarizationMap(atom_lattice, image):
    sublattice_A = atom_lattice.sublattice_list[0]
    sublattice_B = atom_lattice.sublattice_list[1]
    
    za0, za1 = sublattice_A.zones_axis_average_distances[0:2]
    polarization = sublattice_A.get_polarization_from_second_sublattice(za0, za1, sublattice_B, color = 'cyan')
    vector_list = polarization.metadata.vector_list

    x, y = [i[0] for i in vector_list], [i[1] for i in vector_list]
    u, v = [i[2] for i in vector_list], [i[3] for i in vector_list]
    x, y, u, v = np.asarray(x), np.asarray(y), np.asarray(u), np.asarray(v)
    
    sampling = image.axes_manager[-1].scale
    units = image.axes_manager[-1].units
    save = None
    plot_style = 'colormap'
    overlay = True  # the vectors will be plotted on the image
    unit_vector = False  # formerly called normalise
    vector_rep = 'angle'  # 'magnitude' or 'angle'
    degrees = True  # Set to True for degrees, False for radians
    angle_offset = None
    title = ""
    # color = 'yellow'  # may be ignored depending on the plot_style
    cmap = 'hsv'  # may be ignored depending on the plot_style
    alpha = 1.0  # transparency of image or vectors, depending on plot_style
    image_cmap = 'gray'
    ticks = None
    scalebar = False
    monitor_dpi = 50  # set to ~200 to make too-large images a bit smaller
    no_axis_info = True  
    invert_y_axis = True
    antialiased = False  # relevant for the contour mapping
    levels = 20  # relevant for the contour mapping
    remove_vectors = False
    scale = 0.030  # set to 0.001-0.01 to change arrow size
    width = 0.008  # set to ~0.005 for chunky (thicker) arrows
    minshaft = 2
    minlength = 3
    headwidth = 5
    headlength = 3
    headaxislength = 3
    quiver_units = 'width'
    pivot = 'middle'
    angles = 'xy'
    scale_units = 'xy'


    fig, ax = plot_polarisation_vectors(
        x, y, u, v, image.data, sampling=sampling, units=units,
        plot_style=plot_style, overlay=overlay, unit_vector=unit_vector,
        vector_rep=vector_rep, degrees=degrees, angle_offset=angle_offset,
        save=save, title=title  , cmap=cmap, alpha=alpha,
        image_cmap=image_cmap, monitor_dpi=monitor_dpi,
        no_axis_info=no_axis_info, invert_y_axis=invert_y_axis, ticks=ticks,
        scalebar=scalebar, antialiased=antialiased, levels=levels,
        remove_vectors=remove_vectors, quiver_units=quiver_units, pivot=pivot,
        angles=angles, scale_units=scale_units, scale=scale, headwidth=headwidth,
        headlength=headlength, headaxislength=headaxislength, width=width,
        minshaft=minshaft, minlength=minlength)
    
    return fig, ax



# peaksA = am.get_feature_separation(image, separation_range=(32,35))
# peaksA.plot()

# separation = 20 # set based on lower limit of good fit from peaks^

# sublattice_A, image_noA = GetLatticeA(image, separation)
# # sublattice_A.plot()
# # # # # image_noA.plot()

# # peaksB = am.get_feature_separation(image_noA, separation_range=(32,35))
# # peaksB.plot()

# sublattice_B = GetLatticeB(image_noA, separation = 24)

# atom_lattice = am.Atom_Lattice(image = image, name = 'test', sublattice_list=[sublattice_A, sublattice_B])
# atom_lattice.plot()

# polarization, ax = PolarizationMap(atom_lattice, image)
# polarization.plot()




