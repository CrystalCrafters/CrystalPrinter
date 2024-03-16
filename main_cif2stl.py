# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 13:15:08 2023

@author: PhysicsGavin
"""

import cif2stl as cs
import trimesh

filepath = 'Yb2Si2O7.cif'

desired_atom = 'Yb'

atom_data = cs.LoadCif(filepath)

spheres, des_atom_coords = cs.DrawSpheres(atom_data, desired_atom, 2, 1, 2)


cylinders = cs.CreateBond(des_atom_coords, num_nn = 4)
# # Combine all sphere meshes
combined = trimesh.util.concatenate(spheres, cylinders)

# # Save the mesh to a .stl file
combined.export('cif2stl_ybsio.stl')