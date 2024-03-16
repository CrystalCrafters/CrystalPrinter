import trimesh
import numpy as np
from ase import io


def LoadCif(filepath):
    '''This function reads the .cif file and holds the data in an ase "Atoms" data type.
    '''
    atoms_data = io.read(filepath) # Open the file and read the contents
    
    return atoms_data



def draw_sphere_at_position(x, y, z, radius=0.5, color=[0, 255, 0, 255]):
    """
    Create a sphere mesh at position (x, y, z).
    
    Parameters:
    x, y, z (float): 3D coordinates where the sphere is to be drawn.
    radius (float): The radius of the sphere.
    color (list): The RGBA color of the sphere.
    
    Returns:
    trimesh.Trimesh: Sphere mesh at specified position.
    """
    # create a sphere mesh
    sphere = trimesh.creation.uv_sphere(radius=radius)
    sphere.visual.vertex_colors = color
    
    # translate the sphere to desired location
    translation = np.array([x, y, z])
    sphere.apply_translation(translation)
    
    return sphere




def BondMeasure(coordinates):
    '''Parameters:
        coordinates: the x,y,z coordinates of atoms
    '''
    
    
    distance = []
    i_vals = []
    j_vals = []

    for i in range(len(coordinates)):
        for j in range(i+1,len(coordinates)):
            if i != j:
                calc_distance = np.sqrt((coordinates[i][0] - coordinates[j][0])**2 + (coordinates[i][1] - coordinates[j][1])**2 + (coordinates[i][2] -  coordinates[j][2])**2)
            
                i_vals.append(i)
                j_vals.append(j)
                distance.append(calc_distance)
            else:
                continue

            
    sort_order = np.argsort(distance)

    distance_sorted = []
    i_sorted = []
    j_sorted = []


    for i in range(len(distance)):
        distance_sorted.append(distance[sort_order[i]])
        i_sorted.append(i_vals[sort_order[i]])
        j_sorted.append(j_vals[sort_order[i]])

    sorted_coords = np.array([i_sorted, j_sorted, distance_sorted])
    return sorted_coords




def CreateBond(atom_coords, num_nn):
    '''parameters:
        atoms_coords: A list of 3 element lists for the atomic coordinates. This is one of the outputs of DrawSpheres
        num_nn: The number of nearest neighbors to create bonds for. This may 
        be larger then the nominal value one would expect due to rounding errors in the calculation of the distances.
        Just the nearest neighbor bond would be 0.
    '''
    
    lengths = BondMeasure(atom_coords)
    NN_lengths = np.unique(lengths[2,:])
    print('Number of nearest neighbor sets: ', len(NN_lengths))
    cylinders = []
        
        
    for v in range(num_nn):
        if num_nn > len(NN_lengths): #THis needs to turn into an actual error
            print('ERROR: You have requested a larger number of nearest neighbors then was detected. No cylinders drawn.')
            break
        current_length = NN_lengths[v]
        
        for i in range(len(atom_coords)):
            for j in range(len(atom_coords)):
                if i != j:
                
                    pos1 = atom_coords[i]
                    pos2 = atom_coords[j]
                    
                    
                    
                    vec = [pos2[0] - pos1[0], pos2[1] - pos1[1], pos2[2] - pos1[2]]
                    
                    len_vec = np.linalg.norm(vec) #getting the length of the vector between the spheres,
                    #the distance between the spheres.
                    
                    if len_vec <= 1.01*current_length and len_vec >= 0.99*current_length:
                        #The above if statement makes sure that the bond is the correct length to be drawn. 
                        #It allows for a 1% difference between the locally calculated value and the value caluculated by BondMeasure function
                        
                        
                        midpoint_vec = [(pos2[0] - pos1[0])/2, (pos2[1] - pos1[1])/2, (pos2[2] - pos1[2])/2]
                        midpoint_vec = [midpoint_vec[0] + pos1[0], midpoint_vec[1] + pos1[1], midpoint_vec[2] + pos1[2]]
                        #The above step is necessary because the translation function wants to translate relative to the origin,
                        #but the first line of calculation for the midpoint gives the midpoint relative to the first atom.
                        #We add the position of atom 1 to make the midpoint relative to the origin.
                
                    
                
                        vec = vec / len_vec #making the vector unitary in length
                
                        cyl_direc = [0, 0, 1] #The cylinders are generically drawn vertically
                
                
                        cylinder = trimesh.creation.cylinder(radius=0.1, height= len_vec) #Creating the cylinder.
                
                
                        #Finding a vector perpendicular to the cylinder and the vector connecting the atoms. 
                        #This will be the axis to rotate about to get the cylinder orientation correct.
                        perp_vec = np.cross(vec, cyl_direc)
                        
                        if np.linalg.norm(perp_vec) != 0: #This covers the limited times that the desired bond is parallel to the default cylinder direction. 
                            perp_vec = perp_vec / np.linalg.norm(perp_vec)
                            full_ang = -1 * np.arccos(np.dot(vec, cyl_direc) / np.linalg.norm(cyl_direc))
                
                            rotat_mat = trimesh.transformations.rotation_matrix(full_ang, perp_vec, point=None) #first argument is angle in radians, second argument is the direction to rotate around, third is what point to rotate around
                            cylinder.apply_transform(rotat_mat) 
                
                        #The cylinder should now be in the right place, we just need to move the cylinder to the midpoint of the bond.
                
                        cylinder.apply_translation(midpoint_vec) #moving the cylinder to the midpoint.
                        
                        cylinders.append(cylinder)
                    
    
    return cylinders




def DrawSpheres(atoms_data, des_Atom, num_x_ucells, num_y_ucells, num_z_ucells):
    ''' parameters:
            atoms_data: the ase object containing all the atom information
            des_atom: the two letter symbol for the atoms to draw. Currently only one can be provided.
            num_x_ucells: the number of unit cells in the x direction (normally a-axis)
            num_y_ucells: the number of unit cells in the y direction.
            num_z_ucells: the number of unit cells in the z direction.
        
        improvements:
            - make number of unit cells optional.
            - expand to allow for multiple atom types.
    '''

    atom_pos = atoms_data.get_positions()
    atom_symbol = atoms_data.get_chemical_symbols()
    cell_params = atoms_data.get_cell_lengths_and_angles()
   
    atom_coords = []
    spheres = [] #preallocate space
    
    for jj in range(len(atom_pos)): #looping through all possible atoms
        if atom_symbol[jj] == des_Atom: #Checking if the atom is the one we want to draw.
            for nn in range(0,num_x_ucells):
                for mm in range(0,num_y_ucells):
                    for oo in range(0,num_z_ucells):
                        #Drawing the spheres at each position.
                        spheres.append(draw_sphere_at_position(atom_pos[jj,0] + nn*cell_params[0] , atom_pos[jj,1] + mm*cell_params[1], atom_pos[jj,2] + oo*cell_params[2]))
                        atom_coords.append([atom_pos[jj,0] + nn*cell_params[0] , atom_pos[jj,1] + mm*cell_params[1], atom_pos[jj,2] + oo*cell_params[2]])
        elif des_Atom == '':
            for nn in range(0,num_x_ucells):
                for mm in range(0,num_y_ucells):
                    for oo in range(0,num_z_ucells):
                        #Drawing the spheres at each position.
                        spheres.append(draw_sphere_at_position(atom_pos[jj,0] + nn*cell_params[0] , atom_pos[jj,1] + mm*cell_params[1], atom_pos[jj,2] + oo*cell_params[2]))
            
    return spheres, atom_coords

