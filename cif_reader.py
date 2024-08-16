from pymatgen.io.cif import CifParser
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.core.lattice import Lattice
import numpy as np
# from scipy.spatial import Voronoi


def calculate_reciprocal_lattice_vectors(real_lattice_vectors):
    # Assuming real_lattice_vectors is a 3x3 matrix where rows are a, b, c vectors
    # Calculate the volume of the unit cell
    volume = np.dot(real_lattice_vectors[0], np.cross(real_lattice_vectors[1], real_lattice_vectors[2]))

    # Calculate the reciprocal lattice vectors
    b1 = np.cross(real_lattice_vectors[1], real_lattice_vectors[2]) / volume
    b2 = np.cross(real_lattice_vectors[2], real_lattice_vectors[0]) / volume
    b3 = np.cross(real_lattice_vectors[0], real_lattice_vectors[1]) / volume

    return np.array([b1, b2, b3])


def generate_reciprocal_lattice_data(input_data, real_lattice_vectors):
    reciprocal_lattice_vectors = calculate_reciprocal_lattice_vectors(real_lattice_vectors)
    reciprocal_data = []
    input_data = [reciprocal_lattice_vectors[0], reciprocal_lattice_vectors[1], reciprocal_lattice_vectors[2]]

    for position in input_data:
        # Convert fractional position to reciprocal lattice Cartesian position
        # fractional_position = np.array(atom['cartesian_position'])
        # cartesian_position_reciprocal = np.dot(fractional_position, reciprocal_lattice_vectors)

        new_atom = {
            'atom_label': 'H',
            'oxi_atom_label': 'H',
            'fractional_position': 0,
            'cartesian_position': position, # New reciprocal lattice Cartesian position
            'connected_atoms': [],
            'magnetic_spin': {},
            'site_index': 0
        }
        reciprocal_data.append(new_atom)

    return reciprocal_data

# Example usage within the existing function
def get_structure_with_cif(file_path, num_unit_cells=None, is_primitive=False, target_atoms=None,
                           magnetic_spin_atoms=None, site_index_spin=None, use_reciprocal=False):
    if num_unit_cells is None:
        num_unit_cells = [1, 1, 1]

    parser = CifParser(file_path)
    structure = parser.parse_structures(primitive=is_primitive)[0]

    x_unit_cell, y_unit_cell, z_unit_cell = map(int, np.ceil(num_unit_cells))

    unique_atoms = []
    lattice = structure.lattice.matrix

    print(lattice)

    for nx in range(x_unit_cell):
        for ny in range(y_unit_cell):
            for nz in range(z_unit_cell):
                for idx, site in enumerate(structure):
                    fractional_coords = site.frac_coords + np.array([nx, ny, nz])
                    if all(fractional_coords[i] < num_unit_cells[i] for i in range(3)):
                        atom_label = site.species_string[:-2]
                        if target_atoms is None or atom_label in target_atoms:
                            atom_info = {
                                "atom_label": atom_label,
                                "oxi_atom_label": site.species_string,
                                "fractional_position": fractional_coords.tolist(),
                                "cartesian_position": (site.coords + np.dot([nx, ny, nz], lattice)).tolist(),
                                "connected_atoms": [],
                                "magnetic_spin": {},
                                "site_index": idx,
                            }
                            unique_atoms.append(atom_info)

    if magnetic_spin_atoms or site_index_spin:
        unique_atoms = add_magnetic_spin_info(unique_atoms, magnetic_spin_atoms, site_index_spin)

    # Convert positions to reciprocal lattice if use_reciprocal is True
    if use_reciprocal:
        unique_atoms = generate_reciprocal_lattice_data(unique_atoms, lattice)

    print(unique_atoms)

    return unique_atoms


def add_magnetic_spin_info(unique_atoms, magnetic_spin_atoms=None, site_index_spin=None):
    for atom in unique_atoms:
        atom_label = atom['atom_label']
        site_index = atom['site_index']
        if site_index_spin and site_index in site_index_spin:
            atom['magnetic_spin'] = {"direction": site_index_spin[site_index]}
        elif magnetic_spin_atoms and atom_label in magnetic_spin_atoms:
            atom['magnetic_spin'] = {"direction": magnetic_spin_atoms[atom_label]}
        else:
            atom['magnetic_spin'] = {"direction": [0, 0, 0]}  # No spin
    return unique_atoms

async def fetch_materials(**kwargs):
    api_key = "cSFVj0Awg9nQlro7yWhYacD4TRst78YZ"

    try:
        with MPRester(api_key) as mpr:
            results = mpr.summary.search(**kwargs)
            if not results:
                return "No data found with the given search parameters."
            return results
    except Exception as e:
        return f"Failed to fetch data: {str(e)}"

async def get_structure_with_api(structure, num_unit_cells=None, target_atoms=None):
    if num_unit_cells is None:
        num_unit_cells = [1, 1, 1]

    x_unit_cell, y_unit_cell, z_unit_cell = map(int, np.ceil(num_unit_cells))

    unique_atoms = []
    structure = structure[0].structure
    lattice = structure.lattice.matrix

    nn = CrystalNN()
    graph = StructureGraph.from_local_env_strategy(structure, nn)

    for nx in range(x_unit_cell):
        for ny in range(y_unit_cell):
            for nz in range(z_unit_cell):
                for idx, site in enumerate(structure):
                    fractional_coords = site.frac_coords + np.array([nx, ny, nz])
                    if all(fractional_coords[i] < num_unit_cells[i] for i in range(3)):
                        atom_label = site.species_string
                        if target_atoms is None or atom_label in target_atoms:
                            atom_info = {
                                "atom_label": atom_label,
                                "fractional_position": fractional_coords.tolist(),
                                "cartesian_position": (site.coords + np.dot([nx, ny, nz], lattice)).tolist(),
                                "connected_atoms": []
                            }
                            connected_sites = graph.get_connected_sites(idx)
                            for connected_site in connected_sites:
                                connected_fractional_coords = connected_site.site.frac_coords + np.array([nx, ny, nz])
                                connection = {
                                    "connected_to": connected_site.site.species_string,
                                    "bond_length": connected_site.weight,
                                    "connected_fractional_position": connected_fractional_coords.tolist(),
                                    "connected_cartesian_position": (
                                            connected_site.site.coords + np.dot([nx, ny, nz], lattice)).tolist(),
                                    "site_index": connected_site.index
                                }
                                atom_info["connected_atoms"].append(connection)

                            unique_atoms.append(atom_info)

    return unique_atoms

def bond_by_proximity(data, tolerance=0.1):
    positions = np.array([atom['cartesian_position'] for atom in data])
    distance_matrix = np.linalg.norm(positions[:, np.newaxis] - positions, axis=2)
    np.fill_diagonal(distance_matrix, np.inf)
    min_distances = np.min(distance_matrix, axis=1)
    nearest_neighbors_indices = np.argmin(distance_matrix, axis=1)

    for atom in data:
        atom['connected_atoms'] = []

    for index, atom in enumerate(data):
        nearest_distance = min_distances[index]
        close_indices = np.where((distance_matrix[index] <= nearest_distance * (1 + tolerance)) &
                                 (distance_matrix[index] >= nearest_distance * (1 - tolerance)))[0]

        for close_index in close_indices:
            if close_index != index:
                connection = {
                    'nearest_neighbor_index': close_index,
                    'connected_cartesian_position': data[close_index]['cartesian_position']
                }
                atom['connected_atoms'].append(connection)

    return data


import re

def bond_by_crystalNN(data, tolerance=0.1, weighting_factor=0.5):
    positions = np.array([atom['cartesian_position'] for atom in data])
    distance_matrix = np.linalg.norm(positions[:, np.newaxis] - positions, axis=2)
    np.fill_diagonal(distance_matrix, np.inf)

    for atom in data:
        atom['connected_atoms'] = []
        atom['selected_oxidation_state'] = None

    for index, atom in enumerate(data):
        element = atom['atom_label']
        oxi_atom_label = atom.get('oxi_atom_label', None)
        if oxi_atom_label:
            # Extract the oxidation state from the oxi_atom_label, e.g., 'Yb3+' -> 3
            oxidation_state_match = re.search(r'([+-]?\d+)', oxi_atom_label)
            if oxidation_state_match:
                oxidation_state = int(oxidation_state_match.group(1))
            else:
                oxidation_state = None
        else:
            oxidation_state = None

        if oxidation_state is None:
            # If oxidation state is not provided, you might skip or assign a default state
            continue  # Or assign a default oxidation state if appropriate

        min_distances = np.min(distance_matrix, axis=1)
        nearest_distance = min_distances[index]
        close_indices = np.where((distance_matrix[index] <= nearest_distance * (1 + tolerance)) &
                                 (distance_matrix[index] >= nearest_distance * (1 - tolerance)))[0]

        best_score = float('inf')
        best_connections = []

        weighted_distances = []
        current_connections = []

        for close_index in close_indices:
            if close_index == index:
                continue

            neighbor_oxi_atom_label = data[close_index].get('oxi_atom_label', None)
            if neighbor_oxi_atom_label:
                neighbor_oxidation_state_match = re.search(r'([+-]?\d+)', neighbor_oxi_atom_label)
                if neighbor_oxidation_state_match:
                    neighbor_oxidation_state = int(neighbor_oxidation_state_match.group(1))
                else:
                    neighbor_oxidation_state = None
            else:
                neighbor_oxidation_state = None

            if neighbor_oxidation_state is None:
                continue  # Skip neighbors with unknown oxidation states

            # Apply a more complex weighting based on the current oxidation state
            expected_coordination = get_expected_coordination(element, oxidation_state)
            actual_distance = distance_matrix[index, close_index]
            weighted_distance = actual_distance * weighting_factor / expected_coordination

            connection = {
                'nearest_neighbor_index': close_index,
                'connected_cartesian_position': data[close_index]['cartesian_position'],
                'weighted_distance': weighted_distance
            }
            current_connections.append(connection)

        # Calculate the score for this oxidation state (e.g., sum of weighted distances)
        score = sum(conn['weighted_distance'] for conn in current_connections)
        if score < best_score:
            best_score = score
            best_connections = current_connections

        # Assign the best connections to the atom
        atom['connected_atoms'] = best_connections
        atom['selected_oxidation_state'] = oxidation_state

    return data

def get_expected_coordination(element, oxidation_state):
    # A dictionary mapping (element, oxidation_state) to their typical coordination numbers
    coordination_dict = {
        # Hydrogen
        ('H', 1): 1,
        ('H', -1): 1,

        # Oxygen
        ('O', -2): 2,
        ('O', -1): 1,
        ('O', 0): 2,

        # Nitrogen
        ('N', -3): 4,
        ('N', 0): 3,
        ('N', 1): 4,
        ('N', 2): 3,
        ('N', 3): 4,
        ('N', 4): 6,
        ('N', 5): 3,

        # Carbon
        ('C', -4): 4,
        ('C', -3): 4,
        ('C', -2): 4,
        ('C', -1): 3,
        ('C', 0): 4,
        ('C', 1): 3,
        ('C', 2): 4,
        ('C', 3): 4,
        ('C', 4): 4,

        # Sulfur
        ('S', -2): 2,
        ('S', 0): 2,
        ('S', 2): 4,
        ('S', 4): 4,
        ('S', 6): 4,

        # Chlorine
        ('Cl', -1): 6,
        ('Cl', 0): 2,
        ('Cl', 1): 6,
        ('Cl', 3): 4,
        ('Cl', 5): 4,
        ('Cl', 7): 4,

        # Sodium
        ('Na', 1): 6,

        # Magnesium
        ('Mg', 2): 6,

        # Aluminum
        ('Al', 3): 6,

        # Silicon
        ('Si', 4): 4,

        # Phosphorus
        ('P', -3): 3,
        ('P', 0): 3,
        ('P', 3): 4,
        ('P', 5): 4,

        # Potassium
        ('K', 1): 6,

        # Calcium
        ('Ca', 2): 6,

        # Titanium
        ('Ti', 3): 6,
        ('Ti', 4): 6,

        # Chromium
        ('Cr', 2): 6,
        ('Cr', 3): 6,
        ('Cr', 6): 4,

        # Iron
        ('Fe', 2): 6,
        ('Fe', 3): 6,

        # Nickel
        ('Ni', 2): 6,
        ('Ni', 3): 6,

        # Copper
        ('Cu', 1): 2,
        ('Cu', 2): 4,

        # Zinc
        ('Zn', 2): 4,

        # Silver
        ('Ag', 1): 2,

        # Tin
        ('Sn', 2): 6,
        ('Sn', 4): 6,

        # Lead
        ('Pb', 2): 6,
        ('Pb', 4): 6,

        # Gold
        ('Au', 1): 2,
        ('Au', 3): 4,

        # Mercury
        ('Hg', 1): 2,
        ('Hg', 2): 2,

        # Platinum
        ('Pt', 2): 4,
        ('Pt', 4): 6,

        # Additional elements and their oxidation states could be added here
    }
    return coordination_dict.get((element, oxidation_state), 1)  # Default to 1 if unknown

