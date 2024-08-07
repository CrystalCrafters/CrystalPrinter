from pymatgen.io.cif import CifParser
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.graphs import StructureGraph
from mp_api.client import MPRester
import numpy as np

def get_structure_with_cif(file_path, num_unit_cells=None, is_primitive=False, target_atoms=None,
                           magnetic_spin_atoms=None, site_index_spin=None):
    if num_unit_cells is None:
        num_unit_cells = [1, 1, 1]

    parser = CifParser(file_path)
    structure = parser.parse_structures(primitive=is_primitive)[0]

    x_unit_cell, y_unit_cell, z_unit_cell = map(int, np.ceil(num_unit_cells))

    unique_atoms = []
    lattice = structure.lattice.matrix
    nn = CrystalNN()
    graph = StructureGraph.from_local_env_strategy(structure, nn)

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
                                "site_index": idx
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

    if magnetic_spin_atoms or site_index_spin:
        unique_atoms = add_magnetic_spin_info(unique_atoms, magnetic_spin_atoms, site_index_spin)

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

def bond_by_nearest_neighbors(data, tolerance=0.1):
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
