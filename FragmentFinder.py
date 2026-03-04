"""
Fragment Finder

This module provides utilities to identify and extract a common fragment from a set of
molecules.  From `.xyz` file and a specificity level, it identifies
atom indices belonging to a fragment that is shared among molecules in the same
directory.
The workflow broadly consists of:

1. Reading molecules from `.xyz` files and computing their chemical connectivity.
2. Selecting a fragment interactively from a reference molecule using a 3D viewer.
3. Finding matches of this fragment in other molecules by graph isomorphism.
4. Reporting atoms of interest and their neighbors for downstream calculations.
"""

import os
import csv
from pathlib import Path
import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism
from ase.io import read
from vedo import Sphere, Tube, Plotter, Text3D, Assembly, Text2D
from collections import Counter
from ase.neighborlist import NeighborList, natural_cutoffs

###############################################################################
# Section 1: Connectivity and Graph
###############################################################################
import numpy as np
from ase.data import covalent_radii, atomic_numbers

# Default maximum valences (conservative for standard organic chemistry)
MAX_VALENCE_DEF = {
    'H': 1, 'C': 4, 'N': 3, 'O': 2,
    'F': 1, 'Cl': 1, 'Br': 1, 'I': 1,
    'P': 5, 'S': 6,
}

# Pair-specific scales
PAIR_SCALE_DEF = {
    ('H','O'): 1.05,
    ('O','H'): 1.05,
    ('H','N'): 1.05,
    ('N','H'): 1.05,
    ('H','C'): 1.08,
    ('C','H'): 1.08,
    ('C','O'): 1.08,
    ('O','C'): 1.08,
    ('C','C'): 1.10,
    ('O','O'): 1.06,
    ('N','O'): 1.08,
    ('O','N'): 1.08,
    ('N','N'): 1.08,
}

def calculate_connectivity_matrix(
    mol,
    base_scale: float = 1.10,
    pair_scale: dict | None = None,
    max_valence: dict | None = None,
    allow_HH: bool = False,
    iter_max: int = 4,
    debug: bool = False
):

    if pair_scale is None:
        pair_scale = PAIR_SCALE_DEF
    if max_valence is None:
        max_valence = MAX_VALENCE_DEF

    Z = mol.get_atomic_numbers()
    S = mol.get_chemical_symbols()
    n = len(mol)

    D = mol.get_all_distances(mic=True)
    np.fill_diagonal(D, np.inf)

    R = np.array([covalent_radii[z] for z in Z])

    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        ri = R[i]
        si = S[i]
        for j in range(i+1, n):
            sj = S[j]
            # skip H-H contacts unless allowed
            if not allow_HH and si == 'H' and sj == 'H':
                continue
            rj = R[j]
            scale = pair_scale.get((si, sj), base_scale)
            threshold = scale * (ri + rj)
            if D[i, j] <= threshold:
                A[i, j] = 1
                A[j, i] = 1

    # iterative pruning to enforce maximum valence
    for _ in range(iter_max):
        changes = 0
        for i in range(n):
            vmax = max_valence.get(S[i], 4)
            neighbors = np.where(A[i] == 1)[0]
            deg = len(neighbors)
            if deg > vmax:
                # sort neighbors by increasing distance and keep the nearest ones
                order = neighbors[np.argsort(D[i, neighbors])]
                to_remove = order[vmax:]
                for j in to_remove:
                    if A[i, j] == 1:
                        A[i, j] = 0
                        A[j, i] = 0
                        changes += 1
                if debug and len(to_remove) > 0:
                    kept = order[:vmax]
                    print(f"[VALENCE] {S[i]}{i+1}: deg={deg}>vmax={vmax} -> "
                          f"keep {list(kept+1)}, drop {list(to_remove+1)}")
        if changes == 0:
            break

    return A

def matrix_to_graph(matrix, symbols):
    """Convert a connectivity matrix and a list of symbols into a NetworkX graph."""
    G = nx.Graph()
    n = len(matrix)
    for i in range(n):
        G.add_node(i, label=symbols[i])
    for i in range(n):
        for j in range(i + 1, n):
            if matrix[i, j] == 1:
                G.add_edge(i, j)
    return G


###############################################################################
# Section 2: Fragmentation and Matching
###############################################################################
def remove_duplicate_matches(matches: list[list[int]]) -> list[list[int]]:
    """
    Remove duplicate fragment matches regardless of ordering.
    """
    unique_matches: list[list[int]] = []
    seen: set[frozenset[int]] = set()
    for match in matches:
        key = frozenset(match)
        if key not in seen:
            seen.add(key)
            unique_matches.append(match)
    return unique_matches


def match_fragment(molecule_matrix: np.ndarray, fragment_matrix: np.ndarray,
                   molecule_symbols: list[str], fragment_symbols: list[str]) -> list[list[int]]:
    """
    Find all occurrences of a fragment within a molecule using graph isomorphism.
    """
    G_molecule = matrix_to_graph(molecule_matrix, molecule_symbols)
    G_fragment = matrix_to_graph(fragment_matrix, fragment_symbols)
    GM = isomorphism.GraphMatcher(
        G_molecule, G_fragment,
        node_match=lambda n1, n2: n1['label'] == n2['label']
    )
    all_matches: list[list[int]] = []
    for mapping in GM.subgraph_isomorphisms_iter():
        ordered = [mol_index for mol_index, frag_index in sorted(mapping.items(), key=lambda kv: kv[1])]
        all_matches.append(ordered)
    return remove_duplicate_matches(all_matches)


###############################################################################
# Section 3: Atom and Neighbor Information
###############################################################################
def print_unique_atoms_with_neighbors(unique_atoms: list[tuple[int, str]], connectivity_matrix: np.ndarray,
                                      molecule_symbols: list[str]) -> tuple[dict, dict]:
    """
    Print the list of unique atoms along with their neighbors and return dictionaries
    with this information.
    """
    num_neighbors_dict: dict = {}
    neighbor_dict: dict = {}
    print("Unique atoms:", unique_atoms)
    for idx, (atom, symbol) in enumerate(unique_atoms):
        neighbors = [n for n in range(len(connectivity_matrix)) if connectivity_matrix[atom - 1, n] == 1]
        neighbor_atoms = ", ".join(f"{n + 1}({molecule_symbols[n]})" for n in neighbors)
        print(f"{atom}({symbol}) \t {neighbor_atoms}")
        num_neighbors_dict[f"{idx}{molecule_symbols[atom - 1]}"] = [len(neighbors)]
        neighbor_dict[f"{atom}({molecule_symbols[atom - 1]})"] = [f"{n + 1}({molecule_symbols[n]})" for n in neighbors]
    return num_neighbors_dict, neighbor_dict


def include_neighbors(unique_atoms: list[tuple[int, str]], connectivity_matrix: np.ndarray,
                      molecule_symbols: list[str]) -> list[tuple[int, str]]:
    """
    Given a list of unique atoms, include all of their neighbors and return the updated list.
    """
    new_atoms = set(unique_atoms)
    for atom, _ in unique_atoms:
        neighbors = [n for n in range(len(connectivity_matrix)) if connectivity_matrix[atom - 1, n] == 1]
        for neighbor in neighbors:
            new_atoms.add((neighbor + 1, molecule_symbols[neighbor]))
    return sorted(new_atoms)


def calculate_fragment_connectivity_matrix(selected_atoms: list[tuple[int, str]], molecule_matrix: np.ndarray):
    """
    Compute the connectivity matrix for the selected fragment.
    """
    indices = [atom - 1 for atom, _ in selected_atoms]
    fragment_matrix = molecule_matrix[np.ix_(indices, indices)]
    return fragment_matrix, indices


def calculate_neighbor_counts(molecule_matrix: np.ndarray, indices: list[int], fragment_symbols: list[str],
                              molecule_symbols: list[str]) -> tuple[dict, dict]:
    """
    Calculate the number of neighbors for each atom in a fragment match.
    """
    adjusted_indices = [i + 1 for i in indices]
    num_neighbors_dict: dict = {}
    neighbor_dict: dict = {}
    for idx, atom_index in enumerate(adjusted_indices):
        neighbors = [n for n in range(len(molecule_matrix)) if molecule_matrix[atom_index - 1, n] == 1]
        num_neighbors_dict[f"{idx}{fragment_symbols[idx]}"] = [len(neighbors)]
        neighbor_dict[f"{atom_index}({fragment_symbols[idx]})"] = [f"{n + 1}({molecule_symbols[n]})" for n in neighbors]
    return num_neighbors_dict, neighbor_dict


###############################################################################
# Section 4: Reading and Searching Molecules
###############################################################################
def read_molecules_from_xyz_folder(folder: str, mol: str):
    """
    Read XYZ molecules from the specified folder.
    """
    molecules: list[tuple[str, np.ndarray, list[str], any]] = []
    if mol == 'none':
        for filename in os.listdir(folder):
            if filename.endswith('.xyz'):
                path = os.path.join(folder, filename)
                molecule = read(path)
                matrix = calculate_connectivity_matrix(molecule)
                symbols = molecule.get_chemical_symbols()
                molecules.append((filename, matrix, symbols, molecule))
    else:
        molecule = read(mol)
        matrix = calculate_connectivity_matrix(molecule)
        symbols = molecule.get_chemical_symbols()
        molecules.append((mol, matrix, symbols, molecule))
    return molecules


def search_fragment_in_molecules(molecules: list[tuple[str, np.ndarray, list[str], any]],
                                fragment_matrix: np.ndarray,
                                fragment_symbols: list[str]) -> tuple[list, list, list]:
    """
    Search for a fragment in each molecule and return matches.
    """
    results = []
    found = []
    not_found = []
    for name, molecule_matrix, molecule_symbols, _ in molecules:
        matches = match_fragment(molecule_matrix, fragment_matrix, molecule_symbols, fragment_symbols)
        if matches:
            for match in matches:
                fragment_in_molecule = [molecule_symbols[idx] for idx in match]
                num_dict, neighbor_dict = calculate_neighbor_counts(molecule_matrix, match, fragment_in_molecule,
                                                                    molecule_symbols)
                results.append((name, match, fragment_in_molecule, num_dict, neighbor_dict))
            found.append(name)
        else:
            not_found.append(name)
    return results, found, not_found


###############################################################################
# Section 5: 3D Molecular Graphics
###############################################################################
def get_element_color(symbol: str) -> str:
    """
    Return a display color for a given chemical element symbol.  Default is ochre.
    """
    colors = {
        'H': '#FFFFFF',
        'C': '#B0B0B0',
        'O': 'red',
        'N': 'navy',
        'Cl': 'limegreen',
        'Br': 'darkorange',
        'P': '#FFA500',
        'F': '#DDA0DD',
        'S': '#CCCC00',
        'I': 'purple'
    }
    return colors.get(symbol, '#CC7722')


def get_element_radius(symbol: str) -> float:
    """
    Return a display radius for a given chemical element symbol.
    Radii are scaled relative to hydrogen as the smallest and heavier elements are larger.
    """
    radii = {
        'H': 0.3,
        'O': 0.35,
        'C': 0.4,
        'N': 0.4,
        'S': 0.4,
        'F': 0.4,
        'Cl': 0.5,
        'Br': 0.6,
        'P': 0.6,
        'I': 0.6
    }
    return radii.get(symbol, 0.4)


def select_atoms_interactive(molecule):
    """
    Visualize a molecule in 3D and allow interactive atom selection.

    Mouse click toggles selection (highlighted in pink).  Keyboard shortcuts:

    * `e` – toggle atom labels
    * `n` – include neighbors of the current selection based on the connectivity matrix
    * `m` – clear the selection
    * `q` – close the window

    """
    positions = molecule.get_positions()
    symbols = molecule.get_chemical_symbols()
    n_atoms = len(symbols)

    # Compute chemical connectivity matrix
    A = calculate_connectivity_matrix(molecule)

    # Bonds (using A)
    bonds = []
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            if A[i, j] == 1:
                bonds.append(Tube([positions[i], positions[j]], r=0.08, c='white'))

    # Sphere + label
    atom_assemblies = []
    for i, (pos, sim) in enumerate(zip(positions, symbols)):
        radius = get_element_radius(sim)
        z_offset = radius + 0.1
        color = get_element_color(sim)

        sp = Sphere(pos=pos, r=radius, c=color).lighting('glossy')
        sp.pickable(True); sp.idx = i

        label_str = f"{sim}{i+1}"
        txt = Text3D(label_str, pos=(pos[0], pos[1], pos[2]+z_offset), s=0.2, c='black', justify='center')
        txt.follow_camera(); txt.lighting('off'); txt.pickable(False); txt.alpha(0)

        assembly = Assembly(sp, txt)
        assembly.pickable(True); assembly.idx = i
        atom_assemblies.append(assembly)

    # Scene
    plt = Plotter(axes=0, title="FragmentFinder")
    info_text = Text2D("Shortcuts: e=labels  n=neighbors  m=clear  q=exit",
                       pos="top-left", c='white', bg='black', alpha=0.7)
    texto_info = Text2D("", pos="bottom-left", c='white', bg='black', alpha=0.7)
    plt.add(info_text)

    selected = []
    labels_visible = False

    # Click
    def callback_click(evt):
        if not evt.actor or not hasattr(evt.actor, 'idx'):
            return
        idx = evt.actor.idx
        sphere = atom_assemblies[idx].unpack(0)
        if idx in selected:
            # deselect: restore the original color
            sphere.color(get_element_color(symbols[idx]))
            selected.remove(idx)
        else:
            # select: highlight in pink
            sphere.color('hotpink')
            selected.append(idx)
        texto_info.text(f"Selected atoms: {[k+1 for k in selected]}")
        plt.render()

    plt.add_callback("mouse click", callback_click)

    # Keyboard
    def key_pressed(evt):
        nonlocal labels_visible, selected
        k = (evt.keypress or "").lower()
        if k == "e":
            labels_visible = not labels_visible
            for assembly in atom_assemblies:
                lab = assembly.unpack(1)
                lab.alpha(1 if labels_visible else 0)
            plt.render()

        elif k == "n":
            if selected:
                # include neighbors of the current selection
                current = set(selected)
                for i in list(current):
                    neighbors = np.where(A[i] == 1)[0]
                    current.update(neighbors.tolist())
                selected[:] = sorted(current)
                for a in atom_assemblies:
                    sphere = a.unpack(0)
                    sphere.color('hotpink' if a.idx in selected else get_element_color(symbols[a.idx]))
                texto_info.text(f"Selected atoms: {[k+1 for k in selected]}")
                plt.render()

        elif k == "m":
            selected[:] = []
            for a in atom_assemblies:
                a.unpack(0).color(get_element_color(symbols[a.idx]))
            texto_info.text("Selected atoms: []")
            plt.render()

        elif k == "q":
            plt.close()

    plt.add_callback("key press", key_pressed)

    # Add and show
    plt.add(bonds); plt.add(atom_assemblies); plt.add(texto_info)
    plt.show(resetcam=True, interactive=True)
    return selected


def select_interest_fragment(molecule, fragment_indices):
    """
    Interface to choose atoms of interest within a fragment (0-based indices).

    Interactively select atoms within the previously selected fragment.  The
    controls are similar to :func:`select_atoms_interactive` but only operate
    within the fragment:

    * Mouse click – toggle selection (pink)
    * `e` – toggle labels
    * `n` – include neighbors within the fragment
    * `m` – clear selection
    * `q` – close the window

    """
    positions = molecule.get_positions()
    symbols = molecule.get_chemical_symbols()

    # Compute global connectivity and restrict to the fragment
    A = calculate_connectivity_matrix(molecule)
    frag_set = set(fragment_indices)

    # Bonds only within the fragment
    bonds = []
    for i in fragment_indices:
        for j in fragment_indices:
            if j > i and A[i, j] == 1:
                bonds.append(Tube([positions[i], positions[j]], r=0.05, c='gray'))

    atom_assemblies = []
    for i in fragment_indices:
        sim = symbols[i]
        radius = get_element_radius(sim)
        z_offset = radius + 0.1
        color0 = get_element_color(sim)

        sp = Sphere(pos=(0,0,0), r=radius, c=color0).lighting('glossy')
        sp.pickable(True); sp.idx = i

        txt = Text3D(f"{sim}{i+1}", pos=(0,0,z_offset), s=0.2, c='black', justify='center')
        txt.follow_camera(); txt.lighting('off'); txt.pickable(False); txt.alpha(0)

        ass = Assembly(sp, txt)
        ass.pickable(True)
        ass.idx = i
        ass.original_color = color0
        ass.pos(positions[i])
        atom_assemblies.append(ass)

    # Scene
    plt = Plotter(axes=0, title="Select atoms of interest (click). 'q' to exit")
    info_text = Text2D("Shortcuts: e=labels  n=neighbors  m=clear  q=exit",
                       pos="top-left", c='white', bg='black', alpha=0.7)
    texto_info = Text2D("", pos="bottom-left", c='white', bg='black', alpha=0.7)
    plt.add(info_text)

    selected = []
    labels_visible = False

    def callback_click(evt):
        if not evt.actor or not hasattr(evt.actor, 'idx'):
            return
        idx = evt.actor.idx
        # toggle selection only within the fragment
        for ass in atom_assemblies:
            if ass.idx == idx:
                sp = ass.unpack(0)
                if idx in selected:
                    sp.color(ass.original_color)
                    selected.remove(idx)
                else:
                    sp.color('hotpink')
                    selected.append(idx)
                break
        texto_info.text(f"Selected atoms: {[i+1 for i in selected]}")
        plt.render()

    plt.add_callback("mouse click", callback_click)

    def key_pressed(evt):
        nonlocal labels_visible, selected
        k = (evt.keypress or "").lower()
        if k == 'e':
            labels_visible = not labels_visible
            for ass in atom_assemblies:
                lab = ass.unpack(1)
                lab.alpha(1 if labels_visible else 0)
            plt.render()
        elif k == 'n':
            if selected:
                current = set(selected)
                for i in list(current):
                    neighbors = np.where(A[i] == 1)[0]
                    # only include neighbors that are also in the fragment
                    current.update([v for v in neighbors if v in frag_set])
                selected[:] = sorted(current)
                for ass in atom_assemblies:
                    ass.unpack(0).color('hotpink' if ass.idx in selected else ass.original_color)
                texto_info.text(f"Selected atoms: {[i+1 for i in selected]}")
                plt.render()
        elif k == 'm':
            selected[:] = []
            for ass in atom_assemblies:
                ass.unpack(0).color(ass.original_color)
            texto_info.text("Selected atoms: []")
            plt.render()
        elif k == 'q':
            plt.close()

    plt.add_callback("key press", key_pressed)

    plt.add(bonds)
    plt.add(atom_assemblies)
    plt.add(texto_info)
    plt.show(resetcam=True, interactive=True)
    return selected


def neighbor_count_signature(dic_num: dict) -> tuple:
    """
    Convert a dictionary like ``{'0C':[3], '1N':[2], ...}`` into a permutation-invariant signature.

    Each key encodes the match index and the element symbol (e.g. ``'0C'`` for the
    first carbon atom) and the value is a list containing the number of neighbors.
    """
    def symbol_from_key(k: str) -> str:
        i = 0
        while i < len(k) and k[i].isdigit():
            i += 1
        return k[i:]

    signature = Counter()
    for k, v in dic_num.items():
        symbol = symbol_from_key(k)
        degree = v[0] if isinstance(v, (list, tuple)) else int(v)
        signature[(symbol, degree)] += 1
    # Return as a sorted tuple to make it comparable
    return tuple(sorted(signature.items()))


def main(fragment_matrix: np.ndarray, fragment_symbols: list[str], directory: str,
         req: str) -> tuple[list, list, list]:
    """
    Search for a fragment within molecules located in ``directory``.
    """
    if req == 'all':
        mols = read_molecules_from_xyz_folder(directory, mol='none')
    else:
        mols = read_molecules_from_xyz_folder(directory, req)
    res, found, not_found = search_fragment_in_molecules(mols, fragment_matrix, fragment_symbols)
    return res, found, not_found


def generate_csv_report(molecule_counts: dict, directory: str, filename: str = "fragment_counts.csv"):
    """
    Generate a CSV report listing all processed molecules and the number of
    times the fragment was found in each.
    """
    filepath = Path(directory) / filename
    try:
        with open(filepath, mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerow(["Molecule", "Fragment_Count"])
            
            # Sort for better readability
            for mol in sorted(molecule_counts.keys()):
                writer.writerow([mol, molecule_counts[mol]])
                
        print(f"\n--> CSV report generated: {filepath.name}")
    except Exception as e:
        print(f"\nError writing CSV report: {e}")

#############################################
# Function start()
#############################################
def start(file_path: str, specificity: str, req: str = 'all') -> tuple[dict, list, dict]:
    """
    Interactively select a fragment from a reference molecule and search for it in other molecules.

    This function presents a 3D viewer for the reference molecule located at ``file_path``.
    The user selects atoms to define a fragment and optionally expands the selection to include
    neighbors.  The fragment connectivity is extracted and used to find occurrences of the
    fragment in other molecules within the same directory.  Atoms of interest within the
    fragment can be selected for subsequent property calculations.

    Parameters
    ----------
    file_path : str
        Path to the ``.xyz`` file of the reference molecule.
    specificity : str
        Either ``'0'`` or ``'1'``.  If ``'1'``, fragment matches must have the same degree
        distribution signature as the reference fragment.  Use ``'0'`` for a more relaxed match.
    req : str, optional
        Either the name of a specific ``.xyz`` file to search or ``'all'`` to search all
        molecules in the directory.  Defaults to ``'all'``.

    Returns
    -------
    tuple
        ``(results_dict, atoms_of_interest, neighbor_dict)`` where ``results_dict`` maps
        molecule names to match information, ``atoms_of_interest`` is a list of labels for
        selected atoms within the fragment, and ``neighbor_dict`` maps selected atom labels to
        their neighbors in the reference molecule.
    """
    # Validate the input file and get base data
    directory = Path(file_path).parent
    reference_molecule = read(file_path)
    molecule_matrix = calculate_connectivity_matrix(reference_molecule)
    molecule_symbols = reference_molecule.get_chemical_symbols()

    print("Connectivity matrix of the reference molecule:")
    print(molecule_matrix)
    print("Atom symbols of the reference molecule:")
    print(molecule_symbols)

    # Define the base fragment
    while True:
        print("\nSelect atoms in the 3D view (press 'n' to include neighbors):")
        selected = select_atoms_interactive(reference_molecule)
        if not selected:
            print("You must select at least one atom of the fragment. Please try again.")
            continue
        # unique_atoms: list of tuples (1-based index, symbol)
        unique_atoms = [(i + 1, molecule_symbols[i]) for i in selected]
        print("\nSelected fragment:")
        for idx, sym in unique_atoms:
            print(f"{idx}: {sym}")
        break

    # Display information about the base fragment
    neighbor_counts_dict, neighbor_dict_interest = print_unique_atoms_with_neighbors(
        unique_atoms, molecule_matrix, molecule_symbols)
    base_signature = neighbor_count_signature(neighbor_counts_dict)
    # Compute the connectivity matrix of the selected fragment
    new_fragment_matrix, fragment_indices = calculate_fragment_connectivity_matrix(unique_atoms, molecule_matrix)
    print("\nConnectivity matrix of the selected fragment:")
    print(new_fragment_matrix)
    print("\nList of atoms in the fragment:")
    for i, (atom, sym) in enumerate(unique_atoms):
        print(f"{i + 1}: {atom}({sym})")

    fragment_symbols = [sym for _, sym in unique_atoms]
    # Find the fragment in the same molecule to get a canonical ordering
    initial_matches = match_fragment(
        molecule_matrix,
        new_fragment_matrix,
        molecule_symbols,
        fragment_symbols
    )
    if not initial_matches:
        raise RuntimeError("Could not find the fragment in the original molecule.")
    # Take the first match
    fragment_indices = initial_matches[0]
    fragment_labels = [
        f"{idx + 1}({molecule_symbols[idx]})"
        for idx in fragment_indices
    ]
    print("\nBase fragment defined in the first molecule (graph order):")
    print(fragment_labels)

    # Selection of atoms of interest within the base fragment
    while True:
        print("\nSelect atoms of interest within the base fragment (highlighted in pink).")
        interest_indices = select_interest_fragment(reference_molecule, fragment_indices)
        if not interest_indices:
            print("You must select at least one atom of interest. Please try again.")
            continue
        break
    try:
        interest_rel = [fragment_indices.index(i) for i in interest_indices if i in fragment_indices]
    except ValueError:
        print("Error: Some selected atoms are not in the base fragment. Check your selection.")
        return None

    # Build the list of labels for atoms of interest using the fragment order
    atoms_of_interest = [fragment_labels[j] for j in interest_rel]
    print("\nSelected atoms of interest:")
    print(atoms_of_interest)

    # Results container
    results_dict: dict = {}

    # Search only in the same XYZ file if req == 'none'
    _req = file_path if str(req).lower() == 'none' else req

    results, found, not_found = main(
        new_fragment_matrix,
        [molecule_symbols[idx] for idx in fragment_indices],
        str(directory),
        req=_req
    )

    def record_match(name: str, indices: list[int], fragment_in_molecule: list[str],
                     dic_num: dict, dic_vec: dict, interest_rel_idx: list[int]):
        real_indices = [i + 1 for i in indices]
        full_fragment_labels = [f"{real_indices[i]}({fragment_in_molecule[i]})" for i in range(len(real_indices))]
        selected_fragment_labels = [full_fragment_labels[i] for i in interest_rel_idx if i < len(full_fragment_labels)]
        selected_fragment_indices = [real_indices[i] for i in interest_rel_idx if i < len(real_indices)]
        ordered_neighbors = {
            label: dic_vec.get(label, [])
            for label in selected_fragment_labels
        }
        key = Path(name).stem
        results_dict.setdefault(key, []).append({
            'fragment_indices': real_indices,
            'fragment_atoms': fragment_in_molecule,
            'selected_atoms': selected_fragment_labels,
            'neighbor_dict': ordered_neighbors,
            'interest_atom_indices': selected_fragment_indices
        })
        print(f"\n---> Molecule: {key}")
        print("Fragment atoms (full):", full_fragment_labels)
        print("Selected fragment atoms:", selected_fragment_labels)
        print(f"Neighbor dictionary: {ordered_neighbors}")

    print("\nResults of the fragment search:")
    matched_files: set[str] = set()
    for name, indices, frag_in_molecule, dic_num, dic_vec in results:
        if specificity == '1':
            if neighbor_count_signature(dic_num) != base_signature:
                continue  # discard if the signature (symbol, degree) does not match
        record_match(name, indices, frag_in_molecule, dic_num, dic_vec, interest_rel)
        matched_files.add(name)

    print(f"\nFragment found in {len(matched_files)} file(s).")
    if not_found:
        print("Not found in:", not_found)

    # Prepare data for CSV export
    molecule_counts = {}
    for name in found + not_found:
        molecule_counts[name] = 0
    for name, _, _, _, _ in results:
        molecule_counts[name] += 1

    generate_csv_report(molecule_counts, str(directory))

    return results_dict, atoms_of_interest, neighbor_dict_interest


if __name__ == '__main__':
    # Entry point when running this module directly.  Prompt the user for a reference `.xyz` file and a specificity
    # level, then initiate the fragment search.
    while True:
        file_path = input("Enter the path to the .xyz file of the reference molecule: ").strip()
        if not os.path.isfile(file_path):
            print("The file does not exist. Please enter a valid path.")
            continue
        break

    # Prompt for specificity (only 0 or 1)
    while True:
        specificity = input("Enter the specificity level (0 or 1): ").strip()
        if specificity not in ["0", "1"]:
            print("Only 0 or 1 are accepted. Please try again.")
            continue
        break

    start(file_path, specificity, req='all')
