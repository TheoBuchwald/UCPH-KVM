import argparse
from operators import E, BRA, t, amplitude
from math import factorial

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""A script designed to quickly get the code necessary to calculate a Coupled Cluster matrix element.""", epilog="""For help contact
    Theo Juncker von Buchwald
    tjvbu@kemi.dtu.dk""")

    parser.add_argument("bra", type=str, nargs=1, help="""The excitation level of the bra.
    For a left excitation vector write LI where I is the excitation level.
    For a left transformation write EI where I is the excitation level.
    Use HF for a Hartree-Fock state.""")
    parser.add_argument("commutator", type=str, nargs=1, help="""The commutator to expand...Ex. [[XT3]E2]
    For one-electron operators use F or X.
    For all types of amplitudes use TI where I is the excitation level.
    For a right transformation use EI where I is the excitation level.""")
    parser.add_argument('--no_t1', action='store_false', help='Include if t1-transformed.', dest='t1_transformed')
    parser.add_argument('--unrestricted', action='store_false', help='Include to use unrestricted box.', dest='restricted')
    parser.add_argument('--no-perm', action='store_false', help='Include to skip permutation check.', dest='perm_check')

    return parser.parse_args()

def zip_merge_arrays(arr1: list, arr2: list) -> list:
    """Zip merge two arrays."""
    merged_arr = []
    for i, j in zip(arr1, arr2):
        merged_arr.append(i)
        merged_arr.append(j)
    return merged_arr

def commutator_indexing(bra: str, commutator: str, ket: str) -> tuple[dict, int, int]:
    """Assign indices to components in the bra, commutator, and ket.

    In the following I will be used to indicate the position of an integer.
    args:
        bra: String representation of bra.
                LI is a left excitation vector of order I.
                EI is an excited bra of order I.
                HF is the Hartree Fock state.
        commutator: String representation of commutator.
                    H and P are accepted two-electron operators.
                    F and X are accepted one-electron operators.
                    EI is an excitation operator of order I.
                    TI is an amplitude of order I.
        ket: String representation of ket.
                Only HF is accepted.
        restricted: Whether the commutator should be indexed as if restricted or unrestricted.

    return:
        A dictionary containing:
            A factor.
            An indexed bra.
            The indexed amplitudes.
            The indexed excitation operators.
            The indexed summation.
            The ket.
        The number of virtual indices.
        The number of occupied indices.

    Examples of bra, commutator, and ket:
        - bra='HF', commutator='[[HT2]E1]', ket='HF'
        - bra='3', commutator='[FT3]', ket='HF'
        - bra='L2', commutator='[HT3]', ket='HF'
        - bra='2', commutator='[XT2]', ket='HF'
    """
    assert ket == "HF", "The matrix element script is not set up to have any other ket than the HF state."
    start_labelling_E = "E" in commutator
    start_labelling_bra = "E" in bra
    assert not (start_labelling_E and start_labelling_bra), "The bra and commutator cannot both contain an E operator."
    # Remove "[" from commutator string
    commutator = commutator.replace("[","")
    # Seperate nested commutators
    commutator_terms = commutator.split("]")
    if commutator_terms[-1] == "":
        commutator_terms.pop(-1)
    # If an E is present in the communtator
    # use this as stating point for index numbering
    virtual_index_counter = 0
    occupied_index_counter = 0
    summation = []
    factor = 1
    if start_labelling_E:
        for c, component in enumerate(commutator_terms):
            if "E" not in component:
                continue
            order = int(component.split("E")[-1])
            virtual_indices = [f"v{i}" for i in range(virtual_index_counter, virtual_index_counter + order)]
            occupied_indices = [f"o{i}" for i in range(occupied_index_counter, occupied_index_counter + order)]
            indices = zip_merge_arrays(virtual_indices, occupied_indices)
            indexed_E = E(indices)
            virtual_index_counter += order
            occupied_index_counter += order
            commutator_terms.pop(c)
            break
    else:
        indexed_E = E([])
    # If we start labelling by the bra
    if start_labelling_bra or "L" in bra:
        if "L" in bra:
            order = int(bra.replace("L",""))
        else:
            order = int(bra.replace("E",""))
        virtual_indices = [f"v{i}" for i in range(virtual_index_counter, virtual_index_counter + order)]
        occupied_indices = [f"o{i}" for i in range(occupied_index_counter, occupied_index_counter + order)]
        if "L" in bra:
            summation += virtual_indices
            summation += occupied_indices
        indices = zip_merge_arrays(virtual_indices, occupied_indices)
        indexed_bra = BRA(indices, "L" in bra)
        virtual_index_counter += order
        occupied_index_counter += order
    else:
        indexed_bra = BRA([], False)
    # Now we label all amplitudes
    t_indices = []
    for c, component in enumerate(commutator_terms):
        order = int(component.split("T")[-1])
        factor *= 1 / factorial(order)
        virtual_indices = [f"v{i}" for i in range(virtual_index_counter, virtual_index_counter + order)]
        occupied_indices = [f"o{i}" for i in range(occupied_index_counter, occupied_index_counter + order)]
        summation += virtual_indices
        summation += occupied_indices
        indices = zip_merge_arrays(virtual_indices, occupied_indices)
        t_indices.append(indices)
        virtual_index_counter += order
        occupied_index_counter += order
    indexed_t = t(t_indices)
    matrix_element = {
        'factor': factor,
        'summation': summation,
        'bra': indexed_bra,
        't': indexed_t,
        'E': indexed_E,
        'ket': ket,
    }
    return matrix_element, virtual_index_counter, occupied_index_counter

def main():
    """Main function."""
    arguments = parse_arguments()
    bra = arguments.bra[0]
    commutator = arguments.commutator[0]
    t1_transformed = arguments.t1_transformed
    restricted = arguments.restricted
    perm_check = arguments.perm_check
    one_electron = False
    one_electron_type = "F"
    if "F" in commutator:
        one_electron = True
    if "X" in commutator:
        one_electron = True
        one_electron_type = "X"
    matrix_element, virtual_index_counter, occupied_index_counter = commutator_indexing(bra, commutator, "HF")

if __name__ == "__main__":
    main()