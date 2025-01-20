import argparse
import box_13_2
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

def commutator_expansion(matrix_element: dict[str: t | E | BRA | str | list[str]]) -> list[dict]:
    """Expand commutator after indexing.

    args:
        matrix_element: Dictionary containing the bra, ket, summation and any amplitudes, and excitation operators.

    return:
        All commutators resulting from the commutator expanded matrix_element. Each commutator contains:
            The summation.
            The bra.
            The amplitudes.
            The excitation operators.
            The ket.
            The commutator_factor, pre_string, and commutator_string defined as seen below.

                        2 * E_ai [H, E_bj]
                        ^   ^    ^commutator_string
                        ^   ^pre_string
                        ^commutator_factor
        """
    assert "t" in matrix_element.keys() or "E" in matrix_element.keys(), "Commutator expansion needs either a T or an E in the commutator."
    commutator_string: list[list] = []
    commutator_factor = [matrix_element["factor"]]
    pre_string = [E([])]

    amplitudes: list[amplitude] = matrix_element["t"].amplitudes
    E_indices: list[str] = matrix_element["E"].indices
    for amp in amplitudes:
        commutator_string.append(amp.indices)
    if E_indices != []:
        commutator_string.append(E_indices)
    commutator_string = [commutator_string]
    for nr in range(4):
        new_commutator_string: list[list] = []
        new_commutator_factor: list[int] = []
        new_pre_string: list[list] = []
        for commutator, factor, string in zip(commutator_string, commutator_factor, pre_string):
            # If we have run over all components of the commutator it is done
            if len(commutator) < nr+1:
                new_commutator_string.append(commutator)
                new_commutator_factor.append(factor)
                new_pre_string.append(string)
                continue
            component = commutator[nr]
            # If there is only one excitation operator in the current component it is (so far) done
            if len(component) == 2:
                new_commutator_string.append(commutator)
                new_commutator_factor.append(factor)
                new_pre_string.append(string)
                continue
            # Start with the expansion of an amplitude
            # This commutator expansion is the one where excitation operators are added in front of the commutator (to the pre_string)
            if component[0] not in E_indices:
                new_commutator_string.append([comp if i != nr else comp[:2] for i, comp in enumerate(commutator)])
                new_pre_string.append(string + E(component[2:]))
                new_commutator_factor.append(factor * (len(component) // 2))
                # This commutator expansion is the one where the commutator is increased in size
                if len(commutator) < 4:
                    new_commutator_string.append([comp if i != nr else comp[:2] for i, comp in enumerate(commutator)])
                    new_commutator_string[-1].append(component[2:])
                    new_pre_string.append(string)
                    new_commutator_factor.append(factor)
                continue
            # Then the expansion of an excitation operator
            # This commutator expansion is the one where excitation operators are added in front of the commutator (to the pre_string)
            for pair_index in range(0, len(component)-1, 2):
                new_commutator_string.append([comp if i != nr else comp[pair_index:pair_index+2] for i, comp in enumerate(commutator)])
                if pair_index == 0:
                    new_pre_string.append(string + E(component[2:]))
                elif pair_index == len(component)-2:
                    new_pre_string.append(string + E(component[:-2]))
                else:
                    new_pre_string.append(string + E(component[:pair_index]) + E(component[pair_index+2:]))
                new_commutator_factor.append(factor)
                # This commutator expansion is the one where the commutator is increased in size
                if len(commutator) < 4:
                    new_commutator_string.append([comp if i != nr else comp[pair_index:pair_index+2] for i, comp in enumerate(commutator)])
                    if pair_index == 0:
                        new_commutator_string[-1].append(component[2:])
                        new_pre_string.append(string)
                    elif pair_index == len(component)-2:
                        new_commutator_string.pop(-1)
                        continue
                    else:
                        new_commutator_string[-1].append(component[pair_index+2:])
                        new_pre_string.append(E(component[:pair_index]))
                    new_commutator_factor.append(factor)
        commutator_string = new_commutator_string
        commutator_factor = new_commutator_factor
        pre_string = new_pre_string

    new_matrix_elements = []
    for commutator, factor, string in zip(commutator_string, commutator_factor, pre_string):
        new_matrix_elements.append({
            'summation': matrix_element["summation"],
            'bra': matrix_element["bra"],
            't': matrix_element["t"],
            'E': matrix_element["E"],
            'ket': matrix_element["ket"],
            'factor': factor,
            'pre_string': string,
            'commutator': commutator,
        })

    return new_matrix_elements

def commutator_box(matrix_elements: list[dict], virtual_index_counter: int, occupied_index_counter: int, restricted: bool, t1_transformed: bool, one_electron: bool) -> list[dict]:
    """Translate an expanded commutator using Box 13.2 and the (unreleased) Unrestricted Box.

    args:
        matrix_elements: The commutator expanded matrix_elements.
        virtual_index_counter: The value of the next virtual index should one be needed.
        occupied_index_counter: The value of the next occupied index should one be needed.
        restricted: Whether the calculation is restricted or not.
        t1_transformed: Whether the one- and two-electron operators are T1-transformed or not.
        one_electron: Whether the first operator in the commutator is a one-electron operator.

    return:
        Non reduced mathematical expressions for each matrix element by using Box 13.2 and the (unreleased) unrestricted box:
            factor: The prefactor for the expression.
            summation: The indices that are summed over.
            bra: The bra.
            pre_string: The excitation operators left over from the commutator expansion.
            permutation: The permutation operator.
            box_E: The excitation operators resulting from using Box 13.2.
            ket: The ket.
            integrals: The integrals resulting from using Box 13.2.
            t: The amplitudes.
            E: The original excitation operator.
    """
    mathematical_expressions: list[dict] = []
    for element in matrix_elements:
        commutator = element["commutator"]
        pre_string = element["pre_string"]
        bra: BRA = element["bra"]
        nesting = len(commutator)
        minus_shift = len(pre_string)
        bra_rank = len(bra)
        order = bra_rank - minus_shift
        # Skip element if there are to few or to many excitation operators
        if order > 2:
            continue
        if order < 0:
            continue
        if abs(nesting - order) > 2:
            continue
        translater = box_13_2.get_translater(nesting, order, t1_transformed, one_electron, restricted)
        mathematical_expressions += translater(element, virtual_index_counter, occupied_index_counter)
    return mathematical_expressions

def perform_permutations(mathematical_expressions: list[dict]) -> list[list[dict]]:
    """Perform the permutation operators in the mathematical expressions.

    args:
        mathematical_expressions: Non reduced mathematical expressions.

    return:
        All expressions after the permutation operators have been applied:
            factor: The prefactor for the expression.
            summation: The indices that are summed over.
            bra: The bra.
            pre_string: The excitation operators left over from the commutator expansion.
            ket: The ket.
            integrals: The integrals resulting from using Box 13.2.
            t: The amplitudes.
            E: The original excitation operator.
    """
    permuted_expressions: list[list[dict]] = []
    for expression in mathematical_expressions:
        # Skip elements without permutation operators
        if "permutation" not in expression.keys():
            permuted_expressions.append([expression])
            continue
        permutation_list = []
        # Do permutations
        for old_indices, new_indices in expression["permutation"]:
            perm_express = {
                "bra": expression["bra"],
                "ket": expression["ket"],
                "factor": expression["factor"],
                "summation": expression["summation"],
            }
            perm_express["t"] = expression["t"]
            perm_express["pre_string"] = expression["pre_string"] + expression["box_E"].update_indices(old_indices, new_indices)
            perm_express["E"] = expression["E"]
            perm_express["integrals"] = expression["integrals"].update_indices(old_indices, new_indices)
            permutation_list.append(perm_express)
        permuted_expressions.append(permutation_list)
    return permuted_expressions

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
    expanded_matrix_element = commutator_expansion(matrix_element)
    mathematical_expressions = commutator_box(expanded_matrix_element, virtual_index_counter, occupied_index_counter, restricted, t1_transformed, one_electron)
    permuted_expressions = perform_permutations(mathematical_expressions)

if __name__ == "__main__":
    main()