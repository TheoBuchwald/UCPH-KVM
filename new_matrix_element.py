import argparse
import box_13_2

from copy import deepcopy
from operators import E, BRA, t, amplitude, P
from math import factorial
from itertools import permutations

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""A script designed to quickly get the code necessary to calculate a Coupled Cluster matrix element.
The input must contain a bra and a commutator. Examples being:
    L2 [[HT2]T2]
    E3 [FT3]
    HF [HE2]""", epilog="""For help contact
    Theo Juncker von Buchwald
    tjvbu@kemi.dtu.dk"""
    )

    parser.add_argument("bra", type=str, nargs=1, help="""The excitation level of the bra. Example: L4
    For a left excitation vector write LI where I is the excitation level.
    For a left transformation write EI where I is the excitation level.
    Use HF for a Hartree-Fock state.""")
    parser.add_argument("commutator", type=str, nargs=1, help="""The commutator to expand. Example: [[XT3]E2]
    For one-electron operators use F or X.
    For two-electron operators use P.
    For both one- and two-electron operators use H.
    For all types of amplitudes use TI where I is the excitation level.
    For a right transformation use EI where I is the excitation level.""")
    parser.add_argument('--no_t1', action='store_false', help='Disable t1-transformation.', dest='t1_transformed')
    parser.add_argument('--unrestricted', action='store_false', help='Include to use unrestricted box.', dest='restricted')
    parser.add_argument('--no-perm', action='store_false', help='Disable permutation check.', dest='perm_check')
    parser.add_argument('--explicit-sym', action='store_true', help='Do explicit symmetrization.', dest='explicit_sym')

    return parser.parse_args()

def zip_merge_arrays(arr1: list, arr2: list) -> list:
    """Zip merge two arrays."""
    merged_arr = []
    for i, j in zip(arr1, arr2):
        merged_arr.append(i)
        merged_arr.append(j)
    return merged_arr

def progressbar(iterable, prefix):
    """Progressbar taken from https://stackoverflow.com/questions/3160699/python-progress-bar.

    args:
        iterable: Iterable.
        prefix: Text before the progressbar.

    return:
        A generator over iterable
    """
    count = len(iterable)
    size = 60
    def show(j):
        if j == 0:
            x = 0
        else:
            x = int(size*j/count)
        print(f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count}", end='\r', flush=True)
    show(0) # avoid div/0
    for i, item in enumerate(iterable):
        yield item
        show(i+1)
    print("", flush=True)

def commutator_indexing(bra: str, commutator: str, ket: str) -> tuple[dict, int, int]:
    """Assign indices to components in the bra, commutator, and ket.

    In the following I will be used to indicate the position of an integer.
    args:
        bra: String representation of bra.
                LI is a left excitation vector of order I.
                EI is an excited bra of order I.
                HF is the Hartree Fock state.
        commutator: String representation of commutator.
                    H is the accepted one- + two-electron operator.
                    P is the accepted two-electron operator.
                    F and X are the accepted one-electron operators.
                    EI is an excitation operator of order I.
                    TI is an amplitude of order I.
        ket: String representation of ket.
                Only HF is accepted.

    return:
        A dictionary containing:
            A symmetrization operator.
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
        - bra='E3', commutator='[FT3]', ket='HF'
        - bra='L2', commutator='[HT3]', ket='HF'
        - bra='E2', commutator='[XT2]', ket='HF'
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
    indexed_E = E([])
    symmetrization_operator = P([])
    if start_labelling_E:
        for c, component in enumerate(commutator_terms):
            if "E" not in component:
                continue
            order = int(component.split("E")[-1])
            virtual_indices = [f"v{i}" for i in range(virtual_index_counter, virtual_index_counter + order)]
            occupied_indices = [f"o{i}" for i in range(occupied_index_counter, occupied_index_counter + order)]
            indices = zip_merge_arrays(virtual_indices, occupied_indices)
            indexed_E = E(indices)
            symmetrization_operator = P(indices)
            virtual_index_counter += order
            occupied_index_counter += order
            commutator_terms.pop(c)
            break
    # If we start labelling by the bra
    if start_labelling_bra or "L" in bra:
        if "L" in bra:
            order = int(bra.replace("L",""))
            if start_labelling_E:
                factor *= 1 / factorial(order)
        else:
            order = int(bra.replace("E",""))
        virtual_indices = [f"v{i}" for i in range(virtual_index_counter, virtual_index_counter + order)]
        occupied_indices = [f"o{i}" for i in range(occupied_index_counter, occupied_index_counter + order)]
        if "L" in bra:
            summation += virtual_indices
            summation += occupied_indices
        indices = zip_merge_arrays(virtual_indices, occupied_indices)
        if "E" in bra:
            symmetrization_operator = P(indices)
        indexed_bra = BRA(indices, "L" in bra)
        virtual_index_counter += order
        occupied_index_counter += order
    else:
        indexed_bra = BRA([], False)
    # Now we label all amplitudes
    t_indices = []
    for c, component in enumerate(commutator_terms):
        if "T" in component:
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
        'symmetry_operator': symmetrization_operator,
        'factor': factor,
        'summation': summation,
        'bra': indexed_bra,
        't': indexed_t,
        'E': indexed_E,
        'ket': ket,
    }
    return matrix_element, virtual_index_counter, occupied_index_counter

def print_matrix_element(matrix_element: dict[str: t | E | BRA | str | list[str]], operator: str) -> None:
    """Print the matrix element after indexing.

    args:
        matrix_element: Dictionary containing the bra, ket, summation and any amplitudes, and excitation operators.
        operator: String representation of the operator used in the commutator.
    """
    indexed_matrix_element = translate_to_normal_indices([matrix_element])[0]
    commutator = "[" * len(indexed_matrix_element["t"])
    if indexed_matrix_element["E"] != E([]):
        commutator += "["
    commutator += operator
    if indexed_matrix_element["t"] != t([]):
        for amp in indexed_matrix_element["t"].amplitudes:
            commutator += f", {str(amp)[:-1]}]"
    if indexed_matrix_element["E"] != E([]):
        commutator += f", {str(indexed_matrix_element['E'])[:-1]}]"
    left_vector = ""
    if indexed_matrix_element["bra"].left_excitation_vector:
        left_vector = f" \\bar{{t}}^{{{''.join(indexed_matrix_element['bra'].indices[::2])}}}_{{{''.join(indexed_matrix_element['bra'].indices[1::2])}}}"
    print("Working on the matrix element")
    print(f"\sum_{{{''.join(indexed_matrix_element['summation'])}}}{left_vector} {indexed_matrix_element['bra']}{commutator}|HF>")
    print("")

def commutator_expansion(matrix_element: dict[str: t | E | BRA | str | list[str]]) -> list[dict]:
    """Expand commutator after indexing.

    args:
        matrix_element: Dictionary containing the bra, ket, summation and any amplitudes, and excitation operators.

    return:
        All commutators resulting from the commutator expanded matrix_element. Each commutator contains:
            The symmetrization operator.
            The factor.
            The summation.
            The bra.
            The amplitudes.
            The excitation operators.
            The ket.

    The commutator_factor, pre_string, and commutator_string defined as below.

                            2 * E_ai [H, E_bj]
                            ^   ^    ^commutator_string
                            ^   ^pre_string
                            ^commutator_factor
        """
    assert "t" in matrix_element.keys() or "E" in matrix_element.keys(), "Commutator expansion needs either a T or an E in the commutator."
    commutator_string: list[list] = []
    commutator_factor = [matrix_element["factor"]]
    pre_string = [E([])]

    print("Expanding commutator expression.")
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
            component_length = len(component) // 2
            # If there is only one excitation operator in the current component it is (so far) done
            if component_length == 1:
                new_commutator_string.append(commutator)
                new_commutator_factor.append(factor)
                new_pre_string.append(string)
                continue
            # This commutator expansion is the one where excitation operators are added in front of the commutator (to the pre_string)
            new_commutator_string.append([comp if i != nr else comp[:2] for i, comp in enumerate(commutator)])
            new_pre_string.append(string + E(component[2:]))
            new_commutator_factor.append(factor * component_length)
            # This commutator expansion is the one where the commutator is increased in size
            if len(commutator) < 4:
                new_commutator_string.append([comp if i != nr else comp[:2] for i, comp in enumerate(commutator)])
                new_commutator_string[-1].append(component[2:])
                new_pre_string.append(string)
                new_commutator_factor.append(factor)
                if component_length == 2:
                    continue
                for expansion in range(4, 2*component_length, 2):
                    new_commutator_string.append([comp if i != nr else comp[:2] for i, comp in enumerate(commutator)])
                    new_commutator_string[-1].append(component[2:expansion])
                    new_pre_string.append(string + E(component[expansion:]))
                    new_commutator_factor.append(factor)
        commutator_string = new_commutator_string
        commutator_factor = new_commutator_factor
        pre_string = new_pre_string

    new_matrix_elements = []
    for commutator, factor, string in zip(commutator_string, commutator_factor, pre_string):
        new_matrix_elements.append({
            'symmetry_operator': matrix_element["symmetry_operator"],
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

def collect_commutator_terms(matrix_elements: list[dict]) -> list[dict]:
    """Collect equal commutator terms from the commutator expansion.

    args:
        matrix_elements: The commutator expanded matrix_elements.

    return:
        Unique commutator expanded matrix elements.
    """
    collected_elements = []
    for e, element in enumerate(progressbar(matrix_elements[::-1], f"{'Collecting commutator terms: ':<50}"), 1):
        found_match = False
        for comparison in matrix_elements[:-e]:
            if element["pre_string"] != comparison["pre_string"]:
                continue
            if element["commutator"] != comparison["commutator"]:
                continue
            found_match = True
            comparison["factor"] += element["factor"]
        if not found_match:
            collected_elements.append(element)
    return collected_elements


def commutator_box(matrix_elements: list[dict], virtual_index_counter: int, occupied_index_counter: int, restricted: bool, t1_transformed: bool, one_electron: bool, two_electron: bool) -> list[dict]:
    """Translate an expanded commutator using Box 13.2 and the (unreleased) Unrestricted Box.

    args:
        matrix_elements: The commutator expanded matrix_elements.
        virtual_index_counter: The value of the next virtual index should one be needed.
        occupied_index_counter: The value of the next occupied index should one be needed.
        restricted: Whether the calculation is restricted or not.
        t1_transformed: Whether the one- and two-electron operators are T1-transformed or not.
        one_electron: Whether the first operator in the commutator is a one-electron operator.
        two_eelectron: Whether the first operator in the commutator is a two-electron operator.

    return:
        Non reduced mathematical expressions for each matrix element by using Box 13.2 and the (unreleased) unrestricted box:
            symmetry_operator: The symmetrization operator.
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
    for element in progressbar(matrix_elements, f"{'Translating matrix elements using Box 13.2: ':<50}"):
        commutator = element["commutator"]
        pre_string = element["pre_string"]
        bra: BRA = element["bra"]
        nesting = len(commutator)
        minus_shift = len(pre_string)
        bra_rank = len(bra)
        order = bra_rank - minus_shift
        # Skip element if there are too few or too many excitation operators
        if order > 2:
            continue
        if order < 0:
            continue
        if abs(nesting - order) > 2:
            continue
        translater = box_13_2.get_translater(nesting, order, t1_transformed, one_electron, two_electron, restricted)
        mathematical_expressions += translater(element, virtual_index_counter, occupied_index_counter)
    return mathematical_expressions

def perform_permutations(mathematical_expressions: list[dict]) -> list[list[dict]]:
    """Perform the permutation operators in the mathematical expressions.

    args:
        mathematical_expressions: Non reduced mathematical expressions.

    return:
        All expressions after the permutation operators have been applied:
            symmetry_operator: The symmetrization operator.
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
    for expression in progressbar(mathematical_expressions, f"{'Permuting elements: ':<50}"):
        # Skip elements without permutation operators
        if "permutation" not in expression.keys():
            perm_express = deepcopy(expression)
            perm_express["pre_string"] += expression.pop("box_E")
            permuted_expressions.append([perm_express])
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
            perm_express["symmetry_operator"] = expression["symmetry_operator"]
            perm_express["t"] = expression["t"]
            perm_express["pre_string"] = expression["pre_string"] + expression["box_E"].update_indices(old_indices, new_indices)
            perm_express["E"] = expression["E"]
            perm_express["integrals"] = expression["integrals"].update_indices(old_indices, new_indices)
            permutation_list.append(perm_express)
        permuted_expressions.append(permutation_list)
    return permuted_expressions

def match_reduce_indices(permuted_expressions: list[list[dict]]) -> list[list[dict]]:
    """Match and reduce indices in each mathematical expression separately.

    args:
        List of permuted mathematical expressions.

    return:
        List of expressions where indicies have been matched and reduced:
            symmetry_operator: The symmetrization operator.
            factor: The prefactor for the expression.
            summation: The indices that are summed over.
            bra: The bra.
            ket: The ket.
            integrals: The integrals resulting from using Box 13.2.
            t: The amplitudes.
            E: The original excitation operator.
    """
    # Start by matching indices in the excitation operators to the bra
    matched_permuted_expressions = []
    for permutation_list in progressbar(permuted_expressions, f"{'Matching indices: ':<50}"):
        matched_permutation_list = []
        for expression in permutation_list:
            summation: list[str] = deepcopy(expression["summation"])
            E_indices = expression["E"].indices
            bra_indices = expression["bra"].indices
            bra_virtual = [i for i in bra_indices if i[0] == "v"]
            bra_occupied = [i for i in bra_indices if i[0] == "o"]
            match_from = []
            match_to = []
            for a, i in expression["pre_string"]:
                if a in E_indices:
                    old_index = bra_virtual.pop(0)
                    match_from.append(old_index)
                    match_to.append(a)
                    summation.remove(old_index)
                else:
                    new_index = bra_virtual.pop(0)
                    match_from.append(a)
                    match_to.append(new_index)
                    summation.remove(a)
                if i in E_indices:
                    old_index = bra_occupied.pop(0)
                    match_from.append(old_index)
                    match_to.append(i)
                    summation.remove(old_index)
                else:
                    new_index = bra_occupied.pop(0)
                    match_from.append(i)
                    match_to.append(new_index)
                    summation.remove(i)
            matched_permutation_list.append({
                "symmetry_operator": expression["symmetry_operator"],
                "bra": expression["bra"].update_indices(match_from, match_to),
                "ket": expression["ket"],
                "factor": expression["factor"],
                "summation": summation,
                "t": expression["t"].update_indices(match_from, match_to),
                "integrals": expression["integrals"].update_indices(match_from, match_to),
                "E": expression["E"],
            })
        matched_permuted_expressions.append(matched_permutation_list)
    # Then reduce the indices
    reduced_permuted_expressions = []
    for permutation_list in progressbar(matched_permuted_expressions, f"{'Reducing indices: ':<50}"):
        reduced_permutation_list = []
        for expression in permutation_list:
            summation: list[str] = deepcopy(expression["summation"])
            used_indices = []
            used_indices += expression["t"].indices
            used_indices += expression["integrals"].indices
            used_indices += expression["bra"].indices
            used_virtual_indices = list({i for i in used_indices if i[0] == "v"})
            used_occupied_indices = list({i for i in used_indices if i[0] == "o"})
            used_virtual_indices.sort(key=lambda numerical: int(numerical[1:]))
            used_occupied_indices.sort(key=lambda numerical: int(numerical[1:]))
            number_virtual = len(used_virtual_indices)
            number_occupied = len(used_occupied_indices)
            reduce_from = used_virtual_indices + used_occupied_indices
            reduce_to = [f"v{i}" for i in range(number_virtual)] + [f"o{i}" for i in range(number_occupied)]
            for old, new in zip(reduce_from, reduce_to):
                if old not in summation:
                    continue
                summation.remove(old)
                summation.append(new)
            reduced_permutation_list.append({
                "symmetry_operator": expression["symmetry_operator"],
                "bra": expression["bra"].update_indices(reduce_from, reduce_to),
                "ket": expression["ket"],
                "factor": expression["factor"],
                "summation": summation,
                "t": expression["t"].update_indices(reduce_from, reduce_to),
                "integrals": expression["integrals"].update_indices(reduce_from, reduce_to),
                "E": expression["E"],
            })
        reduced_permuted_expressions.append(reduced_permutation_list)
    return reduced_permuted_expressions

def permutation_check(reduced_permuted_expressions: list[list[dict]]) -> list[dict]:
    """Perform permutation check on the expressions resulting from each individual permutation operator.

    args:
        reduced_permuted_expressions: A list containg groups of permuted expressions.

    return:
        A list of all permutationally unique terms:
            symmetry_operator: The symmetrization operator.
            factor: The prefactor for the expression.
            summation: The indices that are summed over.
            bra: The bra.
            ket: The ket.
            integrals: The integrals resulting from using Box 13.2.
            t: The amplitudes.
            E: The original excitation operator.
    """
    # Define a function that check permutations
    def check(permutation_list: list[dict], pg: int, max_pg: int) -> list[dict]:
        check_list = []
        for e, element in enumerate(progressbar(permutation_list[::-1], f"{f'Performing permutation check on group {pg} of {max_pg}: ':<50}"), 1):
            summation = element["summation"]
            permutable_virtual = [i for i in summation if i[0] == "v"]
            permutable_occupied = [i for i in summation if i[0] == "o"]
            found_match = False
            for old, new in element["symmetry_operator"]:
                if found_match:
                    break
                permuted_element = deepcopy(element)
                permuted_element["bra"] = permuted_element["bra"].update_indices(old, new)
                permuted_element["t"] = permuted_element["t"].update_indices(old, new)
                permuted_element["integrals"] = permuted_element["integrals"].update_indices(old, new)
                for comparison in permutation_list[:-e]:
                    if type(element["integrals"]) != type(comparison["integrals"]):
                        continue
                    if found_match:
                        break
                    for virtual_permutation in permutations(permutable_virtual):
                        if found_match:
                            break
                        for occupied_permutation in permutations(permutable_occupied):
                            old_indices = permutable_virtual + permutable_occupied
                            new_indices = list(virtual_permutation) + list(occupied_permutation)
                            if permuted_element["bra"].update_indices(old_indices, new_indices) != comparison["bra"]:
                                continue
                            if permuted_element["t"].update_indices(old_indices, new_indices) != comparison["t"]:
                                continue
                            if permuted_element["integrals"].update_indices(old_indices, new_indices) != comparison["integrals"]:
                                continue
                            found_match = True
                            comparison["factor"] += element["factor"]
                            break
            if not found_match:
                check_list.append(element)
        return check_list
    # Do a check over permutations from each permutation operator seperately
    perm_check = []
    for pg, permutation_group in enumerate(reduced_permuted_expressions):
        if len(permutation_group) == 1:
            print(f"Only one element in group {pg}, skipping permutation check")
            perm_check.append(permutation_group[0])
            continue
        perm_check += check(permutation_group, pg, len(reduced_permuted_expressions))
    print("")
    return perm_check

def perform_explicit_symmetrization(terms: list[dict]) -> list[dict]:
    """Perform explicit symmetrization.

    args:
        List of terms with arbitrary indices.

    return:
        List of terms with arbitrary indices.
    """
    symmetrized_elements = []
    for term in terms:
        symmetrization_operator: P = term["symmetry_operator"]
        assert isinstance(symmetrization_operator, P)
        for old_indices, new_indices in symmetrization_operator:
            symmetrized_elements.append({
                "symmetry_operator": P([]),
                "factor": term["factor"],
                "summation": term["summation"],
                "bra": term["bra"].update_indices(old_indices, new_indices),
                "t": term["t"].update_indices(old_indices, new_indices),
                "integrals": term["integrals"].update_indices(old_indices, new_indices),
                "ket": term["ket"],
                "E": term["E"].update_indices(old_indices, new_indices),
            })
    return symmetrized_elements


def translate_to_normal_indices(terms: list[dict]) -> list[dict]:
    """Translate the v0... and o0... indices to a, i, b, j, and so on.

    args:
        List of terms with arbitrary indices.

    return:
        List of terms with normal indices.
    """
    translation_dict = {
        "v0": "a", "v1": "b", "v2": "c", "v3": "d", "v4": "e", "v5": "f", "v6": "g", "v7": "h", "v8": "q", "v9": "r", "v10": "t", "v11": "w", "v12": "y",
        "o0": "i", "o1": "j", "o2": "k", "o3": "l", "o4": "m", "o5": "n", "o6": "o", "o7": "p", "o8": "s", "o9": "u", "o10": "v", "o11": "x", "o12": "z",
        "v13": "A", "v14": "B", "v15": "C", "v16": "D", "v17": "E", "v18": "F", "v19": "G", "v20": "H", "v21": "Q", "v22": "R", "v23": "T", "v24": "W", "v25": "Y",
        "o13": "I", "o14": "J", "o15": "K", "o16": "L", "o17": "M", "o18": "N", "o19": "O", "o20": "P", "o21": "S", "o22": "U", "o23": "V", "o24": "X", "o25": "Z",
    }
    translated_terms = []
    for term in terms:
        used_indices = term["t"].indices + term.get("integrals", E([])).indices + term["E"].indices + term["bra"].indices
        new_indices = [translation_dict[i] for i in used_indices]
        summation = [translation_dict[i] for i in term["summation"]]
        translated_terms.append({
            "symmetry_operator": term["symmetry_operator"].update_indices(used_indices, new_indices, translate=True),
            "factor": term["factor"],
            "summation": summation,
            "bra": term["bra"].update_indices(used_indices, new_indices, translate=True),
            "t": term["t"].update_indices(used_indices, new_indices, translate=True),
            "integrals": term.get("integrals", E([])).update_indices(used_indices, new_indices, translate=True),
            "ket": term["ket"],
            "E": term["E"].update_indices(used_indices, new_indices, translate=True),
        })
    return translated_terms

def print_expression(terms: list[dict], one_electron_type: str) -> None:
    """Print all terms in latex format.

    args:
        List of terms.
        The one electron integral string representation.
        The permutation operator that symmetrizes the E indices.
    """
    for term in terms:
        summation = "".join(term["summation"])
        left_vector = ''
        if term["bra"].left_excitation_vector:
            left_vector += f'\\bar{{t}}_{{{"".join(term["bra"].indices)}}} '
        expression = f'{term["factor"]} \sum_{{{summation}}} {term["symmetry_operator"]} {str(term["integrals"]).replace("F", one_electron_type)}{left_vector}{term["t"]}'
        print(expression)
    print("")

def print_python_code(terms: list[dict], one_electron_type: str) -> None:
    """Print all terms as Python code.

    args:
        List of terms.
        The one electron integral string representation.
        The permutation operator that symmetrizes the E indices.
    """
    if terms == []:
        return
    for t, term in enumerate(terms):
        # Start with the indices
        if term["bra"].is_HF:
            left_side = "".join(term["bra"].indices)
            right_side = "".join(term["E"].indices)
        if term["bra"].left_excitation_vector:
            left_side = "".join(term["bra"].indices) + ","
            right_side = "".join(term["E"].indices)
        else:
            left_side = ""
            right_side = "".join(term["bra"].indices)
        for amp in term["t"].amplitudes:
            left_side += "".join(amp.indices) + ","
        left_side += "".join(term["integrals"].indices)
        # Then each component
        component = ""
        if term["bra"].left_excitation_vector:
            component += f"l{len(term['bra'])},"
        t_counter = dict()
        for amp in term["t"].amplitudes:
            nr = t_counter.get(len(amp), 1)
            component += f"t{len(amp)}_{nr},"
            if nr == 1:
                t_counter[len(amp)] = 2
            else:
                t_counter[len(amp)] += 1
        component += f"{term['integrals'].type.replace('F', one_electron_type)}_{''.join(term['integrals'].dims)}"
        sign = "+"
        if term['factor'] < 0:
            sign = "-"
        if t == 0:
            print(f"output = {term['factor']} * mem.contract(\"{left_side}->{right_side}\",{component})")
        else:
            print(f"output {sign}= {abs(term['factor'])} * mem.contract(\"{left_side}->{right_side}\",{component})")
    print("")
    if len(terms[0]["symmetry_operator"]) != 1:
        print("return symmetrize_index(output)")
    else:
        print("return output")

def main():
    """Main function."""
    arguments = parse_arguments()
    bra = arguments.bra[0]
    commutator = arguments.commutator[0]
    t1_transformed = arguments.t1_transformed
    restricted = arguments.restricted
    perm_check = arguments.perm_check
    explicit_sym = arguments.explicit_sym
    one_electron = False
    two_electron = False
    one_electron_type = "F"
    if "F" in commutator:
        one_electron = True
        operator = "F"
    if "X" in commutator:
        one_electron = True
        one_electron_type = "X"
        operator = "X"
    if "P" in commutator:
        two_electron = True
        operator = "P"
    if "H" in commutator:
        one_electron = True
        two_electron = True
        operator = "H"
    matrix_element, virtual_index_counter, occupied_index_counter = commutator_indexing(bra, commutator, "HF")
    print_matrix_element(matrix_element, operator)
    expanded_matrix_element = commutator_expansion(matrix_element)
    collected_matrix_elements = collect_commutator_terms(expanded_matrix_element)
    mathematical_expressions = commutator_box(collected_matrix_elements, virtual_index_counter, occupied_index_counter, restricted, t1_transformed, one_electron, two_electron)
    permuted_expressions = perform_permutations(mathematical_expressions)
    reduced_permuted_expressions = match_reduce_indices(permuted_expressions)
    if perm_check:
        permutation_checked = permutation_check(reduced_permuted_expressions)
    else:
        permutation_checked = sum(reduced_permuted_expressions, [])
    if explicit_sym:
        permutation_checked = perform_explicit_symmetrization(permutation_checked)
    normal_indices = translate_to_normal_indices(permutation_checked)
    print_expression(normal_indices, one_electron_type)
    print_python_code(normal_indices, one_electron_type)

if __name__ == "__main__":
    main()