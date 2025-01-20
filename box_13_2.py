from collections.abc import Callable
from operators import F, g, L, E, P
from copy import deepcopy

def n1o0(t1_transformed: bool, restricted: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, _virtual_index_counter: int, _occupied_index_counter: int) -> list[dict]:
        commutator = matrix_element["commutator"]
        a = commutator[0][0]
        i = commutator[0][1]
        new_matrix_element = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"]
        }
        if restricted:
            new_matrix_element["factor"] = 2 * matrix_element["factor"]
        else:
            new_matrix_element["factor"] = matrix_element["factor"]
        new_matrix_element["summation"] = matrix_element["summation"]
        new_matrix_element["integrals"] = F([i,a], t1_transformed)
        new_matrix_element["box_E"] = E([])
        return [new_matrix_element]
    return translater

def n1o1(t1_transformed: bool, one_electron: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, virtual_index_counter: int, occupied_index_counter: int) -> list[dict]:
        commutator = matrix_element["commutator"]
        a = commutator[0][0]
        i = commutator[0][1]
        b = f"v{virtual_index_counter}"
        j = f"o{occupied_index_counter}"
        new_matrix_element_1 = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"]
        }
        new_matrix_element_1["factor"] = matrix_element["factor"]
        new_matrix_element_1["summation"] = matrix_element["summation"] + [b]
        new_matrix_element_1["integrals"] = F([b,a], t1_transformed)
        new_matrix_element_1["box_E"] = E([b,i])
        new_matrix_element_2 = deepcopy(new_matrix_element_1)
        new_matrix_element_2["factor"] = -matrix_element["factor"]
        new_matrix_element_2["summation"] = matrix_element["summation"] + [j]
        new_matrix_element_2["integrals"] = F([i,j], t1_transformed)
        new_matrix_element_2["box_E"] = E([a,j])
        if not one_electron:
            new_matrix_element_3 = deepcopy(new_matrix_element_1)
            new_matrix_element_3["factor"] = matrix_element["factor"]
            new_matrix_element_3["summation"] = matrix_element["summation"] + [b, j]
            new_matrix_element_3["integrals"] = L([b,j,i,a], t1_transformed)
            new_matrix_element_3["box_E"] = E([b,j])
            return [new_matrix_element_1, new_matrix_element_2, new_matrix_element_3]
        return [new_matrix_element_1, new_matrix_element_2]
    return translater

def n1o2(t1_transformed: bool, one_electron: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, virtual_index_counter: int, occupied_index_counter: int) -> list[dict]:
        if one_electron:
            return []
        commutator = matrix_element["commutator"]
        a = commutator[0][0]
        i = commutator[0][1]
        b = f"v{virtual_index_counter}"
        j = f"o{occupied_index_counter}"
        c = f"v{virtual_index_counter + 1}"
        k = f"o{occupied_index_counter + 1}"
        new_matrix_element_1 = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"]
        }
        new_matrix_element_1["factor"] = matrix_element["factor"]
        new_matrix_element_1["summation"] = matrix_element["summation"] + [b,j,c]
        new_matrix_element_1["integrals"] = g([b,j,c,a], t1_transformed)
        new_matrix_element_1["box_E"] = E([b,j,c,i])
        new_matrix_element_2 = deepcopy(new_matrix_element_1)
        new_matrix_element_2["factor"] = -matrix_element["factor"]
        new_matrix_element_2["summation"] = matrix_element["summation"] + [b,j,k]
        new_matrix_element_2["integrals"] = g([b,j,i,k], t1_transformed)
        new_matrix_element_2["box_E"] = E([b,j,a,k])
        return [new_matrix_element_1, new_matrix_element_2]
    return translater

