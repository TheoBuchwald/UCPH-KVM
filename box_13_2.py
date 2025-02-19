from collections.abc import Callable
from operators import F, g, L, E, P
from copy import deepcopy

def n0o0(one_electron: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(_matrix_element: dict, _virtual_index_counter: int, occupied_index_counter):
        output = []
        if one_electron:
            print("This term is very strange. If this is correct, it would be best to do this by hand.")
            quit()
        return output
    return translater

def n0o1(t1_transformed: bool, one_electron: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, virtual_index_counter: int, occupied_index_counter):
        output = []
        a = f"v{virtual_index_counter}"
        i = f"o{occupied_index_counter}"
        new_matrix_element = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"]
        }
        if one_electron:
            new_matrix_element["factor"] = matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [a, i]
            new_matrix_element["integrals"] = F([a, i], t1_transformed)
            new_matrix_element["box_E"] = E([a, i])
            output.append(new_matrix_element)
        return output
    return translater

def n0o2(t1_transformed: bool, two_electron: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, virtual_index_counter: int, occupied_index_counter):
        output = []
        a = f"v{virtual_index_counter}"
        i = f"o{occupied_index_counter}"
        b = f"v{virtual_index_counter + 1}"
        j = f"o{occupied_index_counter + 1}"
        new_matrix_element = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"]
        }
        if two_electron:
            new_matrix_element["factor"] = 0.5 * matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [a, i, b, j]
            new_matrix_element["integrals"] = g([a, i, b, j], t1_transformed)
            new_matrix_element["box_E"] = E([a, i, b, j])
            output.append(new_matrix_element)
        return output
    return translater

def n1o0(t1_transformed: bool, one_electron: bool, restricted: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, _virtual_index_counter: int, _occupied_index_counter: int) -> list[dict]:
        output = []
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
        if one_electron:
            if restricted:
                new_matrix_element["factor"] = 2 * matrix_element["factor"]
            else:
                new_matrix_element["factor"] = matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"]
            new_matrix_element["integrals"] = F([i,a], t1_transformed)
            new_matrix_element["box_E"] = E([])
            output.append(new_matrix_element)
        return output
    return translater

def n1o1(t1_transformed: bool, one_electron: bool, two_electron) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, virtual_index_counter: int, occupied_index_counter: int) -> list[dict]:
        output = []
        commutator = matrix_element["commutator"]
        a = commutator[0][0]
        i = commutator[0][1]
        b = f"v{virtual_index_counter}"
        j = f"o{occupied_index_counter}"
        new_matrix_element = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"]
        }
        if one_electron:
            new_matrix_element["factor"] = matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [b]
            new_matrix_element["integrals"] = F([b,a], t1_transformed)
            new_matrix_element["box_E"] = E([b,i])
            output.append(new_matrix_element)
            new_matrix_element = deepcopy(new_matrix_element)
            new_matrix_element["factor"] = -matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [j]
            new_matrix_element["integrals"] = F([i,j], t1_transformed)
            new_matrix_element["box_E"] = E([a,j])
            output.append(new_matrix_element)
        if two_electron:
            new_matrix_element = deepcopy(new_matrix_element)
            new_matrix_element["factor"] = matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [b, j]
            new_matrix_element["integrals"] = L([b,j,i,a], t1_transformed)
            new_matrix_element["box_E"] = E([b,j])
            output.append(new_matrix_element)
        return output
    return translater

def n1o2(t1_transformed: bool, two_electron: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, virtual_index_counter: int, occupied_index_counter: int) -> list[dict]:
        output = []
        commutator = matrix_element["commutator"]
        a = commutator[0][0]
        i = commutator[0][1]
        b = f"v{virtual_index_counter}"
        j = f"o{occupied_index_counter}"
        c = f"v{virtual_index_counter + 1}"
        k = f"o{occupied_index_counter + 1}"
        new_matrix_element = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"]
        }
        if two_electron:
            new_matrix_element["factor"] = matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [b,j,c]
            new_matrix_element["integrals"] = g([b,j,c,a], t1_transformed)
            new_matrix_element["box_E"] = E([b,j,c,i])
            output.append(new_matrix_element)
            new_matrix_element = deepcopy(new_matrix_element)
            new_matrix_element["factor"] = -matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [b,j,k]
            new_matrix_element["integrals"] = g([b,j,i,k], t1_transformed)
            new_matrix_element["box_E"] = E([b,j,a,k])
            output.append(new_matrix_element)
        return output
    return translater

def n2o0(t1_transformed: bool, two_electron: bool, restricted: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, _virtual_index_counter: int, _occupied_index_counter: int) -> list[dict]:
        output = []
        commutator = matrix_element["commutator"]
        a = commutator[0][0]
        i = commutator[0][1]
        b = commutator[1][0]
        j = commutator[1][1]
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
        if two_electron:
            new_matrix_element["summation"] = matrix_element["summation"]
            new_matrix_element["integrals"] = L([i,a,j,b], t1_transformed)
            new_matrix_element["box_E"] = E([])
            output.append(new_matrix_element)
        return output
    return translater

def n2o1(t1_transformed: bool, one_electron: bool, two_electron: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, virtual_index_counter: int, occupied_index_counter: int) -> list[dict]:
        output = []
        commutator = matrix_element["commutator"]
        a = commutator[0][0]
        i = commutator[0][1]
        b = commutator[1][0]
        j = commutator[1][1]
        c = f"v{virtual_index_counter}"
        k = f"o{occupied_index_counter}"
        new_matrix_element = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"]
        }
        if one_electron:
            new_matrix_element["factor"] = -matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"]
            new_matrix_element["permutation"] = P([a,i,b,j])
            new_matrix_element["integrals"] = F([i,b], t1_transformed)
            new_matrix_element["box_E"] = E([a,j])
            output.append(new_matrix_element)
        if two_electron:
            new_matrix_element = deepcopy(new_matrix_element)
            new_matrix_element["factor"] = -matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [k]
            new_matrix_element["permutation"] = P([a,i,b,j])
            new_matrix_element["integrals"] = L([i,k,j,b], t1_transformed)
            new_matrix_element["box_E"] = E([a,k])
            output.append(new_matrix_element)
            new_matrix_element = deepcopy(new_matrix_element)
            new_matrix_element["factor"] = matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [c]
            new_matrix_element["permutation"] = P([a,i,b,j])
            new_matrix_element["integrals"] = L([c,a,j,b], t1_transformed)
            new_matrix_element["box_E"] = E([c,i])
            output.append(new_matrix_element)
        return output
    return translater

def n2o2(t1_transformed: bool, two_electron: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, virtual_index_counter: int, occupied_index_counter: int) -> list[dict]:
        output = []
        commutator = matrix_element["commutator"]
        a = commutator[0][0]
        i = commutator[0][1]
        b = commutator[1][0]
        j = commutator[1][1]
        c = f"v{virtual_index_counter}"
        k = f"o{occupied_index_counter}"
        d = f"v{virtual_index_counter + 1}"
        l = f"o{occupied_index_counter + 1}"
        new_matrix_element = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"]
        }
        if two_electron:
            new_matrix_element["factor"] = -matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [c,k]
            new_matrix_element["permutation"] = P([a,i,b,j])
            new_matrix_element["integrals"] = g([i,b,c,k], t1_transformed)
            new_matrix_element["box_E"] = E([a,j,c,k])
            output.append(new_matrix_element)
            new_matrix_element = deepcopy(new_matrix_element)
            new_matrix_element["factor"] = -matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [c,k]
            new_matrix_element["permutation"] = P([a,i,b,j])
            new_matrix_element["integrals"] = g([i,k,c,b], t1_transformed)
            new_matrix_element["box_E"] = E([a,k,c,j])
            output.append(new_matrix_element)
            new_matrix_element = deepcopy(new_matrix_element)
            new_matrix_element["factor"] = matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [k,l]
            new_matrix_element["integrals"] = g([i,k,j,l], t1_transformed)
            new_matrix_element["box_E"] = E([a,k,b,l])
            output.append(new_matrix_element)
            new_matrix_element = deepcopy(new_matrix_element)
            new_matrix_element["factor"] = matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [c,d]
            new_matrix_element["integrals"] = g([c,a,d,b], t1_transformed)
            new_matrix_element["box_E"] = E([c,i,d,j])
            output.append(new_matrix_element)
        return output
    return translater

def n3o1(t1_transformed: bool, two_electron: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, _virtual_index_counter: int, _occupied_index_counter: int) -> list[dict]:
        output = []
        commutator = matrix_element["commutator"]
        a = commutator[0][0]
        i = commutator[0][1]
        b = commutator[1][0]
        j = commutator[1][1]
        c = commutator[2][0]
        k = commutator[2][1]
        new_matrix_element = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"]
        }
        if two_electron:
            new_matrix_element["factor"] = -matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"]
            new_matrix_element["permutation"] = P([a,i,b,j,c,k])
            new_matrix_element["integrals"] = L([j,b,i,c], t1_transformed)
            new_matrix_element["box_E"] = E([a,k])
            output.append(new_matrix_element)
        return output
    return translater

def n3o2(t1_transformed: bool, two_electron: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, virtual_index_counter: int, occupied_index_counter: int) -> list[dict]:
        output = []
        commutator = matrix_element["commutator"]
        a = commutator[0][0]
        i = commutator[0][1]
        b = commutator[1][0]
        j = commutator[1][1]
        c = commutator[2][0]
        k = commutator[2][1]
        d = f"v{virtual_index_counter}"
        l = f"o{occupied_index_counter}"
        new_matrix_element = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"]
        }
        if two_electron:
            new_matrix_element["factor"] = matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [l]
            new_matrix_element["permutation"] = P([a,i,b,j,c,k])
            new_matrix_element["integrals"] = g([i,l,j,c], t1_transformed)
            new_matrix_element["box_E"] = E([a,l,b,k])
            output.append(new_matrix_element)
            new_matrix_element = deepcopy(new_matrix_element)
            new_matrix_element["factor"] = -matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"] + [d]
            new_matrix_element["permutation"] = P([a,i,b,j,c,k])
            new_matrix_element["integrals"] = g([i,b,d,c], t1_transformed)
            new_matrix_element["box_E"] = E([a,j,d,k])
            output.append(new_matrix_element)
        return output
    return translater

def n4o2(t1_transformed: bool, two_electron: bool) -> Callable[[dict, int, int], list[dict]]:
    def translater(matrix_element: dict, _virtual_index_counter: int, _occupied_index_counter: int) -> list[dict]:
        output = []
        commutator = matrix_element["commutator"]
        a = commutator[0][0]
        i = commutator[0][1]
        b = commutator[1][0]
        j = commutator[1][1]
        c = commutator[2][0]
        k = commutator[2][1]
        d = commutator[3][0]
        l = commutator[3][1]
        new_matrix_element = {
            "bra": matrix_element["bra"],
            "t": matrix_element["t"],
            "E": matrix_element["E"],
            "ket": matrix_element["ket"],
            "pre_string": matrix_element["pre_string"],
        }
        if two_electron:
            new_matrix_element["factor"] = 0.5 * matrix_element["factor"]
            new_matrix_element["summation"] = matrix_element["summation"]
            new_matrix_element["permutation"] = P([a,i,b,j,c,k,d,l])
            new_matrix_element["integrals"] = g([i,d,j,c], t1_transformed)
            new_matrix_element["box_E"] = E([a,l,b,k])
            output.append(new_matrix_element)
        return output
    return translater

def get_translater(nesting: int, order: int, t1_transformed: bool, one_electron: bool, two_electron: bool, restricted: bool) -> Callable[[dict, int, int], list[dict]]:
    """Find the correct translater in the restricted case."""
    if nesting == 0:
        if order == 0:
            return n0o0(one_electron)
        if order == 1:
            return n0o1(t1_transformed, one_electron)
        if order == 2:
            return n0o2(t1_transformed, two_electron)
    if nesting == 1:
        # [H,E_ai]|HF>
        if order == 0:
            return n1o0(t1_transformed, one_electron, restricted)
        elif order == 1:
            return n1o1(t1_transformed, one_electron, two_electron)
        elif order == 2:
            return n1o2(t1_transformed, two_electron)
    elif nesting == 2:
        # [[H,E_ai],E_bj]|HF>
        if order == 0:
            return n2o0(t1_transformed, two_electron, restricted)
        elif order == 1:
            return n2o1(t1_transformed, one_electron, two_electron)
        elif order == 2:
            return n2o2(t1_transformed, two_electron)
    elif nesting == 3:
        # [[[H,E_ai],E_bj],E_ck]|HF>
        if order == 1:
            return n3o1(t1_transformed, two_electron)
        if order == 2:
            return n3o2(t1_transformed, two_electron)
    elif nesting == 4:
        # [[[[H,E_ai],E_bj],E_ck],E_dl]|HF>
        if order == 2:
            return n4o2(t1_transformed, two_electron)
    raise ValueError(f"get_translater recieved a nesting ({nesting}) and order ({order}) input that was not recognized.")
