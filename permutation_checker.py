from itertools import permutations
from typing import List, Tuple, Union, Dict
from copy import deepcopy
import numpy as np

def g_checker(g1: list, g2: list) -> bool:
    try:
        assert len(g1) == 4 and len(g2) == 4, f"The length of g1 and g2 has to be the same and 4 - here they were {len(g1)} and {len(g2)}"
    except AssertionError as err:
        print(err)
        exit()
    
    if g1 == g2:
        return True
    
    p,q,r,s = g1

    if [q,p,r,s] == g2:
        return True
    elif [p,q,s,r] == g2:
        return True
    elif [q,p,s,r] == g2:
        return True
    elif [r,s,p,q] == g2:
        return True
    elif [s,r,p,q] == g2:
        return True
    elif [r,s,q,p] == g2:
        return True
    elif [s,r,q,p] == g2:
        return True
    else:
        return False

def L_checker(L1: list, L2: list) -> bool:
    try:
        assert len(L1) == 4 and len(L2) == 4, f"The length of L1 and L2 has to be the same and 4 - here they were {len(L1)} and {len(L2)}"
    except AssertionError as err:
        print(err)
        exit()
    
    if L1 == L2:
        return True
    
    p1,q1,r1,s1 = L1
    p2,q2,r2,s2 = L2
    
    return g_checker([p1,q1,r1,s1],[p2,q2,r2,s2]) and g_checker([p1,s1,r1,q1],[p2,s2,r2,q2])

def t_checker(t1: list, t2: list) -> bool:
    try:
        assert len(t1) == 4 and len(t2) == 4, f"The length of t1 and t2 has to be the same and 4 - here they were {len(t1)} and {len(t2)}"
    except AssertionError as err:
        print(err)
        exit()

    if t1 == t2:
        return True

    p,q,r,s = t1

    if [r,s,p,q] == t2:
        return True
    else:
        return False

def RV_checker(RV1: list, RV2: list) -> bool:
    try:
        assert (len(RV1) == 4 and len(RV2) == 4) or (len(RV1) == 2 and len(RV2) == 2), f"The length of RV1 and RV2 has to be the same and either 2 or 4 - here they were {len(RV1)} and {len(RV2)}"
    except AssertionError as err:
        print(err)
        exit()

    if len(RV1) == 2:
        return F_checker(RV1, RV2)
    elif len(RV2) == 4:
        return t_checker(RV1, RV2)

def LV_checker(LV1: list, LV2: list) -> bool:
    try:
        assert (len(LV1) == 4 and len(LV2) == 4) or (len(LV1) == 2 and len(LV2) == 2), f"The length of LV1 and LV2 has to be the same and either 2 or 4 - here they were {len(LV1)} and {len(LV2)}"
    except AssertionError as err:
        print(err)
        exit()

    if len(LV1) == 2:
        return F_checker(LV1, LV2)
    elif len(LV2) == 4:
        return t_checker(LV1, LV2)

def F_checker(F1: list, F2: list) -> bool:
    assert len(F1) == 2 and len(F2) == 2, f"The length of F1 and F2 has to be the same and 2 - here they were {len(F1)} and {len(F2)}"

    return F1 == F2

def permutationChecker(*, P: List[str], F: str = None, L: str = None, g: str = None, RV: str = None, LV: str = None, t: Union[List[str],str] = None, bra: Union[List[str],str] = 'ai', E: Union[List[str],str] = None) -> List[Dict[str, List[str]]]:
    """Generates permutation and automatically elimineates indicies afterwards

    Args:
        P (List[str]): Contains the indicies that are permuted
        F (str, optional): The indicies of the F component. Defaults to None.
        L (str, optional): The indicies of the L component. Defaults to None.
        g (str, optional): The indicies of the g component. Defaults to None.
        RV (str, optional): The indicies of the right excitaion vector. Defaults to None.
        LV (str, optional): The indicies of the left excitaion vector. Defaults to None.
        t (Union[List[str],str], optional): The indicies of the t component - if there are multiple t components they should be written in a list. Defaults to None.
        bra (Union[List[str],str], optional): The indicies in the bra - if the bra is more than singly excited they should be written in a list. Defaults to 'ai'.
        E (Union[List[str],str], optional): The indicies of the E operators - if there are multiple E operators they should be written in a list. Defaults to None.

    Returns:
        List[Dict[str, List[str]]]: Contains a list of all permutations with the keys of the dictionary being F, L, g, and R depending on what has been input
    """
    # Checking that inputs are correct
    try:
        assert isinstance(P, list), f"P must be of type string; here it was {type(P)}"
        assert isinstance(F, str) or F == None, f"F must be of type string; here it was {type(F)}"
        assert isinstance(L, str) or L == None, f"L must be of type string; here it was {type(L)}"
        assert isinstance(g, str) or g == None, f"g must be of type string; here it was {type(g)}"
        assert isinstance(RV, str) or RV == None, f"RV must be of type string; here it was {type(RV)}"
        assert isinstance(LV, str) or LV == None, f"RV must be of type string; here it was {type(LV)}"
        assert isinstance(t, list) or isinstance(t, str) or t == None, f"t must be of type string or list; here it was {type(t)}"
        assert isinstance(E, list) or isinstance(E, str) or E == None, f"E must be of type string or list; here it was {type(E)}"
        assert isinstance(bra, list) or isinstance(bra, str), f"bra must be of type string or list; here it was {type(bra)}"
    except AssertionError as err:
        print(err)
        exit()

    # Changing t, E, and bra to lists if they arent already
    if isinstance(E, str):
        E = [E]
    if isinstance(t, str):
        t = [t]
    if isinstance(bra, str):
        bra = [bra]

    # Indicies recognised by function
    vir = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
    occ = ('i', 'j', 'k', 'l', 'm', 'n', 'o', 'p')

    # Collecting all indicies used which are not reserved - (a, i, b, j) being reserved
    vir_occ_used = set()
    reserved = {'a', 'i', 'b', 'j'}
    # reserved = set()
    
    # Finds vir and occ indicies to permutate
    P_vir = [idx for script in P for idx in script if idx in vir]
    P_occ = [idx for script in P for idx in script if idx in occ]

    # Permutating vir and occ indicies
    perm1 = permutations(P_vir)
    perm2 = permutations(P_occ)

    # Storing index of indicies for use when permutating
    if F:
        F_indexing = indexing(F, vir, occ, P_vir, P_occ)
        vir_occ_used |= set(i for i in F)
    
    if L:
        L_indexing = indexing(L, vir, occ, P_vir, P_occ)
        vir_occ_used |= set(i for i in L)

    if g:
        g_indexing = indexing(g, vir, occ, P_vir, P_occ)
        vir_occ_used |= set(i for i in g)

    if RV:
        vir_occ_used |= set(i for i in RV)

    if LV:
        vir_occ_used |= set(i for i in LV)

    if t:
        vir_occ_used |= set(i for j in t for i in j)

    bra_vir = []
    bra_occ = []
    if bra:
        for excitation in bra:
            for idx in excitation:
                if idx in vir:
                    bra_vir += idx
                elif idx in occ:
                    bra_occ += idx
                else:
                    print("Something went wrong")
                    print(f"Index {idx} was not found in {occ=} or {vir=}")
                    exit()
        vir_occ_used |= set(i for j in bra for i in j)

    if E:
        E_indexing = []
        reserved_vir_in_E = 0
        reserved_occ_in_E = 0
        for operator in E:
            operator_indexing = []
            for idx in operator:
                try:
                    index = P_vir.index(idx)
                    operator_indexing += [f'v{index}']
                    continue
                except ValueError: ...
                try:
                    index = P_occ.index(idx)
                    operator_indexing += [f'o{index}']
                    continue
                except ValueError: ...
                try:
                    index = vir.index(idx)
                    operator_indexing += [f'x{idx}']
                    if vir[index] in bra_vir:
                        reserved_vir_in_E += 1  
                    continue
                except ValueError: ...
                try:
                    index = occ.index(idx)
                    operator_indexing += [f'x{idx}']
                    if occ[index] in bra_occ:
                        reserved_occ_in_E += 1
                    continue
                except ValueError: ...
                print("Something went wrong")
                print(f"Index {idx} was not found in {P_occ=}, {P_vir=}, {occ=} or {vir=}")
                exit()
            E_indexing += [operator_indexing]

        len_E = len(E)
        vir_occ_used |= set(i for j in E for i in j)
    
        vir_avail = sorted([i for i in vir_occ_used if i in vir])[:-len_E + reserved_vir_in_E]
        occ_avail = sorted([i for i in vir_occ_used if i in occ])[:-len_E + reserved_occ_in_E]

    vir_avail = sorted(list(set(vir_avail) - reserved))
    occ_avail = sorted(list(set(occ_avail) - reserved))
    
    permutation_collector = []

    # Looping over permutations
    for vir_perm, occ_perm in zip(perm1,perm2):
        # R and t are not permuted, so they are copied
        RV_perm = []
        if RV:
            RV_perm = [idx for idx in RV]
        
        LV_perm = []
        if LV:
            LV_perm = [idx for idx in LV]
        
        t_perm = []
        if t:
            for correction in t:
                t_perm += [[idx for idx in correction]]

        # F, L, g, and E are permuted, so the current permutation is done here
        F_perm = []
        if F:
            F_perm = permute_indicies(F_indexing, vir_perm, occ_perm)

        L_perm = []
        if L:
            L_perm = permute_indicies(L_indexing, vir_perm, occ_perm)
        
        g_perm = []
        if g:
            g_perm = permute_indicies(g_indexing, vir_perm, occ_perm)
        
        E_perm = []
        if E:
            for operator in E_indexing:
                E_perm += [permute_indicies(operator, vir_perm, occ_perm)]

        for n, (v, o) in enumerate(E_perm):
            if F:
                F_perm = eliminate_E_component(reserved, bra_vir, bra_occ, F_perm, n, v, o)
            if L:
                L_perm = eliminate_E_component(reserved, bra_vir, bra_occ, L_perm, n, v, o)
            if g:
                g_perm = eliminate_E_component(reserved, bra_vir, bra_occ, g_perm, n, v, o)
            if RV:
                RV_perm = eliminate_E_component(reserved, bra_vir, bra_occ, RV_perm, n, v, o)
            if LV:
                LV_perm = eliminate_E_component(reserved, bra_vir, bra_occ, LV_perm, n, v, o)
            if t:
                for nn, correction in enumerate(t_perm):
                    t_perm[nn] = eliminate_E_component(reserved, bra_vir, bra_occ, correction, n, v, o)

        indicies_used = set(F_perm) | set(L_perm) | set(g_perm) | set(RV_perm) | set(LV_perm) | set(i for j in t_perm for i in j)

        vir_used = sorted(list(indicies_used - reserved & set(vir)))
        vir_unused = sorted([i for i in set(vir) - reserved if i < max(vir_used) and i not in vir_used])
        vir_used.reverse()

        for unused, used in zip(vir_unused, vir_used[:len(vir_unused)]):
            F_perm = renameIndicies(used, unused, F_perm)
            L_perm = renameIndicies(used, unused, L_perm)
            g_perm = renameIndicies(used, unused, g_perm)
            RV_perm = renameIndicies(used, unused, RV_perm)
            LV_perm = renameIndicies(used, unused, LV_perm)
            for i, correction in enumerate(t_perm):
                t_perm[i] = renameIndicies(used, unused, correction)

        occ_used = sorted(list(indicies_used - reserved & set(occ)))
        occ_unused = sorted([i for i in set(occ) - reserved if i < max(occ_used) and i not in occ_used])
        occ_used.reverse()

        for unused, used in zip(occ_unused, occ_used[:len(occ_unused)]):
            F_perm = renameIndicies(used, unused, F_perm)
            L_perm = renameIndicies(used, unused, L_perm)
            g_perm = renameIndicies(used, unused, g_perm)
            RV_perm = renameIndicies(used, unused, RV_perm)
            LV_perm = renameIndicies(used, unused, LV_perm)
            for i, correction in enumerate(t_perm):
                t_perm[i] = renameIndicies(used, unused, correction)

        current_perm = dict()
        if LV:
            current_perm['LV'] = LV_perm
        if F:
            current_perm['F'] = F_perm
        if L:
            current_perm['L'] = L_perm
        if g:
            current_perm['g'] = g_perm
        if t:
            current_perm['t'] = t_perm
        if RV:
            current_perm['RV'] = RV_perm

        permutation_collector.append(current_perm)
    return permutation_collector

def eliminate_E_component(reserved: set, bra_vir: List[str], bra_occ: List[str], permutaton_indicies: List[str], n: int, v: str, o: str) -> List[str]:
    perm = deepcopy(permutaton_indicies)
    try:
        if v in reserved:
            iv = perm.index(bra_vir[n])
            perm[iv] = v
        else:
            iv = perm.index(v)
            perm[iv] = bra_vir[n]
    except ValueError: ...
    try:
        if o in reserved:
            iv = perm.index(bra_occ[n])
            perm[iv] = o
        else:
            io = perm.index(o)
            perm[io] = bra_occ[n]
    except ValueError: ...

    return perm

def indexing(component_indicies: str, vir: Tuple[str], occ: Tuple[str], P_vir: List[str], P_occ: List[str]) -> List[str]:
    indexing = []
    for idx in component_indicies:
        try:
            index = P_vir.index(idx)
            indexing += [f'v{index}']
            continue
        except ValueError: ...
        try:
            index = P_occ.index(idx)
            indexing += [f'o{index}']
            continue
        except ValueError: ...
        try:
            vir.index(idx)
            indexing += [f'x{idx}']
            continue
        except ValueError: ...
        try:
            occ.index(idx)
            indexing += [f'x{idx}']
            continue
        except ValueError: ...
        print("Something went wrong")
        print(f"Index {idx} was not found in {P_occ=}, {P_vir=}, {occ=} or {vir=}")
        exit()
    return indexing

def permute_indicies(indexing: List[str], vir_perm: Tuple[str], occ_perm: Tuple[str]) -> List[str]:
    permutation = []
    for idx in indexing:
        if idx[0] == 'o':
            permutation += occ_perm[int(idx[1])]
        elif idx[0] == 'v':
            permutation += vir_perm[int(idx[1])]
        elif idx[0] == 'x':
            permutation += idx[1]
        else:
            print("Something went wrong in permuted_indicies")
            print(f"{idx[0]=} did not match o, v, or x")
            exit()

    return permutation

def renameIndicies(letter_to_replace: str, letter_to_replace_with: str, lst: List[str]) -> List[str]:
    return [i if i != letter_to_replace else letter_to_replace_with for i in lst]

def permutationComparison(permutations: List[Dict[str, List[Union[List[str],str]]]], summation: str) -> np.array:
    try:
        assert isinstance(permutations, list), f"permutations must be of type string or list; here it was {type(permutations)}"
        assert isinstance(summation, str), f"summation must be of type string; here it was {type(summation)}"
    except AssertionError as err:
        print(err)
        exit()

    permutation_keys = permutations[0].keys()
    F = 'F' in permutation_keys
    L = 'L' in permutation_keys
    g = 'g' in permutation_keys
    t = 't' in permutation_keys
    RV = 'RV' in permutation_keys
    LV = 'LV' in permutation_keys

    vir = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
    occ = ('i', 'j', 'k', 'l', 'm', 'n', 'o', 'p')

    summation = [i for i in summation]

    vir_idx = set(vir) & set(summation)
    occ_idx = set(occ) & set(summation)

    perm_copy = deepcopy(permutations)
    for i, perm in enumerate(permutations):
        for key, lst in perm.items():
            if key == 't':
                for n, correction in enumerate(lst):
                    for nn, idx in enumerate(correction):
                        if idx in vir_idx:
                            perm_copy[i][key][n][nn] = 'vir'
                        elif idx in occ_idx:
                            perm_copy[i][key][n][nn] = 'occ'
            else:
                for n, idx in enumerate(lst):
                    if idx in vir_idx:
                        perm_copy[i][key][n] = 'vir'
                    elif idx in occ_idx:
                        perm_copy[i][key][n] = 'occ'

    permutations_compared = np.zeros((len(permutations),len(permutations)))

    for n, perm1 in enumerate(perm_copy):
        for nn, perm2 in enumerate(perm_copy[n:]):
            same = 1
            if F and same:
                same *= F_checker(perm1['F'], perm2['F'])
            if L and same:
                same *= L_checker(perm1['L'], perm2['L'])
            if g and same:
                same *= g_checker(perm1['g'], perm2['g'])
            if t and same:
                for i, _ in enumerate(perm1['t']):
                    same *= t_checker(perm1['t'][i], perm2['t'][i])
            if RV and same:
                same *= RV_checker(perm1['RV'], perm2['RV'])
            if LV and same:
                same *= LV_checker(perm1['LV'], perm2['LV'])
            permutations_compared[n,n+nn] = same

    idx_same = set()
    remaining_permutations = []
    for i, line in enumerate(permutations_compared):
        if i in idx_same:
            continue
        idx_same |= set(np.where(line == 1)[0])
        remaining_permutations += [permutations[i]]

    return remaining_permutations


if __name__ == '__main__':
    perms = permutationChecker(P=['cde','klm'], g='knle', E=['cn','dm'], RV='em', bra=['ai','bj'], t='ckdl')
    
    print("All permutations")
    for i in perms:
        print(i)

    print("\nUnique permutations")

    perms_compared = permutationComparison(perms,'cdeklm')
    for i in perms_compared:
        print(i)
