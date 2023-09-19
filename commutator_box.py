import argparse

def get_occ_and_vir_indices(term : str,bra : list):

    before = [*term.split("[")[0].split(",")[:-1]]

    if len(before) == 0:
        before = []
    elif len(before[0]) == 1 or len(before[0]) == 0:
        before = []

    vir_indices_inside = []
    occ_indices_inside = []

    for com in term.split("[")[-1].split("]"):
        vir_indices_inside += [i[0] for i in com.split(',')[1:]]
        occ_indices_inside += [i[1] for i in com.split(',')[1:]]
   
    occ_indices = occ_indices_inside.copy()
    vir_indices = vir_indices_inside.copy()

    for pair in bra:
    
        vir_indices += [pair[0]]
        occ_indices += [pair[1]]
    
    if len(before) > 0:
        for pair in before:
            vir_indices += [pair[0]]
            occ_indices += [pair[1]]
 
    return before, vir_indices_inside, occ_indices_inside, vir_indices, occ_indices


def get_contribution_from_box(commutator : int, bra_minus_before : int):

    result = None

    # Based on the excitation level of the commutator
    # and the difference between the bra and potential
    # prefixed excitation operators, we find the specific
    # term from Box 13.2

    if commutator == 0:

        # First line in box: H
        if bra_minus_before == 0:
            # sum_i (h_ii + F_ii)
            print("This term is very strange. If this is correct, it would be best to do this by hand.")
            quit()
        elif bra_minus_before == 1:
            # sum_ai F_ai E_ai
            result = " plus -O1 v1o1 -E v1o1"
        elif bra_minus_before == 2:
            # 1/2 sum_aibj g_aibj E_ai E_bj
            result = "plus 1/2 -O2 v1o1v2o2 -E v1o1 v2o2"

    elif commutator == 1:

        # Second line in box [H,ai]
        if bra_minus_before == 0:
            # 2 F_ia
            result = "plus 2 -O1 x1y1"
        elif bra_minus_before == 1:
            # sum_b F_ba E_bi - sum_j F_ij E_aj + sum_bj L_bjia E_bj
            result = """plus -O1 v1y1 -E v1x1
minus -O1 x1o1 -E y1o1
plus  -O3 v1o1x1y1 -E v1o1"""
        elif bra_minus_before == 2:
            # sum_bjc g_bjca E_bj E_ci - sum_bjk g_bjik E_bj E_ak
            result = """plus -O2 v1o1v2y1 -E v1o1 v2x1
minus -O2 v1o1x1o2 -E v1o1 y1o2"""

    elif commutator == 2: 
        
        # Third line in box [[H,ai],bj]
        if bra_minus_before == 0:
            # 2 L_iajb
            result = "plus 2 -O3 x1y1x2y2"
        elif bra_minus_before == 1:
            # - P_ij^ab (F_ib E_aj + sum_k L_ikjb E_ak - L_cajb E_ci)
            result="""minus -P y1y2 x1x2 -O1 x1y2 -E y1x2
minus -P y1y2 x1x2 -O3 x1o1x2y2 -E y1o1
plus -P y1y2 x1x2 -O3 v1y1x2y2 -E v1x1"""
        elif bra_minus_before == 2:
            # -P_ij^ab (sum_ck g_ibck E_aj E_ck + sum_ck g_ikcb E_ak E_cj) + sum_kl g_ikjl E_ak E_bl + sum_cd g_cabd E_ci E_dj
            result="""minus -P y1y2 x1x2 -O2 x1y2v1o1 -E y1x2 v1o1
minus -P y1y2 x1x2 -O2 x1o1v1y2 -E y1o1 v1x2
plus -O2 x1o1x2o2 -E y1o1 y2o2
plus -O2 v1y1v2y2 -E v1x1 v2x2"""
    
    elif commutator == 3:
    
        # Fourth line in box [[[H,ai],bj],ck]
        if bra_minus_before == 1:
            # - P_ijk^abc L_jbic E_ak
            result="minus -P y1y2y3 x1x2x3 -O3 x2y2x1y3 -E y1x3"
        elif bra_minus_before == 2:
            # P_ijk^abc (sum_l g_iljc E_ak E_bk - sum_d g_ibdc E_aj E_dk)
            result="""plus -P y1y2y3 x1x2x3 -O2 x1o1x2y3 -E y1o1 y2x3
minus -P y1y2y3 x1x2x3 -O2 x1y2v1y3 -E y1x2 v1x3"""
    elif commutator == 4:
        # Fifth line in box [[[[H,ai],bj],ck],dl]
        if bra_minus_before == 2:
            # 1/2 P_ijkl^abcd g_idjc E_al E_bk
            result = "plus 1/2 -P y1y2y3y4 x1x2x3x4 -O2 x1y4x2y3 -E y1x4 y2x3"
    
    return result

def check_output(output : str, operator : str, transformed : bool):

    output = output.splitlines()

    terms_to_remove = []

    for i in range(len(output)):

        # For one-electron operators, remove two-electron terms
        if operator in ['F','X']:
            if 'O2' in output[i] or 'O3' in output[i]:
                terms_to_remove.append(i)

        # For canonical H and F, remove the off-diagonal Fock matrices
        if operator in ['H','F'] and not transformed:
            if '-O1 x1y1' in output[i] or '-O1 x1y2' in output[i] or '-O1 v1o1' in output[i]:
                terms_to_remove.append(i)
        
        if operator == 'P':
            if not transformed:
                # For canonical P, remove all Fock terms
                if '-O1' in output[i]:
                    terms_to_remove.append(i)
            elif transformed:
                # For transformed P, remove the oo and vv Fock matrices
                if '-O1 v1y1' in output[i] or '-O1 x1o1' in output[i]:
                    #terms_to_remove.append(i)
                    pass
    
    terms_to_remove = sorted(terms_to_remove,reverse=True)

    for elem in terms_to_remove:
        output.pop(elem)

    if output == []:
        output = None
    else:
        output = '\n'.join(output)

        if operator == 'X':
            singles_op = 'X'
        elif operator in ['H','P','F']:
            singles_op = 'F'
        else:
            print(f"ERROR! The operator {operator} was not recognized! Please use one of H, P, X and F!")
            quit()

        output = output.replace('O1',singles_op)

    return output

def get_box_terms(term : str, bra : list, transform : bool) -> list:

    init_global_variables()
    
    # Determine the operator
    operator = term.split("[")[-1][0]

    # Determine the indices used in the commutator and in total
    before, v_in, o_in, v_total, o_total = get_occ_and_vir_indices(term,bra)

    # Determine the indices to use for summations
    occ_not_used = [i for i in OCC if i not in o_total]
    vir_not_used = [i for i in VIR if i not in v_total]

    # Check if term gives something from Slater-Condon
    # If so, calculate the corresponding term from box 13.2
    if len(bra) < len(before) or len(bra) > len(before) + 2:
        output = None
    else:
        output = get_contribution_from_box(len(v_in),len(bra)-len(before))

    # Determine whether a term should be kept/removed due to
    # the operator or a one-index transformation
    if output is not None:
        output = check_output(output,operator, transform)

    # Have to check again as all the terms may have been removed
    if output is not None:

        # Insert two-electron terms
        # If a one-electron operator is used, the two-electron
        # terms were removed in check_output
        output = output.replace('O2','g')
        output = output.replace('O3','L')
        

        # Insert the summation indices
        if 'o' in output or 'v' in output:
            output = output.replace('o1',occ_not_used[0])
            output = output.replace('v1',vir_not_used[0])
            if 'o2' in output or 'v2' in output:
                output = output.replace('o2',occ_not_used[1])
                output = output.replace('v2',vir_not_used[1])

        # Insert the indices from the parent commutator
        if o_in:
            output = output.replace('x1',o_in[0])
            output = output.replace('y1',v_in[0])
            if len(o_in) > 1:
                output = output.replace('x2',o_in[1])
                output = output.replace('y2',v_in[1])
            if len(o_in) > 2:
                output = output.replace('x3',o_in[2])
                output = output.replace('y3',v_in[2])
            if len(o_in) > 3:
                output = output.replace('x4',o_in[3])
                output = output.replace('y4',v_in[3])
        
        output = output.splitlines()
      
        # Add the prefixed excitation operators, if any
        if before and not any('E' in string for string in output):
            for i in range(len(output)):
               output[i] += ' -E ' + ' '.join(before)
        elif before:
            for i in range(len(output)):
                output[i] += ' '+' '.join(before)

    return output

def init_global_variables() -> None:

    global VIR, OCC

    VIR = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
    OCC = ('i', 'j', 'k', 'l', 'm', 'n', 'o', 'p')

def main() -> None:

    init_global_variables()

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""Takes a commutator of an arbitrary operator with excitation operators to be applied to the HF state and returns the
       conversion from Box 13.2 in "Molecular Electronic Structure Theory". Keeps track of the indices used in the excitation vectors.""", epilog="""For help contact
       Magnus Bukhave Johansen
       qhw298@alumni.ku.dk""")

    parser.add_argument('term', type=str, nargs=1, help='The central part of the matrix element to calculate...Ex. ck,dl,[[H,ai],bj] - Accepted operators are H and P for two-electron operators and F and X for one-electron operators')
    parser.add_argument('-bra', type=str, nargs='+',                                help='Excitations in the bra........................Ex. -bra ai bj. If not included the HF state is assumed.')
    parser.add_argument('-transformed', action='store_true', help='Whether the operator is one-index transformed, e.g T1 or trial-vector', dest='transform')

    # Parses the arguments
    args = parser.parse_args()

    term = args.term[0]
    bra = args.bra
    transform = args.transform
    if bra is None:
        bra = []

    output = get_box_terms(term,bra,transform)

    if output is not None:   
        # Print final result
        for line in output:
            print(line)
    else:
        print("This term does not give anything")

if __name__ == '__main__':
    main()
