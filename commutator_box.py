import argparse

def get_occ_and_vir_indices(term,bra):

    before = [*term.split("[")[0].split(",")]

    if before[0] == '':
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

    for pair in before:
        vir_indices += [pair[0]]
        occ_indices += [pair[1]]

        

    # if len(before) == 1: before = []
    return before, vir_indices_inside, occ_indices_inside, vir_indices, occ_indices


def get_contribution_from_box(virtual_indices,occupied_indices,bra_minus_before):

    commutator = len(virtual_indices)

    result = None

    if commutator == 0:
        # First line in box: H
        if bra_minus_before == 0:
            # sum_i (h_ii + F_ii)
            print("This term is very strange")
            quit()
        elif bra_minus_before == 1:
            # sum_ai F_ai E_ai
            result = " plus -O1 v1o1 -E v1o1"
        elif bra_minus_before == 2:
            #1/2 sum_aibj g_aibj E_ai E_bj
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
minus -O2 v1o1x1o2 -E v1o1 y1o2
            """
    elif commutator == 2: 
        # Third line in box [[H,ai],bj]
        if bra_minus_before == 0:
            # 2 L_iajb
            result = "plus 2 -O3 x1y1x2y2"
        elif bra_minus_before == 1:
            # - P_ij^ab (F_ib E_aj + sum_k L_ikjb E_ak - L_cajb E_ci)
            result="""minus -P y1y2 x1x2 -O1 x1y2 -E y1x2
minus -P y1y2 x1x2 -O3 x1o1x2y2 -E y1o1
plus -P y1y2 x1x2 -O3 v1y1x2y2 -E v1x1
            """
        elif bra_minus_before == 2:
            # -P_ij^ab (sum_ck g_ibck E_aj E_ck + sum_ck g_ikcb E_ak E_cj) + sum_kl g_ikjl E_ak E_bl + sum_cd g_cabd E_ci E_dj
            result="""minus -P y1y2 x1x2 -O2 x1y2v1o1 -E y1x2 v1o1
minus -P y1y2 x1x2 -O2 x1o1v1y2 -E y1o1 v1x2
plus -O2 x1o1x2o2 -E y1o1 y2o2
plus -O2 v1y1v2y2 -E v1x1 v2x2
            """
    elif commutator == 3:
        # Fourth line in box [[[H,ai],bj],ck]
        if bra_minus_before == 1:
            # - P_ijk^abc L_jbic E_ak
            result="minus -P y1y2y3 x1x2x3 -O3 x2y2x1y3 -E y1x3"
        elif bra_minus_before == 2:
            # P_ijk^abc (sum_l g_iljc E_ak E_bk - sum_d g_ibdc E_aj E_dk)
            result="""plus -P y1y2y3 x1x2x3 -O2 x1o1x2y3 -E y1o1 y2x3
minus -P y1x1 y2x2 y3x3 -O2 x1y2v1y3 -E y1x2 v1x3
            """
    elif commutator == 4:
        # Fifth line in box [[[[H,ai],bj],ck],dl]
        if bra_minus_before == 2:
            # 1/2 P_ijkl^abcd g_idjc E_al E_bk
            result = "plus 1/2 -P y1y2y3y4 x1x2x3x4 -O2 x1y4x2y3 -E y1x4 y2x3"
    
    return result


def init_global_variables() -> None:

    global VIR, OCC

    VIR = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
    OCC = ('i', 'j', 'k', 'l', 'm', 'n', 'o', 'p')

def main() -> None:

    init_global_variables()

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""A script designed to be used to quickly expand a commutator using commutator relations""", epilog="""For help contact
        Theo Juncker von Buchwald
        fnc970@alumni.ku.dk""")

    parser.add_argument('term', type=str, nargs=1, help='The central part of the matrix element to calculate...Ex. ck,dl[[H,ai],bj] - Accepted operators are H and P for two-electron operators and F and X for one-electron operators')
    parser.add_argument('-bra', type=str, nargs='+',                                help='Excitations in the bra........................Ex. -bra ai bj. If not included the HF state is assumed.')
    parser.add_argument('-transformed', action='store_false', help='Whether the operator is one-index transformed, e.g T1 or trial-vector', dest='transform')

    # Parses the arguments
    args = parser.parse_args()

    term = args.term[0]
    bra = args.bra
    if bra is None:
        bra = []


    before, v_in, o_in, v_total, o_total = get_occ_and_vir_indices(term,bra)

    occ_not_used = [i for i in OCC if i not in o_total]
    vir_not_used = [i for i in VIR if i not in v_total]

    # check if term gives something from Slater-Condon
    if len(bra) < len(before) or len(bra) > len(before) + 2:
        output = None
    else:
        output = get_contribution_from_box(v_in,o_in,len(bra)-len(before))

    if output is not None:

        if before and 'E' not in output:
            output += ' -E ' + ' '.join(before)
        elif before:
            output += ' '+' '.join(before)

        output = output.replace('O1','F')
        output = output.replace('O2','g')
        output = output.replace('O3','L')

        if 'o' in output or 'v' in output:
            output = output.replace('o1',occ_not_used[0])
            output = output.replace('v1',vir_not_used[0])
            if 'o2' in output or 'v2' in output:
                output = output.replace('o2',occ_not_used[1])
                output = output.replace('v2',vir_not_used[1])

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
        
        print(output)
    else:
        print("This term does not give anything")

    



    
if __name__ == '__main__':
    main()
