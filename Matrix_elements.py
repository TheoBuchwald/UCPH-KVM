import argparse
from typing import List, Set, Tuple, Union, Dict
from permutation_checker import permutationChecker, print_check, print_compared,permutationComparison
from commutator_relations import commutator_expansion
from commutator_box import get_box_terms


def convert_to_float(frac_str):
    """
    Stolen from https://stackoverflow.com/questions/1806278/convert-fraction-to-float
    Used for converting fractions (1/2) and integers to floats
    """

    try:
        return float(frac_str)
    except ValueError:
        try:
            num, denom = frac_str.split('/')
        except ValueError:
            return None
        try:
            leading, num = num.split(' ')
        except ValueError:
            return float(num) / float(denom)        
        if float(leading) < 0:
            sign_mult = -1
        else:
            sign_mult = 1
        return float(leading) + sign_mult * (float(num) / float(denom))

def init_global_variables() -> None:

    global VIR, OCC

    VIR = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
    OCC = ('i', 'j', 'k', 'l', 'm', 'n', 'o', 'p')

def commutator_relations(commutator, operator, start_brackets, end_brackets):

    if operator == 'H' or operator == 'P':
        max_com = 4
    elif operator == 'F' or operator == 'X':
        max_com = 2 
    else:
        print(f"ERROR! Operator {operator} not recognized! Please put either H, P, F, X")
        quit()

    if start_brackets != end_brackets:
        print("Check your brackets, idiot!")
        quit()
    elif start_brackets > max_com:
        print("This is zero, you nimrod!")
        quit()

    # Expands commutator
    commutator_expanded_terms = commutator_expansion(commutator,operator)
    # Sort the terms here, mostly for debugging
    commutator_expanded_terms = sorted(commutator_expanded_terms,key=len,reverse=True)

    for term in sorted(commutator_expanded_terms,key=len,reverse=True):
        # Only prints unique terms with the associated weight
        print(term)

    return commutator_expanded_terms


def commutator_box(term,bra,transform=False):
    """
    Only does canonical integrals for now
    """

    if bra is None:
        bra = []

    if transform:
        print("ERROR! No transformed integrals, as of yet!")
        quit()

    output = get_box_terms(term,bra,transform)

    # if output is not None:   
    #     # Print final result
    #     for line in output:
    #         print(line)
    # else:
    #     print("This term does not give anything")

    return output

def main():

    init_global_variables()

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""Takes a commutator of an arbitrary operator with excitation operators to be applied to the HF state and returns the
       conversion from Box 13.2 in "Molecular Electronic Structure Theory". Keeps track of the indices used in the excitation vectors.""", epilog="""For help contact
       Magnus Bukhave Johansen
       qhw298@alumni.ku.dk""")

    parser.add_argument('commutator', type=str, nargs=1, help='The commutator to expand...Ex. [[X,ai,bj,ck],dl,em] - Accepted operators are H and P for two-electron operators and F and X for one-electron operators')
    parser.add_argument('-bra', type=str, nargs='+',                                help='Excitations in the bra........................Ex. -bra ai bj. If not included the HF state is assumed.')
    parser.add_argument('-F', default=None, nargs=1, type=str,                      help='The indicies of the F component...............Ex. -F ci')
    parser.add_argument('-L', default=None, nargs=1, type=str,                      help='The indicies of the L component...............Ex. -L cile')
    parser.add_argument('-g', default=None, nargs=1, type=str,                      help='The indicies of the g component...............Ex. -g cile')
    parser.add_argument('-t', default=None, nargs='+', type=str,                    help='The indicies of each individual t component...Ex. -t cile dlem')
    parser.add_argument('-LV', default=None, nargs=1, type=str,                     help='The indicies of the left excitaiton vector....Ex. -LV ci')
    parser.add_argument('-RV', default=None, nargs=1, type=str,                     help='The indicies of the right excitaiton vector...Ex. -RV ck')
    parser.add_argument('-res', '--reserved', nargs=1, type=str,  help='The indicies reserved for the outer sum.......Ex. -res aibj')
    parser.add_argument('-sum', '--summation', default='', nargs=1, type=str,       help='The indicies of the summation.................Ex. -sum cdeklm')
    parser.add_argument('-no-reduce', action='store_false',                         help='Whether to reduce indicies after excitations', dest='reduce')
    #parser.add_argument('-transformed', action='store_true', help='Whether the operator is one-index transformed, e.g T1 or trial-vector', dest='transform')

    # Parses the arguments
    args = parser.parse_args()

    commutator = args.commutator[0]
    bra = args.bra

    t: List[str] = args.t

    if isinstance(args.LV,list):
        LV: List[str] = args.LV[0]
    else:
        LV: str = args.LV

    if isinstance(args.RV,list):
        RV: List[str] = args.RV[0]
    else:
        RV: str = args.RV

    if isinstance(args.summation,list):
        summation: List[str] = args.summation[0]
    else:
        summation: str = args.summation

    reserved: List[str] = [i for i in args.reserved[0]]
    reserved = set(reserved).difference(set(summation))
    print(bra)

    arguments: Dict[str, Union(List[str],str)] = {
        'bra': args.bra,
        't': t,
        'LV': LV,
        'RV': RV,
        'reserved': reserved,
        'summation': summation,
        'reduce': args.reduce
    }

    operator = commutator.split("[")[-1][0]
    start_brackets = commutator.count("[")
    end_brackets = commutator.count("]")

    output_1 = commutator_relations(commutator,operator,start_brackets,end_brackets)

    for term in output_1:
        output2 = commutator_box(term,bra)
        if output2 is not None:
            for term2 in output2:
                # Get sign and prefactors, if any
                total_pre = term2.split('-')[0]
                prefactor = 1
                if len(total_pre.split()) == 2: prefactor *= convert_to_float(total_pre.split()[1])
                if "minus" in total_pre :
                    prefactor *= -1

                # Get the excitation operators, if any
                if "E" in term2:
                    E: List[str] = term2.split("E")[1].split("-")[0].split()
                else:
                    E: List[str] = []
                
                arguments["E"] = E 

                # Get the permutation operator, if any
                if "P" in term2:
                    P: List[str] = term2.split("P")[1].split("-")[0].split()
                else:
                    P: List[str] = ["ai"] # PATCHWORK, PLEASE CHANGE
                
                arguments["P"] = P

                # Get the Fock matrix terms
                if "F" in term2:
                    F: str = term2.split("F")[1].split("-")[0].strip()
                else:
                    F: str = None
                
                arguments["F"] = F 

                # Get the two-electron integrals
                if "g" in term2:
                    g: str = term2.split("g")[1].split("-")[0].strip()
                else:
                    g: str = None
                
                arguments["g"] = g

                # Get the L terms integrals
                if "L" in term2:
                    L: str = term2.split("L")[1].split("-")[0].strip()
                else:
                    L: str = None
                
                arguments["L"] = L

                perms, idx_used = permutationChecker(**arguments)

                res_used = idx_used & reserved
                vir_res = ''.join(i for i in sorted(res_used) if i in VIR)
                occ_res = ''.join(i for i in sorted(res_used) if i in OCC)

                bra_out = [i.pop("bra") for i in perms]

                print_check(bra_out,perms)

                perms_compared = permutationComparison(perms, summation, idx_used, vir_res, occ_res)

                print_compared(bra_out,perms_compared, summation)

if __name__ == '__main__':
    main()