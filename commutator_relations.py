
import argparse
from typing import Union
from itertools import permutations

def relation_unpack(E: list) -> tuple[str]:

    # Pure magic
    com1 = ",".join(E[:-1]) + f"],{E[-1]}]"
    com2 = ",".join(E[:-1]) + "]"
    bef2 = f"{E[-1]}"
    com3 = f"{E[-1]}" + "]"
    bef3 = ",".join(E[:-1])

    return com1, com2, com3, bef2, bef3

def relation(after: list[str], before: str, operator: str) -> Union[list[str],bool]:

    # Loop over nested commutators
    for n, i in enumerate(after):

        # E operators in the nested level of the commutator
        E = i.split(",")[1:]

        # Handle inner commutator
        inner = after[:n] or after[0][0]
        if len(inner) > 1:
            inner = f"{inner[0]}]" + "]".join(inner[1:]) + "]"
        elif type(inner) is list:
            inner = f"{inner[0].split(',')[0]},{inner[0].split(',')[1]}]"

        # Only do any relations when there are more than one E operator
        if len(E) > 1:

            # Unpack the commutator using the commutator relation
            com1, com2, com3, bef2, bef3 = relation_unpack(E)

            # Finds needed amount of [ in the commutator
            bef_com1 = (com1.count("]")+len(after[n+1:])+inner.count("]"))*"["
            bef_com2 = (com2.count("]")+len(after[n+1:])+inner.count("]"))*"["
            bef_com3 = (com3.count("]")+len(after[n+1:])+inner.count("]"))*"["

            # Collect terms
            term1 = before + bef_com1 + inner + "," + com1 + "]".join(after[n+1:]) + (len(after[n+1:]) > 0)*"]"
            term2 = bef2 + "," + before + bef_com2 + inner + "," + com2 + "]".join(after[n+1:]) + (len(after[n+1:]) > 0)*"]"
            term3 = bef3 + "," + before + bef_com3 + inner + "," + com3 + "]".join(after[n+1:]) + (len(after[n+1:]) > 0)*"]"

            # Add terms to output only if the commutator is non-zero due to particle rank

            if operator == 'H' or operator == 'P':
                max_com = 4
            elif operator == 'F' or operator == 'X':
                max_com = 2

            output = []
            if len(bef_com1) <= max_com:
                output.append(term1)
            if len(bef_com2) <= max_com:
                output.append(term2)
            if len(bef_com3) <= max_com:
                output.append(term3)

            return output

    # If no terms can be expanded return False
    return False

def commutator_expansion(commutator: str, operator : str) -> list[str]:

    # Put the input commutator in a list for future processing
    terms = [commutator]

    # Initialize iteration number
    iter = 0

    # Number of excitation operators
    # Used to break while loop, as the amount of E operators equal the max amount of iterations needed for the worst case
    nr_E = commutator.count(",")

    while True:

        # Initialize collection variable
        new_terms = []

        # Loop over terms
        for term in terms:

            # Split commutator based on [ and take everything after last occurence
            after = term.split("[")[-1].split("]")[:-1]

            # Split commutator based on [ and take everything before first occurence
            # Needed for iteration two and onwards
            before = term.split("[")[0]

            # Handle commutator based on the commutator relation [X,AB] = [[X,A],B] + B[X,A] + A[X,B]
            relation_out = relation(after,before,operator)

            # If new terms are returned add these to the new terms
            if relation_out:
                new_terms += relation_out
            # Else add the old term
            else:
                new_terms.append(term)

        # Prepare next iteration
        terms = new_terms

        # Break if the amount of iterations have reached the amount of E operators
        if iter == nr_E:
            break

        iter += 1

    return terms

def get_occ_and_vir_indices(term):

    vir_indices = [*term.split("[")[0][::3]]
    occ_indices = [*term.split("[")[0][1::3]]
   
    for com in term.split("[")[-1].split("]"):
        vir_indices += [i[0] for i in com.split(',')[1:]]
        occ_indices += [i[1] for i in com.split(',')[1:]]
    
    return vir_indices, occ_indices

def check_symmetry_weights(expanded_terms):

    weights = [1 for i in expanded_terms]

    no_of_terms = len(expanded_terms)

    for i in range(no_of_terms):
        term1 = expanded_terms[i]
        vir_indices_1, occ_indices_1 = get_occ_and_vir_indices(term1)

        #Check if this has been shown to be the same as a previous term, if so skip
        if weights[i] == 0:
                continue

        for j in range(i+1,no_of_terms):
            term2 = expanded_terms[j]

            #Check if different type of commutator, if so skip
            if term1.count("[") != term2.count("["):
                continue
            #Check if this has been shown to be the same as a previous term, if so skip
            if weights[j] == 0:
                continue

            # print(term1,term2)
            vir_indices_2, occ_indices_2 = get_occ_and_vir_indices(term2)
            # print(vir_indices_1,vir_indices_2)
            # for v1,v2 in zip(vir_indices_1,vir_indices_2):
            #     if v1 != v2:
            #         print(f"element at index {vir_indices_1.index(v1)} is at {vir_indices_2.index(v1)}")


            perm_occ = permutations(occ_indices_2)
            perm_vir = permutations(vir_indices_2)
            for o2, v2 in zip(perm_occ,perm_vir):
                if list(o2) == occ_indices_1 and list(v2) == vir_indices_1:
                    weights[i] += 1
                    weights[j] = 0
                    break
    
    return weights


def main() -> None:

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""A script designed to be used to quickly expand a commutator using commutator relations""", epilog="""For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk""")

    parser.add_argument('commutator', type=str, nargs=1, help='The commutator to expand...Ex. [[X,ai,bj,ck],dl,em] - Accepted operators are H and P for two-electron operators and F and X for one-electron operators')
    parser.add_argument('-no-reduce', action='store_false',                         help='Whether to reduce indices in the final terms', dest='reduce')

    # Parses the arguments
    args = parser.parse_args()

    commutator = args.commutator[0]
    operator = commutator.split("[")[-1][0]
    start_brackets = commutator.count("[")
    end_brackets = commutator.count("]")
    reduce = args.reduce

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
    
    # Symmetry checker
    # Assumes that all indices are interchangeable and does not yet reduce indices
    if reduce:
        weights = check_symmetry_weights(commutator_expanded_terms)
    else:
        weights = [1 for i in commutator_expanded_terms]

    # Prints terms to terminal
    for term, weight in zip(sorted(commutator_expanded_terms,key=len,reverse=True),weights):
        # Only prints unique terms with the associated weight
        if weight == 0:
            continue
        elif weight > 1:
            print(f"{weight} {term}")
        else:
            print(term)

if __name__ == '__main__':
    main()
