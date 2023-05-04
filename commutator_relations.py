
import argparse
from typing import Union

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

def main() -> None:

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="""A script designed to be used to quickly expand a commutator using commutator relations""", epilog="""For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk""")

    parser.add_argument('commutator', type=str, nargs=1, help='The commutator to expand...Ex. [[X,ai,bj,ck],dl,em] - Accepted operators are H and P for two-electron operators and F and X for one-electron operators')

    # Parses the arguments
    args = parser.parse_args()

    commutator = args.commutator[0]
    operator = commutator.split("[")[-1][0]
    start_brackets = commutator.count("[")
    end_brackets = commutator.count("]")
    
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

    # Prints terms to terminal
    for term in sorted(commutator_expanded_terms,key=len,reverse=True):
        print(term)

if __name__ == '__main__':
    main()
