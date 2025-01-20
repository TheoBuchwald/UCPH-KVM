import argparse

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

if __name__ == "__main__":
    main()