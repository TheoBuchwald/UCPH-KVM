import argparse
import re
from itertools import combinations
from collections import Counter, defaultdict
from fractions import Fraction
import csv
from io import StringIO
from typing import TextIO

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""A script designed to convert opt.einsum contractions into standalone ScaTeLib-depending subroutines for the code necessary to calculate a Coupled Cluster matrix element.
The input must contain a a read file and a write file. Examples being:
    path/to/read_file.txt path/to/write_file.txt
    """, epilog="""For help contact
    Phillip Gustav Iuel Lunøe Dünweber
    pgd@chem.ku.dk"""
    )

    parser.add_argument('Subroutine_name', type=str, nargs=1, help="Name of the generated subroutine")
    parser.add_argument('Contraction_file', type=str, nargs=1, help="Read the opt.einsum contractions from this file")
    parser.add_argument('ScaTeLib_file', type=str, nargs=1, help="Write the subroutine to this file")
    parser.add_argument('--indent', type=int, help="Choose the indentation in spaces, default is 2.", dest="indent")

    return parser.parse_args()

def read_contractions_from_file(file_path: str) -> list[dict]:
    """Read opt_einsum-style contraction lines from file and parse them into a list of dicts."""
    patterns = []
    contraction_str_re = re.compile(r'mem\.contract\(["\'](.*?)["\']')
    args_re = re.compile(r'mem\.contract\(["\'].*?["\']\s*,\s*(.*)\)')

    with open(file_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        pattern = {}
        line = line.strip()

        # Determine sign and factor
        if '-=' in line or '= -' in line:
            pattern['sign'] = '-'
        else:
            pattern['sign'] = ''

        factor_match = re.search(r'([+-]?=)\s*(-?\d+(/\d+)?|\d*\.\d+)?\s*\*\s*mem\.contract', line)
        if factor_match and factor_match.group(2):
            pattern['factor'] = float(Fraction(factor_match.group(2)))
        else:
            pattern['factor'] = 1.0

        # Extract contraction string
        contraction_match = contraction_str_re.search(line)
        if not contraction_match:
            continue  # skip if malformed
        input_indices_str, output_indices_str = contraction_match.group(1).split('->')
        input_indices = input_indices_str.split(',')

        # Extract tensor variable names
        args_match = args_re.search(line)
        if not args_match:
            continue  # skip if malformed
        tensor_args = [arg.strip() for arg in args_match.group(1).split(',')]

        for name, indices in zip(tensor_args, input_indices):
            pattern[name] = indices.strip()

        patterns.append(pattern)

    return patterns

def find_optimal_contraction_order(patterns: list[dict]) -> list[list[str]]:
    """Find the optimal contractions from opt.einsum patterns."""
    def compute_vo_label(string: str) -> str:
        v = sum(c in virt_string for c in string)
        o = sum(c in occ_string for c in string)
        return f'v{v}o{o}'

    def build_tensor_contract(surv: str, out_indices: str, comb: str, keys: list[str], values: list[str], sign: str, factor: float, pattern_index: int, first: bool, last: bool) -> str:
        a_idx, b_idx = comb
        a_key, b_key = keys[a_idx], keys[b_idx]
        a_val, b_val = values[a_idx], values[b_idx]
        expr = f"'C({out_indices}) = beta * A({a_val}) B({b_val})"
        if last and pattern_index!=0:
            expr += f" + alpha * C({out_indices})"
        expr += "'"
        args = f"{a_key},{b_key}"
        beta = f"{sign}{factor}D0" if first else "1.0D0"
        alpha = "alpha=1.0D0," if last and pattern_index!=0 else ""
        if last:
            return f"call tensor_contract(r{surv},{expr},{args},{alpha}beta={beta},buf=buf)"
        else:
            return f"call tensor_contract(w{surv},{expr},{args},{alpha}beta={beta},buf=buf)"

    def canonicalize_C_indices(line: str) -> str:
        canonical_string = 'aibjckdlemfngo'
        order_map = {c: i for i, c in enumerate(canonical_string)}

        def repl(match: re.match) -> str:
            inside = match.group(1)
            # Sort letters according to canonical_string order
            sorted_letters = sorted(inside, key=lambda c: order_map.get(c, 100))
            # Re-group into pairs of two letters
            pairs = [''.join(sorted_letters[i:i+2]) for i in range(0, len(sorted_letters), 2)]
            # Join pairs back into string
            return f"C({''.join(pairs)})"

        return re.sub(r"C\((.*?)\)", repl, line)

    def rename_tensor_contracts(code_lines: list[str]) -> list[str]:
        buffer_counts = defaultdict(int)
        latest_buffer_name = {}  # maps original -> latest suffixed name
        updated_lines = []

        for line in code_lines:
            # Match beginning of line: call tensor_contract(wv1o3,...
            match = re.match(r"(call tensor_contract\()(\w+)(,)", line)
            if not match:
                # Replace inputs with latest names
                for orig, latest in latest_buffer_name.items():
                    line = re.sub(rf'\b{orig}\b', latest, line)
                updated_lines.append(line)
                continue

            call_prefix, out_buffer, call_suffix = match.groups()

            # Strip the output portion temporarily
            remainder = line[len(call_prefix + out_buffer + call_suffix):]

            # Replace only inputs in the remainder
            for orig, latest in latest_buffer_name.items():
                remainder = re.sub(rf'\b{orig}\b', latest, remainder)

            # Rename the output buffer if it's a temp
            if out_buffer.startswith('w'):
                buffer_counts[out_buffer] += 1
                new_out = f"{out_buffer}_{buffer_counts[out_buffer]}"
                latest_buffer_name[out_buffer] = new_out
            else:
                new_out = out_buffer  # keep as-is

            # Reassemble the line
            new_line = f"{call_prefix}{new_out}{call_suffix}{remainder}"
            new_line = canonicalize_C_indices(new_line)
            updated_lines.append(new_line)

        return updated_lines

    contractions = []
    surv_order = ['v4o4','v3o3','v4o0', 'v3o1', 'v2o2', 'v1o3', 'v0o4', 'v3o0', 'v2o1', 'v1o2', 'v2o0', 'v1o1', 'v0o2', 'v1o0', 'v0o1', 'v0o0']
    flop_order = ['v4o4', 'v4o3', 'v3o4', 'v4o2', 'v3o3', 'v2o4', 'v3o2', 'v2o3', 'v3o1', 'v2o2', 'v1o3','v1o2']
    virt_string = 'abcdefgh'
    occ_string = 'ijklmnop'

    for pattern_index, pattern in enumerate(patterns):
        sign = pattern.get('sign', '')
        factor = pattern.get('factor', 1.0)
        keys = [k for k in pattern if k not in {'sign', 'factor'}]
        values = [pattern[k] for k in keys]
        scatelib_pattern = []
        first_iter = True
        final_iter = False

        while len(values) > 1:
            final_iter = len(values) == 2
            combos = list(combinations(range(len(values)), 2))
            results = []

            for i, j in combos:
                a, b = values[i], values[j]
                all_chars = a + b
                counts = Counter(all_chars)
                surv = ''.join([c for c in counts if counts[c] == 1])
                flop = ''.join(sorted(set(all_chars)))
                results.append((surv, flop, (i, j)))

            # Sort by shortest surviving indices
            results.sort(key=lambda x: len(x[0]))

            best = None

            for surv, flop, comb in results:
                flop_label = compute_vo_label(flop)
                surv_label = compute_vo_label(surv)
                if surv_label in surv_order:
                    if best is None or flop_order.index(flop_label) >= flop_order.index(best[1]):
                        if best is not None:
                            if flop_order.index(flop_label) == flop_order.index(best[1]):
                                if surv_order.index(surv_label) > surv_order.index(best[0]):
                                    continue
                                else:
                                    best = (surv, flop_label, comb, surv_label)
                                    continue
                        best = (surv, flop_label, comb, surv_label)

            if not best:
                raise ValueError("No valid contraction found.")

            out_surv, _, best_comb, best_surv = best

            scatelib_pattern.append(build_tensor_contract(best_surv, out_surv, best_comb, keys, values, sign, factor, pattern_index, first_iter, final_iter))

            # Remove used values and keys safely (higher index first to avoid shifting)
            i1, i2 = sorted(best_comb, reverse=True)
            for i in [i1, i2]:
                del values[i]
                del keys[i]

            # Add intermediate
            keys.append(f'w{best_surv}')
            values.append(out_surv)
            first_iter = False

        scatelib_pattern = rename_tensor_contracts(scatelib_pattern)
        contractions.append(scatelib_pattern)

    return contractions

def write_subroutine_to_file(write_file: str, contractions: list[list[str]], subroutine_name: str, indent: str) -> None:
    """Write a Fortran subroutine for the contraction"""
    def find_unique_inputs(contractions: list[list[str]]) -> dict[str,str]:
        variable_dict = {}

        for contraction in contractions:
            for call in contraction:
                match = re.search(r'call tensor_contract\((.*)\)', call)
                content_inside = match.group(1)
                reader = csv.reader(StringIO(content_inside), skipinitialspace=True)
                split_values = next(reader)

                if split_values[0] not in variable_dict:
                    if split_values[0][0] == "w" or split_values[0][0] == "r":
                        variable_dict[f'{split_values[0]}'] = f'type(tensor), intent(inout) :: {split_values[0]}'
                    else:
                        variable_dict[f'{split_values[0]}'] = f'type(tensor), intent(in) :: {split_values[0]}'

                if split_values[2] not in variable_dict:
                    if split_values[2][0] == "w":
                        variable_dict[f'{split_values[2]}'] = f'type(tensor), intent(inout) :: {split_values[2]}'
                    else:
                        variable_dict[f'{split_values[2]}'] = f'type(tensor), intent(in) :: {split_values[2]}'

                if split_values[3] not in variable_dict:
                    if split_values[3][0] == "w":
                        variable_dict[f'{split_values[3]}'] = f'type(tensor), intent(inout) :: {split_values[3]}'
                    else:
                        variable_dict[f'{split_values[3]}'] = f'type(tensor), intent(in) :: {split_values[3]}'

        return variable_dict

    def write_subroutine_header(writefile_object: TextIO, variable_dict: dict[str,str], subroutine_name: str, indent: str) -> None:
        header_indent = len('subroutine') + len(f'{subroutine_name}') + len('(&')
        writefile_object.write(f'subroutine {subroutine_name}(   buf &\n')
        for key in variable_dict.keys():
            writefile_object.write(f"{' '*header_indent}&, {key} &\n")
        writefile_object.write(f"{' '*header_indent}&)\n")
        writefile_object.write('\n')
        writefile_object.write('  implicit none\n')
        writefile_object.write('\n')
        writefile_object.write(f"{' '*indent}real(tensor_dp), pointer, intent(inout) :: buf\n")
        for value in variable_dict.values():
            writefile_object.write(f"{' '*indent}{value}\n")
        writefile_object.write('\n')

    def write_subroutine_body(writefile_object: TextIO, contractions: list[list[str]], indent: str) -> None:
        max_line_length = 130
        base_indent = ' ' * indent
        cont_indent = base_indent + '  & '

        for contraction in contractions:
            writefile_object.write('\n')
            for call in contraction:
                full_line = f"{base_indent}{call}\n"
                if len(full_line) <= max_line_length:
                    writefile_object.write(full_line)
                    continue

                # Extract contents inside parentheses
                match = re.match(r'call tensor_contract\((.*)\)', call)
                if not match:
                    writefile_object.write(full_line)  # fallback
                    continue

                inside = match.group(1)
                parts = [p.strip() for p in re.split(r',(?![^()]*\))', inside)]  # split on commas not in parentheses

                # First line with trailing &
                first_line = f"{base_indent}call tensor_contract({parts[0]}, &\n"
                writefile_object.write(first_line)

                # All middle parts with both leading and trailing &
                for part in parts[1:-1]:
                    writefile_object.write(f"{cont_indent}{part}, &\n")

                # Final part without trailing comma, closing parenthesis
                writefile_object.write(f"{cont_indent}{parts[-1]})\n")

    variable_dict = find_unique_inputs(contractions)

    with open(write_file, 'a') as writefile_object:
        write_subroutine_header(writefile_object, variable_dict, subroutine_name, indent)
        write_subroutine_body(writefile_object, contractions, indent)
        writefile_object.write(f'end subroutine {subroutine_name}')
        writefile_object.write('\n')

def main() -> None:
    """Main function."""
    arguments = parse_arguments()
    subroutine_name = arguments.Subroutine_name[0]
    read_file = arguments.Contraction_file[0]
    write_file = arguments.ScaTeLib_file[0]
    indent = arguments.indent if arguments.indent else 2

    patterns=read_contractions_from_file(read_file)
    contractions=find_optimal_contraction_order(patterns)
    write_subroutine_to_file(write_file, contractions, subroutine_name, indent)

if __name__ == "__main__":
    main()