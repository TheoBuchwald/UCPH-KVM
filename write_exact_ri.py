import re
import sys
from typing import Counter, Dict, Generator, List, Tuple

def main() -> None:
    filename = sys.argv[1]

    with open(filename,'r') as f:
        lines = f.readlines()
        lines = (i for i in lines)

    quant_now = -1
    labels = ["S","P","D","F","G","H","I"]

    lines_to_write = []

    for line in lines:
        if "TYPE FUNCTIONS" in line:

            exp = []
            quant_num = []

            exp, quant_num, next_block = read_basis(lines, exp, quant_now, quant_num)

            RI_quant_contraction = RI_basis_set_contraction(quant_num)

            newline = lines_to_write[-1]
            lines_to_write[-1] = f"{' '.join(newline.split()[:2]):<15}" + '(' + ','.join([f'{RI_quant_contraction[i]}{labels[i].lower()}' for i in RI_quant_contraction]) + ')' + '\n'

            lines_to_write = basis_set_to_ri(exp, quant_num, lines_to_write)
            lines_to_write.append(''.join(next_block))
        else:
            lines_to_write.append(''.join(line))

    if "exact" in filename:
        outfile = f"{filename}-RI"
    else:
        outfile = f"{filename}-exact-RI"

    lines_to_write[0] = lines_to_write[0].replace(filename, outfile)

    with open(outfile,'w') as f:
        f.writelines(lines_to_write)

def RI_basis_set_contraction(quant_num: List[int]) -> Dict[int, int]:
    quant_contraction = Counter(quant_num)
    new_quant_contraction = {}
    for key1, val1 in quant_contraction.items():
        for key2, val2 in quant_contraction.items():
            if key1 == key2:
                for i in range(val1):
                    for j in range(i,val2):
                        try:
                            new_quant_contraction[key1+key2] += 1
                        except KeyError:
                            new_quant_contraction[key1+key2] = 1
            elif key1 < key2:
                for i in range(val1):
                    for j in range (val2):
                        try:
                            new_quant_contraction[key1+key2] += 1
                        except KeyError:
                            new_quant_contraction[key1+key2] = 1
    return new_quant_contraction


def read_basis(lines: Generator[str, None, None], exp: List[float], quant_now: int, quant_num: List[int]) -> Tuple[List[float], List[int], str]:
    quant_now += 1
    amount_of_functions = int(next(lines).split()[0])
    quant_num += amount_of_functions * [quant_now]
    for _ in range(amount_of_functions):
        basis_set_function = next(lines)
        if '1.00000' in basis_set_function:
            exp.append(float(basis_set_function.split()[0]))
        else:
            exit(f"""Encountered a line that did not have a 1.00000 as a contraction coefficient:
    {basis_set_function}
This is a problem for the script. It may arise due to the basis set being contracted, for which this script is not applicable,
or the basis set may contain more than 6 contraction coefficients, in which case the 1.00000 may be on the next line of the
basis set file. For this script to work, all the contraction coeffecients need to be on the same line.""")
    
    next_block = next(lines)
    if "TYPE FUNCTIONS" in next_block:
        exp, quant_num, next_block = read_basis(lines, exp, quant_now, quant_num)

    return exp, quant_num, next_block


def basis_set_to_ri(exp: List[float], quant_num: List[int], lines_to_write: List[str]) -> List[str]:
    s_func = []
    p_func = []
    d_func = []
    f_func = []
    g_func = []
    h_func = []
    i_func = []

    outputs = [s_func,p_func,d_func,f_func,g_func,h_func,i_func]
    labels = ["S","P","D","F","G","H","I"]

    n = len(exp)

    for i in range(n):
        for j in range(i,n):
            new_q = quant_num[i] + quant_num[j]
            outputs[new_q].append(exp[i] + exp[j])

    for x in outputs:
        x.sort(reverse=True)

    const_1 = 0.0
    const_2 = 1.0

    for y,x in enumerate(outputs[:new_q+1]):
        num_bas = len(x)
        lines_to_write.append(f"$ {labels[y]}-TYPE FUNCTIONS\n")
        lines_to_write.append(f"{num_bas}    {num_bas}    0\n".rjust(17))
        for i, j in enumerate(x):

            amount_of_lines, _ = divmod(num_bas+1, 7)

            line_to_add = [f"{j:>15.7f}"] + [f"{const_1:>12.8f}" for _ in range(min(6, num_bas))]

            if i < 6:
                line_to_add[i+1] = f"{const_2:>12.8f}"
            lines_to_write.append(''.join(line_to_add)+'\n')

            for k in range(amount_of_lines):
                line_to_add = [f"{const_1:>15.7f}"] + [f"{const_1:>12.8f}" for _ in range(min(6, num_bas - 7*(k+1)))]
                char_to_replace = i - 6 - 7*k
                if char_to_replace == 0:
                    line_to_add[char_to_replace] = f"{const_2:>15.7f}"
                elif 1 <= char_to_replace <= 6:
                    line_to_add[char_to_replace] = f"{const_2:>12.8f}"
                lines_to_write.append(''.join(line_to_add)+'\n')

    return lines_to_write

if __name__ == '__main__':
    main()