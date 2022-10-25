import sys

def main() -> None:
    filename = sys.argv[1]

    with open(filename,'r') as f:
        lines = f.readlines()

    exp = []
    quant_num = []

    quant_now = -1

    for i in lines:
        if i.startswith('$'):
            quant_now += 1
            continue
        if '1.00000' in i:
            quant_num.append(quant_now)
            exp.append(float(i.split()[0]))

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

    lines_to_write = []

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

    outfile = 'ribas.txt'

    with open(outfile,'w') as f:
        f.writelines(lines_to_write)

if __name__ == '__main__':
    main()