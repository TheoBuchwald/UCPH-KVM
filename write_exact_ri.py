import sys

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


outputs = [s_func,p_func,d_func,f_func,g_func]
labels = ["S","P","D","F","G"]


n=len(exp)

for i in range(n):
    for j in range(i,n):
        new_q = quant_num[i] + quant_num[j]
        outputs[new_q].append(exp[i] + exp[j])


for x in outputs:
    x.sort(reverse=True)

char_const_1 = "  0.00000000"
char_const_2 = "  1.00000000"

print(len(s_func))

lines_to_write = []

for y,x in enumerate(outputs):
    num_bas = len(x)
    lines_to_write.append(f"$ {labels[y]}-TYPE FUNCTIONS\n")
    lines_to_write.append(f"{num_bas}    {num_bas}    0\n".rjust(17))
    for i, j in enumerate(x):
        line_to_add = [char_const_1 for k in range(num_bas)]
        line_to_add[i] = char_const_2
        line_to_add.insert(0,f"{j:.7f}".rjust(15))
        lines_to_write.append(''.join(line_to_add)+'\n')

outfile = 'ribas.txt'

with open(outfile,'w') as f:
    f.writelines(lines_to_write)
