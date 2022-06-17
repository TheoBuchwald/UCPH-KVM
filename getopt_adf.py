# Loads the sys to additional input in the line
# Loads the re search function

import sys

# Set global variables
start = 0
end = 0


#INPUTS
#-----------------------------------------------------
filename = sys.argv[1]
#-----------------------------------------------------

newfile = str(filename)[:-4] + "_opt.xyz"

with open(filename,'r') as f:
  rline=f.readlines()

for i in range (len(rline)):
    if "Optimized geometry:" in rline[i]:
        start = i
        break

if start == 0:
   print(f"Error, optimized geometry not found in file {filename} !")
   exit()

for m in range (start + 8, len(rline)):
    if "Total System Charge" in rline[m]:
        end = m-1
        break

# Conversion section

lines_to_add = []

lines_to_add.append(str(end-(start+8))+ '\n')
lines_to_add.append('\n')

for line in rline[start+8 : end] :
    words = line.split()
    label = words[1]
    lines_to_add.append(label + line[15:-1]+'\n')

with open(newfile,'w') as f:
  f.writelines(lines_to_add)

