# Loads the sys to additional input in the line
# Loads the re search function

import sys

# Set global variables
start = 0
end = 0

# This comes from the input
# Gives the name to a possible new file

filename = sys.argv[1]
newfile = str(filename)[:-4] + "_opt.xyz"

# Read the entire original file
with open(filename,'r') as f:
    rline = f.readlines()

#Find all the times coordinates are listed in the output; the last are retained
for i in range (len(rline)):
    if "CARTESIAN COORDINATES (ANGSTROEM)" in rline[i]:
        start = i+2
        for m in range (start, len(rline)):
            if "CARTESIAN COORDINATES (A.U.)" in rline[m]:
                end = m-2
                break

# Write lines in output file
lines_to_add = []

lines_to_add.append(f"{end-(start)}\n")
lines_to_add.append('\n')

for line in rline[start : end] :
    words = line.split()
    lines_to_add.append(f"{words[0]}  {words[1]}  {words[2]}  {words[3]}\n")

with open(newfile,'w') as f:
    f.writelines(lines_to_add)

