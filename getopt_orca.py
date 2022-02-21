
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''A script to find the optimized structure from an orca output file''', epilog='''For help contact
    Theo Juncker von Buchwald
    fnc970@alumni.ku.dk''')

    parser.add_argument('infile', type=str, nargs='+', help='The file(s) to extract data from', metavar='.out orca file')


    args = parser.parse_args()
    input_files = args.infile

    for filename in input_files:
        newfile = str(filename)[:-4] + "_opt.xyz"

        start = 0
        end = 0

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

