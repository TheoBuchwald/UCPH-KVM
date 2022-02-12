
import numpy as np
import sys
import os
import subprocess


POV_Atom_Parameters = {'H': [1, 0.324],
					   'C': [12, 0.417],
					   'N': [14, 0.417],
					   'O': [16, 0.378],
					   'S': [32, 0.501],
                       'Ti': [44, 0.8],
                       'Cu': [58, 0.65],
                       'Ag': [104, 0.7],
                       'Au': [158, 0.765]}

BondLengths = {'CC': 1.6,
               'CH': 1.6,
               'HC': 1.6,
               'CO': 1.3,
               'OC': 1.3,
               'CS': 1.9,
               'SC': 1.9,
               'CN': 1.5,
               'NC': 1.5,
               'NH': 1.4,
               'HN': 1.4,
               'SH': 1.4,
               'HS': 1.4,
               'TiO': 2.,
               'OTi': 2.,
               'CuCu': 2.5,
               'AgAg': 2.9,
               'AuAu': 2.9}

BondTypes = [key for key in BondLengths.keys()]


class Atom:
    def __init__(self, species, position=[0, 0, 0]):
        self.species = species
        self.position = np.array(position)

        self.mass, self.rad = POV_Atom_Parameters[self.species]

    def toPOV(self):
        return f'Atom( <{self.position[0]}, {self.position[1]}, {self.position[2]}>, {self.rad}, {self.species}_Texture )\n'

class Bond():
    def __init__(self, startAtom, endAtom):
        self.startAtom = startAtom
        self.endAtom = endAtom

    def toPOV(self):
        halfway_point = (self.startAtom.position - self.endAtom.position) / 2 + self.endAtom.position

        startAtom_cylinder = f'Bond( <{self.startAtom.position[0]}, {self.startAtom.position[1]}, {self.startAtom.position[2]}>, \
    <{halfway_point[0]}, {halfway_point[1]}, {halfway_point[2]}>, 0.25, {self.startAtom.species}_Texture )\n'

        endAtom_cylinder = f'Bond( <{self.endAtom.position[0]}, {self.endAtom.position[1]}, {self.endAtom.position[2]}>, \
    <{halfway_point[0]}, {halfway_point[1]}, {halfway_point[2]}>, 0.25, {self.endAtom.species}_Texture )\n'
        return startAtom_cylinder + endAtom_cylinder


def Get_Structure(file):
    atoms = np.array([])
    with open(file) as data:
        for linenumber, line in enumerate(data):
            text = line.split()
            if linenumber >= 2:
                atoms = np.append(atoms, Atom(text[0], [float(text[1]),	float(text[2]),	float(text[3])]))
    return atoms

def Get_Camera(file):
    text = []
    Camera = []
    with open(file) as data:
        for line in data:
            if 'light_source {' in line:
                break
            text.append(line)

    End = len(text)
    for linenumber, line in enumerate(text):
        if 'camera {' in line:
            Start = linenumber
            break

    for line in text[Start:End]:
        if 'location' in line:
            Camera.append(line)
        elif 'angle' in line:
            Camera.append(line)
        elif 'direction' in line:
            Camera.append(line)
    return Camera




if __name__ == '__main__':
    files = [sys.argv[1], sys.argv[2]]

    for i in files:
        if '.xyz' in i:
            xyzfile = i
        elif '.pov' in i:
            povfile = i
        else:
            print(f'{i} was not recognised as either an XYZ or POV file. Exciting programme')
            exit()

    molecule = Get_Structure(xyzfile)

    Camera = Get_Camera(povfile)

    newpovfile = povfile.replace('.pov','') + '_new.pov'

    pos1 = float(Camera[0].split()[1][1:-1])
    pos2 = float(Camera[0].split()[2][:-1])
    pos3 = float(Camera[0].split()[3][:-1])

    default_settings = f'''#version 3.7;
global_settings {{ assumed_gamma 1.8 }}
background {{color rgb <0.744, 0.784, 0.896>}}

#declare camera_location = <{pos1},{pos2},{pos3}> * 0.7;

camera {{
    location camera_location
    {Camera[1][1:-1]}
    {Camera[2][1:-2]}
	up <0, 1, 0>
	right <1, 0, 0> * 1.33333
}}

light_source {{
    camera_location * 1.05
    color rgb <1, 1, 1>
}}

#macro Metal(BaseColor)
    #local CVect3 = BaseColor - <0.00, 0.10, 0.20>;
    #local CVect5 = BaseColor - <0.00, 0.00, 0.00>;
    // Cast CVect as an rgb vector
    #local P3 = rgb CVect3;
    #local P5 = rgb CVect5;
    // Reflection colors, derived from pigment color, 'grayed down' a touch.
    #local RE = P5 * 0.50 + 0.25;
    // Ambient colors, derived from base color
    #local AE = P5 * 0.02 + 0.1;
    // Diffuse values
    // Calculated as 1 - (ambient+reflective+specular values)
    #local DE = max(1-(((RE.red+RE.green+RE.blue)/3) + ((AE.red+AE.green+AE.blue)/3)),0);
    #local F_Metal = finish {{
        brilliance 6
        diffuse DE
        ambient AE
        reflection RE
        metallic 0.9
        specular 0.80
        roughness 1/120
    }}
    texture{{
        pigment{{
            P3
        }}
        finish{{
            F_Metal
        }}
    }}
#end

#macro Non_Metal(BaseColor)
    texture {{

    pigment {{
        color rgbt BaseColor
    }}
    finish {{
        ambient 0.2
        diffuse 0.8
        specular 0.8
    }}
}}
#end

#macro Atom( position, radii, tex )
    sphere{{ position, radii
    texture {{ tex }}
}}
#end

#macro Bond( pos1, pos2, radii, tex )
    cylinder{{ pos1
    pos2
    radii
    texture {{ tex }}
}}
#end

#declare Au_Texture = Metal(<1.00, 0.875, 0.575>)
#declare C_Texture = Non_Metal(<0.15, 0.15, 0.15, 0.1>)
#declare S_Texture = Non_Metal(<1, 0.749, 0, 0.1>)
#declare H_Texture = Non_Metal(<1, 1, 1, 0.1>)
#declare N_Texture = Non_Metal(<0, 0, 1, 0.1>)
#declare O_Texture = Non_Metal(<1, 0, 0, 0.1>)

union {{
'''

    with open(newpovfile, "w") as newpov:
        newpov.write(default_settings)

        for atom in molecule:
            newpov.write(atom.toPOV())

        for i, atom1 in enumerate(molecule):
            for atom2 in molecule[i+1:]:
                Bond_Type = atom1.species + atom2.species
                Bond_Length = np.linalg.norm(atom1.position - atom2.position)
                if (Bond_Type in BondTypes) and (abs(Bond_Length) <= BondLengths[Bond_Type]):
                    bond = Bond(atom1, atom2)
                    newpov.write(bond.toPOV())

        newpov.write('\n}')

    # Check to see if POVRAY is installed
    Check_for_povray = subprocess.run(['dpkg', '-s', 'povray'], capture_output=True, text=True)
    Check_for_povray = Check_for_povray.stdout
    # If so generate the picture
    if 'Status: install ok installed' in Check_for_povray[1]:
        subprocess.run(['povray', f'{newpovfile}', '+A0.1', '+AM2', '+AG0', '+R5', '-J'])
        # Runs with the settings, +A0.1: Antialliasing set to 0.1 threshold, +AM2: Antialiasing method 2,
        #   +AG0: Gamma set to 0, +R5: Depth set to 5, -J: Jitter set to off
