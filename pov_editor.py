
import sys
import subprocess
import re

if __name__ == '__main__':
    for file in sys.argv[:]:
        if '.pov' in file:
            povfile = file
        else:
                print(f'{file} was not recognised as either a POV file. Moving on to next file')
                continue

        newpovfile = povfile.replace('.pov','') + '_new.pov'

        global_settings = f'''#version 3.7;
global_settings {{ assumed_gamma 1.8 }}
background {{color rgb <0.744, 0.784, 0.896>}}
    '''

        color_settings = f'''
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

#declare Au_Texture = Metal(<1.00, 0.875, 0.575>)
#declare C_Texture = Non_Metal(<0.15, 0.15, 0.15, 0.1>)
#declare S_Texture = Non_Metal(<1, 0.749, 0, 0.1>)
#declare H_Texture = Non_Metal(<1, 1, 1, 0.1>)
#declare N_Texture = Non_Metal(<0, 0, 1, 0.1>)
#declare O_Texture = Non_Metal(<1, 0, 0, 0.1>)

    '''

        with open(newpovfile, "w") as newpov:
            newpov.write(global_settings)
            newpov.write(color_settings)

            with open(povfile, "r") as oldpov:
                lines_to_skip = 0
                for line in oldpov:
                    if lines_to_skip > 0:
                        lines_to_skip -= 1
                        continue
                    elif 'sky_sphere' in line:
                        lines_to_skip = 10
                        continue
                    elif 'background' in line:
                        continue
                    elif 'global_settings { assumed_gamma' in line:
                        continue
                    elif '#version' in line:
                        continue
                    elif re.compile(r'#declare color[A-O]').search(line):
                        lines_to_skip = 4
                        continue
                    elif 'texture { colorB }' in line:
                        newpov.write('texture { Au_Texture }')
                    elif 'texture { colorA }' in line:
                        newpov.write('texture { O_Texture }')
                    elif 'texture { colorO }' in line:
                        newpov.write('texture { H_Texture }')
                    elif 'texture { colorN }' in line:
                        newpov.write('texture { C_Texture }')
                    elif 'texture { colorF }' in line:
                        newpov.write('texture { S_Texture }')
                    else:
                        newpov.write(line)

        try:
            subprocess.run(['povray', f'{newpovfile}', '+A0.1', '+AM2', '+AG0', '+R5', '-J'])
            # +A0.1 +AM2 +AG0 +R5 -J
            # Runs with the settings, +A0.1: Antialliasing set to 0.1 threshold, +AM2: Antialiasing method 2,
            #   +AG0: Gamma set to 0, +R5: Depth set to 5, -J: Jitter set to off
        except subprocess.CalledProcessError or OSError or FileNotFoundError:
            print('Could not run POVRAY for some reason')
            print('Maybe it is not installed')
