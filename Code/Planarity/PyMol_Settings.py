from pymol import cmd

# (1) Remove all other unnecessary stuffs
# (2) color the ligand different color
# (3) run “set stick_h_scale, 1”
# (4) hydrogen→remove nonpolar hydrogen
# (5) remove valences
# (6) Setting→rending→shadows→non
# (7) display → background → white
# (8) remove depth cue (under display)
# (9) find the interactions
# (10) color the interactions
# (11) ray 1000,1000
# (12) file → export image as → save

def DFHBI():
    cmd.show('sticks')
    cmd.color('magenta', 'resn 38E')
    cmd.remove('hydro and (elem H) and (neighbor elem C)')
    cmd.remove('resn HOH')
    cmd.set('valence', 0)
    cmd.bg_color('white')
    cmd.set('depth_cue', 0)
