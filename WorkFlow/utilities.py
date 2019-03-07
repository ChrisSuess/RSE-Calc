# Utilities.py: utility functions for the Notebook
from __future__ import print_function

def download(pdb_id):
    """
    Downloads a structure from the Protein Data Bank.
    
    Args:
        pdb_id (str): The four-letter code for the structure
        
    Returns:
        mdtraj.Trajectory: The structure.
    """
    import sys
    
    if sys.version_info[0] < 3: # Python 2.x
        from urllib import urlretrieve
    else: # Python 3.x
        from urllib.request import urlretrieve
    import mdtraj as mdt
    import os
    
    pdb_file = pdb_id + '.pdb'
    if not os.path.exists(pdb_file):
        try:
            path = urlretrieve('http://files.rcsb.org/download/{}'.format(pdb_file), pdb_file)
        except:
            print('Error - can not find PDB code {}'.format(pdb_id))
            return None
    t = mdt.load(pdb_file)
    return t

def select(traj, selection):
    """
    Given an mdtraj trajectory, return a new one with just the selected atoms.
    
    Args:
        traj (mdtraj.Trajectory): The full trajectory.
        selection (str): The atom selection string.
        
    Returns:
        mdtraj.Trajectory
    """
    import mdtraj as mdt
    
    selection = traj.topology.select(selection)
    newtop = traj.topology.subset(selection)
    crds = traj.xyz[:, selection]
    return mdt.Trajectory(crds, newtop)

def clean_up(structure):
    """
    Cleans up a PDB file to make it ready for use by Amber tools.
    
    Args:
        structure (object): Something with a .save() method -e.g. mdtraj.Trajectory.
        
    Returns:
        xflow.CompressedFileObject: cleaned-up fata in PDB format.
    """
    from xbowflow import xflowlib
    
    pdb4amber = xflowlib.SubprocessKernel('amber-shell pdb4amber -y -i bad.pdb -o fixed.pdb')
    pdb4amber.set_inputs(['bad.pdb'])
    pdb4amber.set_outputs(['fixed.pdb'])
    fixed = pdb4amber.run(structure)
    return fixed

def reduction(structure):

    from xbowflow import xflowlib

    reduce = xflowlib.SubprocessKernel('amber-shell reduce start.pdb > reduced.pdb')
    reduce.set_inputs(['start.pdb'])
    reduce.set_outputs(['reduced.pdb'])
    reduced = reduce.run(structure)
    return reduced

def params(mol2file):

    from xbowflow import xflowlib
    mol2 = xflowlib.load(mol2file)
 
    antechamber = xflowlib.SubprocessKernel('amber-shell antechamber -i chimeraOut.mol2 -fi mol2 -at gaff -an y -du y -o antechOut.prepc -fo prepc -c gas')    
    antechamber.set_inputs(['chimeraOut.mol2'])
    antechamber.set_outputs(['antechOut.prepc'])
    prepfile = antechamber.run(mol2)
    return prepfile   

def parmcheck(prepc):

    from xbowflow import xflowlib
    
    parmchk = xflowlib.SubprocessKernel('amber-shell parmchk2 -i ligandHchimera.prepc -a Y -f prepc -o ligandHchimera.frcmod')
    parmchk.set_inputs(['ligandHchimera.prepc'])
    parmchk.set_outputs(['ligandHchimera.frcmod'])
    frcmod = parmchk.run(prepc)
    return frcmod

def run_leap(parameters, script, structure):
    """
    Run the Amber leap command on the given structure, using the given script.
    
    Args:
        parameters (str or list of strs): Names of parameter files.
        script (str): The leap input script.
        structure (object): SOmething with a .save() method that can produce a pdb format file.
        
    Returns:
        topology: Amber topology
        coordinates: Amber coordinates
        
    """
    from xbowflow import xflowlib
    import os

    if isinstance(parameters, str):
        params = [parameters]
    else:
        params = parameters
    with open('leap.in', 'w') as f:
        for p in params:
            f.write('source {} \n'.format(p))
        f.write('x = loadpdb x.pdb\n')
        f.write(script)
        f.write('saveamberparm x x.prmtop x.rst7\nquit\n')
    leapin = xflowlib.load('leap.in')
    os.remove('leap.in')
    leap = xflowlib.SubprocessKernel('tleap -f leap.in')
    leap.set_inputs(['leap.in', 'x.pdb'])
    leap.set_outputs(['x.prmtop', 'x.rst7'])
    topology, coordinates = leap.run(leapin, structure)
    print(leap.STDOUT)
    return topology, coordinates


def run_leap_mol2(parameters, script, mol2, frcmod, prepc):
    """
    Run the Amber leap command on the given structure, using the given script.
    
    Args:
        parameters (str or list of strs): Names of parameter files.
        script (str): The leap input script.
        structure (object): SOmething with a .save() method that can produce a pdb format file.
        
    Returns:
        topology: Amber topology
        coordinates: Amber coordinates
        
    """
    from xbowflow import xflowlib
    import os
    
    if isinstance(parameters, str):
        params = [parameters]
    else:
        params = parameters
    with open('leap.in', 'w') as f:
        for p in params:
            f.write('source {} \n'.format(p))
        f.write(script)
        f.write('saveamberparm x x.prmtop x.rst7\nquit\n')
    leapin = xflowlib.load('leap.in')
    os.remove('leap.in')
    leap = xflowlib.SubprocessKernel('tleap -f leap.in')
    leap.set_inputs(['leap.in', 'ligandHchimera.mol2', 'ligandHchimera.frcmod', 'ligandHchimera.prepc'])
    leap.set_outputs(['x.prmtop', 'x.rst7'])
    topology, coordinates = leap.run(leapin, mol2, frcmod, prepc)
    print(leap.STDOUT)
    return topology, coordinates

def display(systems, topology=None):
    """
    Produce an ngl wiget for a molecular system.
    
    Args:
        system (object): An mdtraj.Trajectory or something that can be made into one
        
    Returns:
        ngl wiget for the system.
    """
    import nglview
    import mdtraj as mdt
    
    bigsys = None
    if not isinstance(systems, list):
        lsys = [systems]
    else:
        lsys = systems

    for system in lsys:
        if not isinstance(system, mdt.Trajectory):
            if topology is None:
                sys = mdt.load(system.as_file())
            else:
                sys = mdt.load(system.as_file(), top=topology.as_file())
        else:
            sys = system

        if bigsys is None:
            bigsys = sys
        else:
            bigsys = bigsys + sys

    print(len(bigsys))
    view = nglview.show_mdtraj(bigsys)
    view.add_representation('line', selection='water')
    view.add_representation('licorice', selection='protein')
    return view

def minimize(script, coordinates, topology):
    """
    Run an energy minimisation job.
    
    Args:
        script (str): The script of control parameters.
        coordinates (Amber .crd format): starting coordinates.
        topology (Amber .prmtop format): topology and forcefield information.
        
    Returns:
        final_coordinates (Amber .ncrst format)
        log_file (text file)
    """
    from xbowflow import xflowlib
    import os
    
    with open('script_file', 'w') as f:
        f.write(script)
    script_file = xflowlib.load('script_file')
    os.remove('script_file')
    mini = xflowlib.SubprocessKernel('pmemd -i x.in -c x.crd -p x.prmtop -r x.ncrst -o x.log')
    mini.set_inputs(['x.in', 'x.crd', 'x.prmtop'])
    mini.set_outputs(['x.ncrst', 'x.log'])
    return mini.run(script_file, coordinates, topology)


def run_md(script, coordinates, topology):
    """
    Run an MD job.
    
    Args:
        script (str): The script of control parameters.
        coordinates (Amber .crd format): starting coordinates.
        topology (Amber .prmtop format): topology and forcefield information.
        
    Returns:
        final_coordinates (Amber .ncrst format)
        trajectory (Amber .nc format)
        log_file (text file)
    """
    from xbowflow import xflowlib
    import os
    
    with open('script_file', 'w') as f:
        f.write(script)
    script_file = xflowlib.load('script_file')
    os.remove('script_file')
    md = xflowlib.SubprocessKernel('pmemd -i x.in -c x.crd -p x.prmtop -x x.nc -r x.ncrst -o x.log')
    md.set_inputs(['x.in', 'x.crd', 'x.prmtop'])
    md.set_outputs(['x.ncrst', 'x.nc', 'x.log'])
    return md.run(script_file, coordinates, topology)

def plot_energies(logfile):
    """
    Extract total energy data from a minimization log file
    
    Args:
        logfile: Amber pmemd/sander log file
        
    Returns:
        nsteps (list of ints): step numbers.
        enetgies (list of floats): total energies.
        
    """
    from matplotlib import pyplot as plt
    
    nsteps = []
    energies = []
    with open(logfile.as_file()) as f:
        keep = False
        for line in f:
            if keep:
                nsteps.append(int(line.split()[0]))
                energies.append(float(line.split()[1]))
                keep = False
            else:
                words = line.split()
                if len(words) > 0:
                    keep = line.split()[0] == 'NSTEP'
    plt.plot(nsteps, energies, '-xg')
    plt.xlabel('step number')
    plt.ylabel('total energy (kCal/mol)')
