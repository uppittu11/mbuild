from collections import OrderedDict

from warnings import warn

import mbuild as mb
from mbuild.compound import Compound
from mbuild.compound import clone
from mbuild.coordinate_transform import force_overlap
from copy import deepcopy
import multiprocessing as mp

import pdb

__all__ = ['reverse_map']

def reverse_map(coarse_grained, mapping_moieties, target=None, parallel=True):
    """ Reverse map an mb.Compound

    Parameters
    ---------
    coarse_grained : mb.Compound
        original structure. Generated from a MDTraj trajectory
    mapping_moieties : dictionary
        Relate CG molecule names to a list of finer-detailed mbuild
        Compounds. Care must be taken that bead indices match with
        list indices.
    target_structure : dictionary
        mb.Compound, optional, default=False
        A target atomistic structure which can be used to reconstruct
        bonding.
        Bond network in the reverse-mapped structure will be completely
        overridden by the bond network from the target atomistic structure
        Care must be taken that atom indices match perfectly

    """

    aa_system = Compound()

    # For each bead, replace it with the appropriate mb compound
    # Iterate through each molecule (set of particles that are bonded together)
    if parallel:
        pool = mp.Pool(processes=mp.cpu_count())
        inp = zip(coarse_grained.children,
                [target]*len(coarse_grained.children),
                [mapping_moieties]*len(coarse_grained.children))
        molecule_list = pool.starmap(reverse_map_molecule, inp)
        for molecule in molecule_list:
            aa_system.add(molecule)
    else:
        [aa_system.add(reverse_map_molecule(molecule, target, mapping_moieties)) for molecule in coarse_grained.children]

    return aa_system

def reverse_map_molecule(molecule, target, mapping_moieties):
    cg_molecule = clone(molecule) # CG molecule
    aa_template = target[molecule.name] # full aa Compound for molecule
    aa_moieties = mapping_moieties[molecule.name] # list of lists of indices for each bead
    aa_molecule = Compound() # this will have the final aa molecule
    cg_to_aa = [] # list of tuples containing (real index, aa atom)

    # now cycle through beads
    for index, bead in enumerate(cg_molecule.particles()):
        aa_atoms = Compound() # placeholder for molecule atoms
        [aa_atoms.add(clone(aa_template.children[i])) for
            i in aa_moieties[index]]
        aa_atoms.translate_to(bead.pos) # shift to cg_bead position
        print(aa_atoms.pos[2])
        cg_to_aa += list(zip(aa_moieties[index], aa_atoms.children))

    # sort atoms in cg_to_aa and add them to the aa_molecule
    cg_to_aa = sorted(cg_to_aa)
    for atom in cg_to_aa:
        aa_molecule.add(clone(atom[1]))

    # add bonds from the template
    aa_template = aa_template.to_trajectory()
    for i,j in aa_template.top.bonds:
        aa_molecule.add_bond([aa_molecule[i.index], aa_molecule[j.index]])

    # equilibrate molecule and shift back to center
    # if the atom names match OpenBabel forcefield naming convention:
    try:
        aa_molecule.energy_minimization(steps=2500)

    # otherwise rename with just element names:
    except:
        atomnames = [i.name for i in aa_molecule] # get the atomnames
        for atom in aa_molecule: # make the atomnames elements
            atom.name=atom.name[0]

        aa_molecule.energy_minimization(steps=2500)

        for i, atom in enumerate(atomnames):
            aa_molecule[i].name = atomnames[i]

    #aa_molecule.translate(center)
    aa_molecule.name = molecule.name
    return aa_molecule




"""

        new_molecule = Compound()
        # Rather than sort through the molecule, which may be unsorted
        # Look at the parent's particles, which will be sorted
        for bead in molecule[0].parent.particles():
            new_atom = clone(mapping_moieties[bead.name])
            cg_to_aa[bead] = new_atom
            new_atom.translate(bead.pos)
            new_molecule.add(new_atom)
        aa_system.add(new_molecule)

    # Go back and include bonds
    if target_structure:
        # If a target atomistic structure is provided, just its bond graph
        # to the reverse-mapped structure
        aa_system.root.bond_graph = None
        target_traj = target_structure.to_trajectory()

        for (i,j) in target_traj.topology.bonds:
            aa_system.add_bond([aa_system[i.index], aa_system[j.index]])

    else:
        # If no target atomistic structure is provided, look at each molecule,
        # working inwards from the ends of the molecule

        cg_bonds = list(coarse_grained.bonds())
        # Repeatedly iterate through the coarse grained bonds, but only bond
        # particles that have a certain number of available ports
        while len(cg_bonds) > 0:
            for p_i, p_j in cg_bonds:
                if 0 < len(cg_to_aa[p_i].available_ports()) <= 1 or \
                    0 < len(cg_to_aa[p_j].available_ports()) <= 1:
                            p_i_port, p_j_port = _find_matching_ports(cg_to_aa[p_i],
                                cg_to_aa[p_j])
                            force_overlap(cg_to_aa[p_i], from_positions=p_i_port,
                                to_positions=p_j_port)
                            cg_bonds.remove((p_i, p_j))


    # Put molecules back after energy minimization
    for cg_particle, aa_particles in cg_to_aa.items():
        aa_particles.translate_to(cg_particle.pos)
"""




def _find_matching_ports(i, j):
    """ Find corresponding ports on two mBuild compounds"""

    i_ports = i.available_ports()
    j_ports = j.available_ports()
    i_port_names = [p.name for p in i.available_ports()]
    j_port_names = [p.name for p in j.available_ports()]
    common_name = list(set(i_port_names).intersection(j_port_names))
    if len(common_name) != 1:
        warn("{} ports were found with corresponding names for"
                " particles {} and {}".format(len(common_name), i,j))
    i_port = [p for p in i.available_ports() if p.name == common_name[0]]
    j_port = [p for p in j.available_ports() if p.name == common_name[0]]
    #for j_port in j_ports:
        #if j_port.name == i_port.name:
            #return i_port, j_port
    return i_port[0], j_port[0]


