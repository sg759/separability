#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script serves as an example to run a LAMMPS simulation and directly supply
it's output to the SIF module to infer crack tiop location and it's associated K-field
using one or all 3 approaches.

INPUTS:
    infile:      a LAMMPS input file pr string of commands
    
    r2:          outer radius of annulus
    
    alpha:       scaling factor that determines the size of the inner radius of the annulus
                 
    guess:       initial guess of crack top location. Should be an array of two
                 coordinates, eg, [x,y]
    
    timestep:    timestep at which the LAMMPS simulation should be run and the output
                 should use used to infer the crack tip location
    
    end:         timestep at which the LAMMPS simulation should stop
    
    constants:   values of elastic constants, if known. The default values are (in order)
                 poisson =  0.25          # poisson's ratio
                 E = 70                   # Young's modulus
                 shear_modulus = 0.25  
                 
                 These are used to calculate the associated K-field values.
    
    Control variables: 
                geo: use Geometrical method (geo = 1) or not (geo = 0)
                Sep: use Separability method (Sep = 1) or not (Sep = 0)
                DC: use Displacement Corr method (DC = 1) or not (DC = 0)

OUTPUTS:
   cracktip1:    crack tip location found using geometrical method
   K_field1:     K field associated with cracktip found using geometrical method

   cracktip2:   crack tip location found using seperability 
   K_field2:    K field associated with the cracktip found using seperability 
   
   cracktip3:   crack tip location found using DC 
   K_field3:    K field associated with the cracktip found using DC 

@author: Swati Gupta

"""

from lammps import lammps
from mpi4py import MPI
import numpy as np
import sys
import SIF_final as SIF
import pdb


def inferCrack(infile='', R2 = 20, alpha = 2, timestep = 500, end = 5000, guess = [10,40], 
               sep = 1, geo = 0, dc = 0, constants = []):
        
    lmp = lammps()
    if len(infile) > 1:
        # infile = sys.argv[1]
        try:
            lmp.file(infile) #run simulation
        except Exception as e:
            print("An error occured: \n", e)
            return
    else: 
        # if no infile or input command string, run the following set of commands
        
        #the following string is a copy of in.crack fropm lammps examples
        # https://www.lammps.org/inputs/in.crack.txt
        crack2d = '''# 2d LJ crack simulation
        dimension	2
        boundary	s s p
        atom_style	atomic
        neighbor	0.3 bin
        neigh_modify	delay 5
        # create geometry
        lattice		hex 0.93
        region		box block 0 100 0 40 -0.25 0.25
        create_box	5 box
        create_atoms	1 box
        mass		1 1.0
        mass		2 1.0
        mass		3 1.0
        mass		4 1.0
        mass		5 1.0
        # LJ potentials
        pair_style	lj/cut 2.5
        pair_coeff	* * 1.0 1.0 2.5
        # define groups
        region	        1 block INF INF INF 1.25 INF INF
        group		lower region 1
        region		2 block INF INF 38.75 INF INF INF
        group		upper region 2
        group		boundary union lower upper
        group		mobile subtract all boundary
        region		leftupper block INF 20 20 INF INF INF
        region		leftlower block INF 20 INF 20 INF INF
        group		leftupper region leftupper
        group		leftlower region leftlower
        set		group leftupper type 2
        set		group leftlower type 3
        set		group lower type 4
        set		group upper type 5
        # Define Settings 
        #compute         disp all displace/atom refresh check
        #compute         disp all displace/atom
        # initial velocities
        compute	  	new mobile temp
        velocity	mobile create 0.01 887723 temp new
        velocity	upper set 0.0 0.3 0.0
        velocity	mobile ramp vy 0.0 0.3 y 1.25 38.75 sum yes
        # fixes
        fix		1 all nve
        fix		2 boundary setforce NULL 0.0 0.0
        #fix 		disp all store/state vx vy vz
        # run
        #timestep	0.003
        #thermo		200
        thermo_modify	temp new
        neigh_modify	exclude type 2 3
        '''
        
        lmp = lammps(cmdargs=["-nocite"])
        lmp.commands_string(crack2d)
    
    # The simulations is now alive and we can perform computations
    start = 0
    step= 0
    cmd = "run " + str(step)
    lmp.commands_list([
        cmd
     ])
    
    # get the atomic positions at the first run
    coords = lmp.numpy.extract_atom("x");
    # assuming 2D crack
    ref_coords = np.vstack([coords[:,0], coords[:,1]]).T
    i = 0
    trials = int(np.floor(end/timestep))
    cracktip1 = np.zeros((trials,2)) # Geometrical (geo) method
    cracktip2 = np.zeros((trials,2)) # Separability (Sep)
    cracktip3 = np.zeros((trials,2)) # Displacement Correlation (DC)
    K_field1 = np.zeros((trials,6))  # Geo
    K_field2 = np.zeros((trials,6))  # Sep
    K_field3 = np.zeros((trials,6))  # Dc
    prev = guess
                
    while step <= (end - timestep):
        cmd = "run " + str(timestep)
        lmp.commands_list([
           cmd
        ])
        step = step + timestep
        coords = lmp.numpy.extract_atom("x")
        coords2d =  np.vstack([coords[:,0], coords[:,1]]).T
        # dx = np.array(coords2d - ref_coords)
        cracktip1[i],cracktip2[i],cracktip3[i],K_field1[i], K_field2[i], K_field3[i] = \
        SIF.SIF_projection(synData = 2, coords = coords2d, coords_ref = ref_coords,
                           guess = prev, geo = geo, DC = dc, Sep = sep) 
        prev = cracktip2[i]
        i = i+1
        
    
    print("\n \n PYTHON OUTPUT:")
   
    natoms = lmp.extract_global("natoms")
    mass = lmp.extract_atom("mass")
    print("\n #atoms, mass =",natoms,mass[1])
    
    print("\n cracktips using separability =",cracktip2)
    if dc:
        print("\n cracktips using displaceemnt correlation =",cracktip3)
    
   
   