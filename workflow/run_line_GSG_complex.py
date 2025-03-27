import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'modules')))

import modules.util_stackup_reader as stackup_reader
import modules.util_gds_reader as gds_reader
import modules.util_utilities as utilities
import modules.util_simulation_setup as simulation_setup
import modules.util_meshlines as util_meshlines

from openEMS import openEMS
import numpy as np
import matplotlib.pyplot as plt

# Model comments
# This model has one challenge: 
# Ports are composite, with each 50 Ohm GSG port in the S-params created from two EM ports of 100 Ohm between signal line and each side ground 
# Layout from GDSII requires preprocessing, see switch below
# S2P output is created from "fake" reverse path data, assuming symmetry

# ======================== workflow settings ================================

# preview model/mesh only?
# postprocess existing data without re-running simulation?
preview_only = True
postprocess_only = False

# ===================== input files and path settings =======================

gds_filename = "gsg_through_50ohm.gds"   # geometries
XML_filename = "SG13G2.xml"               # stackup

# preprocess GDSII for safe handling of cutouts/holes?
preprocess_gds = True

# merge via polygons with distance less than .. microns, set to 0 to disable via merging.
merge_polygon_size = 0

# get path for this simulation file
script_path = utilities.get_script_path(__file__)

# use script filename as model basename
model_basename = utilities.get_basename(__file__)

# set and create directory for simulation output
sim_path = utilities.create_sim_path (script_path,model_basename)
print('Simulation data directory: ', sim_path)


# ======================== simulation settings ================================

unit   = 1e-6   # geometry is in microns
margin = 100    # distance in microns from GDSII geometry boundary to simulation boundary 

fstart = 0
fstop  = 350e9
numfreq = 401

refined_cellsize = 1.5  # mesh cell size in conductor region

# choices for boundary: 
# 'PEC' : perfect electric conductor (default)
# 'PMC' : perfect magnetic conductor, useful for symmetries
# 'MUR' : simple MUR absorbing boundary conditions
# 'PML_8' : PML absorbing boundary conditions
Boundaries = ['PEC', 'PEC', 'PEC', 'PEC', 'PEC', 'MUR']

cells_per_wavelength = 20   # how many mesh cells per wavelength, must be 10 or more
energy_limit = -40          # end criteria for residual energy (dB)

# ports from GDSII Data, polygon geometry from specified special layer
# note that for multiport simulation, excitations are switched on/off in simulation_setup.createSimulation below
# CPW port consists of two CSX ports, which are in parallel. One of them requires opposite direction, so that both have 'plus" terminal on the center line
Z0 = 50
simulation_ports = simulation_setup.all_simulation_ports()
simulation_ports.add_port(simulation_setup.simulation_port(portnumber=1, voltage=1, port_Z0=2*Z0, source_layernum=201, target_layername='TopMetal2', direction= 'y'))
simulation_ports.add_port(simulation_setup.simulation_port(portnumber=2, voltage=1, port_Z0=2*Z0, source_layernum=202, target_layername='TopMetal2', direction='-y'))
simulation_ports.add_port(simulation_setup.simulation_port(portnumber=3, voltage=1, port_Z0=2*Z0, source_layernum=203, target_layername='TopMetal2', direction= 'y'))
simulation_ports.add_port(simulation_setup.simulation_port(portnumber=4, voltage=1, port_Z0=2*Z0, source_layernum=204, target_layername='TopMetal2', direction='-y'))

# ======================== simulation ================================

# get technology stackup data
materials_list, dielectrics_list, metals_list = stackup_reader.read_substrate (XML_filename)
# get list of layers from technology
layernumbers = metals_list.getlayernumbers()
layernumbers.extend(simulation_ports.portlayers)

# read geometries from GDSII, only purpose 0
allpolygons = gds_reader.read_gds(gds_filename, layernumbers, purposelist=[0], metals_list=metals_list, preprocess=preprocess_gds, merge_polygon_size=merge_polygon_size)


# calculate maximum cellsize from wavelength in dielectric
wavelength_air = 3e8/fstop / unit
max_cellsize = (wavelength_air)/(np.sqrt(materials_list.eps_max)*cells_per_wavelength) 


########### create model, run and post-process ###########

# Create simulation for port 1+2 excitation
# Excited GSG port on left side is composite from CSX ports 1+2, opposite polarity, so we excite 1+2 simultaneously
excite_ports = [1,2]  # list of ports that are excited for this simulation run

FDTD = openEMS(EndCriteria=np.exp(energy_limit/10 * np.log(10)))
FDTD.SetGaussExcite( (fstart+fstop)/2, (fstop-fstart)/2 )
FDTD.SetBoundaryCond( Boundaries )

# use polygon-based meshing
FDTD = simulation_setup.setupSimulation (excite_ports, 
                                         simulation_ports, 
                                         FDTD, 
                                         materials_list, 
                                         dielectrics_list, 
                                         metals_list, 
                                         allpolygons, 
                                         max_cellsize, 
                                         refined_cellsize, 
                                         margin, 
                                         unit, 
                                         xy_mesh_function=util_meshlines.create_xy_mesh_from_polygons)

# run simulation
sub1_data_path = simulation_setup.runSimulation (excite_ports, FDTD, sim_path, model_basename, preview_only, postprocess_only)


########## evaluation of results with composite GSG ports ###########

if not preview_only:

    # get results, CSX port definition is read from simulation ports object
    # left side CPW port consists of CSX ports 1 and 2
    CSX_port1 = simulation_ports.get_port_by_number(1).CSXport
    CSX_port2 = simulation_ports.get_port_by_number(2).CSXport
    # right side CPW port consists of CSX ports 3 and 4
    CSX_port3 = simulation_ports.get_port_by_number(3).CSXport
    CSX_port4 = simulation_ports.get_port_by_number(4).CSXport

    # S-parameters must combine results from both 100 Ohm CSX ports into one 50 Ohm GSG port
    f = np.linspace(fstart,fstop,numfreq)
    CSX_port1.CalcPort(sub1_data_path, f, 2*Z0)
    CSX_port2.CalcPort(sub1_data_path, f, 2*Z0)
    CSX_port3.CalcPort(sub1_data_path, f, 2*Z0)
    CSX_port4.CalcPort(sub1_data_path, f, 2*Z0)

    s11 = (CSX_port1.uf_ref + CSX_port2.uf_ref) / (CSX_port1.uf_inc + CSX_port2.uf_inc)
    s21 = (CSX_port3.uf_ref + CSX_port4.uf_ref) / (CSX_port1.uf_inc + CSX_port2.uf_inc)

    Zin = 0.5 * (CSX_port1.uf_tot + CSX_port2.uf_tot)  / (CSX_port1.if_tot + CSX_port2.if_tot)
    
    # S12, S22 is NOT available because we have NOT simulated port2 excitation
    # fake it by assuming symmetry
    s22 = s11
    s12 = s21

    # write Touchstone S2P file
    s2p_name = os.path.join(sim_path, model_basename + '.s2p')
    utilities.write_snp (np.array([[s11, s21],[s12,s22]]),f, s2p_name)


    # define dB function for S-parameters
    def dB(value):
        return 20.0*np.log10(np.abs(value))        

    # define phase function for S-parameters
    def phase(value):
        return np.angle(value, deg=True) 

      
    fig, axis = plt.subplots(num='Return Loss', tight_layout=True)
    axis.plot(f/1e9, dB(s11), 'k-',  linewidth=2, label='S11 (dB)')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_title('Return Loss')
    axis.legend()
    
    fig, axis = plt.subplots(num='Insertion Loss', tight_layout=True)
    axis.plot(f/1e9, dB(s21), 'k-',  linewidth=2, label='S21 (dB)')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_title('Insertion Loss')
    axis.legend()

    fig, axis = plt.subplots(num='Transmission Phase', tight_layout=True)
    axis.plot(f/1e9, phase(s21), 'k-',  linewidth=2, label='S21 (dB)')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_title('Transmission Phase')
    axis.legend()

    fig, axis = plt.subplots(num='Input Impedance', tight_layout=True)
    axis.plot(f/1e9, np.real(Zin), 'k-',  linewidth=2, label='Real(Zin)')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_title('Input Impedance')
    axis.legend()

    #portVariable.u_data.ui_val[0] = U over time
    #portVariable.i_data.ui_val[0] = I over time
    #portVariable.u_data.ui_time[0] = timesteps
    #portVariable.i_data.ui_time[0] = timesteps (should be same as for voltage)

    u1 = CSX_port1.u_data.ui_val[0]
    u2 = CSX_port2.u_data.ui_val[0]
    u3 = CSX_port3.u_data.ui_val[0]
    u4 = CSX_port4.u_data.ui_val[0]
    t  = CSX_port1.u_data.ui_time[0]

    fig, axis = plt.subplots(num='Port voltages (V)', tight_layout=True)
    axis.plot(t*1e9, u1, 'r-',  linewidth=2, label='u1')
    axis.plot(t*1e9, u2, 'k:',  linewidth=2, label='u2')
    axis.plot(t*1e9, u3, 'b-',  linewidth=2, label='u3')
    axis.plot(t*1e9, u4, 'y:',  linewidth=2, label='u4')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Time (ns)')
    axis.set_title('Input Impedance')
    axis.legend()

   # show all plots
    plt.show()
