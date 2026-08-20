import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'modules')))
from modules import *

from openEMS import openEMS
import numpy as np
import matplotlib.pyplot as plt

# Model comments
# Test derived layers used to model IHP resistors in SG13G2_resistors_200um.xml
# The resistor connected to ports on layers 201,202 is Rsil with nominal value of 50 Ohm.


# ======================== workflow settings ================================

preview_only = False  # preview model/mesh only?
postprocess_only = False # postprocess existing data without re-running simulation?
no_gui = False # if set to True, there is no mesh/model preview in AppCSXCAD, and simulation starts immediately. 

# ===================== input files and path settings =======================

gds_filename = "resistors_with_ports.gds"   # geometries
XML_filename = "SG13G2_with_resistors_200um.xml"  # stackup

# merge via polygons with distance less than .. microns, set to 0 to disable via merging.
merge_polygon_size = 1.0

# get path for this simulation file
script_path = utilities.get_script_path(__file__)

# use script filename as model basename
model_basename = utilities.get_basename(__file__)

# set and create directory for simulation output
sim_path = utilities.create_sim_path (script_path,model_basename)
print('Simulation data directory: ', sim_path)

# change current path to model script path
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# ======================== simulation settings ================================

unit   = 1e-6   # geometry is in microns
margin = 10    # distance in microns from GDSII geometry boundary to simulation boundary 

fstart = 0
fstop  = 100e9
numfreq = 401

refined_cellsize = 0.2  # mesh cell size in conductor region

numThreads = 0 # set 0 for automatic thread selection by openEMS

# choices for boundary: 
# 'PEC' : perfect electric conductor (default)
# 'PMC' : perfect magnetic conductor, useful for symmetries
# 'MUR' : simple MUR absorbing boundary conditions
# 'PML_8' : PML absorbing boundary conditions
Boundaries = ['PEC', 'PEC', 'PEC', 'PEC', 'PEC', 'PEC']

cells_per_wavelength = 20   # how many mesh cells per wavelength, must be 10 or more
energy_limit = -40          # end criteria for residual energy (dB)

# ports from GDSII Data, polygon geometry from specified special layer
# note that for multiport simulation, excitations are switched on/off in simulation_setup.createSimulation below
simulation_ports = simulation_setup.all_simulation_ports()
# instead of in-plane port specified with target_layername, we here use via port specified with from_layername and to_layername
simulation_ports.add_port(simulation_setup.simulation_port(portnumber=1, voltage=1, port_Z0=50, source_layernum=201, from_layername='Metal1', to_layername='Metal2', direction='z'))
simulation_ports.add_port(simulation_setup.simulation_port(portnumber=2, voltage=1, port_Z0=50, source_layernum=202, from_layername='Metal1', to_layername='Metal2', direction='z'))

# ======================== simulation ================================

# get technology stackup data
materials_list, dielectrics_list, metals_list = stackup_reader.read_substrate (XML_filename)
# get list of layers from technology
layernumbers = metals_list.getlayernumbers()
layernumbers.extend(simulation_ports.portlayers)

# read geometries from GDSII, only purpose 0
allpolygons = gds_reader.read_gds(gds_filename, layernumbers, purposelist=[0], metals_list=metals_list, merge_polygon_size=merge_polygon_size)


# calculate maximum cellsize from wavelength in dielectric
wavelength_air = 3e8/fstop / unit
max_cellsize = (wavelength_air)/(np.sqrt(materials_list.eps_max)*cells_per_wavelength) 



########### create model, run and post-process ###########

FDTD = openEMS(EndCriteria=np.exp(energy_limit/10 * np.log(10)))
FDTD.SetGaussExcite( (fstart+fstop)/2, (fstop-fstart)/2 )
FDTD.SetBoundaryCond( Boundaries )

if preview_only==False:

    # run all port excitations, one after another
    for port in simulation_ports.ports:
        simulation_setup.setupSimulation   ([port.portnumber], 
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
        
        simulation_setup.runSimulation  ([port.portnumber], 
                                            FDTD, 
                                            sim_path, 
                                            model_basename, 
                                            preview_only, 
                                            postprocess_only,
                                            no_gui=no_gui)


    # Initialize an empty matrix for S-parameters
    num_ports = simulation_ports.portcount
    s_params = np.empty((num_ports, num_ports, numfreq), dtype=object)

    # Define frequency resolution. Due to FFT from Empire time domain results, 
    # this is postprocessing and we can change it again at any time.
    f = np.linspace(fstart,fstop,numfreq)

    # Populate the S-parameter matrix with simulation results
    for i in range(1, num_ports + 1):
        for j in range(1, num_ports + 1):
            s_params[i-1, j-1] = utilities.calculate_Sij(i, j, f, sim_path, simulation_ports)

    # Write to Touchstone *.snp file
    snp_name = os.path.join(sim_path, model_basename + '.s' + str(num_ports) + 'p')
    utilities.write_snp(s_params, f, snp_name, z0=simulation_ports.get_reference_impedance())

    print('Created S-parameter output file at ', snp_name)