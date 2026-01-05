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


# Model comments
# 
# This is a generic model running port excitation for all ports defined below, 
# to get full [S] matrix data.
# Output is stored to Touchstone S-parameter file. 
# No data plots are created by this script.


# ======================== workflow settings ================================

# preview model/mesh only?
# postprocess existing data without re-running simulation?
preview_only = False
postprocess_only = False
simulate_dumps = False          #Turn this on to allow creation of dumps !!High storage consumption!!
preview_first_excitation = True   #Turn this on to preview the model for first port excitation

# ===================== input files and path settings =======================

gds_filename = "line_simple_viaport.gds"   # geometries
XML_filename = "SG13G2_nosub.xml"          # stackup

# preprocess GDSII for safe handling of cutouts/holes?
preprocess_gds = False

# merge via polygons with distance less than .. microns, set to 0 to disable via merging.
merge_polygon_size = 0


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

unit   = 1e-6  # geometry is in microns
margin = 50    # distance in microns from GDSII geometry boundary to simulation boundary 

fstart =  0e9
fstop  = 110e9
numfreq = 401

refined_cellsize = 1  # mesh cell size in conductor region

# choices for boundary: 
# 'PEC' : perfect electric conductor (default)
# 'PMC' : perfect magnetic conductor, useful for symmetries
# 'MUR' : simple MUR absorbing boundary conditions
# 'PML_8' : PML absorbing boundary conditions
Boundaries = ['PEC', 'PEC', 'PEC', 'PEC', 'PEC', 'PEC']

cells_per_wavelength = 20   # how many mesh cells per wavelength, must be 10 or more
energy_limit = -40          # end criteria for residual energy (dB)

# ports from GDSII Data, polygon geometry from specified special layer
# note that for multiport simulation, active port excitation are defined in the loop calling simulation_setup.setupSimulation
simulation_ports = simulation_setup.all_simulation_ports()
# in-plane port is specified with target_layername=
# via port is specified with from_layername= and to_layername=
simulation_ports.add_port(simulation_setup.simulation_port(portnumber=1, voltage=1, port_Z0=50, source_layernum=201, from_layername='Metal1', to_layername='TopMetal2', direction='z'))
simulation_ports.add_port(simulation_setup.simulation_port(portnumber=2, voltage=1, port_Z0=50, source_layernum=202, from_layername='Metal1', to_layername='TopMetal2', direction='z'))

# ======================== simulation ================================

# get technology stackup data
materials_list, dielectrics_list, metals_list = stackup_reader.read_substrate (XML_filename)
# get list of layers from technology
layernumbers = metals_list.getlayernumbers()
# we must also read the layers where we added ports, these are not included in technology layers
layernumbers.extend(simulation_ports.portlayers)

# read geometries from GDSII, only purpose 0
allpolygons = gds_reader.read_gds(gds_filename, layernumbers, purposelist=[0], metals_list=metals_list, preprocess=preprocess_gds, merge_polygon_size=merge_polygon_size)

# calculate maximum cellsize from wavelength in dielectric
wavelength_air = 3e8/fstop / unit
max_cellsize = (wavelength_air)/(np.sqrt(materials_list.eps_max)*cells_per_wavelength) 

# define excitation and stop criteria and boundaries
FDTD = openEMS(EndCriteria=np.exp(energy_limit/10 * np.log(10)))
FDTD.SetGaussExcite( (fstart+fstop)/2, (fstop-fstart)/2 )
FDTD.SetBoundaryCond( Boundaries )


########### create model, run and post-process ###########

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
                                        xy_mesh_function=util_meshlines.create_xy_mesh_from_polygons,
                                        field_dumps=simulate_dumps)
    
    simulation_setup.runSimulation  ([port.portnumber], 
                                        FDTD, 
                                        sim_path, 
                                        model_basename, 
                                        preview_only, 
                                        postprocess_only,
                                        preview_first_excitation)


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
utilities.write_snp(s_params, f, snp_name)

print('Created S-parameter output file at ', snp_name)
