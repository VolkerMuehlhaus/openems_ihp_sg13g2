import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'modules')))

import modules.util_stackup_reader as stackup_reader
import modules.util_gds_reader as gds_reader
import modules.util_utilities as utilities
import modules.util_simulation_setup as simulation_setup
import modules.util_meshlines as util_meshlines

from pylab import *
from CSXCAD import ContinuousStructure
from CSXCAD import AppCSXCAD_BIN
from openEMS import openEMS
from openEMS.physical_constants import *

# Model comments
# 
# Ports have different impedance level, numbers of turns = 1:2 (would expect 1:4 impedance ratio, but actual value closer to 1:8)
# GDSII file has small vias, enabled via array merging in gds_reader.read_gds()
# Meshing: using default mesh, uniform mesh size within drawing region


# ======================== workflow settings ================================

# preview model/mesh only?
# postprocess existing data without re-running simulation?
preview_only = True
postprocess_only = False

# ===================== input files and path settings =======================

gds_filename = "transformer.gds"   # geometries
XML_filename = "SG13G2.xml"               # stackup

# preprocess GDSII for safe handling of cutouts/holes?
preprocess_gds = False

# merge via polygons with distance less than .. microns, set to 0 to disable via merging.
merge_polygon_size = 1.5

# get path for this simulation file
script_path = utilities.get_script_path(__file__)

# use script filename as model basename
model_basename = utilities.get_basename(__file__)

# set and create directory for simulation output
sim_path = utilities.create_sim_path (script_path,model_basename)
print('Simulation data directory: ', sim_path)


# ======================== simulation settings ================================

unit   = 1e-6  # geometry is in microns
margin = 80    # distance in microns from GDSII geometry boundary to simulation boundary 

fstart =  20e9
fstop  = 100e9
numfreq = 401

refined_cellsize = 0.7  # mesh cell size in conductor region

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
simulation_ports.add_port(simulation_setup.simulation_port(portnumber=1, voltage=1, port_Z0=50,  source_layernum=201, target_layername='TopMetal2', direction='x'))
simulation_ports.add_port(simulation_setup.simulation_port(portnumber=2, voltage=1, port_Z0=200, source_layernum=202, target_layername='TopMetal2', direction='x'))

# ======================== simulation ================================

# get technology stackup data
materials_list, dielectrics_list, metals_list = stackup_reader.read_substrate (XML_filename)
# get list of layers from technology
layernumbers = metals_list.getlayernumbers()
layernumbers.extend(simulation_ports.portlayers)

# read geometries from GDSII, only purpose 0, merge vias with distance < 1.5 um
allpolygons = gds_reader.read_gds(gds_filename, layernumbers, purposelist=[0], metals_list=metals_list, preprocess=preprocess_gds, merge_polygon_size=merge_polygon_size)


# calculate maximum cellsize from wavelength in diecletric
wavelength_air = 3e8/fstop
max_cellsize = (wavelength_air/unit)/(sqrt(materials_list.eps_max)*cells_per_wavelength) 


########### create model, run and post-process ###########

# Create simulation for port 1 excitation
excite_ports = [1]  # list of ports that are excited for this simulation run

FDTD = openEMS(EndCriteria=exp(energy_limit/10 * log(10)))
FDTD.SetGaussExcite( (fstart+fstop)/2, (fstop-fstart)/2 )
FDTD.SetBoundaryCond( Boundaries )
FDTD = simulation_setup.setupSimulation (excite_ports, simulation_ports, FDTD, materials_list, dielectrics_list, metals_list, allpolygons, max_cellsize, refined_cellsize, margin, unit)
sub1_data_path = simulation_setup.runSimulation (excite_ports, FDTD, sim_path, model_basename, preview_only, postprocess_only)



########## evaluation of results with composite GSG ports ###########

if not preview_only:

    # get results, CSX port definition is read from simulation ports object
    CSX_port1 = simulation_ports.get_port_by_number(1).CSXport
    CSX_port2 = simulation_ports.get_port_by_number(2).CSXport
    
    f = np.linspace(fstart,fstop,numfreq)
    CSX_port1.CalcPort(sub1_data_path, f, simulation_ports.get_port_by_number(1).port_Z0)
    CSX_port2.CalcPort(sub1_data_path, f, simulation_ports.get_port_by_number(2).port_Z0)

    # note the mixed impedances, we need to adjust this in S-param below!
    s11 = CSX_port1.uf_ref  / CSX_port1.uf_inc
    s21 = (CSX_port2.uf_ref  / CSX_port1.uf_inc)  * sqrt((simulation_ports.get_port_by_number(1).port_Z0/simulation_ports.get_port_by_number(2).port_Z0))
    Zin = 0.5 * CSX_port1.uf_tot  / CSX_port1.if_tot 

    s11_dB = 20.0*np.log10(np.abs(s11))
    s21_dB = 20.0*np.log10(np.abs(s21))

    figure()
    plot(f/1e9, s11_dB, 'k-', linewidth=2, label='S11 [dB]')
    plot(f/1e9, s21_dB, 'r-', linewidth=2, label='S11 [dB]')
    grid()
    legend()
    xlabel('Frequency (GHz)')

    figure()
    plot(f/1e9, real(Zin), 'k-', linewidth=2, label='Real(Zin)')
    grid()
    legend()
    xlabel('Frequency (GHz)')

    u1 = CSX_port1.u_data.ui_val[0]
    u2 = CSX_port2.u_data.ui_val[0]
    t  = CSX_port1.u_data.ui_time[0]

    figure()
    plot(t*1e9,u1, 'r-', label='u1')
    plot(t*1e9,u2, 'g-', label='u2')
    grid()
    ylabel('Port voltages (V)')
    xlabel('Time (ns)')
    legend()


    # Show all plots
    show()

