
import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from deaerator import Deaerator
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
    AmountBasis,
    )

from pyomo.environ import ConcreteModel, SolverFactory, SolverStatus, TerminationCondition, Block, TransformationFactory, units, Objective, value, Constraint, Var, Param

# This defines the base case for all tests
m = pyo.ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.water = HelmholtzParameterBlock(
                    pure_component="h2o",
                    phase_presentation=PhaseType.LG,
                    state_vars=StateVars.PH,
                    amount_basis=AmountBasis.MASS,
                    )
# Create deaerator with 1 inlets
m.fs.deaerator = Deaerator(property_package=m.fs.water, num_inlets=1)

# Set water inlets
T1_in = 85 # C
P1_in = 1.0 # bar g
m.fs.deaerator.inlet_1.flow_mass[0].fix(80/3.6)  # Mass flow rate in kg/s
m.fs.deaerator.inlet_1.enth_mass[0].fix(value(m.fs.water.htpx(T=(T1_in+273)*units.K, p=(P1_in+1)*units.bar)))
m.fs.deaerator.inlet_1.pressure[0].fix((P1_in+1)*units.bar) 

# set steam inlet
T1_in = 160 # C
P1_in = 4.5 # bar g
m.fs.deaerator.inlet_steam.flow_mass[0].fix(5/3.6)  # Mass flow rate in kg/s
m.fs.deaerator.inlet_steam.enth_mass[0].fix(value(m.fs.water.htpx(T=(T1_in+273)*units.K, p=(P1_in+1)*units.bar)))
m.fs.deaerator.inlet_steam.pressure[0].fix((P1_in+1)*units.bar) 

# # Set pump properties - need to move these to Deaerator config so thats its m.fs.deaerator.deaerator_pressure.fix()
m.fs.deaerator.transfer_pump.outlet.pressure[0].fix(10e5)  # Fix pump outlet pressure
m.fs.deaerator.transfer_pump.efficiency_pump[0].fix(0.8)  # Fix pump outlet pressure

# # use flash to flash down so that 0.2% of the the total flow is vapor which is then split
# Target vent fraction (0.2% of total inlet molar flow)
m.fs.vent_frac = Param(initialize=0.002, mutable=True)
m.fs.deaerator.flash.heat_duty.fix(0)
m.fs.deaerator.flash.control_volume.deltaP.fix(-2*units.bar)

# SPEC: vapor outlet molar flow = vent_frac * inlet molar flow
# m.fs.vent_spec = Constraint(expr=
#     m.fs.deaerator.flash.vap_outlet.flow_mass[0] == m.fs.vent_frac * m.fs.deaerator.flash.inlet.flow_mass[0]
# )

# m.fs.deaerator.flash.control_volume.deltaP.fix(-7*units.bar)
# m.fs.deaerator.flash.control_volume.deltaP.unfix()
# probably needs initialisation here

#Optional) sensible bounds/initial guess on pressure drop (negative = drop)
# m.fs.deaerator.flash.control_volume.deltaP.setlb(-2e6*units.Pa)
# m.fs.deaerator.flash.control_volume.deltaP.setub(-1e3*units.Pa)
# m.fs.deaerator.flash.control_volume.deltaP.value = -2e5*units.Pa


# add stream table

m.fs.deaerator.initialize()
solver = SolverFactory("ipopt")
solver.options = {"tol": 1e-3, "max_iter": 5000}
print(degrees_of_freedom(m.fs))
result = solver.solve(m, tee=True)

m.fs.deaerator.report()
m.fs.deaerator.mixer.report()
m.fs.deaerator.transfer_pump.report()
m.fs.deaerator.flash.report()

m.fs.deaerator.flash.vap_outlet.enth_mass.display() 
assert result.solver.termination_condition == TerminationCondition.optimal


