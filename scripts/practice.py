
from pyomo.environ import ConcreteModel, SolverFactory, SolverStatus, TerminationCondition, Block, TransformationFactory, units, Objective, value
from pyomo.network import SequentialDecomposition, Port, Arc
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import report_statistics, degrees_of_freedom
from idaes.models.unit_models.heater import Heater
from idaes.models.unit_models import Compressor, PressureChanger
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
    AmountBasis,
)
import idaes.logger as idaeslog
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption


### Build Model
m = ConcreteModel()
m.fs1 = FlowsheetBlock(dynamic= False)

# Property Package (need to specificy property package for each component in each flowsheet)
m.fs1.water = HelmholtzParameterBlock(
  pure_component="h2o",
  phase_presentation=PhaseType.LG,
  state_vars=StateVars.PH,
  amount_basis=AmountBasis.MASS,
)

# Unit ops 
m.fs1.heater1 = Heater(property_package= m.fs1.water, has_pressure_change=False) # but if no pressure change then leave out as it defaults to 0
m.fs1.comp1 = Compressor(property_package= m.fs1.water)
m.fs1.expv1 = PressureChanger(property_package= m.fs1.water, compressor=False, thermodynamic_assumption=ThermodynamicAssumption.adiabatic)

# Arcs (these can be more than streams but we mostly use them for streams)
m.fs1.stream1 = Arc(source=m.fs1.heater1.outlet, destination=m.fs1.comp1.inlet)
m.fs1.stream2 = Arc(source=m.fs1.comp1.outlet, destination=m.fs1.expv1.inlet)
TransformationFactory("network.expand_arcs").apply_to(m) # this is telling Pyomo to treat the arcs as a connection/stream between the two unit ops
print(degrees_of_freedom(m.fs1))

#fIX vARIABLES
m.fs1.heater1.inlet.flow_mass.fix(1)
m.fs1.heater1.inlet.pressure.fix(101.3e3)
m.fs1.heater1.inlet.enth_mass.fix(m.fs1.water.htpx(T=300*units.K, p=101.3*units.kPa)) #m.fs1.heater1.inlet.enthalpy.fix(m.fs1.water.htpx(T=200-273.15, P=101.3e3))
#m.fs1.heater1.outlet.enth_mass.fix(m.fs1.water.htpx(x=0.5, p=101.3*units.kPa))
m.fs1.heater1.heat_duty.fix(2200e3)

# Compressor
m.fs1.comp1.outlet.pressure.fix(200e3)
m.fs1.comp1.efficiency_isentropic.fix(0.8)

# Expander
m.fs1.expv1.outlet.pressure.fix(150e3)

print(degrees_of_freedom(m.fs1))


#m.fs1.heater1.initialize(outlvl=idaeslog.INFO) # can specifiy initialisation method but by default it has one, this line is just showing the output log of the default initialisation
solver = SolverFactory("ipopt")
result = solver.solve(m, tee=True)

# Show stream data
m.fs1.heater1.report()
m.fs1.comp1.report()
m.fs1.expv1.report()
print('heater outlet temp',value(m.fs1.heater1.control_volume.properties_out[0].temperature))

m.fs1.visualize('flowsheet1', loop_forever=True)





'''
# 'Optimise' the heat duty to attain a desired quality of steam. (this is more like an excel goalseek or solver than)
m.fs1.objfn = Objective(expr=(m.fs1.water.htpx(x=0.5, p=101.3*units.kPa) - m.fs1.heater1.outlet.enth_mass[0])**2) 
m.fs1.heater1.heat_duty.unfix()
m.fs1.heater1.heat_duty.setlb(0)
m.fs1.heater1.heat_duty.setub(1000e5)
result = solver.solve(m, tee=True)

m.fs1.heater1.report()
print(m.fs1.objfn())
'''
