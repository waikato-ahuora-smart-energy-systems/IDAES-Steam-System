# Import Pyomo libraries
from pyomo.environ import ConcreteModel, SolverFactory, SolverStatus, TerminationCondition, Block, TransformationFactory, units, Objective, value, Constraint, Var
from pyomo.network import SequentialDecomposition, Port, Arc


# Import IDAES libraries
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import report_statistics, degrees_of_freedom
import idaes.logger as idaeslog
from idaes.core.util.tables import _get_state_from_port  

# Import required models
from idaes.models.unit_models import (
    Feed,
    Mixer,
    Heater,
    Compressor,
    Product,
    MomentumMixingType,
)
from idaes.models.unit_models import Separator as Splitter
from idaes.models.unit_models import Compressor, PressureChanger
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
    AmountBasis,
    )
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption, Turbine
from turbine_base_model import TurbineBase
from idaes.core.util import DiagnosticsToolbox


def build_model(m):
    # Define model components and blocks
    m.fs1 = FlowsheetBlock(dynamic=False)
    m.fs1.water = HelmholtzParameterBlock(
                    pure_component="h2o",
                    phase_presentation=PhaseType.LG,
                    state_vars=StateVars.PH,
                    amount_basis=AmountBasis.MASS,
                    )
    
    #m.fs1.turbine = Turbine(property_package=m.fs1.water)
    calculation_method = "part_load_willans" #"part_load_willans" # "isentropic"  # or "simple willans" "Tsat_willans"
    m.fs1.turbine = TurbineBase(property_package=m.fs1.water, calculation_method=calculation_method)
    
    

def set_inputs(m):
    m_in = 187.0 # t/h
    P_in = 41.3 # bar g
    P_out = 10.4 # bar g
    T_in = 381 # C

    m.fs1.turbine.inlet.flow_mass.fix(m_in/3.6) # why not [0] on mass
    m.fs1.turbine.inlet.enth_mass[0].fix(value(m.fs1.water.htpx(T=(T_in+273)*units.K, p=(P_in+1)*units.bar)))
    m.fs1.turbine.inlet.pressure[0].fix((P_in+1)*units.bar) # why the [0]
    m.fs1.turbine.outlet.pressure[0].fix((P_out+1)*units.bar)
    
    if  m.fs1.turbine.config.calculation_method == "isentropic":
        m.fs1.turbine.efficiency_isentropic.fix(0.75)
    
    elif  m.fs1.turbine.config.calculation_method == "simple willans":
        m.fs1.turbine.willans_slope.fix(3.24*1000*units.J/units.mol)  # Willans slope 
        m.fs1.turbine.willans_intercept.fix(500*1000*units.W)  # Willans intercept
        m.fs1.turbine.willans_max_mol.fix(100*15.4*units.mol / units.s)  # Willans intercept
    
    elif m.fs1.turbine.config.calculation_method == "part_load_willans":
        m.fs1.turbine.efficiency_motor.fix(1.0)  
        m.fs1.turbine.willans_max_mol.fix(217.4*15.4)  # Willans intercept
        m.fs1.turbine.willans_a.fix(1.5435)  # Willans slope 
        m.fs1.turbine.willans_b.fix(0.2*units.kW)  # Willans intercept
        m.fs1.turbine.willans_c.fix(0.3759)  # Willans intercept
    
    elif m.fs1.turbine.config.calculation_method == "Tsat_willans":
        m.fs1.turbine.efficiency_motor.fix(1.0) 
        m.fs1.turbine.willans_max_mol.fix(217.4*15.4)  # Willans intercept
        # m.fs1.turbine.willans_a.fix(1.5435)  # Willans slope 
        # m.fs1.turbine.willans_b.fix(0.2*units.kW)  # Willans intercept
        # m.fs1.turbine.willans_c.fix(0.3759)  # Willans intercept
       
       



m = ConcreteModel()
build_model(m)  # build flowsheet
set_inputs(m)

solver = SolverFactory("ipopt")
solver.options = {"tol": 1e-3, "max_iter": 5000}
print(degrees_of_freedom(m.fs1))
result = solver.solve(m, tee=False)



# dt = DiagnosticsToolbox(m)
# dt.report_structural_issues()
# dt.display_underconstrained_set()
# dt.display_overconstrained_set()
# dt.display_components_with_inconsistent_units()

from pyomo.util.check_units import assert_units_consistent
assert_units_consistent(m.fs1)

m.fs1.turbine.report()
m.fs1.turbine.control_volume.properties_in[0].temperature_sat.display()
m.fs1.turbine.control_volume.properties_out[0].temperature_sat.display()
print(m.fs1.turbine.willans_a[0].value)
print(m.fs1.turbine.willans_b[0].value)
print(m.fs1.turbine.willans_c[0].value)
print(m.fs1.turbine.willans_slope[0].value)
print(m.fs1.turbine.willans_intercept[0].value)

'''
m.fs1.turbine.willans_slope.display()
m.fs1.turbine.willans_intercept.display()
m.fs1.turbine.control_volume.properties_in[0].enth_mol.display()
m.fs1.turbine.properties_isentropic[0].enth_mol.display()
 '''




