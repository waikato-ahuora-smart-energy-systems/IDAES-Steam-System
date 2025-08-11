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
    m.fs1.turbine = TurbineBase(property_package=m.fs1.water, calculation_method="simple_willans")
    
    

def set_inputs(m):
    m.fs1.turbine.inlet.flow_mass.fix(2/3.6) # why not [0] on mass
    m.fs1.turbine.inlet.enth_mass[0].fix(value(m.fs1.water.htpx(T=(400+273)*units.K, p=45*units.bar)))
    m.fs1.turbine.inlet.pressure[0].fix(45*units.bar) # why the [0]
    m.fs1.turbine.outlet.pressure[0].fix(12.5*units.bar)
    #m.fs1.turbine.efficiency_isentropic.fix(0.75)
    m.fs1.turbine.willans_slope.fix(3.24*1000*units.J/units.mol)  # Willans slope 
    m.fs1.turbine.willans_intercept.fix(500*1000*units.W)  # Willans intercept
    m.fs1.turbine.willans_max_mol.fix(100*15.4*units.mol / units.s)  # Willans intercept



m = ConcreteModel()
build_model(m)  # build flowsheet
set_inputs(m)
solver = SolverFactory("ipopt")
print(degrees_of_freedom(m.fs1))
result = solver.solve(m, tee=False)


m.fs1.turbine.outlet.pressure.display()
m.fs1.turbine.report()

''' 
#QA porosity process at axiam 4
#testing cage - chain length fence? 2 how can we set up quick
#aurel from Rio Tinto props incl 1
#glue options & testing for foam-flange 3
how to use 2" elf sensors

all done some work within 1-2 weeks 
'''


