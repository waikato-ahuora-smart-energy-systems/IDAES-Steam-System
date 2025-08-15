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
from splitter_new_model import Separator_new
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
    # Define unit ops
    #m.fs1.turbine = Turbine(property_package=m.fs1.water)
    calculation_method = "CT_willans" #"part_load_willans" # "isentropic"  # or "simple_willans" "Tsat_willans" " BPST_willans" "CT_willans"
    m.fs1.splitter = Separator_new(property_package=m.fs1.water, outlet_list=["outlet1", "outlet2"])
    m.fs1.mixer = Mixer(
        property_package=m.fs1.water,
        inlet_list=["supply1", "supply2"],
    )

    m.fs1.stream1 = Arc(
        source=m.fs1.mixer.outlet,
        destination=m.fs1.splitter.inlet,
    )

    # Expand arcs
    TransformationFactory("network.expand_arcs").apply_to(m)


    
def set_inputs(m):
    m_supply_1 = 1.81 * 3.6 # t/h
    T_supply_1 = 250 # C
    P_supply_1 = 7 # bar a

    m_supply_2 = 1.81 * 3.6 # t/h
    P_supply_2 = 8 # bar a
    X_supply_2 = 0.8 # mass fraction of vapor in steam

    # Mixer
    m.fs1.mixer.supply1.flow_mass[0].fix(m_supply_1 / 3.6)  # Convert t/h to kg/s
    m.fs1.mixer.supply1.enth_mass[0].fix(value(m.fs1.water.htpx(T=(T_supply_1 + 273) * units.K, p=(P_supply_1) * units.bar)))
    m.fs1.mixer.supply1.pressure[0].fix((P_supply_1) * units.bar)

    m.fs1.mixer.supply2.flow_mass[0].fix(m_supply_2 / 3.6)  # Convert t/h to kg/s
    m.fs1.mixer.supply2.enth_mass[0].fix(value(m.fs1.water.htpx(x=X_supply_2, p=(P_supply_2) * units.bar)))
    m.fs1.mixer.supply2.pressure[0].fix((P_supply_2) * units.bar)

    # Splitter
    m.fs1.splitter.split_fraction[0, 'outlet1'].fix(0.5)

def initialise(m):
    # Initialize the ops
    m.fs1.mixer.initialize()
    m.fs1.splitter.initialize()
   
   
m = ConcreteModel()
build_model(m)  # build flowsheet
set_inputs(m)
initialise(m)

solver = SolverFactory("ipopt")
solver.options = {"tol": 1e-3, "max_iter": 5000}
print(degrees_of_freedom(m.fs1))
result = solver.solve(m)

# dt = DiagnosticsToolbox(m)
# dt.report_structural_issues()
# dt.display_underconstrained_set()
# dt.display_overconstrained_set()
# dt.display_components_with_inconsistent_units()

# from pyomo.util.check_units import assert_units_consistent
# assert_units_consistent(m.fs1)

m.fs1.mixer.report()
m.fs1.splitter.report()
assert result.solver.termination_condition == TerminationCondition.optimal