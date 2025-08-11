# Import Pyomo libraries
from pyomo.environ import ConcreteModel, SolverFactory, SolverStatus, TerminationCondition, Block, TransformationFactory, units, Objective, value, Constraint, Var, maximize
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
from idaes.models.unit_models.separator import SplittingType
from idaes.models.unit_models import Compressor, PressureChanger
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    PhaseType,
    StateVars,
    AmountBasis,
    )
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption, Turbine
from .turbine_base_model import TurbineBase


def build_model(m):
    # Define model components and blocks
    m.fs1 = FlowsheetBlock(dynamic=False)

    # Attach property package to flowsheet
    m.fs1.water = HelmholtzParameterBlock(
                    pure_component="h2o",
                    phase_presentation=PhaseType.LG,
                    state_vars=StateVars.PH,
                    amount_basis=AmountBasis.MASS,
                    )
    
    ## HP stage
    m.fs1.HP_stage = TurbineBase(property_package=m.fs1.water)
    
    # Add a splitter after the HP stage
    m.fs1.MP_splitter = Splitter(
        property_package=m.fs1.water,
        #MomentumMixingType=MomentumMixingType.none,  # No momentum mixing
        outlet_list=["MP_passout", "MP_next_stage"],
        split_basis =SplittingType.totalFlow,  # Split based on mass flow
    )

    m.fs1.MP_header_splitter = Splitter(
        property_package=m.fs1.water,
        #MomentumMixingType=MomentumMixingType.none,  # No momentum mixing
        outlet_list=["MP_demand", "MP_to_letdown"],
        split_basis =SplittingType.totalFlow,  # Split based on mass flow
    )


    ## LP stage
    m.fs1.LP_stage = TurbineBase(property_package=m.fs1.water)


    # Connect HP stage outlet to splitter inlet
    m.fs1.HP_passout = Arc(
        source=m.fs1.HP_stage.outlet,
        destination=m.fs1.MP_splitter.inlet
    )


    # Connect MP passout spliiter to header
    m.fs1.MP_passout = Arc(
        source=m.fs1.MP_splitter.MP_passout,
        destination=m.fs1.MP_header_splitter.inlet,
    )

    # Connect MP spliiter outlet to LP stage inlet
    m.fs1.LP_stageing = Arc(
        source=m.fs1.MP_splitter.MP_next_stage,
        destination=m.fs1.LP_stage.inlet
    )


    # Expand arcs
    TransformationFactory("network.expand_arcs").apply_to(m)
    
    
def set_inputs(m, params):
    # Unpack params
    HP_inlet_flow = params['HP_inlet_flow'] * (1000 / 3600)  # Convert t/h to kg/s
    LP_passout_limit = params['LP_passout_limit'] * (1000 / 3600)  # Convert t/h to kg/s
    MP_demand_flow = params['MP_demand_flow'] * (1000 / 3600)  # Convert t/h to kg/s
    LP_demand_flow = params['LP_demand_flow'] * (1000 / 3600)  # Convert t/h to kg/s
    HP_pressure = params['HP_pressure'] * units.bar
    MP_pressure = params['MP_pressure'] * units.bar
    LP_pressure = params['LP_pressure'] * units.bar
    HP_temp = (params['HP_temperature'] + 273.15) * units.K  # Convert C to K

    # HP turbine
    m.fs1.HP_stage.inlet.flow_mass.fix(HP_inlet_flow) # why not [0] on mass
    m.fs1.HP_stage.inlet.enth_mass[0].fix(value(m.fs1.water.htpx(T=HP_temp, p=HP_pressure)))
    m.fs1.HP_stage.inlet.pressure[0].fix(HP_pressure) # why the [0]
    m.fs1.HP_stage.outlet.pressure[0].fix(MP_pressure)
    m.fs1.HP_stage.efficiency_isentropic.fix(0.75)

    # MP splitter
    #m.fs1.MP_splitter.split_fraction[0, 'MP_next_stage'].fix(LP_passout_limit/HP_inlet_flow)
    #m.fs1.MP_header_splitter.split_fraction[0, 'MP_demand'].fix(MP_demand_flow/m.fs1.MP_header_splitter.inlet.flow_mass[0])  
    
    
    m.fs1.cons1 = Constraint(expr=(m.fs1.MP_splitter.MP_next_stage.flow_mass[0] <=  LP_passout_limit))
    #m.fs1.cons2 = Constraint(expr=(m.fs1.MP_splitter.MP_next_stage.flow_mass[0] + m.fs1.MP_header_splitter.MP_to_letdown.flow_mass[0] == LP_demand_flow))
    m.fs1.cons3 = Constraint(expr=(m.fs1.MP_header_splitter.MP_demand.flow_mass[0] == MP_demand_flow))
    m.fs1.objfn = Objective(expr=(m.fs1.LP_stage.work_mechanical[0]))

    # LP stage
    m.fs1.LP_stage.outlet.pressure[0].fix(LP_pressure)
    m.fs1.LP_stage.efficiency_isentropic.fix(0.65)

    
def initialise(m):
    # Initialize the HP stage
    m.fs1.HP_stage.initialize()
    m.fs1.MP_splitter.initialize()
    m.fs1.MP_header_splitter.initialize()

    # Initialize the LP stage
    #m.fs1.LP_stage.initialize()


def report(m):
    m.fs1.HP_stage.report()
    m.fs1.MP_splitter.report()
    m.fs1.MP_header_splitter.report()
    m.fs1.LP_stage.report()

def series_tubine(m, params):
    solver = SolverFactory("ipopt")
    solver.options = {"tol": 1e-3, "max_iter": 1000}


    build_model(m)  # build flowsheet
    set_inputs(m, params)
    initialise(m)  # initialize model

    print("DOF before solve: ", degrees_of_freedom(m))
    
    result = solver.solve(m, tee=False)
    report(m)
   
    assert result.solver.termination_condition == TerminationCondition.optimal
    #m.fs1.visualize('flowsheet1', loop_forever=True)



''' 

'''


