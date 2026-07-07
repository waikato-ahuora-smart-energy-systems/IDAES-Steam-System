from pyomo.environ import Suffix, Var
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.network import Arc, SequentialDecomposition
from pyomo.core.base.reference import Reference
import pyomo.environ as pyo

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.math import smooth_max
from idaes.models.unit_models.separator import SplittingType, EnergySplittingType, Separator
from idaes.models.unit_models.mixer import Mixer
from idaes.models.unit_models.heater import Heater
from scripts.deaerator.custom_separator import CustomSeparator
from simple_separator import SimpleSeparator

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("SteamHeader")
class SteamHeaderData(UnitModelBlockData):
    """
    Steam Header unit operation:
    Mixer -> Cooler -> Phase Separator -> Splitter
    Separates 100% liquid to liquid_outlet and 100% vapor to splitter outlets.
    Uses Sequential Decomposition for optimal initialization order.
    """

    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
        ),
    )
    CONFIG.declare(
        "num_inlets",
        ConfigValue(
            default=2,
            domain=int,
            description="Number of inlets to add" \
            "Index [-1]: Steam makeup",
        ),
    )
    CONFIG.declare(
        "num_outlets",
        ConfigValue(
            default=2,
            domain=int,
            description="Number of outlets to add" \
            "Index [-1]: Let-down" \
            "Index [-2]: Condensate",
        ),
    )
    
    def build(self):
        super().build()
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # Create internal units
        self.mixer = Mixer(
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            num_inlets=self.config.num_inlets
        )
        
        self.cooler = Heater(
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            has_pressure_change=True,
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup
        )
        
        self.phase_separator = CustomSeparator(
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            num_outlets=2,
            split_basis=SplittingType.phaseFlow,
            energy_split_basis=EnergySplittingType.enthalpy_split
        )
        self.phase_separator.split_fraction[:, "outlet_1", "Vap"].fix(1.0)
        self.phase_separator.split_fraction[:, "outlet_1", "Liq"].fix(0.0)

        self.splitter = SimpleSeparator(
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            num_outlets=self.config.num_outlets
        )
        
        # Updated internal arcs
        self.mixer_to_cooler_arc = Arc(
            source=self.mixer.outlet,
            destination=self.cooler.inlet
        )
        
        self.cooler_to_separator_arc = Arc(
            source=self.cooler.outlet,
            destination=self.phase_separator.inlet
        )
        
        self.separator_to_splitter_arc = Arc(
            source=self.phase_separator.outlet_1,
            destination=self.splitter.inlet
        )

        # Create external ports and references to map them to internal units
        # Prepare tracking lists
        self.inlet_list = []
        self.outlet_list = []

        self.inlet_blocks = []
        self.outlet_blocks = []
        
        # Add external inlet ports (connected to mixer inlets)
        for i in range(1, self.config.num_inlets + 1):
            # Get the internal state block from the mixer
            inlet_name = f"inlet_{i}"
            mixer_inlet_i_state = getattr(self.mixer, f"inlet_{i}_state")
            
            # Create port referencing the mixer state block
            self.add_port(name=inlet_name, block=mixer_inlet_i_state)
            
            # Create reference to easily access variables 
            # This allows accessing like: header.inlet_1_vars[time].flow_mol
            setattr(self, f"{inlet_name}_vars", Reference(mixer_inlet_i_state[:]))
            
            # Add to tracking lists
            self.inlet_list.append(inlet_name)
            self.inlet_blocks.append(mixer_inlet_i_state)
        
        # Add external vapor outlet ports (connected to splitter outlets)
        for i in range(1, self.config.num_outlets + 1):
            # Get the internal state block
            outlet_name = f"outlet_{i}"
            splitter_outlet_i_state = getattr(self.splitter, f"outlet_{i}_state")
            
            # Create port referencing the splitter state block
            self.add_port(name=outlet_name, block=splitter_outlet_i_state)
            
            # Create reference for easy access
            setattr(self, f"{outlet_name}_vars", Reference(splitter_outlet_i_state[:]))
            
            # Add to tracking lists
            self.outlet_list.append(outlet_name)
            self.outlet_blocks.append(splitter_outlet_i_state)
        
        # Add liquid outlet from phase separator
        outlet_name = "liquid_outlet"
        liquid_outlet_state = getattr(self.phase_separator, "outlet_2_state")
        self.add_port(name=outlet_name, block=liquid_outlet_state)
        self.liquid_outlet_vars = Reference(liquid_outlet_state[:])
        
        setattr(self, f"{outlet_name}_vars", Reference(self.liquid_outlet_vars))

        # Add to tracking lists
        self.outlet_list.append(outlet_name)
        self.outlet_blocks.append(liquid_outlet_state)

        self.balance_flow_mol = Var(self.flowsheet().time, initialize=0.0, doc="Balance molar flow")
        self.partial_inlet_flow_mol = Var(self.flowsheet().time, initialize=0.0, doc="Partial total inlet molar flow")

        # Additional bounds and constraints
        self._enforce_nonnegative_port_flows(self.inlet_blocks, self.outlet_blocks)
        self._add_overall_material_balance(self.inlet_blocks, self.outlet_blocks, self.balance_flow_mol)
        self._assign_balance_flow(self.balance_flow_mol, self.inlet_blocks[-1], self.outlet_blocks[-2], self.partial_inlet_flow_mol, self.inlet_blocks)
        
        # Expand arcs
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)

        self.split_fraction = Reference(self.splitter.split_flow)

    @property
    def heat_duty(self):
        return self.cooler.heat_duty
        
    @property
    def deltaP(self):
        return self.cooler.deltaP
    
    @property
    def total_flow_mol(self):
        return Reference(self.cooler.control_volume.properties_out[:].flow_mol)

    @property
    def total_flow_mass(self):
        return Reference(self.cooler.control_volume.properties_out[:].flow_mass)
    
    @property
    def pressure(self):
        return Reference(self.cooler.control_volume.properties_out[:].pressure)

    @property
    def temperature(self):
        return Reference(self.cooler.control_volume.properties_out[:].temperature)

    @property
    def enth_mol(self):
        return Reference(self.cooler.control_volume.properties_out[:].enth_mol)  

    @property
    def vapor_frac(self):
        return Reference(self.cooler.control_volume.properties_out[:].vapor_frac)            

    def _enforce_nonnegative_port_flows(self, inlet_blocks, outlet_blocks):
        """
        Set lower bounds on flow variables for all external ports.
        """
        [state_block[t].flow_mol.setlb(0.0)
            for state_block in (inlet_blocks + outlet_blocks)
            for t in state_block]

    def _add_overall_material_balance(self, inlet_blocks, outlet_blocks, balance_flow_mol):
        """
        Add overall material balance equation.
        """
        # Get phase component list(s)
        pc_set = outlet_blocks[0].phase_component_set

        # Write phase-component balances
        @self.Constraint(self.flowsheet().time, doc="Material balance equation")
        def material_balance_equation(b, t):
            return 0 == sum(
                sum(
                    sum(
                        o[t].get_material_flow_terms(p, j)
                        for o in inlet_blocks[:-1]
                    )                    
                    - 
                    sum(
                        o[t].get_material_flow_terms(p, j)
                        for o in outlet_blocks[:-2] + outlet_blocks[-1:]
                    )                    
                    for j in outlet_blocks[0].component_list
                    if (p, j) in pc_set
                )
                for p in outlet_blocks[0].phase_list
            ) - balance_flow_mol[t]

    def _assign_balance_flow(self, balance_flow_mol, makeup, letdown, partial_inlet_flow_mol, inlet_blocks):
        pc_set = inlet_blocks[0].phase_component_set
        @self.Constraint(self.flowsheet().time, doc="Material balance of known inlet flows. " \
                                                    "Used for scaling the steam makeup and letdown flow balance.")
        def partial_material_balance_equation(b, t):
            return 0 == sum(
                sum(
                    sum(
                        o[t].get_material_flow_terms(p, j)
                        for o in inlet_blocks[:-1]
                    )                    
              
                    for j in inlet_blocks[0].component_list
                    if (p, j) in pc_set
                )
                for p in inlet_blocks[0].phase_list
            ) - partial_inlet_flow_mol[t]
        
        eps = 1e-5 # smoothing parameter; smaller = closer to exact max, larger = smoother
        @self.Constraint(self.flowsheet().time, doc="Steam letdown flow balance.")
        def letdown_flow_balance(b, t):
            return 0 == (
                smooth_max( balance_flow_mol[t] / (partial_inlet_flow_mol[t] + 1e-6), 0.0, eps) * (partial_inlet_flow_mol[t] + 1e-6)
                - letdown[t].flow_mol
            )         

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        self.mixer.calculate_scaling_factors()
        self.cooler.calculate_scaling_factors()
        self.phase_separator.calculate_scaling_factors()
        self.splitter.calculate_scaling_factors()

    def initialize_build(self, outlvl=idaeslog.NOTSET, **kwargs):
        """
        Initialize the Header unit using Sequential Decomposition to determine optimal order
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        
        init_log.info("Starting CustomHeader initialization using Sequential Decomposition")
        
        # create Sequential Decomposition object
        seq = SequentialDecomposition()
        seq.options.select_tear_method = "heuristic"
        seq.options.tear_method = "Wegstein"
        seq.options.iterLim = 1

        # create computation graph
        G = seq.create_graph(self)
        heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
        # get calculation order
        order = seq.calculation_order(G)

        for o in heuristic_tear_set:
            print(o.name)
        
        for o in order:
            print(o[0].name)
        
        # define unit initialisation function
        def init_unit(unit):
            unit.initialize(outlvl=outlvl, **kwargs)
        
        # run sequential decomposition
        seq.run(self, init_unit)

    def _get_stream_table_contents(self, time_point=0):
        """
        Create stream table showing all inlets, outlets, and liquid outlet
        """
        io_dict = {}
        
        for inlet_name in self.inlet_list:
            io_dict[inlet_name] = getattr(self, inlet_name)
        
        for outlet_name in self.outlet_list:
            io_dict[outlet_name] = getattr(self, outlet_name)
            
        io_dict["liquid_outlet"] = self.liquid_outlet
            
        return create_stream_table_dataframe(io_dict, time_point=time_point)