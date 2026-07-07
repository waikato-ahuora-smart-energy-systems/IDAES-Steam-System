from pyomo.environ import Suffix, Var, Param
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
from idaes.models.unit_models import Pump, Flash

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("Deaerator")
class DeaeratorData(UnitModelBlockData):
    """
    Deaerator unit operation:
    Mixer -> Pump -> Flash ->
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
            "Steam line is automatically created seperate to this init",
        ),
    )
   
   
    def build(self):
        super().build()
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        # User parameters
    
        # Create internal units
        self.transfer_pump = Pump(
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.mixer = Mixer(
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            num_inlets=self.config.num_inlets + 1 # additonal inlet for steam line
        )

        self.flash = Flash(
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            has_heat_transfer=True,
            has_pressure_change=True,
            ideal_separation=False,        
            energy_split_basis=EnergySplittingType.equal_temperature,
        )
        
        # Updated internal arcs
        self.mixer_to_pump_arc = Arc(
            source=self.mixer.outlet,
            destination=self.transfer_pump.inlet
        )

        self.pump_to_flash_arc = Arc(
            source=self.transfer_pump.outlet,
            destination=self.flash.inlet
        )

        #self.transfer_pump.outlet.pressure.fix(self.deaerator_pressure)

        # Expand arcs
        pyo.TransformationFactory("network.expand_arcs").apply_to(self)


        # Create external ports and references to map them to internal units
        # Prepare tracking lists
        self.inlet_list = []
        self.outlet_list = []

        self.inlet_blocks = []
        self.outlet_blocks = []

        # Add external inlet ports (connected to mixer inlets)
        for i in range(1, self.config.num_inlets + 2):
            # Get the internal state block from the mixer
            inlet_name = f"inlet_steam" if i == self.config.num_inlets + 1 else f"inlet_{i}"
            mixer_inlet_i_state = getattr(self.mixer, f"inlet_{i}_state")
            
            # Create port referencing the mixer state block
            self.add_port(name=inlet_name, block=mixer_inlet_i_state)
            
            # Create reference to easily access variables 
            # This allows accessing like: header.inlet_1_vars[time].flow_mol
            setattr(self, f"{inlet_name}_vars", Reference(mixer_inlet_i_state[:]))
            
            # Add to tracking lists
            self.inlet_list.append(inlet_name)
            self.inlet_blocks.append(mixer_inlet_i_state)

        
      
        ''' '''
        # Add external outlet ports (connected to pump outlets)
        for i in range(1, len(self.transfer_pump.outlet) + 1):
            # Get the internal state block
            outlet_name = f"outlet"
            pump_outlet_i_state = getattr(self.transfer_pump.control_volume, f"properties_out")
            
            # Create port referencing the pump state block
            self.add_port(name=outlet_name, block=pump_outlet_i_state)
            
            # Create reference for easy access
            setattr(self, f"{outlet_name}_vars", Reference(pump_outlet_i_state[:]))
            
            # Add to tracking lists
            self.outlet_list.append(outlet_name)
            self.outlet_blocks.append(pump_outlet_i_state)
        

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        self.mixer.calculate_scaling_factors()
        self.transfer_pump.calculate_scaling_factors()

    def initialize_build(self, outlvl=idaeslog.NOTSET, **kwargs):
        """
        Initialize the deaerator unit using Sequential Decomposition to determine optimal order
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        
        init_log.info("Starting deaerator initialization using Sequential Decomposition")
        
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
            
        return create_stream_table_dataframe(io_dict, time_point=time_point)


         