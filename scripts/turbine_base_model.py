#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Standard IDAES pressure changer model.
"""
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

# Changing existing config block attributes
# pylint: disable=protected-access

# TODO: Keegan
# Add different mech work calculation modes in turbine_willins.py, also need to figure out how 
# test on series_turbine.py 
# test different units for willans params
# bring in other willans models

# Import Python libraries
from enum import Enum

# Import Pyomo libraries
from pyomo.environ import (
    Block,
    value,
    Var,
    Expression,
    Constraint,
    Reference,
    check_optimal_termination,
    Reals,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    EnergyBalanceType,
    MomentumBalanceType,
    MaterialBalanceType,
    ProcessBlockData,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.exceptions import PropertyNotSupportedError, InitializationError
from idaes.core.util.config import is_physical_parameter_block
import idaes.logger as idaeslog
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.initialization import SingleControlVolumeUnitInitializer
from idaes.core.util import to_json, from_json, StoreSpec
from idaes.core.util.math import smooth_max, safe_sqrt, sqrt, smooth_min
from pyomo.environ import units as pyunits


__author__ = "Emmanuel Ogbe, Andrew Lee"
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("TurbineBase")
class TurbineBaseData(UnitModelBlockData):
    """
    Standard Compressor/Expander Unit Model Class
    """

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.useDefault.
**Valid values:** {
**MaterialBalanceType.useDefault - refer to property package for default
balance type
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.useDefault.
**Valid values:** {
**EnergyBalanceType.useDefault - refer to property package for default
balance type
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
            default=MomentumBalanceType.pressureTotal,
            domain=In(MomentumBalanceType),
            description="Momentum balance construction flag",
            doc="""Indicates what type of momentum balance should be
constructed, **default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Phase equilibrium construction flag",
            doc="""Indicates whether terms for phase equilibrium should be
constructed, **default** = False.
**Valid values:** {
**True** - include phase equilibrium terms
**False** - exclude phase equilibrium terms.}""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property
calculations, **default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property
block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "calculation_method",
        ConfigValue(
            default='isentropic',
            domain=str,
            description="Calculation method used to model mechanical work",
            doc="""Property parameter object used to define property
calculations, **default** - 'isentropic'.
**Valid values:** {
**isentropic** - default method, uses isentropic efficiency to determine work
**simple_willans** - simple willans line requring slope and intercept.}""",
        ),
    )

    def build(self):
        """

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build
        super().build()

        # Add a control volume to the unit including setting up dynamics.
        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        # Add geometry variables to control volume
        if self.config.has_holdup:
            self.control_volume.add_geometry()

        # Add inlet and outlet state blocks to control volume
        self.control_volume.add_state_blocks(
            has_phase_equilibrium=self.config.has_phase_equilibrium
        )

        # Add mass balance
        # Set has_equilibrium is False for now
        # TO DO; set has_equilibrium to True
        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=self.config.has_phase_equilibrium,
        )

        # Add energy balance
        eb = self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type, has_work_transfer=True
        )

        # add momentum balance
        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type, has_pressure_change=True
        )

        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()

        # Set Unit Geometry and holdup Volume
        if self.config.has_holdup is True:
            self.volume = Reference(self.control_volume.volume[:])

        self.work_mechanical = Reference(self.control_volume.work[:])

        # Add Momentum balance variable 'deltaP'
        self.deltaP = Reference(self.control_volume.deltaP[:])

        # Performance Variables
        self.ratioP = Var(self.flowsheet().time, initialize=1.0, doc="Pressure Ratio")

        # Pressure Ratio
        @self.Constraint(self.flowsheet().time, doc="Pressure ratio constraint")
        def ratioP_calculation(self, t):
            return (
                self.ratioP[t] * self.control_volume.properties_in[t].pressure
                == self.control_volume.properties_out[t].pressure
            )

        units_meta = self.control_volume.config.property_package.get_metadata()

        # Get indexing sets from control volume
        # Add isentropic variables
        self.efficiency_isentropic = Var(
            self.flowsheet().time,
            initialize=0.8,
            doc="Efficiency with respect to an isentropic process [-]",
        )
        self.work_isentropic = Var(
            self.flowsheet().time,
            initialize=0.0,
            doc="Work input to unit if isentropic process",
            units=units_meta.get_derived_units("power"),
        )

        # Add willans line parameters
        if self.config.calculation_method == "simple_willans":
            self.willans_slope = Var(
                self.flowsheet().time,
                initialize=1.0,
                doc="Slope of willans line",
                units=units_meta.get_derived_units("energy") / units_meta.get_derived_units("amount"),
            )

            self.willans_intercept = Var(
                self.flowsheet().time,
                initialize=1.0,
                doc="Intercept of willans line",
                units=units_meta.get_derived_units("power"),
            )

            self.willans_max_mol = Var(
                self.flowsheet().time,
                initialize=1.0,
                doc="Max molar flow of willans line",
                units=units_meta.get_derived_units("amount") / units_meta.get_derived_units("time"),
            )

            
        # Build isentropic state block
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = self.config.has_phase_equilibrium
        tmp_dict["defined_state"] = False

        self.properties_isentropic = self.config.property_package.build_state_block(
            self.flowsheet().time, doc="isentropic properties at outlet", **tmp_dict
        )

        # Connect isentropic state block properties
        @self.Constraint(
            self.flowsheet().time, doc="Pressure for isentropic calculations"
        )
        def isentropic_pressure(self, t):
            return (
                self.properties_isentropic[t].pressure
                == self.control_volume.properties_out[t].pressure
            )

        # This assumes isentropic composition is the same as outlet
        self.add_state_material_balances(
            self.config.material_balance_type,
            self.properties_isentropic,
            self.control_volume.properties_out,
        )

        # This assumes isentropic entropy is the same as inlet
        @self.Constraint(self.flowsheet().time, doc="Isentropic assumption")
        def isentropic(self, t):
            return (
                self.properties_isentropic[t].entr_mol
                == self.control_volume.properties_in[t].entr_mol
            )
        
        @self.Expression(
                self.flowsheet().time,
                doc="calculate ideal amount of work per mole of fluid"
        )
        def work_isentropic_mol(self, t):
            return self.properties_isentropic[t].enth_mol - self.control_volume.properties_in[t].enth_mol
        
        @self.Expression(
                self.flowsheet().time,
                doc="calculate actual amount of work per mole of fluid"
        )
        def work_mechanical_mol(self, t):
            return self.properties_isentropic[t].enth_mol - self.control_volume.properties_in[t].enth_mol

        # Actual work
        @self.Constraint(
            self.flowsheet().time, doc="Actual mechanical work calculation"
        )
        def actual_work(self, t):
            # if config.calc method == isentropic:
            if self.config.calculation_method == "isentropic":
                return self.work_mechanical_mol[t] == (
                    self.work_isentropic_mol[t] * self.efficiency_isentropic[t]
                )
            elif self.config.calculation_method == 'simple_willans':
                eps = 1e-4  # smoothing parameter; smaller = closer to exact max, larger = smoother
                
                return self.work_mechanical[t] == smooth_min(
                    -(self.willans_slope[t] * self.control_volume.properties_in[t].flow_mol - self.willans_intercept[t]) / (self.willans_slope[t] * self.willans_max_mol[t] - self.willans_intercept[t]),
                    0.0 * pyunits.W,
                    eps
                    ) * (self.willans_slope[t] * self.willans_max_mol[t] - self.willans_intercept[t])
             
                
        self.add_mechanical_work_definition()

        # Property packages should define
        # properties_in.enth_mol
        # properties_in.entr_mol
        # properties_out.flow_mol
        # .pressure
        # .temperature
        # 
    
    def add_mechanical_work_definition(self):

        # Isentropic work
        @self.Constraint(
            self.flowsheet().time, doc="Calculate work of isentropic process"
        )
        def isentropic_energy_balance(self, t):
            return self.work_isentropic[t] == ( self.work_isentropic_mol[t] ) * self.control_volume.properties_in[t].flow_mol
        


    def initialize_build(
        blk,
        state_args=None,
        routine=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for pressure changer initialization routines

        Keyword Arguments:
            routine : str stating which initialization routine to execute
                        * None - use routine matching thermodynamic_assumption
                        * 'isentropic' - use isentropic initialization routine
                        * 'isothermal' - use isothermal initialization routine
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None, use default solver)

        Returns:
            None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        cv = blk.control_volume
        t0 = blk.flowsheet().time.first()
        state_args_out = {}

        if state_args is None:
            state_args = {}
            state_dict = cv.properties_in[t0].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        # Get initialisation guesses for outlet and isentropic states
        for k in state_args:
            if k == "pressure" and k not in state_args_out:
                # Work out how to estimate outlet pressure
                if cv.properties_out[t0].pressure.fixed:
                    # Fixed outlet pressure, use this value
                    state_args_out[k] = value(cv.properties_out[t0].pressure)
                elif blk.deltaP[t0].fixed:
                    state_args_out[k] = value(state_args[k] + blk.deltaP[t0])
                elif blk.ratioP[t0].fixed:
                    state_args_out[k] = value(state_args[k] * blk.ratioP[t0])
                else:
                    # Not obvious what to do, use inlet state
                    state_args_out[k] = state_args[k]
            elif k not in state_args_out:
                state_args_out[k] = state_args[k]

        # Initialize state blocks
        flags = cv.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
            state_args=state_args,
        )
        cv.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=False,
            state_args=state_args_out,
        )

        init_log.info_high("Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize Isentropic block

        blk.properties_isentropic.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )

        init_log.info_high("Initialization Step 2 Complete.")

        # Skipping step 3 because Isothermal had problems.

        # ---------------------------------------------------------------------
        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Initialization Step 4 {}.".format(idaeslog.condition(res)))


        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.control_volume.release_state(flags, outlvl)

        if not check_optimal_termination(res):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        init_log.info(f"Initialization Complete: {idaeslog.condition(res)}")

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        if hasattr(self, "deltaP"):
            var_dict["Mechanical Work"] = self.work_mechanical[time_point]
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]
        if hasattr(self, "ratioP"):
            var_dict["Pressure Ratio"] = self.ratioP[time_point]
        if hasattr(self, "efficiency_pump"):
            var_dict["Efficiency"] = self.efficiency_pump[time_point]
        if hasattr(self, "efficiency_isentropic"):
            var_dict["Isentropic Efficiency"] = self.efficiency_isentropic[time_point]

        return {"vars": var_dict}

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if hasattr(self, "work_fluid"):
            for t, v in self.work_fluid.items():
                iscale.set_scaling_factor(
                    v,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                )

        if hasattr(self, "work_mechanical"):
            for t, v in self.work_mechanical.items():
                iscale.set_scaling_factor(
                    v,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                )

        if hasattr(self, "work_isentropic"):
            for t, v in self.work_isentropic.items():
                iscale.set_scaling_factor(
                    v,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                )

        if hasattr(self, "ratioP_calculation"):
            for t, c in self.ratioP_calculation.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].pressure,
                        default=1,
                        warning=True,
                    ),
                    overwrite=False,
                )

        if hasattr(self, "fluid_work_calculation"):
            for t, c in self.fluid_work_calculation.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.deltaP[t], default=1, warning=True
                    ),
                    overwrite=False,
                )

        if hasattr(self, "actual_work"):
            for t, c in self.actual_work.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                    overwrite=False,
                )

        if hasattr(self, "isentropic_pressure"):
            for t, c in self.isentropic_pressure.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].pressure,
                        default=1,
                        warning=True,
                    ),
                    overwrite=False,
                )

        if hasattr(self, "isentropic"):
            for t, c in self.isentropic.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.properties_in[t].entr_mol,
                        default=1,
                        warning=True,
                    ),
                    overwrite=False,
                )

        if hasattr(self, "isentropic_energy_balance"):
            for t, c in self.isentropic_energy_balance.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                    overwrite=False,
                )

        if hasattr(self, "zero_work_equation"):
            for t, c in self.zero_work_equation.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(
                        self.control_volume.work[t], default=1, warning=True
                    ),
                )

        if hasattr(self, "state_material_balances"):
            cvol = self.control_volume
            phase_list = cvol.properties_in.phase_list
            phase_component_set = cvol.properties_in.phase_component_set
            mb_type = cvol._constructed_material_balance_type
            if mb_type == MaterialBalanceType.componentPhase:
                for (t, p, j), c in self.state_material_balances.items():
                    sf = iscale.get_scaling_factor(
                        cvol.properties_in[t].get_material_flow_terms(p, j),
                        default=1,
                        warning=True,
                    )
                    iscale.constraint_scaling_transform(c, sf)
            elif mb_type == MaterialBalanceType.componentTotal:
                for (t, j), c in self.state_material_balances.items():
                    sf = iscale.min_scaling_factor(
                        [
                            cvol.properties_in[t].get_material_flow_terms(p, j)
                            for p in phase_list
                            if (p, j) in phase_component_set
                        ]
                    )
                    iscale.constraint_scaling_transform(c, sf)
            else:
                # There are some other material balance types but they create
                # constraints with different names.
                _log.warning(f"Unknown material balance type {mb_type}")

