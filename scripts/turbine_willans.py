from turbine_base_model import TurbineBaseData
from pyomo.environ import Var, Constraint
from idaes.core.util.math import smooth_min, smooth_max

@declare_process_block_class("TurbineWillans")
class TurbineWillinsData(TurbineBaseData):

    def add_mechanical_work_definition(blk):
        """"
        Define a willins formulation here
        """
        blk.a = Var(
            blk.flowsheet().time,
            initialize=1.0,
            doc="Coefficient for mechanical work calculation",
        )
        @blk.Constraint(blk.flowsheet().time, doc="Mechanical work balance")
        def mechanical_work_eqn(b, t):
            return b.mechanical_work[t] == b.electric_power[t] - b.heat_duty[t]