from pyomo.environ import ConcreteModel
from scripts.series_turbine import series_tubine
import time

if __name__ == "__main__":
    m = ConcreteModel()

    params = {'HP_inlet_flow': 428, # t/h
              'LP_passout_limit':150, # t/h
              'MP_demand_flow':225, # t/h
              'LP_demand_flow': 204.5, # t/h
              "HP_pressure": 45, # bar
              "MP_pressure": 12.5, # bar
              "LP_pressure": 4.5, # bar
              "HP_temperature": 400, # C
    }

    start_time = time.time()
    series_tubine(m, params)
    end_time = time.time()
    print(f"Execution time: {end_time - start_time} seconds")

    #chapter 2 remove ML, copy paste and submit to energies
    #chapter 3 update when feedback comes back
    #chapter 4 submit as is but check with Oji about sensitive parts, need lit review in intro. Organise meeting with them to discuss
    #chapter 5&6 together by cutting some parts 
  