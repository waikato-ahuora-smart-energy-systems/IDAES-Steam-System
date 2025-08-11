# IDAES-Steam-System

Open-source, equation-based synthesis and optimisation of heat exchanger networks 
Models apply the Gekko modelling language (https://github.com/BYU-PRISM/GEKKO) and solve using the open-source APOPT solver (https://apopt.com/) or COIN-OR solvers (https://www.coin-or.org)

## 🚀 Installation (via Conda + setup.py)

Follow these steps to install and run the IDAES-Steam-System using a Conda environment.

### 1. Clone the Repository

```bash
git clone 
cd IDAES-Steam-System
```

---

### 2. Install Miniconda (if not already installed)

Download and install **Miniconda** from:  
👉 https://docs.conda.io/en/latest/miniconda.html

> During setup, check the box to "Add Miniconda to my PATH environment variable" if you want to use it from any terminal.

Once installed, open **Anaconda Prompt** (Windows) or terminal (macOS/Linux). Do not use virtual environments as they dont work with packages outside of Python

---

### 3. Create and Activate a Conda Environment

```bash
conda create -n idaes-env python=3.12
conda activate idaes-env
```

---

### 5. Install the Package (Using setup.py)

From the project root:

```bash
pip install -e .
```

This installs the `idaeswork` package in **editable mode** and uses `requirements.txt` automatically.

---

### 6. Install Solvers (IDAES Extensions)

```bash
idaes get-extensions --extra petsc
```


---

### 7. Optional: Install COIN-OR Solvers

Download IPOPT or other COIN-OR binaries from:  
👉 https://www.jdhp.org/docs/notebook/python_pyomo_getting_started_0_installation_instructions_pyomo_and_solvers.html

Then set the solver in your script like this:

```python
from pyomo.environ import SolverFactory
solver = SolverFactory("ipopt")
solver.set_executable("C:/your/path/to/ipopt.exe")
```

---


## 🧼 Deleting the Conda Environment

To delete the environment:

```bash
conda deactivate
conda remove -n idaes-env --all
```

This will **not affect** any other Conda environments or your base Python install.

---

## 💡 Notes

- If using **VSCode**, make sure to install the **Python extension by Microsoft**, and select the `idaes-env` interpreter.
- If you modify the codebase, the `-e .` install ensures changes are reflected automatically without re-installing.

---


