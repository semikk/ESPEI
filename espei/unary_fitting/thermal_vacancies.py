import operator

# TODO (Brandon): take this as input from a JSON file
phase_models = {
    "components": ["W", "VA"],
    "phases": {
        "BCC_A2": {
            "sublattice_model": [["W", "VA"]],
            "sublattice_site_ratios": [1.0],
        }
    }
}

import numpy as np
from sympy import Symbol, Eq

# TODO (Brandon): these will eventually have to be a part of ESPEI datasets
# the next parts will be a part of datasets that are parsed
T_grid = np.arange(0, 2500, 10)
G_f_grid = 337120-26.605*T_grid
G_s_grid = 776676-44.423*T_grid


datasets = {
    'G_SUBLIMATION': G_s_grid,
    'G_VA_FORMATION': G_f_grid,
    'T': T_grid
}

# TODO (Usman): Make this function real
def phase_fit_thermal_vacancies(dbf, phase_name, subl_model, site_ratios, datasets):
    # this function takes all the input arguments, fits to the datasets, then modifies the database by adding the parameter
    # for the time being, assume that datasets is just a dictrionary of {'G_SUBLIMATION': G_s_grid, ...}, see above defintion
    return
    G_s_grid = datasets['G_SUBLIMATION']

# this part takes the grid of temperature values and Gibbs energy formation and
# sublimation and calculates what the L parameters should be on the T grid
L0 = Symbol('L0')
L1 = Symbol('L1')

L0_grid = []
L1_grid = []
for idx in range(T_grid.shape[0]):
    # solve two equations
    eq_1 = Eq(L0+L1, G_f_grid[idx])
    eq_2 = Eq(L0-L1, G_s_grid[idx])
    result = sympy.solve([L0+L1-G_f_grid[idx], L0-L1-G_s_grid[idx]], L0, L1)
    L0_grid.append(float(result[L0]))
    L1_grid.append(float(result[L1]))
#     result = sympy.solveset([eq_1, eq_2], L0, L1)

L0_grid = np.array(L0_grid)
L1_grid = np.array(L1_grid)

# this part fits the linear equation Ax=b where x is the temperature dependent terms for each interaction,
# L0 first
A = np.stack([np.ones(T_grid.size), T_grid]).T
b = L0_grid.reshape((L0_grid.size, 1))
L0_A, L0_B = np.matmul(np.linalg.pinv(A),b)

# L1
b = L1_grid.reshape((L1_grid.size, 1))
L1_A, L1_B = np.matmul(np.linalg.pinv(A),b)

# TODO (Usman): add code to put these parameters in the database,
# dbf.add_parameter()

#### END thermal vacancies function

# generate_thermal_vacancies_tdb

# this code needs phase_models to run
from pycalphad import Database
dbf = Database()
dbf.elements = set(phase_models['components'])
for phase_name, phase_obj in sorted(phase_models['phases'].items(), key=operator.itemgetter(0)):
    # Perform parameter selection and single-phase fitting based on input
    # TODO: More advanced phase data searching
    site_ratios = phase_obj['sublattice_site_ratios']
    subl_model = phase_obj['sublattice_model']
    dbf.add_phase(phase_name, dict(), site_ratios)
    dbf.add_phase_constituents(phase_name, subl_model)
    dbf.add_structure_entry(phase_name, phase_name)
    # modifies the database in place
    phase_fit_thermal_vacancies(dbf, phase_name, subl_model, site_ratios, datasets)
