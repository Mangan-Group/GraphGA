Supporting code for the unpublished article: Dreyer, KS et al. GCAD: a Computational Framework for Mammalian Genetic Program Computer-Aided Design.

## Summary of README contents

+ [Code overview](#code-overview)
+ [Setup and installation](#setup-and-installation)
+ [Run instructions](#run-instructions)


## Code overview 

This code contains the following key components:
+ Files comprising the GA algorithm setup and functions
+ A file for each test case included in the associated paper (Amplifier, Signal Conditioner, and Pulse Generator), with associated methods to simulate test case ODEs and calculate objectives for optimization
+ A file to run a given test case with specified settings (e.g., GA hyperparameters) and plot optimal results
+ Files to further analyze GA results (e.g., to assess sub-optimal solutions or compare results across GA initialization seeds)
+ Files to setup and run hyperparameter optimization
+ Files to run a combinatorial search for each test case
+ Files to analyze the Pulse Generator experimental time series data
+ Files to simulate circuits with different Z matrix manifestations to calculate a confidence interval for the population model objective functions
+ A file to produce the plots included in the main and supplementary figures of the paper


## Setup and installation

1. To clone the repository, navigate to the desired directory for the repository on your computer and run the following command:

```bash
$ git clone https://github.com/Mangan-Group/GraphGA.git
```


2. We recommend using a conda environment for dependency management. The GCAD_environment.yml file can be used to create an environment with the necessary dependencies (specified versions indicate versions used to produce results reported in the paper). 

    To create an environment with the .yml file, run the following command:

```bash
$ conda env create -f GCAD_environment.yml
```


3. There are some hard-coded paths within the code, in the files listed below. These need to be updated to include the base path to the repository on your computer (<path_to_your_repository>). The following is a list of these instances and how they should be updated:

    - flow_cytometry_calculations.py (top of file):\
        plt.style.use("<path_to_your_repository>/paper.mplstyle.py")

    - plot_search_results.py (top of file):\
        plt.style.use("<path_to_your_repository>/paper.mplstyle.py")

    - load_files_pop.py (top of file):\
        repo_path = "<path_to_your_repository>"

    - load_Z_mat_samples.py (top of file):\
        repo_path = "<path_to_your_repository>"

    - Reference_pop.py (top of file):\
        repo_path = "<path_to_your_repository>"\
        Z_path = "<path_to_your_repository>/Z_matrix_samples/"


## Run instructions

The file run_test_case.py can be used to run the GA for one of the test cases included in the paper, with specified settings (e.g., hyperparameters). The following details how to specify run settings and execute the code:

1. Run settings are set using a settings.json file. Here, we include an example settings file for each test case (Amplifier_settings.json, SignalConditioner_settings.json, and PulseGenerator_settings.json). Detailed descriptions of each settings parameter are as follows:

    - "test_case": a string specifying which test case will be run ("Amplifier", "SignalConditioner", or "PulseGenerator").

    - "promo_node": a string specifying the promoter node name (used to specify the promoter parameters). Set this as "P1" to use the endogenous promoter parameters used in the paper simulations ("P_exp_amp" is used to simulate the "test scenario" Amplifiers in the paper).

    - "dose_specs": a list specifying [minimum_dose, maximum_dose, dose increment] to create a list of possible plasmid doses to allowed in the GA search space. For the paper, we use [5, 75, 5].

    - "max_part": an integer specifying the maximum number of parts allowed in a circuit (i.e., maximum activators and/or inhibitors). Allowable values are 1 or 2.

    - "inhibitor": a boolean specifying whether to allow for inhibitors in the circuit.

    - "DsRed_inhibitor": a boolean specifying whether to use the DsRED inhibitors (as opposed to the weaker synTF inhibitors). For the paper, only DsRED inhibitors are used.

    - "num_dict": a dictionary specifying the number of 1-part circuits and 2-part circuits in the GA search space. The total number of circuits must be an even number. Both dictionary keys must be specified (e.g., if using all 1-part circuits, specify "2": 0 in num_dict). ***Note that either this parameter OR "pop_size" AND "pop_ratio" are set***

    - "pop_size": a float specifying ***1/2*** of the total population size. This number is multiplied by 2 in the code to ensure an even number of circuits is used. ***Note that either this parameter AND "pop_ratio" are set, OR "num_dict is set.***

    - "pop_ratio": a float specifying the ratio of 1- to 2-part circuits (e.g., a pop_ratio of 0.25 and pop_size of 100 means 25 circuits will contain 1 part and 75 will contain 2 parts). ***Note that either this parameter AND "pop_size" are set, OR "num_dict is set.***
 
    - "n_gen": an int specifying the number of generations to run the GA.

    - "probability_mutation": a float specifying the probability of mutation.

    - "probability_crossover: a float specifying the probability of crossover.

    - "mutate_dose": a boolean speficying whether to include dose mutation as a mutation type. For the paper, we set this to true.

    - "pop": a boolean specifying whether to use the population model (default population model is a 20-cell model).
    
    - "num_processes": an integer specifying the number of processes when using parallelization (with Python's multiprocessing package). Only relevant when running with "pop" = true.

    - "obj_labels": a list of strings specifying the names of the objectives. These will remain the values specified in the example settings files unless running the Pulse Generator test case. See pulse_generator_problem.py for other options.

    - "max_time": an integer specifying the maximum simulation time in hours. For the paper, 42 was used except for in the Pulse Generator test case, where 126 was used.

    - "plot": a boolean specifying whether to plot the graph networks for all circuits that lie on the pareto front (only relevant for multi-objective optimization). This will plot the same number of circuits as in the GA population (generally 50+), so use with caution.

    - "repository_path": a string containing the path to the repository directory on your computer.

    - "folder_name": a string specifying the name of the output directory that is created with the GA run, which will contain all GA outputs (data and plots). The current date will be prepended to this folder name.


2. To save results from the GA runs, you must create a GA_results directory within the repository directory on your computer before running the code. A sub-directory will be created each time you run run_test_case.py. 


3. Each time the GA is run, an initialization seed is used for reproducibility. Currently, run_test_case.py is set up to run the GA 10 times, for seed=0-9. Commented code/notes within the file indicate the necessary change to run the GA once, with a specified seed.


4. run_test_case.py can be executed in the command line, after completing the above steps. Simply run this command, after activating your conda environment:

```bash
$ python run_test_case.py <"path_to_settings.json">
```
where "path_to_settings.json" is the path to the settings file.

If setting the seed with modifications from step (3), run:

```bash
$ python run_test_case.py <"path_to_settings.json"> <seed>
```
where seed is the integer seed value for the run.
