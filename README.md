# StringNet Herding: Adversarial Swarm Herding

This repository contains MATLAB simulation code for **herding an adversarial swarm** using a StringNet-based control strategy. The framework supports multiple herding scenarios through separate main scripts, each corresponding to a different simulation or experimental setup.

All simulation, planning and control, and environment parameters are defined in dedicated parameter files, enabling flexible customization without modifying the core simulation logic.

---

## Repository Structure
```
StringNet_Herding/
│
├── mainHerd.m # Core herding simulation
├── mainHerdGazebo.m # Gazebo-compatible simulation setup
├── mainHerdExperiments.m # Experimental demo setup
├── setupPaths.m # Script to setup all paths
│
├── config/
│ ├── AllParameters.m
│ ├── AllParametersGazebo.m
│ └── AllParametersExperiments.m
│
├── methods/ # Herding and control methods
├── utils/ # Utility and helper functions,including plotting tools
├── data/ # Saved simulation and experimental data
├── tests/ # Test scripts and validation cases
│
└── README.md
```

---

## Requirements

- MATLAB (R2020a or newer recommended)
- Gurobi (MATLAB) License
- No additional MATLAB toolboxes required unless explicitly stated in the code
- (Optional) Gazebo installation and configuration for `mainHerdGazebo.m`
- (Optional) ROS for experimental demo codes

---

## Simulation and Demo Setups

The repository provides three main herding setups. Each setup is implemented as a standalone MATLAB script and uses a corresponding parameter file located in the `config` directory.

---

### 1. Standard Herding Simulation

**Main script:**
`mainHerd.m`

**Parameter file:**
`config/AllParameters.m`

**Description:**
Baseline adversarial swarm herding simulation intended for algorithm development, debugging, and visualization. Runs entirely within MATLAB.

---

### 2. Gazebo-Based Herding Simulation

**Main script:**
`mainHerdGazebo.m`

**Parameter file:**
`config/AllParametersGazebo.m`

**Description:**
Interfaces with a Gazebo simulation environment to enable more realistic physics and robotics-oriented testing. Requires Gazebo to be running and properly configured.

---

### 3. Experimental Demo Setup

**Main script:**
`mainHerdExperiments.m`

**Parameter file:**
`config/AllParametersExperiments.m`

**Description:**
Designed for running this herding algorithm on real quadrotor platforms (PX4) using ROS framework.

**CPP Implementation:**
For a C++ implementation of this framework intended for experimental demonstrations, please refer to the following repository:

[https://github.com/vishnuchipade/StringNet_Herding_CPP](https://github.com/vishnuchipade/StringNet_Herding_CPP)

---

## Parameter Configuration

All simulation parameters are defined in the corresponding `AllParameters*.m` files.

Typical parameters include:

* Number of swarm agents and herders
* Adversarial behavior characteristics
* StringNet and control gains
* Initial positions and velocities
* Simulation time and integration step
* Environment constraints and obstacles
* Visualization, logging, and data-saving options

### Modifying Parameters

1. Open the appropriate parameter file (e.g., `config/AllParameters.m`)
2. Modify the desired variables
3. Save the file
4. Run the corresponding main script

Normally, no changes to the main scripts are required to adjust parameters, however, there may be additional parameters or flags that are defined in the main script which can be edited there only.

---

## How to Run

1. Open MATLAB
2. Set the repository root as the current directory:

   ```matlab
   cd path/to/StringNet_Herding
   ```
3. Run one of the following commands:

**Standard simulation**

```matlab
mainHerd
```

**Gazebo simulation**

```matlab
mainHerdGazebo
```

**Experimental setup**

```matlab
mainHerdExperiments
```

Ensure the corresponding parameter file is configured correctly before execution.

---

## Output and Data

Depending on the setup and parameter choices, simulations may produce:

* Real-time visualization of swarm and herding agents
* Agent trajectories and control signals
* Performance metrics and convergence measures
* Saved `.mat` files containing simulation or experimental data or initializtion data for reproducibility

Saved results are typically stored in the `data` directory. Output behavior is controlled through the parameter files.

---

## Testing

The `tests` directory contains scripts for:

* Validating individual methods
* Verifying expected herding behavior

Test scripts can be executed directly from MATLAB.

---

## Customization and Extensions

The framework can be extended by:

* Adding new herding or control strategies in the `methods` directory
* Creating additional parameter files for new scenarios
* Integrating learning-based or adaptive herding approaches

---

## Citation

If you use this code in your research, please cite the associated work:

```
@article{chipade2021multiagent,
  title={Multiagent planning and control for swarm herding in 2-D obstacle environments under bounded inputs},
  author={Chipade, Vishnu S and Panagou, Dimitra},
  journal={IEEE Transactions on Robotics},
  volume={37},
  number={6},
  pages={1956--1972},
  year={2021},
  publisher={IEEE}
}
```

---

## Contact

For questions, issues, or contributions, please open a GitHub issue or contact the repository owner.
