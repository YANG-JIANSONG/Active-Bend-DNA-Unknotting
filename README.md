# Topoisomerase-induced DNA Knot Suppression via Active Bending

This repository contains the LAMMPS simulation scripts and Python analysis tools used in the study:

**"Topoisomerase suppresses DNA knotting via active-bending-induced symmetry breaking"**


---

# 1. System Requirements

## Hardware

* Multi-core CPU recommended (10+ cores recommended for batch simulations)
* 50 GB RAM or higher recommended for large-scale simulations

---

## Software Dependencies

### LAMMPS

Tested version:

```text
LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator)
Version: 3 Mar 2020
```

LAMMPS must be compiled with Python support enabled.

Example compilation options:

```bash
make yes-python
```

or using CMake:

```bash
-D PKG_PYTHON=yes
-D BUILD_SHARED_LIBS=yes
```

---

### Python

Tested with:

```text
Python 3.8+
```

Required Python packages:

```bash
pip install numpy pandas pythonknot
```

---

# 2. Directory Structure

```text
project/
│
├── README.md
│   Documentation and usage instructions
│
├── LICENSE
│   MIT License
│
├── Run.sh
│   Main batch simulation launcher
│
├── Calc_check.sh
│   Batch topology verification script
│
├── Active_bend.py
│   Generates local active-bending geometries
│
├── bend_generator.py
│   Generates initial DNA configurations using Bézier curves
│
├── in.py
│   Main LAMMPS simulation controller
│
├── curve_150.txt
├── wanqu_150.txt
│   Pre-generated active-bending geometry templates
│
├── lammps_31_L300_close.data
│   Example initial LAMMPS data file
│
├── core/
│   Core topology templates and reference structures
│   │
│   ├── pos_knot0.txt
│   ├── pos_knot31.txt
│   ├── pos_knot41.txt
│   ├── pos_knot41_tight.txt
│   ├── pos_knot51.txt
│   ├── pos_knot52.txt
│   ├── knot_core.txt
│   └── ...
│
├── source/
│   Main simulation and topology-analysis utilities
│   │
│   ├── in.py
│   ├── check.py
│   ├── check_hug.py
│   ├── hugging.py
│   ├── xyztolammps.py
│   ├── merge_overall_tri.sh
│   ├── table_knot_Alexander_polynomial.txt
│   ├── lammps*.data
│   └── result/
│
└── demo/
    Minimal runnable example for testing
    │
    ├── Run.sh
    ├── Calc_check.sh
    └── source/
        ├── in.py
        ├── check.py
        ├── hugging.py
        ├── example .data files
        └── example trajectory files
```

---

# 3. Core Algorithmic Workflow

The simulation workflow consists of four major stages:

1. **Local Active Bend Generation (`Active_bend.py`)**

   * Generates V-shaped local DNA bending geometries
   * Defines the enzyme-induced bending angle (`degree`)

2. **Global DNA Closure (`bend_generator.py`)**

   * Uses cubic Bézier curves to smoothly connect the remaining DNA chain
   * Produces topologically controlled closed DNA conformations

3. **Dynamic Simulation (`in.py`)**

   * Runs molecular dynamics simulations through the LAMMPS Python interface
   * Implements strand-passage events triggered by the hugging condition

4. **Automated Topology Verification (`Calc_check.sh`)**

   * Batch-checks all trajectories using topology analysis scripts
   * Determines the final knot state of each replica

---

# 4. File Descriptions

| File                    | Description                                                     |
| ----------------------- | --------------------------------------------------------------- |
| `Active_bend.py`        | Generates local active-bending geometries                       |
| `bend_generator.py`     | Generates initial DNA configurations using Bézier interpolation |
| `in.py`                 | Main simulation controller                                      |
| `Run.sh`                | Batch simulation launcher                                       |
| `Calc_check.sh`         | Batch topology verification script                              |
| `source/check.py`       | Determines final knot topology                                  |
| `source/hugging.py`     | Hugging-condition detection module                              |
| `source/xyztolammps.py` | Converts XYZ coordinates into LAMMPS data format                |

---

# 5. Quick Start

## Step 1: Generate Initial Configurations

```bash
python bend_generator.py
```

This script automatically calls `Active_bend.py`.

The active bending angle can be modified by changing the `degree` parameter within the script.

Example:

```python
degree = 150
```

---

## Step 2: Launch Simulations

```bash
bash Run.sh
```

This script automatically:

* Creates folders for different knot topologies
* Generates independent simulation replicas
* Launches LAMMPS simulations in batch mode

Independent replicas are generated using different random seeds.

---

## Step 3: Verify Topological States

After simulations are completed:

```bash
bash Calc_check.sh
```

This script automatically:

* Traverses all topology directories
* Updates and runs `check.py`
* Outputs topology analysis results

Typical output files include:

```text
output_all_save_*.txt
output_index_*.txt
trj_reco_*.xyz
```

---

# 6. Reproducing Manuscript Results

## Active Bending Effect

Compare simulations with:

```text
degree = 0
```

and

```text
degree = 150
```

to quantify the effect of active bending.

---

## Symmetry-Breaking Analysis

After simulations:

* Analyze `output_all_save_*.txt`
* Compare transition fluxes:

  * Knot → Unknot
  * Unknot → Knot

The active-bending mechanism enhances unknotting transitions and suppresses knot formation.

---

## Steady-State Knot Statistics

Use:

```bash
bash Calc_check.sh
```

to calculate:

* Final knot probabilities
* Suppression factors
* Topological steady-state distributions

---

# 7. Typical Runtime

Typical runtime on a standard workstation:

| Task                             | Runtime               |
| -------------------------------- | --------------------- |
| Initial configuration generation | Seconds to minutes    |
| Single MD trajectory             | ~10 minutes           |
| Full batch simulation            | Several hours to days |

Runtime depends on:

* Chain length
* Number of replicas
* CPU resources

---
# 8. Demo Usage

A minimal runnable example is provided in the `demo/` directory.

Enter the demo folder:

```bash id="5g4d8m"
cd demo
```

Launch the simulations using:

```bash id="a0c7ci"
bash Run.sh
```

The demo workflow automatically:

* Generates topology-dependent simulation folders
* Launches multiple independent replicas
* Performs batch molecular dynamics simulations using LAMMPS

Typical simulation time:

```text id="jlwmzt"
4–6 days
```

depending on:

* CPU performance
* Number of available cores
* System workload

After simulations are completed, run:

```bash id="t2oqc7"
bash Calc_check.sh
```

to verify the final topological states and generate statistical outputs.

---

# 9. Reproducibility Notes

For reproducibility:

* All simulations in the manuscript were generated using the scripts provided in this repository
* Independent replicas use different random seeds
* No manual intervention is required after launching `Run.sh`

Before re-running simulations, remove old output files to avoid mixing previous results with new trajectories.

---

# 10. License

This project is distributed under the MIT License.

---

# 11. Citation

If you use this code in your research, please cite:

> Jiansong Yang and Liang Dai,
> "Topoisomerase suppresses DNA knotting via active-bending-induced symmetry breaking" (2025).
