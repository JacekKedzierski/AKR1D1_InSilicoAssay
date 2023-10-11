# AKR1D1 InSilicoAssay

**Description:**

AKR1D1 InSilicoAssay is a computational tool designed for the identification of inhibitors of the AKR1D1 enzyme through virtual screening. This README provides an overview of how to use the tool, including setup, input requirements, and execution. 

**Usage:**

## Setup

1. First, activate the conda environment dedicated to AKR1D1 InSilicoAssay:

    ```bash
    conda activate AKR1D1_InSilicoAssay
    ```

2. Set the SCHRODINGER environment variable (assuming you have Schrödinger suite installed):

    ```bash
    export SCHRODINGER=/opt/schrodinger2020-2
    ```

3. Ensure that the LM_LICENSE_FILE environment variable is set. This variable should be configured with your Schrödinger license information.

## Input

Molecules for screening should be supplied in the SMILES format and placed in the 'input' folder. The tool expects a CSV file called 'compounds.csv' to be present in the input folder. This CSV file should contain a list of the compounds to be screened, along with their relevant information.

## Execution

To start the virtual screening process, execute the following command:

```bash
nohup snakemake --cores 1 --config ligands=compounds.csv&
```

- `nohup` is used to run the process in the background and ensure it continues running even if you log out.
- `snakemake` is the workflow management system used to execute the virtual screening.
- `--cores 1` specifies the number of CPU cores to use. You can adjust this based on your system's capacity.
- `--config ligands=compounds.csv` configures the ligands to be screened using the 'compounds.csv' file in the input folder.

**Output:**

The virtual screening process will generate results that can be found in the output directory. The results typically include information about potential AKR1D1 inhibitors and their associated properties.

**Note:**

- Make sure to configure your Schrödinger license information properly.
- Ensure that the 'compounds.csv' file contains the necessary compound information.
- The tool can be customized further as needed to meet specific screening requirements.

For more information and troubleshooting, please refer to the documentation or contact the developers of AKR1D1 InSilicoAssay.

**References:**

https://www.sciencedirect.com/science/article/pii/S0378427423002205?via%3Dihub
