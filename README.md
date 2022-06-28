# AKR1D1_InSilicoAssay
This is a molecular model intended to identify inhibitors of the AKR1D1 threw out virtual screening.
Molecules should be supplied in the SMILES format in the input folder

Start:
conda activate AKR1D1_InSilicoAssay
export SCHRODINGER=/opt/schrodinger2020-2; export LM_LICENSE_FILE=27000@131.152.94.126
nohup snakemake --cores 1 --config ligands=compounds.csv&