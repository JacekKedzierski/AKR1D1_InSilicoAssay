import pandas as pd

ligands = pd.read_csv('input/compounds.csv', sep=' ').ligand

print(ligands)

rule target:
    input:
        expand("output/poses/{ligand}.pdb",ligand=ligands),

rule create_smi:
    input:
        'input/compounds.csv'
    output:
        expand("output/ligands/smi/{lig}.smi",lig=ligands)
    shell:
        'src/create_smifiles.sh {input}'

rule ligprep:
    input:
        "output/ligands/smi/{ligand}.smi"
    output:
        "output/ligands/prep/{ligand}.pdb"
    log:
        "output/log/prep/{ligand}.log"
    params:
        schrodinger="$SCHRODINGER"

    shell:
        "{params.schrodinger}/ligprep -ismi {input} -osd {output} -g -ph 7.4 -pht 0.1 -epik -s 1 -t 1 -WAIT -NOJOBID > {log}"

rule smina:
    input:
        rec="input/receptor/AKR1D1.pdb",
        lig="output/ligands/prep/{ligand}.pdb"
    output:
        "output/smina/{ligand}.pdb"
    log:
        "output/log/smina/{ligand}.log"
    params:
        smina="src/smina.static",
        center_x = -7.570,
        center_y = -17.030,
        center_z = 33.68,
        boxsize = 15,
        ex = 2,
        n_poses = 10,
    threads: 4
    shell:
        "{params.smina} -r {input.rec} "
        "-l {input.lig} "
        "--center_x {params.center_x} "
        "--center_y {params.center_y} "
        "--center_z {params.center_z} "
        "--size_x {params.boxsize} "
        "--size_y {params.boxsize} "
        "--size_z {params.boxsize} "
        "--seed 23 "
        "-o {output} "
        "--exhaustiveness {params.ex} "
        "--num_modes {params.n_poses} "
        "--cpu {threads} > {log} "

rule DockingSplit:
    input:
        lig="output/smina/{ligand}.pdb"
    output:
        "output/sep/{ligand}.pdb"
    log:
        "output/log/sep/{ligand}.log"
    shell:
        "obabel {input} -O {output} -m"

rule DockingMerge:
    input:
        rec="input/receptor/AKR1D1.pdb",
        lig="output/sep/{ligand}.pdb"
    output:
        "output/poses/{ligand}.pdb"
    log:
        "output/log/poses/{ligand}.log"
    shell:
        "obabel -j -ipdb {input} -O {output}"