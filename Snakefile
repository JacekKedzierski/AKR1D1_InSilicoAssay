import pandas as pd

ligands = pd.read_csv('input/compounds.csv', sep=' ').ligand

print(ligands)

rule target:
    input:
        expand("output/EvaluatePoses/{ligand}.mae",ligand=ligands)

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
        n_poses = 1,
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

rule DockingMerge:
    input:
        rec="input/receptor/AKR1D1.pdb",
        lig="output/smina/{ligand}.pdb"
    output:
        "output/poses/{ligand}.pdb"
    log:
        "output/log/poses/{ligand}.log"
    shell:
        "obabel -j -ipdb {input} -O {output}"

rule Pdb2Mae:
    input:
        "output/poses/{ligand}.pdb"
    output:
        "output/Pdb2Mae/{ligand}.mae"
    log:
        "output/log/Pdb2Mae/{ligand}.log"
    params:
        schrodinger="$SCHRODINGER"
    shell:
        "{params.schrodinger}/utilities/structconvert -ipdb {input} -omae {output} && chmod 777 {output}"

rule EvaluatePoses:
    input:
        "output/Pdb2Mae/{ligand}.mae"
    output:
        "output/EvaluatePoses/{ligand}.mae"
    log:
        "output/log/EvaluatePoses/{ligand}.log"
    params:
        schrodinger="$SCHRODINGER"
    shell:
        "{params.schrodinger}/run {input} {output} -a 'res.num 58' -hbond 1 -a 'res.num 120' -hbond 2 -m all -lig_asl 'res.num 900' -hbond_dist_max 2.5 -hbond_donor_angle 90.0 -hbond_acceptor_angle 60.0 -contact_dist_max 5.0 -ring_dist_max 5.0 -aromatic_dist_max 5.0 -TMPLAUNCHDIR -WAIT -NOJOBID > {log}"