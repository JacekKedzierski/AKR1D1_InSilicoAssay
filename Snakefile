import pandas as pd
df = pd.read_csv('input/compounds.csv', sep=' ',header = None)
ligands = df[df.columns[1]]

print(ligands)

rule target:
    input:
        expand("output/log/EvaluatePoses/{ligand}.log",ligand=ligands)

rule CreateSmi:
    input:
        'input/compounds.csv'
    output:
        expand("output/CreateSmi/{lig}.smi",lig=ligands)
    shell:
        'src/CreateSmiFiles.sh {input}'

rule LigPrep:
    input:
        "output/CreateSmi/{ligand}.smi"
    output:
        "output/LigPrep/{ligand}.pdb"
    log:
        "output/LigPrep/{ligand}.log"
    params:
        schrodinger="$SCHRODINGER"

    shell:
        "{params.schrodinger}/ligprep -ismi {input} -osd {output} -g -ph 7.4 -pht 0.1 -epik -s 1 -t 1 -WAIT -NOJOBID > {log}"

rule Smina:
    input:
        rec="input/Receptor/AKR1D1.pdb",
        lig="output/LigPrep/{ligand}.pdb"
    output:
        "output/Smina/{ligand}.pdb"
    log:
        "output/log/Smina/{ligand}.log"
    params:
        smina="src/smina.static",
        center_x = 5.59,
        center_y = -18.25,
        center_z = 33.78,
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

# rule Glide:
#     input:
#         grid="input/Receptor/AKR1D1.zip",
#         lig="output/LigPrep/{lig}.mae"
#     output:
#         workdir="output/Glide/{lig}"
#     params:
#         tool="$SCHRODINGER",
#         home = '/mnt/jacek/jkedzierski/Documents/Projects/AKR1D1/AKR1D1_InSilicoAssay/'
#     log: 
#         "output/log/Glide/AKR1D1_{lig}.log"
#     priority: 3
#     threads: 22
#     message: "AKR1D1 and compound {lig} threads: 8"
#     shell:
#         """
#         echo "FORCEFIELD OPLS_2005" >> output/tmp_AKR1D1.in
#         echo "GRIDFILE" {params.home}{input.grid} >> output/tmp_AKR1D1.in
#         echo "LIGANDFILE" {params.home}{input.molecule} >> output/tmp_AKR1D1.in
#         echo "PRECISION SP" >> output/tmp_AKR1D1.in

#         mkdir -p {output.workdir}      
#         cd {output.workdir}   
#         {params.tool}/glide ../../output/tmp_AKR1D1.in -WAIT -HOST "localhost:8" -NJOBS 1          
#         """

rule DockingMerge:
    input:
        rec="input/Receptor/AKR1D1.pdb",
        lig="output/Smina/{ligand}.pdb"
    output:
        "output/DockingMerge/{ligand}.pdb"
    log:
        "output/log/DockingMerge/{ligand}.log"
    shell:
        "obabel -j -ipdb {input} -O {output}"

rule Pdb2Mae:
    input:
        "output/DockingMerge/{ligand}.pdb"
    output:
        "output/Pdb2Mae/{ligand}.maegz"
    log:
        "output/log/Pdb2Mae/{ligand}.log"
    params:
        schrodinger="$SCHRODINGER"
    shell:
        "{params.schrodinger}/utilities/structconvert -ipdb {input} -omae {output} && chmod 777 {output}"

rule EvaluatePoses:
    input:
        "output/Pdb2Mae/{ligand}.maegz"
    output:
        "output/log/EvaluatePoses/{ligand}.log"
    log:
        "output/EvaluatePoses/{ligand}.maegz"
    params:
        schrodinger="$SCHRODINGER"
    shell:
        "{params.schrodinger}/run pose_filter.py {input} {log} -a 'res.num 57' -hbond 1 -a 'res.num 119' -hbond 2 -m all -lig_asl 'res.num 2' -hbond_dist_max 2.5 -hbond_donor_angle 90.0 -hbond_acceptor_angle 60.0 -contact_dist_max 5.0 -ring_dist_max 5.0 -aromatic_dist_max 5.0 -complex -WAIT -NOJOBID > {output}"