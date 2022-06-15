import pandas as pd
df = pd.read_csv('input/compounds.csv', sep=' ',header = None)
ligands = df[df.columns[1]]

print(ligands)

rule target:
    input:
        #expand("output/log/EvaluatePoses/{ligand}.log",ligand=ligands)
        #expand("output/Glide/{ligand}",ligand=ligands)
        #expand("output/log/EvaluatePoses/{ligand}.log",ligand=ligands)
        expand("output/Glide/{ligand}",ligand=ligands)

rule CreateSmi:
    input:
        'input/compounds.csv'
    output:
        expand("output/CreateSmi/{ligand}.smi",ligand=ligands)
    shell:
        'src/CreateSmiFiles.sh {input}'

rule LigPrepMae:
    input:
        expand("output/CreateSmi/{ligand}.smi",ligand=ligands)
    output:
        "output/LigPrep/{ligand}.mae"
    log:
        "output/LigPrep/{ligand}.log"
    params:
        schrodinger="$SCHRODINGER"

    shell:
        "{params.schrodinger}/ligprep -ismi {input} -osd {output} -g -ph 7.4 -pht 0.1 -epik -s 1 -t 1 -WAIT -NOJOBID > {log}"

rule LigPrepPdb:
    input:
        expand("output/CreateSmi/{ligand}.smi",ligand=ligands)
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
        lig=expand("output/LigPrep/{ligand}.pdb",ligand=ligands)
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

rule Glide:
    input:
        grid="input/Receptor/AKR1D1.zip",
        ligandMae=expand("output/LigPrep/{ligand}.mae",ligand=ligands)
    output:
        workdir="output/Glide/{ligand}"
    params:
        tool="$SCHRODINGER",
        home="/mnt/jacek/jkedzierski/Documents/Projects/AKR1D1/AKR1D1_InSilicoAssay/"
    log: 
        "output/log/Glide/AKR1D1_{ligand}.log"
    # priority: 3
    # threads: 22
    # message: "AKR1D1 and compound {liga} threads: 8"
    shell:
        """
        echo "FORCEFIELD OPLS_2005" >> output/AKR1D1.in
        echo "GRIDFILE" {params.home}{input.grid} >> output/AKR1D1.in
        echo "LIGANDFILE" {params.home}{input.ligandMae} >> output/AKR1D1.in
        echo "PRECISION SP" >> output/AKR1D1.in

        mkdir -p {output.workdir}      
        cd {output.workdir}   
        {params.tool}/glide ../../AKR1D1.in -WAIT -HOST "localhost:8" -NJOBS 1          
        """

rule DockingMerge:
    input:
        rec="input/Receptor/AKR1D1.pdb",
        ligandPdb=expand("output/Smina/{ligand}.pdb",ligand=ligands)
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

rule EvaluatePosesPv:
    input:
        "output/Glide/AKR1D1_pv.maegz/AKR1D1_pv.maegz"
    output:
        "output/log/EvaluatePoses/{lig}.log"
    log:
        "output/EvaluatePoses/{lig}.maegz"
    params:
        schrodinger="$SCHRODINGER",
        home = '/mnt/jacek/jkedzierski/Documents/Projects/AKR1D1/AKR1D1_InSilicoAssay/'
    shell:
        "{params.schrodinger}/run pose_filter.py {params.home}{input} {log} -a 'res.num 57' -hbond 1 -a 'res.num 119' -hbond 2 -m all -lig_asl 'res.num 2' -hbond_dist_max 2.5 -hbond_donor_angle 90.0 -hbond_acceptor_angle 60.0 -contact_dist_max 5.0 -ring_dist_max 5.0 -aromatic_dist_max 5.0 -WAIT -NOJOBID > {output}"