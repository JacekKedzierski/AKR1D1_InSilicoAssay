import pandas as pd
import os
cwd = os.getcwd()+"/"
compounds="input/"+config["ligands"]
df = pd.read_csv(compounds, sep=' ',header = None)
ligands = df[df.columns[1]]

wildcard_constraints:
   ligands = '[a-zA-Z]'

rule target:
    input:
        "output/log/HappyHB/Done.txt"

rule LigPrepMae:
    input:
        expand({compounds},compounds=config["ligands"])
    output:
        "output/LigPrep/compounds.mae"
    log:
        "output/LigPrep/compounds.log"
    params:
        schrodinger="$SCHRODINGER"

    shell:
        "{params.schrodinger}/ligprep -icsv {input} -osd {output} -g -ph 7.4 -pht 0.1 -epik -s 1 -t 1 -WAIT -NOJOBID > {log}"

rule Neutral:
    input:
        "output/LigPrep/compounds.mae"
    output:
        "output/Neutral/compounds-out.maegz"
    params:
        schrodinger="$SCHRODINGER"

    shell:
        "{params.schrodinger}/utilities/ligfilter {input} -o {output} -f input/Receptor/rules.lff -WAIT -NOJOBID"      

rule Glide:
    input:
        grid="input/Receptor/AKR1D1.zip",
        ligandMae="output/Neutral/compounds-out.maegz"
    output:
        workdir="output/Glide/",
        infile="output/In/compounds.in",
        pose="output/Glide/compounds_pv.maegz"
    params:
        tool="$SCHRODINGER",
        home=cwd
    log: 
        "output/log/Glide/compounds.log"
    shell:
        """
        echo "FORCEFIELD OPLS_2005" >> {output.infile}
        echo "GRIDFILE" {params.home}{input.grid} >> {output.infile}
        echo "LIGANDFILE" {params.home}{input.ligandMae} >> {output.infile}
        echo "PRECISION SP" >> {output.infile}

        mkdir -p {output.workdir}      
        cd {output.workdir}   
        {params.tool}/glide ../../{output.infile} -WAIT -HOST "localhost:1" -NJOBS 1
        """

rule EvaluatePosesPv:
    input:
        workdir="output/Glide/compounds_pv.maegz"
    output:
        "output/EvaluatePosesPv/compounds.maegz"
    log:
        "output/log/EvaluatePosesPv/compounds.log"
    params:
        schrodinger="$SCHRODINGER",
        home=cwd
    shell:
        "{params.schrodinger}/run pose_filter.py {params.home}{input} {output} -a 'res.num 58' -hbond 1 -a 'res.num 120' -hbond 2 -m all -lig_asl 'SMARTS. [R]=O' -hbond_dist_max 2.5 -hbond_donor_angle 90.0 -hbond_acceptor_angle 60.0 -contact_dist_max 5.0 -ring_dist_max 5.0 -aromatic_dist_max 5.0 -WAIT -NOJOBID > {log}"

rule GlideDockingMergePv:
    input:
        "output/EvaluatePosesPv/compounds.maegz"
    output:
        tmp="output/EvaluatePosesPv/compounds-out_complex.maegz",
        end="output/GlideDockingMerge/compounds-out_complex.maegz"
    log:
        "output/log/GlideDockingMergePv/compounds.log"
    params:
        schrodinger="$SCHRODINGER",
        home=cwd
    shell:
        """
        {params.schrodinger}/run pv_convert.py -m {input} > {log}
        cp {output.tmp} {output.end}
        """
rule MaeGz2Pdb:
    input:
        "output/GlideDockingMerge/compounds-out_complex.maegz"
    log:
        tmp="compounds-out_complex.pdb",
        end="output/MaeGz2Pdb/compounds-out_complex.pdb"
    output:
        "output/log/MaeGz2Pdb/compounds.log"
    params:
        schrodinger="$SCHRODINGER",
        home=cwd
    shell:
        """
        {params.schrodinger}/utilities/structconvert  -imae {input} -opdb {log.tmp} > {output}
        mv *.pdb output/MaeGz2Pdb/.
        """

rule HappyHB:
    input:
        "output/log/MaeGz2Pdb/compounds.log"
    output:
        "output/log/HappyHB/Done.txt"
    params:
        home = '/mnt/jacek/jkedzierski/Documents/Projects/AKR1D1/AKR1D1_InSilicoAssay/'
    shell:
        """
        python3 src/HappyLoop.py --directory output/MaeGz2Pdb
        echo Done > output/log/HappyHB/Done.txt
        """