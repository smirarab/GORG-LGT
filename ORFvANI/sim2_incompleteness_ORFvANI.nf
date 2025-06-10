#!/user/bin/nextflow
nextflow.enable.dsl=2

/*
I. Inputs:

A directory called ./input/sim2_incompleteness/ that contains several tarfiles (e.g. AG-359-G18.tar.gz, AG-390-D15.tar.gz,  etc. ).

Each of these tarfiles is a zipped "experiment" folder.
There are 29 experiment tarfiles in all.

II. What is an experiment?

It contains a set of 31 directories, each represented a simulated genome. All 31 genomes in the experiment are based on the same real source genome for that experiment (e.g. "AG-359-G18")
So how do the 31 directories differ?
They each have a different amount of simulated mutations added to them. (reflected by GND from the source genome)
So in the experiment for source SAG "AG-359-G18", there will be 31 directories named like mutated_BLOSUM62_AG-359-G18_a5_gnd00000, mutated_BLOSUM62_AG-359-G18_a5_gnd001/, etc)
The suffix "gnd00000", "gnd001" etc refer to GND from the source genome, a.k.a. the extent of simulated mutations.
Rather than contain a single fasta, they each each 6 fastas. 
These fastas are Simulated Incomplete Genomes (SIGs).

III. Okay, what and where are the Simulated Incomplete Genomes?

Let me give you an example of SIG fastas:

./input/
    sim2_incompleteness/
        AG-359-G18/
            mutated_BLOSUM62_AG-359-G18_a5_gnd001/
                mutated.fasta
                mutated_p17_s1.fasta
                mutated_p12_s2.fasta
                mutated_p3_s4.fasta
                mutated_p13_s3.fasta
                mutated_p13_s5.fasta


Each ".fasta" above is a SIG, except for "mutated.fasta", which is the original simulated genome (i.e. It has no incompleteness simulated through random genome removal.)
The naming syntax is:
     mutated_p[incompleteness]_s[replicate].fasta

So "mutated_p17_s1.fasta" is the first replicate, with a randomly-drawn incompletness of 17.
Since incompleteness (p) is randomly drawn from a distribution, two different replicates can have the same p value. But those simulated incomplete genomes (SIGs) will still be different since the random seeds (sites for genome removal) are specified by the replicate number.
But if p=0, no regions have been removed from the genome, so the output is identical to the input genome.


IV. So what does this pipeline do?

For a given experiment:
1. Estimates genomic completeness of each genome using checkM1 (and appends completeness value to the fasta name)
2. Tabulates the pairwise ANI (average nucleotide identity); i.e. between any two given genomes in that experiment.
4. Finds ORFs using Prokka
5. Uses BLAST between ORFs to see whether homology is still detectable between mutated genomes, as the mutation levels increase.
6. Integrates this homology info with ANI, resulting in a table that will let you address the above question.
So in total, get 29 tables, for the 29 experiments.

IV. What command did I use to run this nextflow job?
nextflow run sim2_incompleteness_ORFvANI.nf

BTW: ask me (Greg Gavelis, ggavelis@gbigelow.org) for more information if you want to run nextflow yourself.
*/

params.dev = false // Lets user testrun nextflow command (by adding flag '--dev') which will have this pipeline run on JUST ONE SAG (again, as a test)
params.num_inputs = 1

params.outdir = "results"

workflow {
    CH_SAGID_and_TARBALL = Channel.fromPath("./input/sim2_incompleteness/*.tar.gz")
            .take( ( params.dev ? params.num_inputs : -1) )
            .flatten()                                                                                      // Emit each demultiplexed fastq as its own object. E.g. blahblah_R1.fastq.gz
            .map { file -> tuple(file.simpleName, file) }
            //.view()

    // Unzip each experiment (an experiment folder is based on a source genome with 31 permutations (i.e. 31 fastas.) Those permutations are genomes that have been algorithmically "mutated" to have various levels of general nucleotide distance 'gnd' from the source genome.)
    // E.g. for experiment named like "AG-359-G18_a5", there will be fastas like "AG-359-G18_a5_gnd00000", "AG-359-G18_a5_gnd00005", etc., where 'gnd' represents the average general nucleotide distance from the source genome.
    UNZIP_AND_RENAME_FASTAS(CH_SAGID_and_TARBALL)

    CH_eachFasta = UNZIP_AND_RENAME_FASTAS.out.fastas.toSortedList().flatten().take( ( params.dev ? params.num_inputs : -1) )

    CHECKM_v1_1_9(CH_eachFasta)
    APPEND_CHECKM_COMPLETENESS_TO_FASTA_FILENAME(CHECKM_v1_1_9.out)

    CH_ID_AND_FASTA = APPEND_CHECKM_COMPLETENESS_TO_FASTA_FILENAME.out.map{ it -> tuple(it.simpleName, it)} //.view() // e.g. emits [AG-390-D15_gnd007_p0_s1_c89, /path/to/AG-390-D15_gnd007_p0_s1_c89.82.fasta]

    // CH_allFastasSameSAGID = UNZIP_AND_RENAME_FASTAS.out.fastas
    CH_allFastasSameSAGID = APPEND_CHECKM_COMPLETENESS_TO_FASTA_FILENAME.out.map{ it -> tuple(it.simpleName.toString().replaceAll(/_gnd(.+)/,''), it) }.groupTuple() //.view() //derives the experiment ID by dropping filename characters starting at "_gnd" (e.g. "AG-359-G18_a5_gnd00000.fasta" -> "AG-359-G18_a5")

    // Within an experiment, calculate the pairwise ANI between each of the 31 mutated genomes 
    PYANI(CH_allFastasSameSAGID)
    TABULATE_ANI(PYANI.out)

    // Use prokka to find the ORFs in each genome -> .ffn files.
    FIND_ORFS_WITH_PROKKA_v1_14_6(CH_ID_AND_FASTA)
    RENAME_ORFS(FIND_ORFS_WITH_PROKKA_v1_14_6.out)

    // Group the ORF '.fa' files back together by their shared experiment ID (186 fastas per experiment)
    CH_allORFsSameSAGID = RENAME_ORFS.out.map{ it -> tuple(it.simpleName.toString().replaceAll(/_gnd(.+)/,''), it) }.groupTuple() //.view() //derives the experiment ID by dropping filename characters starting at "_gnd" (e.g. "AG-359-G18_a5_gnd00000.fasta" -> "AG-359-G18_a5")

    // In each experiment, count up all 186 genomes' worth of ORFS (To ask: Does number of detected ORFs vary across the 186 genome permutations?)
    COUNT_ORFS(CH_allORFsSameSAGID)

    // In each experiment, BLAST those ORFs against eachother (To ask: Is homology still detectable?)
    BLAST_ORFS_TO_SELF(CH_allORFsSameSAGID)

    CH_FINAL_TABLES = BLAST_ORFS_TO_SELF.out.join(TABULATE_ANI.out).join(COUNT_ORFS.out)
    SUMMARIZE_EXPERIMENT(CH_FINAL_TABLES)
    }


process UNZIP_AND_RENAME_FASTAS {
    tag "Unzipping and organizing all Simulated Incomplete Genomes from experiment ${SAGID}"
    errorStrategy = 'terminate'
    beforeScript 'module load anaconda3/2019.07'
    input: tuple val(SAGID), path(TARBALL)
    output:
        tuple val(SAGID), path("${SAGID}"), emit: groupdir
        path("*fasta"), emit: fastas
    shell:
    '''
    echo "unpacking files"
    tar xvzfk !{TARBALL}

    echo "DONE UNPACKING"

    for SUBDIR in ./!{SAGID}/*/ ; do

        # Get GND, e.g. "gnd00000"
        GND=$(echo "$SUBDIR" | rev | cut -d'_' -f 1 | rev) # Keeping substring after the last  "_" in the subdirectorys name
        GND=${GND::-1} #Remove slash

        # For fasta just called "mutated", append "p0_s0" to name to reflect that it has 0 incompleteness and is the 0th replicate
        mv ${SUBDIR}/mutated.fasta ${SUBDIR}/mutated_p0_s0.fasta

        for FASTA_FULLPATH in ${SUBDIR}/*fasta ; do
            # Get fasta basename, e.g. "mutated_p0_s0.fasta"
            FASTA=$(echo "$FASTA_FULLPATH" | rev | cut -d'/' -f 1 | rev ) # Keeping substring after the last backslash in the fastas path

            # Get incompleteness e.g. "p27"
            INCOMP=$(echo "${FASTA::-6}" | cut -f2 -d'_')

            # Get replicate e.g. "a2"
            ITER=$(echo "${FASTA::-6}" | cut -f3 -d'_')
            echo "${GND}_${INCOMP}_${ITER}"

            # Copy over fasta, giving it more informative name. e.g. "AG-255-A01_gnd00000_p27_a2.fasta"
            cp ${FASTA_FULLPATH} ./!{SAGID}_${GND}_${INCOMP}_${ITER}.fasta
        done
    
    done
    ''' }

process CHECKM_v1_1_9 {
    errorStrategy 'ignore'
    memory = '70.GB'
    container='docker://quay.io/biocontainers/checkm-genome:1.1.9--pyhdfd78af_0'
    input: path(fasta)
    output: tuple path(fasta), path("checkm_results.tsv")
    """
    mkdir tmp_dir; cp ${fasta} ./tmp_dir/input.fasta
    checkm lineage_wf -f checkm_results.tsv --tab_table -q -x fasta -t ${task.cpus} tmp_dir output_dir
    rm -r tmp_dir
    rm -r output_dir
    """ }

process APPEND_CHECKM_COMPLETENESS_TO_FASTA_FILENAME{
    errorStrategy = 'terminate'
    container = 'brwnj/kmernorm:v1.0.0'
    input: tuple val(fasta), path(checkm_tsv)
    output: path("*_c*fasta"), emit: fasta
    """
    #!/usr/bin/env python
    import pandas as pd
    import shutil
    import os

    # Read checkM output table to get completeness value
    STR_completeness = str(pd.read_csv("${checkm_tsv}", sep='\t').loc[0,'Completeness']).replace('.','~')

    STR_fasta_basename = os.path.basename("${fasta}").replace('.fasta','') # e.g. "AG-359-G18_gnd00000_p0_s0"

    STR_new_name = STR_fasta_basename = STR_fasta_basename + "_c" + STR_completeness + ".fasta" # e.g. "AG-359-G18_gnd00000_p0_s0_c88.86.fasta"

    # Copy over fasta with new name. E.g. AG-359-G18_gnd00000_p0_s0.fasta -> AG-359-G18_gnd00000_p0_s0_c88.86.fasta
    shutil.copy("${fasta}", STR_new_name)
    """ }


process PYANI {
    tag "calculate ANI between each genome permutation (all 186 mutated genomes the 31 GND subfolders of experiment ${ID})"
    cpus = 72
    memory = { 50.GB * task.attempt }
    errorStrategy = 'finish'
    container='docker://quay.io/biocontainers/pyani:0.2.9--pyh24bf2e0_0'
    publishDir "results/raw/${ID}/", mode: "copy"
    input: tuple val(ID), path(fastas)
    output: tuple val(ID), path("${ID}_pyani")
    """
    average_nucleotide_identity.py -o "${ID}_pyani" -i ./ -m ANIb --workers "${task.cpus}" --nocompress
    rm -r "${ID}_pyani/blastn_output"
    """ }

process TABULATE_ANI {
    tag "Distill 4 pyANI tables into 1 for experiment ${ID}"
    errorStrategy = 'finish'
    publishDir "results/raw/${ID}/", mode: "copy"
    beforeScript 'module load anaconda3/2019.07' // environment with pandas
    input: tuple val(ID), path("DIR_ani")
    output: tuple val(ID), path("${ID}_ani.csv")
    """
    #!/usr/bin/env python

    ## What it does: Outputs a table with...
    ## * 1 row per genome-pair
    ## * 5 important columns: Genome A, Genome B, ani, ani_aln_cov_ab, ani_aln_cov_ba, comp

    PATH_out = "${ID}_ani.csv"

    PATH_pyani_identity = "${DIR_ani}/ANIb_percentage_identity.tab"
    PATH_pyani_lengths = "${DIR_ani}/ANIb_alignment_lengths.tab"
    PATH_pyani_coverage = "${DIR_ani}/ANIb_alignment_coverage.tab"
    PATH_pyani_simerrors = "${DIR_ani}/ANIb_similarity_errors.tab"
    
    import pandas as pd
    
    min_coverage = 0.03
    
    TAB = "\t"
    
    print("Loading pyani results.")
    load_matrix_table = lambda metric, table: pd.read_csv(table, sep=TAB).melt(id_vars='Unnamed: 0', var_name='Genome B', value_name=metric).rename(columns={'Unnamed: 0':'Genome A'})
    idf = load_matrix_table('ani', PATH_pyani_identity)
    ldf = load_matrix_table('ani_aln_len', PATH_pyani_lengths)
    cdf = load_matrix_table('ani_aln_cov', PATH_pyani_coverage)
    edf = load_matrix_table('ani_similarity_errs', PATH_pyani_simerrors)
    
    bdf = idf.merge(ldf).merge(cdf).merge(edf)
    bdf = bdf[bdf['Genome A'] != bdf['Genome B']]
    
    print(str(len(bdf))+" genome pairs were compared")
    
    ### Make a table "comp" to show what is being compared.
    
    comps = []
    
    for i, l in bdf.iterrows():
        lst = [l['Genome A'], l['Genome B']]
        lst.sort()
        string = "{}_{}".format(lst[0], lst[1])
        comps.append(string)
    
    bdf['comp'] = comps
    bdf = bdf.sort_values(by=['comp','ani'], ascending=False)
    
    dfa = bdf.drop_duplicates(subset=['comp'], keep='first')
    dfb = bdf.drop_duplicates(subset=['comp'], keep='last')
    
    dfa.rename(columns={'ani_aln_cov':'ani_aln_cov_ab', 'ani_aln_len':'ani_aln_len_ab', 
                           'ani_similarity_errs':'ani_similarity_errs_ab'}, inplace=True)
    dfb.rename(columns={'ani_aln_cov':'ani_aln_cov_ba', 'ani_aln_len':'ani_aln_len_ba', 
                           'ani_similarity_errs':'ani_similarity_errs_ba'}, inplace=True)
    
    dfab = dfa.merge(dfb[['comp', 'ani_aln_cov_ba','ani_aln_len_ba','ani_similarity_errs_ba']], on='comp', how='outer')
    print(str(len(bdf) - len(dfab))+" of the comparisons dropped as being redundant (reciprocal)")
    
    DF = dfab[(dfab['ani_aln_cov_ba'] > min_coverage) & (dfab['ani_aln_cov_ab'] > min_coverage)]
    print(str(len(DF)) + ' ANI estimates had over '+str(min_coverage)+' coverage for both compared genomes')

    DF['Genome A'] = DF['Genome A'].str.replace("~",".")
    DF['Genome B'] = DF['Genome B'].str.replace("~",".")

    DF.to_csv(PATH_out, index=False)
    """ }

process FIND_ORFS_WITH_PROKKA_v1_14_6 {
    tag "Output ORFs (as fnn file) of ${ID} (which is a mutated genome from experiment ${ID})"
    cpus = 6
    memory = { 10.GB * task.attempt }
    errorStrategy = 'finish'
    container='docker://quay.io/biocontainers/prokka:1.14.6--pl5262hdfd78af_1'
    input: tuple val(ID), path(fasta)
    output: path("${ID}.ffn")
    """
    prokka --outdir prokka_${ID} --quiet --compliant --force --cpus ${task.cpus} ${fasta}
    # prokka --outdir prokka_${ID} --prefix ${ID} --locustag ${ID} --quiet --compliant --force --cpus ${task.cpus} ${fasta}
    mv prokka_${ID}/*ffn ./${ID}.ffn
    """ }

process RENAME_ORFS {
    errorStrategy = 'terminate'
    container = 'brwnj/kmernorm:v1.0.0'
    input: path(fasta)
    output: path("*.fa")
    script:
    """
    #!/usr/bin/env python

    from Bio import SeqIO
    import os

    PATH_in="${fasta}"

    # Extract ID from basename of .ffn file
    ID = os.path.basename(PATH_in).replace(".ffn","").replace("~",".")
    PATH_out = os.path.basename(PATH_in).replace(".ffn",".fa")
    
    LIST_records = []
    for record in SeqIO.parse(PATH_in, "fasta"):
        OLD_id = str(record.id)
        NUM_contig = int(OLD_id.split('_')[-1])
    
        NEW_id = ID + "_n"+str(NUM_contig)
        record.id = NEW_id
        record.description = ''
    
        LIST_records.append(record)
    
    SeqIO.write(LIST_records, PATH_out, 'fasta')
    """ }

process COUNT_ORFS {
    tag "for experiment ${ID} make a table with the following columns: sag, gene_count, gene_count_500-1500bp (summarizing all 31 mutated genomes)"
    errorStrategy = 'finish'
    beforeScript 'module load anaconda3/2019.07' // python environment
    publishDir "results/raw/${ID}/", mode: "copy"
    input: tuple val(ID), path(FILES)
    output:
        tuple val(ID), path("${ID}_orf_counts.csv")
    """
    #!/usr/bin/env python

    PATH_out = "${ID}_orf_counts.csv"
    
    from glob import glob
    import os
    import itertools
    
    length_threshold = [500, 1500]
    LIST_ffn = glob("*.fa")
    
    def readfa(fh):
        for header, group in itertools.groupby(fh, lambda line: line[0] == '>'):
            if header:
                line = next(group)
                name = line[1:].strip()
            else:
                seq = ''.join(line.strip() for line in group)
                yield name, seq
    
    with open(PATH_out, 'w') as oh:
        print('sag','gene_count','gene_count_{}-{}bp'.format(length_threshold[0], length_threshold[1]), sep = ',', file = oh)
        for ffn in LIST_ffn:
            sag = os.path.basename(ffn).split(".")[0]
            gene_count = 0
            thresh_count = 0
            for name, seq in readfa(open(ffn)):
                gene_count += 1
                gene_length = len(seq)
                if length_threshold[0] <= gene_length and length_threshold[1] >= gene_length:
                    thresh_count += 1
            #print(sag, gene_count, thresh_count)
            print(sag, gene_count, thresh_count, sep = ',', file = oh)
    """ }

process BLAST_ORFS_TO_SELF {
    tag "All v. all BLAST for ORFs from experiment ${ID} (across all 186 mutated genomes)"
    cpus = 72
    memory = { 50.GB * task.attempt }
    errorStrategy = 'finish'
    publishDir "results/", mode: "copy"
    container='docker://quay.io/biocontainers/blast:2.7.1--h4422958_6'
    input: tuple val(ID), path(FILES)
    output:
        tuple val("${ID}"), path("${ID}_blast.tsv")
    
    shell:
    """
    cat *.fa > "!{ID}_orfs.fasta"

    makeblastdb -in "!{ID}_orfs.fasta" -dbtype nucl

    echo "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen" > "!{ID}_blast.tsv"

    blastn -query "!{ID}_orfs.fasta" \
    -db "!{ID}_orfs.fasta" \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
    -num_threads !{task.cpus} -max_target_seqs 2000 \
    -out headerless_blast.tsv
    
    cat headerless_blast.tsv >> "!{ID}_blast.tsv"
    """ }

process SUMMARIZE_EXPERIMENT {
    tag "Collate results into single table for experiment ${ID}"
    errorStrategy = 'finish'
    memory = '100.GB'
    publishDir "results/raw/", mode: "copy"
    beforeScript 'module load anaconda3/2019.07'  // environment with pandas
    input: tuple val(ID), path(BLAST_tsv), path(ANI_csv), path(ORF_csv)
    output:
        tuple val(ID), path("${ID}_results.csv"), emit: tuple
        path("${ID}_results.csv"), emit: csv
    """
    #!/usr/bin/env python

    PATH_blast = "${BLAST_tsv}"
    PATH_ani_table = "${ANI_csv}"
    PATH_orf_table = "${ORF_csv}"
    
    PATH_out = "${ID}_results.csv"
    
    TAB = "\t"
    
    import pandas as pd
    import numpy as np
    
    LIST_pid = [99.9, 99.8, 99.7, 99.6, 99.5, 99, 98.5, 98, 97.5, 97]
    len_threshold = [500, 1500]
    min_coverage = 0.03
    
    ### I. Load and filter BLAST results _______________________
    print("I. Loading BLAST results.")
    DF_blast = pd.read_csv(PATH_blast, sep=TAB)
    DF_blast.insert(0, "qsag", DF_blast["qseqid"].str.rsplit("_",n=1).str.get(0))
    DF_blast.insert(3, "ssag", DF_blast["sseqid"].str.rsplit("_",n=1).str.get(0))
    
    print(TAB+str(len(DF_blast["qsag"].unique())) + " SAGs were BLASTed")
    
    LIST_sags = DF_blast["qsag"].unique().tolist()
    
    DF = DF_blast
    print(TAB+str(len(DF_blast)) + " initial BLAST hits")
    
    # Exclude BLAST self-hits
    DF = DF[DF['qsag'] != DF['ssag']]
    print(TAB+"Dropping "+str(len(DF_blast)-len(DF))+" self-hits")
    
    print(TAB+"Dropping "+str(len(DF) - len(DF[(DF['length'] > (0.8 * DF['qlen'])) | (DF['length'] > 0.8 * DF['slen'])]))+ " hits where less than 80% of the query or 80% of the subject seq are aligned")
    DF = DF[(DF['length'] > (0.8 * DF['qlen'])) | (DF['length'] > 0.8 * DF['slen'])]
    
    ### II. DF_pairwise, a table of BLAST-inferred orthologs shared by pairs of SAGs ___________________________________
    
    DF_pairwise = pd.DataFrame(columns = ['qsag','ssag'])
        
    ni_alpha_mean = DF.groupby(['qsag','ssag'], as_index = False)['pident'].mean().rename(columns = {'pident':'mean_pident'}).merge(
                    DF.groupby(['qsag','ssag'])['pident'].std().reset_index().rename(columns = {'pident':'stdev_pident'})).merge(
                    DF.groupby(['qsag','ssag'], as_index = False)['pident'].count().rename(columns = {'pident':'total_hits'})).merge(
                    DF.groupby(['qsag','ssag'], as_index = False)['pident'].median().rename(columns = {'pident':'median_pident'}))
    
    print("II. Integrating ORF counts per SAG")
    DF_orf_count = pd.read_csv(PATH_orf_table)
    DF_orf_count['sag'] = DF_orf_count['sag'].str.replace("~",".")

    DF_pairwise = DF_pairwise.merge(ni_alpha_mean, how = 'outer').merge(
                DF_orf_count.rename(columns = {i:'s{}'.format(i) for i in DF_orf_count.columns}), how = 'left').merge(
                DF_orf_count.rename(columns = {i:'q{}'.format(i) for i in DF_orf_count.columns}), how = 'left')
    
    print(TAB+"Converting to 'DF_pairwise', which quantifies the orthologs detected across "+str(len(DF_pairwise))+" pairs of SAGs")
    
    ### How many orthologs persist across increasing percent identities? 
    
    for pid in LIST_pid:
            
        len_bound = DF[(DF['pident'] >= pid) & 
                ((DF['qlen'] >= len_threshold[0]) & (DF['qlen'] <= len_threshold[1])) & 
                ((DF['slen'] >= len_threshold[0]) & (DF['slen'] <= len_threshold[1]))]\
                .groupby(['qsag','ssag'], as_index=False)['sseqid']\
                .count().rename(columns = {'sseqid':'{}_pid_{}-{}bp_orthologs'.format(pid, len_threshold[0],len_threshold[1])})
        all_sizes = DF[DF['pident'] >= pid].groupby(['qsag','ssag'], as_index=False)['sseqid']\
                .count().rename(columns = {'sseqid':'{}_pid_orthologs'.format(pid, len_threshold[0],len_threshold[1])})
        
        DF_pairwise = DF_pairwise.merge(
                len_bound, how = 'outer').merge(
                all_sizes, how = 'outer')
    
    print(TAB+"Calculating orthologs across "+str(len(LIST_pid))+" percids ranging from " + str(LIST_pid[0])+ " to " +str(LIST_pid[-1]))
    
    # add a pairs column to have a unique identifier per genome pair:
    pairs = []
    for i, l in DF_pairwise.iterrows():
        lst = [l['qsag'], l['ssag']]
        lst.sort()
        pairs.append("_".join(lst))
    DF_pairwise['pair'] = pairs
    DF_pairwise['pair'] = DF_pairwise['pair'].str.replace("~",".")
    
    ortho_names = [i for i in DF_pairwise.columns if 'orthologs' in i]
    
    # Drop reciprocal ANI hits as they are redundant
    DF_pairwise = DF_pairwise.drop_duplicates(subset = 'pair')
    
    ### III. Combine ortholog info (DF_pairwise) with genome ANI (DF_ani) -> DF_final
    
    print("III. Loading ANI data.")
    DF_ani = pd.read_csv(PATH_ani_table)
    #DF_ani = pd.read_csv(PATH_ani_table, index_col=0)
    DF_ani = DF_ani.rename(columns = {'comp':'pair'})
    DF_ani = DF_ani[['Genome A','Genome B','ani', 'ani_aln_cov_ab','ani_aln_cov_ba', 'pair']].dropna() ### Keep only essential columns
    DF_ani['pair'] = DF_ani['pair'].str.replace("~",".")

    print(TAB+"Incorporating ANI data. Deriving GND.")
    DF_final = DF_pairwise.merge(DF_ani, on = 'pair', how = 'outer')
    DF_final['GND'] = [round(100 - i * 100, 2) for i in DF_final['ani']]
    
    # include all pairs for which ANI was measured in summary calculations:
    DF_final[ortho_names] = DF_final[ortho_names].fillna(0)
    
    tups = [(round(i-0.01, 2), round(i, 2),) for i in np.arange(1, 0.6, -0.01)]
    bins = pd.IntervalIndex.from_tuples(tups)
    DF_final['ani_bin'] = pd.cut(DF_final['ani'], bins)
    print(TAB+"Binning ANI across "+str(len(bins))+" increments.")
    
    # write to csv:
    DF_final.to_csv(PATH_out, index=False)
    print("Done.")
    """ }






