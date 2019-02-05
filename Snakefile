rule all:
    input:
        auspice_tree = "auspice/enterovirus_d68_{seg}_tree.json",
        auspice_meta = "auspice/enterovirus_d68_{seg}_meta.json"


rule files:
    params:
        raw_vipr = "data/genomeEntero-30Jan19.tsv", #raw VIPR download!
        
        #samples sequenced in sweden
        swedish_seqs = "data/ev_d68_genomes_2018.fasta",
        swedish_meta = "data/20190128_Karolinska-region.csv",
        #samples added manually - not from ViPR (maybe found on genbank)
        manual_seqs = "data/manual-seqs.fasta",
        manual_meta = "data/manual-meta.csv",

        dropped_strains = "config/dropped_strains.txt",
        reference = "config/ev_d68_reference_genome.gb",
        colors = "config/colors.tsv",
        clades = "config/clades.tsv",
        auspice_config = "config/auspice_config.json",
        regions = "config/geo_regions.tsv"

files = rules.files.params

import os.path
RERUN = True if os.path.isfile("genbank/current_vipr_download.tsv") else False

#always parse new genbank files - add regions if wanted
#unlike the vp1 script, this doesn't make a new strain name!!
rule parse_vipr_meta:
    input:
        meta = files.raw_vipr,
        regions = ancient(files.regions)
    output:
        "temp/current_vipr_download.tsv"
    message:
        "This rerun will use existing GenBank files! Only new accession numbers will be downloaded" if RERUN else "Starting new run from scratch. All VIPR samples will be downloaded."
    shell:
        "python vipr_parse_genome.py --input {input.meta} --output {output} --regions {input.regions}"


#if new run, ensure to exclude any sequences we already have in 'swedish' or 'manual'
rule remove_dupes:
    input:
        genbank_meta = rules.parse_vipr_meta.output,
        swed_meta = ancient(files.swedish_meta), #do not rerun if other meta changes - won't influence genbank!
        man_meta = ancient(files.manual_meta),
    output:
        "temp/meta_to_download.tsv"
    shell:
        "python find_new.py --input-new {input.genbank_meta} --exclude {input.swed_meta} {input.man_meta} --output {output}"


#if rerun, find only the new ones to download - exclude old, swedish & manual samples we already have
rule find_new:
    input:
        old_meta = ancient("genbank/current_vipr_download.tsv"),
        swed_meta = ancient(files.swedish_meta), #do not rerun if other meta changes - won't influence genbank!
        man_meta = ancient(files.manual_meta),
        new_meta = rules.parse_vipr_meta.output
    output:
        "temp/new_meta_to_download.tsv"
    shell:
        """
        python find_new.py --input-new {input.new_meta} \
            --exclude {input.swed_meta} {input.old_meta} {input.man_meta} \
            --output {output}
        """


#download only new and non-duplicate sequences
rule download_seqs:
    input:
        meta = rules.find_new.output[0] if RERUN else rules.remove_dupes.output[0]
    output:
        sequences = "temp/add_sequences_genome.fasta" if RERUN else "temp/genbank_genome_sequences.fasta", 
        meta = "temp/add_meta_genome.tsv" if RERUN else "temp/genbank_genome_meta.tsv"
    run:
        import pandas as pd
        from Bio import Entrez, SeqIO
        from augur.parse import forbidden_characters
        Entrez.email = "richard.neher@unibas.ch"

        print(input.meta)
        meta = pd.read_csv(input.meta, sep='\t')
        additional_meta = {}
        tooShort = []
        with open(output.sequences, 'w') as fh:
            for ri, row in meta.iterrows():
                try:
                    handle = Entrez.efetch(db="nucleotide", id=row.accession, rettype="gb", retmode="text")
                except:
                    print(row.accession, "did not work")
                    continue
                print(row.strain, row.accession)
                rec = SeqIO.read(handle, 'genbank')
                if len(rec.seq) - rec.seq.count("N") < 6400:
                    print(row.strain, row.accession, "is too short when Ns removed!")
                    tooShort.append(row.strain)
                    continue
                try:
                    authors = rec.annotations['references'][0].authors
                    title = rec.annotations['references'][0].title
                except:
                    authors = ''
                    title = ''

                url = 'https://www.ncbi.nlm.nih.gov/nuccore/'+row.accession
                additional_meta[ri] = {'url':url, 'authors':authors, 'title':title}
                tmp = row.strain
                for c,r in forbidden_characters:
                    tmp=tmp.replace(c,r)
                rec.id = tmp
                rec.name = tmp
                rec.description = ''
                SeqIO.write(rec, fh, 'fasta')

        print(len(tooShort), "sequences were too short after Ns were removed, and were excluded.")
        add_meta = pd.DataFrame(additional_meta).transpose()
        all_meta = pd.concat((meta, add_meta), axis=1)
        all_meta.to_csv(output.meta, sep='\t', index=False)


#concatenate meta and sequences to existing genbank, if a rerun:
#concat meta
rule add_meta:
    input:
        metadata = [ancient("genbank/genbank_genome_meta.tsv"), rules.download_seqs.output.meta]
    output:
        metadata = "temp/genbank_genome_meta.tsv"
    run:
        import pandas as pd
        from augur.parse import fix_dates, forbidden_characters
        md = []
        for fname in input.metadata:
            tmp = pd.read_csv(fname, sep='\t' if fname.endswith('tsv') else ',')
            tmp_name = []
            for x in tmp.strain:
                f = x
                for c,r in forbidden_characters:
                    f=f.replace(c,r)
                tmp_name.append(f)
            tmp.strain = tmp_name
            md.append(tmp)
        all_meta = pd.concat(md)
        all_meta.to_csv(output.metadata, sep='\t', index=False)
#concat sequences
rule add_sequences:
    input:
        ancient("genbank/genbank_genome_sequences.fasta"), rules.download_seqs.output.sequences
    output:
        "temp/genbank_genome_sequences.fasta"
    shell:
        '''
        cat {input} > {output}
        '''

#if all has gone well up to here, copy relevant files to allow rerun next time!
rule make_database:
    input:
        gen_seqs = "temp/genbank_genome_sequences.fasta",
        gen_meta = "temp/genbank_genome_meta.tsv",
        download = "temp/current_vipr_download.tsv",
    output:
        gen_seqs = "genbank/genbank_genome_sequences.fasta",
        gen_meta = "genbank/genbank_genome_meta.tsv",
        download = "genbank/current_vipr_download.tsv",
    message:
        "Genbank files updated with new sequences!" if RERUN else "Genbank files stored. Reruns will only download new accession numbers."
    shell:
        '''
        cp {input.gen_seqs} genbank
        cp {input.gen_meta} genbank
        cp {input.download} genbank
        '''

#concatenate genbank meta with Swedish & manual
#we want this to run if changes to swedish/manual data, even if there aren't new genbank!
rule concat_meta:
    input:
        metadata = [files.swedish_meta, files.manual_meta, "genbank/genbank_genome_meta.tsv"]
    output:
        metadata = "results/metadata.tsv"
    run:
        import pandas as pd
        from augur.parse import fix_dates, forbidden_characters
        md = []
        for fname in input.metadata:
            tmp = pd.read_csv(fname, sep='\t' if fname.endswith('tsv') else ',')
            tmp_name = []
            for x in tmp.strain:
                f = x
                for c,r in forbidden_characters:
                    f=f.replace(c,r)
                tmp_name.append(f)
            tmp.strain = tmp_name
            md.append(tmp)
        all_meta = pd.concat(md)
        all_meta.to_csv(output.metadata, sep='\t', index=False)

#concatenate genbank seqs with Swedish & manual
rule concat_sequences:
    input:
        files.swedish_seqs, files.manual_seqs, "genbank/genbank_genome_sequences.fasta"
    output:
        "results/sequences.fasta"
    shell:
        '''
        cat {input} > {output}
        '''


##############################
# now run usual augur analysis
###############################

#Nextstrain run starts here!

rule filter:
    input:
        sequences = rules.concat_sequences.output,
        metadata = rules.concat_meta.output.metadata,
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered.fasta"
    params:
        sequences_per_category = 200,
        categories = "country year month",
        min_date = 1980
    shell:
        """
        augur filter --sequences {input.sequences} --metadata {input.metadata} \
            --output {output.sequences} \
            --group-by {params.categories} \
            --sequences-per-group {params.sequences_per_category} \
            --exclude {input.exclude}  --min-date {params.min_date}
        """

rule align:
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_genome.fasta"
    shell:
        """
        augur align --sequences {input.sequences} --output {output.alignment} \
            --reference-sequence {input.reference} --fill-gaps
        """

rule sub_alignments:
    input:
        rules.align.output
    output:
        subaln = ["results/aligned_vp4vp1.fasta", "results/aligned_p2p3.fasta"]
    params:
        boundaries = [(732,3315), (3315,7296)]
    run:
        from Bio import AlignIO
        aln = AlignIO.read(input[0], 'fasta')
        print(aln)
        for fname,b in zip(output.subaln, params.boundaries):
            print(b)
            AlignIO.write(aln[:,b[0]:b[1]], fname, 'fasta')


rule tree:
    input:
        alignment = [rules.align.output.alignment] + rules.sub_alignments.output
    output:
        tree = "results/raw_tree_{seg}.nwk"
    shell:
        """
        augur tree --alignment results/aligned_{wildcards.seg}.fasta --output {output.tree}
        """

rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = 'results/aligned_{seg}.fasta',
        metadata = rules.concat_meta.output.metadata,
    output:
        tree = "results/tree_{seg}.nwk",
        node_data = "results/branch_lengths_{seg}.json"
    params:
        clock_filter_iqd = 5
    shell:
        """
        augur refine --tree {input.tree} --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} --output-node-data {output.node_data} \
            --timetree --date-confidence --date-inference marginal --coalescent opt \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = 'results/aligned_{seg}.fasta',
    output:
        nt_data = "results/nt_muts_{seg}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral --tree {input.tree} --alignment {input.alignment} \
            --output {output.nt_data} --inference {params.inference}
        """

rule translate:
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.nt_data,
        reference = 'config/ev_d68_reference_{seg}.gb'
    output:
        aa_data = "results/aa_muts_{seg}.json"
    shell:
        """
        augur translate --tree {input.tree} --ancestral-sequences {input.node_data} \
            --output {output.aa_data} --reference-sequence {input.reference}
        """

rule clades:
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.aa_data,
        nuc_muts = rules.ancestral.output.nt_data,
        clades = files.clades
    output:
        clade_data = "results/clades_{seg}.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output {output.clade_data}
        """

rule traits:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.concat_meta.output.metadata,
    output:
        node_data = "results/traits_{seg}.json",
    params:
        columns = "country"
    shell:
        """
        augur traits --tree {input.tree} --metadata {input.metadata} \
            --output {output.node_data} --confidence --columns {params.columns}
        """

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.concat_meta.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.nt_data,
        aa_muts = rules.translate.output.aa_data,
        colors = files.colors,
        clades = rules.clades.output.clade_data,
        auspice_config = files.auspice_config
    output:
        auspice_tree = "auspice/enterovirus_d68_{seg}_tree.json",
        auspice_meta = "auspice/enterovirus_d68_{seg}_meta.json"
    shell:
        """
        augur export --tree {input.tree} --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.clades}\
            --colors {input.colors} --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} --output-meta {output.auspice_meta}
        """
