rule all:
    input:
        auspice_json = "auspice/rdrp.json",

input_fasta = "data/sequences.fasta",
input_metadata = "data/metadata.tsv",
dropped_strains = "config/dropped_strains.txt",
reference = "config/rdrp_outgroup.gb",
colors = "config/colors.tsv",
lat_longs = "config/nextstrain_emma_lat_longs.tsv",
auspice_config = "config/auspice_config.json"
clades = "config/clades.tsv"

rule filter:
    message:
        """
        Filtering to
          - excluding strains in {input.exclude}
        """
    input:
        sequences = input_fasta,
        metadata = input_metadata,
        exclude = dropped_strains,
        reference = reference
    output:
        sequences = "results/filtered.fasta"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --exclude {input.exclude} \
            --output {output.sequences} \
            --metadata {input.metadata}
        """


rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    params:
        root = "internal-node-root"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """
	
rule prerefine:
    message:
        """
        PreRefining tree
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = input_metadata
    output:
        tree = "results/tree_pre.nwk",
        node_data = "results/branch_lengths_pre.json"
    params:
        root = "internal-node-root"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
        """
		
rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.prerefine.output.tree,
        alignment = rules.align.output,
        metadata = input_metadata
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        root = "internal-node-root",
        clock_rate = 0.0004
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --coalescent {params.coalescent} \
			--keep-root \
			--timetree \
            --date-inference {params.date_inference} \
            --keep-polytomies \
            --precision 0 \
            --clock-rate {params.clock_rate}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = "config/rdrp_outgroup.gb"
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """

rule clades:
    message: " Labeling clades as specified in config/clades.tsv"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = clades
    output:
        node_data = "results/clades.json"
    shell:
        """
        augur clades \
            --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output {output.node_data}
        """

		
rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = input_metadata,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        lat_longs = lat_longs,
        auspice_config = auspice_config,
        clades = rules.clades.output.node_data
    output:
        auspice_json = rules.all.input.auspice_json,
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.clades} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
