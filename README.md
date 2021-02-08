# sarbecovirus-partial-rdrp-phylogeny

Nextstrain tree of Sarbecovirus partial RdRp sequences.

The names of the sequences will eventually be improved, keeping the accession is to make edition of sequences.fasta easier, the character | is chosen to be easily removed.

Once augur is install (pip install augur)

Snamemake --cores 1 

generates the file auspice/rdrp.json

Copying auspice's files plus rdrp.json gives

http://babarlelephant.free-hoster.net/dist/index_rdrp.html?c=clade_membership


The tree is supposed to contain almost all available partial RdRp sequences (many strains have only their partial RdRp sequenced). The main methodology was to blast some partial RdRp sequences to obtain the Sarbecoviruses accessions in the results. A Sarbecovirus is anything with more than 90% alignment in this region to one of the main lineages.

