The file was generated as follows:

On 13 Oct 2015 the website www.ensembl.org/biomart was visited and the following operations were performed:

Select "Ensembl Genes 82" database
Choose dataset "Mus musculus genes (GRCm38.p4)"
Click on "Attributes" section on the left, select "Homologs" tick.
In the "Gene" category, seleect "Ensembl Protein ID", unselect "Ensembl Transcript ID".
In the "ORTHOLOGS (Max select 6 orthologs)" section, tick "Human Ensembl Gene ID" and "Human Ensembl Protein ID" in "Human Orthologs".
Click on "Results" in the upper left corner besides, download results as TSV.


From e-mail Thomas Hume:

One last note, I remember that I actually did this in both ways (selecting either species in the first step) and then merged the results. A few more orthologs would appear that way (probably due to the mess with Gene / Transcripts / Proteins IDs, even though it's supposed to be an all-pairs result...).

======
Networks

from L-GRAAL paper (Noel's website). They are from BioGRID. Reason: largest ones. Identifiers: NCBI Gene IDs.  Would require soem mapping, therefore I don;t do it now.

the data form natalie 2.0: Identifiers: I guess these are ENSEMBL protein IDs. Although they have a number at the beginning, denoting the species. Networks: Natalie2.0 paper: "From the STRING database v8.3 [2], we obtained PPI networks for the following six species:
165 C. elegans (cel), S. cerevisiae (sce), D. melanogaster (dme), R. norvegicus (rno), M. musculus (mmu) and
166 H. sapiens (hsa). We only considered interactions that were experimentally verified."

I removed the species numbers (9606 for hsa and 10090 for mmu) in the networks to make them work with the ENSEMBL orthologs

Also: note that the .gw versions are much smaller.  Reason is that natalie filters the STRING edges when reading in and only considers edges with combined score >= 900.  This has already been done for the .gw graphs




