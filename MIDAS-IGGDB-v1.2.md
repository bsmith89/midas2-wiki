# Reformat MIDAS DB v1.2 into IGGtools

IGGDB is driven by a TOC with four columns

```
genome     species                        representative   genome_is_representative
1190605.3  Enterovibrio_norvegicus_54866  1190605.3        1
349741.6   Akkermansia_muciniphila_55290  349741.6         1
1408424.3  Bacillus_bogoriensis_60417     1408424.3        1
```

## Phyeco Marker Genes Database

Need to create `marker_centroids` under the `marker_genes` folder, inside which is the ID mapping from the SGC marker genes of the representative genomes to the centroid_99 genes of the pan-genome.


## Representative Genomes Database

Rename the FASTA file into `{species_id}/{genome_id}/{genome_id}.fna`

## Pan Genomes Database

We need `centorids.ffn` and `gene_info.txt`. 


