# GIANT
GIANT (Gene-based data Integration and ANalysis Technique) is a method for unified analysis of atlas-level single cell data. GIANT generates a unified gene embedding space across multiple data modalities and tissues. GIANT first constructs gene graphs for cell clusters of each tissue from each data modality. A dendrogram is then built to connect the gene graphs in a hierarchy. GIANT next combines the gene graphs and the dendrogram to recursively embed genes in the graphs to a common latent space. Locations of genes in the space reflect their functions in their cell clusters.

## Dependencies
- python 3.9.7</br>
- anndata 0.7.5</br>
- cython 0.29.30</br>
- goatools 1.1.12</br>
- joblib 1.1.0</br>
- networkx 2.6.3</br>
- numpy 1.20.3</br>
- pandas 1.3.4</br>
- pyensembl 1.9.4</br>
- scanpy 1.8.2</br>
- scikit-learn 1.0.1</br>
- scipy 1.7.2</br>
- setuptools 59.2.0</br>

GIANT supports distributed training based on the [Gensim](https://radimrehurek.com/gensim/apiref.html#api-reference) package implemented in Cython. Following the steps to compile the Cython code:
```
cd src/embedding
python setup.py build_ext --inplace
```

## Usage
Use `python script.py --help` to check the specific usage of each command.

### Build gene graphs
- Run the following commands to build co-expression graphs for scRNA-seq data or Slide-seq data:
```
cd src/build_graphs
python coexpression_knn_graph.py --K 10 --datatype "rna" --datadir "path/to/h5ad" --genescopefile "path/to/gene_scope" --outdir "path/to/out_dir"
```

- Run the following commands to build gene-TF hypergraphs for scATAC-seq data:
```
cd src/build_graphs
python atacseq_hypergraph.py --bpupstream 500 --datadir "path/to/peak_files" --genescopefile "path/to/gene_scope" --generangefile "path/to/gene_ranges" --outdir "path/to/out_dir"
```

- Run the following commands to build spatial co-expression hypergraphs for Slide-seq data:
```
cd src/build_graphs
python spatial_coexpression_knn_graph.py --K 10 --datadir "path/to/h5ad" --genescopefile "path/to/gene_scope" --outdir "path/to/out_dir"
```

### Build dendrogram of graphs
- Run the following commands to build dendrogram:
```
cd src/build_dendrogram
python build_dendrogram.py --K 50 --degenedir "path/to/de_genes" --graphdir "path/to/gene_graphs" --outdir "path/to/save_output" --figdir "path/to/save_output_figure"
```

### Learn gene embeddings
- Run the following commands to learn gene embeddings:
```
cd src/embedding
python main.py --input "path/to/graph.list" --outdir "emb" --hierarchy "path/to/graph.hierarchy" --iter 150 --regstrength 1 --workers 8 --dimension 128
```

### Identify enriched functions or gene regulons in different parts of the embedding space
- Run the following commands to split the embedding space into components (cluster genes in the whole space), then run enrichment analysis for each embedding component:
```
cd src/analysis
python cluster_genes.py --emb "dir/to/gene_embedding" --resolution 50 --nneighbors 30 --annotationfile "dir/to/gene_annotations" --outdir "dir/to/out_dir"
python GO_TF_enrichment.py "dir/to/cluster_0.txt" "dir/to/population.txt" "dir/to/gene_human_graph.gaf" --ev_exc=IEA --pval=0.05 --method=fdr_bh --pval_field=fdr_bh
```

### Infer novel gene functions
- Run the following commands to learn gene embeddings:
```
cd src/analysis
python neighbors_in_graph.py --emb "dir/to/gene_embedding" --outdir "dir/to/out_dir"
python GO_enrichment_in_neighbors.py "dir/to/neighbor_gene_list.txt" "dir/to/population.txt" "dir/to/gene_human_graph.gaf" --ev_exc=IEA --pval=0.05 --method=fdr_bh --pval_field=fdr_bh
```

## Credits
The software is an implementation of the method GIANT, jointly developed by Hao Chen, Nam D. Nguyen, Matthew Ruffalo, and Ziv Bar-Joseph from the [System Biology Group @ Carnegie Mellon University](http://sb.cs.cmu.edu/).

## Contact
Contact us if you have any questions:</br>
Hao Chen: hchen4 at andrew.cmu.edu</br>
Ziv Bar-Joseph: zivbj at andrew.cmu.edu</br>

## License
This project is licensed under the MIT License - see the LICENSE file for details.
