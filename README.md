# GIANT
GIANT (Gene-based data Integration and ANalysis Technique) is a method for unified analysis of atlas-level single cell data. The method generates a unified gene embedding space across multiple data modalities and tissues. GIANT first constructs gene graphs for cell clusters of each tissue from each data modality. A dendrogram is then built to connect the gene graphs in a hierarchy. GIANT next combines the gene graphs and the dendrogram to recursively embed genes in the graphs to a common latent space. Locations of genes in the space reflect their functions in their cell clusters.

## System requirements
### Operating system
The software has been tested on the CentOS Linux 7 system.

### Software requirements
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

### Installation
It is recommended to create a virtual environment using [Conda](https://conda.io/projects/conda/en/latest/index.html). After successfully installing Anaconda/Miniconda, create an environment using the provided `environment.yml` file:
```
conda env create -f environment.yml
conda activate GIANT-env
```

GIANT supports distributed training based on the [Gensim](https://radimrehurek.com/gensim/apiref.html#api-reference) package implemented in Cython. Following the steps to compile the Cython code:
```
cd src/embedding
python setup.py build_ext --inplace
```

## Usage
Below is the description of command usages on example data. Use `python script.py --help` to check the specific usage of each command.

### Build gene graphs
- Run the following commands to build co-expression graphs for scRNA-seq data or Slide-seq data:
```
cd src/build_graphs
python coexpression_knn_graph.py --K 10 --datatype "rna" --tissue 'Heart' --datadir "../../example_data/secondary_analysis.h5ad" --genescopefile "../../example_data/gene_scope.txt" --outdir "../../graphs/edgelists"
```

- Run the following commands to build gene-TF hypergraphs for scATAC-seq data:
```
cd src/build_graphs
python atacseq_hypergraph.py --bpupstream 500 --tissue "Thymus" --datadir "../../example_data/atac_peaks" --genescopefile "../../example_data/gene_scope.txt" --generangefile "../../example_data/GenesRanges.csv" --outdir "../../graphs/edgelists/"
```

- Run the following commands to build spatial co-expression hypergraphs for Slide-seq data:
```
cd src/build_graphs
python spatial_coexpression_knn_graph.py --K 10 --tissue "Kidney" --datadir "../../example_data/spatial_secondary_analysis.h5ad" --genescopefile "../../example_data/gene_scope.txt" --outdir "../../graphs/edgelists/"
```

### Build dendrogram of graphs
- Run the following commands to build dendrogram:
```
cd src/build_dendrogram
python build_dendrogram.py --K 50 --degenedir "../../example_data/DE_genes/" --graphdir "../../graphs/edgelists" --outdir "../../graphs/" --figdir "../../graphs/"
```

### Learn gene embeddings
- Run the following commands to learn gene embeddings:
```
cd src/embedding
python change_id.py --graphdir "../../graphs/edgelists/" #Map gene IDs to numbers for the model input
python main.py --input "../../graphs/graph.list" --outdir "emb" --hierarchy "../../graphs/graph.hierarchy" --iter 20 --regstrength 1 --workers 8 --dimension 128
```

### Identify enriched functions or gene regulons in different parts of the embedding space
- Download gene annotations from the [GOA database](https://www.ebi.ac.uk/GOA/human_release).
- Download the ontology from the [Gene Ontology database](http://geneontology.org/docs/download-ontology/).
- Run the following commands to split the embedding space into components (cluster genes in the whole space), then run enrichment analysis for each embedding component:
```
cd src/analysis
python cluster_genes.py --emb "../../example_data/gene.emb" --resolution 10 --nneighbors 30 --annotationfile "goa_human.gaf" --outdir "outdir"
python GO_TF_enrichment.py "outdir/cluster_0.txt" "outdir/population.txt" "goa_human_graph.gaf" --ev_exc=IEA --pval=0.05 --method=fdr_bh --pval_field=fdr_bh
```

### Infer novel gene functions
- Run the following commands to cluster each gene's embeddings of different graphs, then infer the functions of each gene in each cluster by running enrichment analysis on its neighbor genes in the cluster:
```
cd src/analysis
python neighbors_in_graphs.py --emb "../../example_data/gene.emb" --annotationfile "goa_human.gaf" --outdir "gene_list_outdir"
python GO_enrichment_in_neighbors.py "gene_list_outdir/A2M_nlist_0.txt" "gene_list_outdir/population.txt" "goa_human_name.gaf" --ev_exc=IEA --pval=0.05 --method=fdr_bh --pval_field=fdr_bh
```

## Credits
The software is an implementation of the method GIANT, jointly developed by Hao Chen, Nam D. Nguyen, Matthew Ruffalo, and Ziv Bar-Joseph from the [System Biology Group @ Carnegie Mellon University](http://sb.cs.cmu.edu/).

## Contact
Contact us if you have any questions:</br>
Hao Chen: hchen4 at andrew.cmu.edu</br>
Ziv Bar-Joseph: zivbj at andrew.cmu.edu</br>

## License
This project is licensed under the MIT License - see the LICENSE file for details.
