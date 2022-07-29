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
### Build gene graphs
Run the following commands to build co-expression graphs for scRNA-seq data or Slide-seq data:
```
cd src/build_graphs
python coexpression_knn_graph.py --K 10 --datatype "rna" --datadir "path/to/h5ad" --genescopefile "path/to/gene_scope" --outdir "path/to/out_dir"
```

The usage of this command is listed as follows:
```
usage: coexpression_knn_graph.py [-h] [--K K] [--datatype DATATYPE] [--datadir DATADIR]
                                 [--genescopefile GENESCOPEFILE] [--outdir OUTDIR]

Build co-expression graphs

optional arguments:
  -h, --help            show this help message and exit
  --K K                 Number of nearest neighbors
  --datatype DATATYPE   rna or slide
  --datadir DATADIR     Path to a directory that saves the h5ad files of gene expression
  --genescopefile GENESCOPEFILE
                        Path to the file that lists the genes considered in this study
  --outdir OUTDIR       Path to a directory where the graphs will be saved
```

Run the following commands to build gene-TF hypergraphs for scATAC-seq data:
```
cd src/build_graphs
python
```

The usage of this command is listed as follows:
```
```

Run the following commands to build spatial co-expression hypergraphs for Slide-seq data:
```
cd src/build_graphs
python
```

The usage of this command is listed as follows:
```
```



### Build dendrogram of graphs
Run the following commands to build dendrogram:
```
cd src/build_dendrogram
python build_dendrogram.py --K 50 --degenedir "path/to/de_genes" --graphdir "path/to/gene_graphs" --outdir "path/to/save_output" --figdir "path/to/save_output_figure"
```

The usage of this command is listed as follows:
```
usage: build_dendrogram.py [-h] [--K K] [--degenedir DEGENEDIR] [--graphdir GRAPHDIR] [--outdir OUTDIR]
                           [--figdir FIGDIR]

Build dendrogram of graphs

optional arguments:
  -h, --help            show this help message and exit
  --K K                 The number of top differentailly expressed genes considered to measure cell cluster
                        similarities
  --degenedir DEGENEDIR
                        Path to the lists of differentailly expressed genes of cell clusters
  --graphdir GRAPHDIR   Path to the gene graphs of cell clusters
  --outdir OUTDIR       Path to a directory where the dendrogram is saved for the input of the embedding algorithm
  --figdir FIGDIR       Path to a directory where the figure of the dendrogram is saved
```

### Learn gene embeddings
Run the following commands to learn gene embeddings:
```
cd src/embedding
python main.py --input "path/to/graph.list" --outdir "emb" --hierarchy "path/to/graph.hierarchy" --iter 150 --regstrength 1 --workers 8 --dimension 128
```

The usage of this command is listed as follows:
```
usage: main.py [-h] [--input INPUT] [--outdir OUTDIR] [--hierarchy HIERARCHY] [--dimension DIMENSION]
               [--walk-length WALK_LENGTH] [--num-walks NUM_WALKS] [--window-size WINDOW_SIZE] [--iter ITER]
               [--workers WORKERS] [--regstrength REGSTRENGTH] [--p P] [--q Q] [--l_rate L_RATE] [--weighted]
               [--unweighted] [--directed] [--undirected]

Run encoder

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         Path to a file containing locations of network layers
  --outdir OUTDIR       Path to a directory where results are saved
  --hierarchy HIERARCHY
                        Path to a file containing multi-layer network hierarchy
  --dimension DIMENSION
                        Number of dimensions. Default is 128.
  --walk-length WALK_LENGTH
                        Length of walk per source. Default is 10.
  --num-walks NUM_WALKS
                        Number of walks per source. Default is 20.
  --window-size WINDOW_SIZE
                        Context size for optimization. Default is 5.
  --iter ITER           Number of epochs in SGD
  --workers WORKERS     Number of parallel workers. Default is 8.
  --regstrength REGSTRENGTH
                        Hierarchical regularization strength. Default is 1.
  --p P                 Return hyperparameter. Default is 1.
  --q Q                 Inout hyperparameter. Default is 1.
  --l_rate L_RATE       learning rate. Default is 0.05.
  --weighted            Boolean specifying (un)weighted. Default is unweighted.
  --unweighted
  --directed            Graph is (un)directed. Default is undirected.
  --undirected
```

### Identify enriched functions in different parts of the embedding space
### Infer novel gene functions

## Credits
The software is an implementation of the method GIANT, jointly developed by Hao Chen, Nam D. Nguyen, Matthew Ruffalo, and Ziv Bar-Joseph from the [System Biology Group @ Carnegie Mellon University](http://sb.cs.cmu.edu/).

## Contact
Contact us if you have any questions:</br>
Hao Chen: hchen4 at andrew.cmu.edu</br>
Ziv Bar-Joseph: zivbj at andrew.cmu.edu</br>

## License
This project is licensed under the MIT License - see the LICENSE file for details.
