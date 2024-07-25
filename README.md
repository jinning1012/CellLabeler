# CellLabeler
CellLabeler is implemented as an open source R package for detecting unique marker genes and automatic cell type annotation in  scRNA-seq transcriptomic studies. 

# Installation
```{r}
library(devtools)
install_github("jinning1012/CellLabeler")
```
# Usage
The main function is celllabeler(). You can find the instructions and an example by '?celllabeler.default'. celllabeler() can be performed on both raw counts or celllabeler objects. 

# Example
We subsample a real scRNA-seq data as an example. 
```{r}
data(exampledata)
head(sample.id)
head(cluster.id)

## create object
meta.data = data.frame(sample = sample.id, cluster=cluster.id, row.names = colnames(counts))
object = CreateCellLabelerObject(counts,meta.data)
object

## run celllabeler
object = celllabeler(object, sample.var = "sample", cluster.var = "cluster",
  markers = markers, num.core = 10)
object@prediction
```
The output is a CellLabeler object involoving ude (uniquely differential expressed genes), prediction (a dataframe of cluster and predicted cell type infomation), ModelFits (full statistic results of pairwise Firth LR test) and ModelScores (evaluation score of each potential cell type). Usually we save there results by "res = list(ude = object@ude, prediction = object@prediction, ModelScores = object@ModelScores, ModelFits = object@ModelFits)". Then we offer a function to update the prediction by giving a new marker list.

```{r}
res = AddMarkers(res = res, counts = counts, markers = markers, cluster.id = cluster.id)
```

We also involve a clustering strategy. We adopt the approach following “The molecular cytoarchitecture of the adult mouse brain”. We cluster the cells for multiple rounds and achieve different levels of precision.
```{r}
cluster_df = multiple_clustering(counts, loop.max = 10, k = 50,
  min.cells = 100, num.hvg = 10,num.de = 10,logfc = 0.1, min.umi = 100, verbose = T)
```


# Our group
https://sqsun.github.io/.
