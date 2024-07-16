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
 
meta.data = data.frame(sample = sample.id, cluster=cluster.id, row.names = colnames(counts))
object = CreateCellLabelerObject(counts,meta.data)
object

object = celllabeler(object, sample.var = "sample", cluster.var = "cluster",markers = markers,num.core = 10)
object@prediction
```
                       cluster         prediction score
1     cardiac endothelial cell   Endothelial cell     1
2          cardiac muscle cell      Cardiomyocyte     1
3                  native cell      Cardiomyocyte     1
4           smooth muscle cell Smooth muscle cell     1
5 fibroblast of cardiac tissue         Fibroblast     1
6                   macrophage         Macrophage     1



# Our group
https://sqsun.github.io/.
