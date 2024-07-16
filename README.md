# CellLabeler
CellLabeler is implemented as an open source R package for detecting unique marker genes and automatic cell type annotation in  scRNA-seq transcriptomic studies. 

# Installation

> library(devtools)
> 
> install_github("jinning1012/CellLabeler")

# Usage
The main function is celllabeler(). You can find the instructions and an example by '?celllabeler.default'. celllabeler() can be performed on both raw counts or celllabeler objects. 

> ## create celllabeler object
> data(exampledata)
> head(sample.id)
# [1] "TSP12" "TSP12" "TSP12" "TSP12" "TSP12" "TSP12"
> head(cluster.id)
# [1] "cardiac endothelial cell" "cardiac endothelial cell"
# [3] "cardiac muscle cell"      "cardiac endothelial cell"
# [5] "native cell"              "cardiac muscle cell"     

> meta.data = data.frame(sample = sample.id, cluster=cluster.id, row.names = colnames(counts))
> object = CreateCellLabelerObject(counts,meta.data)
> object
## An object of class CellLabeler 
## Gene number: 1546 
## Cell number: 4285 
## Meta columns: sample, cluster 

> object = celllabeler(object, sample.var = "sample", cluster.var = "cluster",markers = markers,num.core = 10)
Performing log-normalization
0%   10   20   30   40   50   60   70   80   90   100%
[----|----|----|----|----|----|----|----|----|----|
**************************************************|
## Combining similar clusters before DE test.
## ===== CellLabeler INPUT INFORMATION ====## 
## number of total cells:  4285 
## number of total features:  581 
## number of cell clusters:  6 
## number of cores:  10 
## ========== END INFORMATION ============## 

### Logistf on cardiac endothelial cell versus cardiac muscle cell...
### Logistf on cardiac endothelial cell versus native cell...
### Logistf on cardiac endothelial cell versus smooth muscle cell...
### Logistf on cardiac endothelial cell versus fibroblast of cardiac tissue...
### Logistf on cardiac endothelial cell versus macrophage...
### Logistf on cardiac muscle cell versus native cell...
### Logistf on cardiac muscle cell versus smooth muscle cell...
### Logistf on cardiac muscle cell versus fibroblast of cardiac tissue...
### Logistf on cardiac muscle cell versus macrophage...
### Logistf on native cell versus smooth muscle cell...
### Logistf on smooth muscle cell versus macrophage...
### Logistf on fibroblast of cardiac tissue versus macrophage...
### Logistf on smooth muscle cell versus fibroblast of cardiac tissue...
### Logistf on native cell versus macrophage...
### Logistf on native cell versus fibroblast of cardiac tissue...

> 



# Example
Example code is listed in "R/example.R" for reference.

# Our group
https://sqsun.github.io/.
