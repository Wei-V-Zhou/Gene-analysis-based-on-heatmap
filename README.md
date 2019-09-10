# Gene-analysis-based-on-heatmap
Marker genes cluster and new gene expression behavior analysis

## Filter the suitable marker genes
1. remove the genes that expressed in all cells
2. remove the genes that didn't express in any cells
3. remove the genes that expressed in wrong clusters

***Note:*** delete the non-specific genes and wrong genes
***

## Plot the heatmap to decide clusters
```
pheatmap(objective,treeheight_col = 30,treeheight_row = 0,show_colnames = F,clustering_method="ward.D",
           color = colorRampPalette(c("green", "black", "red"))(50),annotation_col = annotation_col)
```
## Add the targets into the heatmap expression
![avatar](C:\Users\asus1\Desktop\360截图178003067579120.png)
