Basic Seurat: <br>
https://satijalab.org/seurat/articles/essential_commands.html

Tutorial: <br>
https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

NormalizeData and ScaleData: <br>
https://github.com/satijalab/seurat/issues/950#issuecomment-444626498 <br>
https://towardsdatascience.com/understand-data-normalization-in-machine-learning-8ff3062101f0

Standarization: <br>
<img src="https://user-images.githubusercontent.com/77600778/120001621-9673b900-bf99-11eb-873d-da6c0e5a9fb9.png" width="400" height="350">

FindVariableFeatures:<br>
1. selection.method = 'vst'<br>
Polynomial fit a curve to predict the log(variance) of each gene as a function of its log(mean). <br>
Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter). <br>
[Detials](https://doi.org/10.1016/j.cell.2019.05.031) <br>
nfeatures = 2000 only used when selection.method is set to 'dispersion' or 'vst' <br>
https://github.com/satijalab/seurat/issues/1410 <br>
<img src="https://user-images.githubusercontent.com/77600778/121581966-149d7a00-c9f4-11eb-918a-bda02186c318.png" width="550" height="350">\
2. selection.method = 'mean.var.plot'<br>
Features are binned by their mean expression, and genes with high standarized dispersions are selected as HVGs in each bin. <br>
y.cutoff = 2: features that are more than 2 SD away from the average dispersion within a bin. <br>
The default X-axis is the mean expression level, and for Y-axis it is the log(Variance/mean). <br>
All mean/variance calculations are not performed in log-space, but the results are reported in log-space. <br>
<img src="https://user-images.githubusercontent.com/77600778/121582022-22eb9600-c9f4-11eb-926d-79e4ad44e07a.png" width="550" height="350">\
3. selection.method = 'dispersion'<br>
Select the genes with the highest dispersion <br>
<img src="https://user-images.githubusercontent.com/77600778/121582036-2717b380-c9f4-11eb-9ac6-1e452a84d7e2.png" width="550" height="350">


