
# Raw Data:

Build from Allen Brain Human Development RNA-seq RPKM Data

Please see "meta.txt" for sample information

http://www.brainspan.org/static/download.html

http://www.brainspan.org/api/v2/well_known_file_download/267666525

# Workflow:

1. Calculate kendall correlation coefficient between development stage and expression value of each gene.

2. Use genes with kendall correlation coefficient smaller than -0.5. 

3. Use the average expression value of each development stage to build this reference.

Detailed scripts are in: https://github.com/jumphone/scRef/tree/master/scripts/ALLEN
 
# Column names:

Column name is the development stage tag. Please see meta.txt for details.

