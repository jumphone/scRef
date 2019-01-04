# Data Type

Human, RPKM

# Data Source:

This reference is built from Allen Human Brain Development RNA-seq Data

Please see "meta.txt" for sample information

http://portal.brain-map.org/

http://www.brainspan.org/static/download.html

http://www.brainspan.org/api/v2/well_known_file_download/267666525

# Workflow:

1. Calculate kendall correlation coefficient between development stage and expression value of each gene.

2. Select genes (2,221 genes) with kendall correlation coefficient smaller than -0.5. 

3. Use the average expression value of each development stage to build this reference.

* Detailed scripts are in: https://github.com/jumphone/scRef/tree/master/scripts/ALLEN

* We have generated a mouse reference by mapping human gene to mouse gene.

# Tips:

Column name is the development stage. Please see meta.txt for details.

