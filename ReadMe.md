This repository contains the analysis code used for the Brain Map Cognitive Ontology paper that is currently under review:

Bolt, Nomi, Aren, Vij, Riedel, Salo, Laird, Eickhoff, and Uddin (Under Review). Towards a Data-Driven Ontology of Cognitive-Neural Mappings

There is current agreement that the way we divide up the mind, our cognitive ontology or taxonomy, doesn't do the job. The idea of this project is to use the relationships among whole-brain BOLD activation patterns as a source of novel data-driven categories. Of course this code can be used for any sort of analysis beyond the goals of this project.

The main features of this analysis pipeline are a fast decomposition of whole-brain activation maps using non-negative matrix factorization, followed by a data-driven clustering analysis. This code borrows a lot from well-developed MATLAB functions and toolboxes so please make sure to have the following code dependenices on your path:

Code dependencies you need to run this code:

1. SVD Initialization provided from the following GitHub page (user: trigeorgis) for semi-NMF: 

https://github.com/trigeorgis/Deep-Semi-NMF/blob/master/matlab/NNDSVD.m

* just the 'NNDSVD.m' function is needed

2. NMF Algorithm provided by Haesun Park on her webpage (and associated functions):

https://www.cc.gatech.edu/~hpark/nmfsoftware.php

3. NMF Instability measure provided by Wu et al. (2016)

https://github.com/bdgp/staNMF

* Just the 'amariMaxerror.m' function is needed

4. Brain Connectivity Toolbox (BCT) for graph clustering,consensus procedure, visualization, etc.

https://sites.google.com/site/bctnet/

5. Adjusted Rand Index to compute stability for consensus clustering algorithm provided on the Network Community Toolbox webpage:

http://commdetect.weebly.com/



