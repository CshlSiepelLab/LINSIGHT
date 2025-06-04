## Synopsis

**LINSIGHT** is a statistical model for estimating negative selection on noncoding sequences in the human genome. The **LINSIGHT** score measures the probability of negative selection on noncoding sites which can be used to prioritize SNVs associated with genetic diseases or quantify evolutionary constraint on regulatory sequences, e.g., enhancers or promoters.

## Motivation

**LINSIGHT** can be viewed as a generalized linear model with an objective function derived from the [INSIGHT evolutionary model](http://mbe.oxfordjournals.org/content/30/5/1159). Given a set of genomic features, **LINSIGHT** seeks a linear combination of these features to best explain the signatures of negative selection in the human lineage. More specifically, if a noncoding site is under negative selection, it will be less likely to have a substitution or SNV in the human lineage. In addition, even if we see a SNV at the site, it will tend to segregate at low frequency because of selection. By combining the idea of generalized linear model and the INSIGHT likelihood function, **LINSIGHT** implicitly aggregates the very weak signatures of natural selection across millions of noncoding sites in the human genome to obtain accurate, high resolution estimate of evolutionary constraint. More detail can be found in the [**LINSIGHT** manuscript](http://biorxiv.org/content/early/2016/08/15/069682).

## Precomputed **LINSIGHT** scores

If you are interested in prioritization of variants or estimation of negative selection in the human genome, you can use the precomputed **LINSIGHT** scores described in our manuscript. You don't need to install and run **LINSIGHT** unless you want to retrain **LINSIGHT** using customized genomic features or extend it to non-human species. The precomputed **LINSIGHT** scores are based on the hg19 assembly and can be downloaded [Download LINSIGHT.bw](http://compgen.cshl.edu/LINSIGHT/LINSIGHT.bw). *(Right-click and choose "Save link as..." if clicking doesn't work.)* Note that precomputed **LINSIGHT** scores currently can only be used for non-commercial purposes because of the licenses associated with some of the genomic features used in training.

## Installation

**LINSIGHT** is implemented in standard C\+\+. A modern complier with C\+\+11 support is required to install **LINSIGHT**. We currently only tested **LINSIGHT** in Linux environment and recommend recent versions of g++ compiler (e.g., version 4.9.3 or higher). To install, please change to the directory of **LINSIGHT** source code and type *make* in Terminal. After compilation, you will obtain three excutable files, *LINSIGHT-prep*, *LINSIGHT-fit* and *LINSIGHT-score*. You can then add the **LINSIGHT** directory to your shell environment variables. *LINSIGHT-prep* is used to convert a genome-wide INSIGHT file to binary format for efficient operations. *LINSIGHT-fit* trains the **LINSIGHT** model and estimates the weights of genomic features. *LINSIGHT-score* computes genome-wide **LINSIGHT** scores.

## Quick guide

**LINSIGHT** requires two input files, a genome-wide INSIGHT file and a genome-wide feature file. The genome-wide INSIGHT input file follows the format described in the [INSIGHT manual](http://compgen.cshl.edu/INSIGHT/downloads/INSIGHT_Manual.pdf). The genome-wide feature file is a tab separated bed file and must be sorted by chromosome names and start positions. The first three columns of the feature file correspond to chromosome name, start position (0-based), and end position (1-based) while the following columns correspond to numerical or binary features associated with each region. Categorical variables can be handled by the standard [dummy coding](https://en.wikipedia.org/wiki/Categorical_variable) in statistics and machine learning. An example of feature file is below.
```
#chr	start	end	feature1	feature2	feature3
chr1	0	9996	0	0	0
chr1	9996	10003	0	1	0.5
chr1    10003   10009	1	0	0.2
chr1    10009   10021	0	0	0.8
...
...
```

Here we give a step-by-step guide to explain how to run **LINSIGHT** for the human genome.

* **Step** 1: Download genome-wide INSIGHT file (hg19 assembly) and genome-wide feature file used in our manuscript.
``` 
$ wget http://compgen.cshl.edu/LINSIGHT/genomic_feature.bed.gz
$ wget http://compgen.cshl.edu/LINSIGHT/INSIGHT_database.ins.gz
$ gunzip genomic_feature.bed.gz
$ gunzip INSIGHT_database.ins.gz
```

* **Step** 2: Convert the genome-wide INSIGHT file to binary files. *-i*: INSIGHT input file; *-o*: output directory of binary files.
```
$ mkdir binary_dir
$ LINSIGHT-prep -i INSIGHT_database.ins -o binary_dir/
```

* **Step** 3: Estimate weights in the **LINSIGHT** model. *-d*: directory of INSIGHT binary files; *-f*: genomic feature file; *-o*: output file of estimated parameters; *-n*: number of epochs. *-b1* and *-b2* corresponds to the genome-wide proportions of low- and intermediate- frequency derived neutral alleles in the **LINSIGHT** model. These two parameters can be calculated using INSIGHT and a file of genome-wide neutral regions (see [INSIGHT paper](http://mbe.oxfordjournals.org/content/30/5/1159) for more detail).
```
$ LINSIGHT-fit -d binary_dir/ -f genomic_feature.bed -o parameter.txt -n 20 -b1 0.7567672441 -b2 0.2053143291 > log.txt
```

* **Step** 4: Compute **LINSIGHT** scores given estimated weights. *-f*: genomic feature file; *-p*: parameter file generated in Step 3; *-o*: **LINSIGHT** score file in bedgraph format.
```
$ LINSIGHT-score -f genomic_feature.bed -p parameter.txt -o LINSIGHT_score.bedgraph
```
