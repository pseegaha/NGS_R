# NGS_R
The code does the following jobs.....Note: Make sure the datasets of one's interest are installed and loaded
## 1) R code for NGS analysis
## 2) Accessing the data stored
## 3) Extracting SNP Information for a set of RS IDS
## 4) Injection in the Reference genome
## 5) Quality Control
## Arguments

```seqname```

The name of the sequence for which to get the SNP locations and alleles.

If as.GRanges is FALSE, only one sequence can be specified (i.e. seqname must be a single string). If as.GRanges is TRUE, an arbitrary number of sequences can be specified (i.e. seqname can be a character vector of arbitrary length).

```as.GRanges```	

TRUE or FALSE. If TRUE, then the SNP locations and alleles are returned in a GRanges object. Otherwise (the default), they are returned in a data frame (see below).

```caching```	

Should the loaded SNPs be cached in memory for faster further retrieval but at the cost of increased memory usage?

```rsids```	

A vector of rs ids. Can be integer or character vector, with or without the "rs" prefix. NAs are not allowed.
