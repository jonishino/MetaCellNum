# MetaCellNum;
### Estimate the founder population size of a metastatic tumor using Exome sequencing data




## Starting with running MetaCellNum

Make sure that ‘doSNOW’ and ’foreach’ packages for parallel processing is installed.


Set the working directory to the directory containing ‘MetaCellNum.R’ and ‘sample.txt’ files (an artificial Exome sequencing data), then read the script and the test data.

```R
source("MetaCellNum.R")
sample <- read.csv("sample.txt")
head(sample.txt)

CHR	LOC		 D1	m1	D2	m2
1	8075288		46	6	35	8
1	11115023	102	13	63	0
1	40705699	66	4	49	5
1	75602914	124	16	152	30
1	84764399	42	8	43	12
1	109268555	168	17	106	21
1	166818507	47	6	63	0
1	186296695	186	14	188	52
1	203140258	157	2	196	0
2	1480860		63	14	54	6
2	33442612	129	13	139	12
2	37349810	43	7	49	0
2	96956476	194	50	207	49
2	112922552	115	7	182	1
2	155099385	183	27	158	13
2	173352986	193	43	242	76
```

Estimate the founder population size of a metastatic tumor using the function ‘mleNb’ assuming 
purity of the primary tumor sample, pur1=0.5, and purity of the metastatic tumor sample, pur2=0.4.

The option disp=T prints the process of the estimation.

```R
pur1=0.5
pur2=0.4

eNb <- mleNb(dsin=sample[c("D1","m1","D2","m2")], disp=T,pur1=pur1,pur2=pur2, min_m1=2) 

```

Show the estimated founder population size.

```R
eNb
```

Get the 100 nonparametric bootstrap samples and calculate the 95% (or 80%) confidence limits specifying two as the number of clusters for parallel computing.

```R
nbs <- mleNb_nbs(dsin=sample[c("D1","m1","D2","m2")],pur1=pur1,pur2=pur2,reps=100,clusters=2,min_m1=2)
ci95 <- quantile(nbs,c(0.025,0.975))
ci80 <- quantile(nbs,c(0.1,0.9))

```

Show the the 95% (or 80%) confidence limits.

```R
ci95
ci80
```
