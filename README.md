# DarkSpace

This workflow post-processes the results from the
[DarkSpaceProject](https://github.com/pporrasebi/darkspaceproject) to try to
stablish a rank for curating PMID articles.

## Status

The current version is the first implementation with little attention paid to
details. 

### Reading the code

There are two relevant files with code:

* lib/rbbt/sources/DarkSpace.rb: where we process slightly the input files from the darkspaceproject repository
* workflow.rb: Where we build the scores in a step-wise process

The most interesting one is the workflow.tb. It's easy to get confused by the
syntax that relays heavily in the Rbbt framework that introduces several
programming paradigms that you might not be familiar with; on the other hand
is very succinct and expressive and gets right to the point. Try to use your
imagination to figure out what is does.

### Results

The results for the different steps are included in the directory 'results'.
The '.info' files contain execution meta-data.

## General approach

The articles that we will rank are those that are found by the different
text-mining resources listed in the DarkSpace Project: EPMC, EVEX, and
STRING_TM. 

For these resources scores are provided at the level of 'triplet'
(protein-protein interaction pair plus PMID). We will use this scores to derive
the 'relevance' score, which measures the chance that the article is relevant
for PPI information. This score will then be weight by 'interest' score that
measures how much the article shines light into the 'Dark Space'

## Computing the Relevance score

To quickly aggregate the several text-mining scores into one we use a
prediction approach. Each article is represented by its vector of scores, which
are the input features. The response variable is the number of 'other'
resources (see below) that cover that PMID which represent the 'Bright Space',
we sum 10 if the 'IMEX' is one of them to give an special weight to these. We
use RandomForest to train a model to predict the response from the features.

The list of 'other' resources are those that are not text mining.

## Computing the interest score

Interest is the sum of 'protein interest' for all proteins described in the
article divided by the coverage of the article in 'other' resources. We sum 1
to avoid infinity.

The 'protein interest score' is calculated as 1 / (N + 1) where N is the number
of PPI-Databases entries. For instance, if PROT1 appears in three PPI described
in IMEX and two PPI described in OmniPath its interest is 1 / (5 + 1). The same
PPI in different resources still counts twice.

For each PMID we get all the PPI pairs that are associated to it and extract
the list of proteins, we then sum up the interest of this proteins to derived
the 'protein interest' score for that article. Proteins in several PPIs are
counted several times.

## Computing the final combination score

The final score is the 'relevance' score divided by the PMID coverage times the
cumulative 'protein interest' score

## Main enhancement points

* Treatment of resources
  * Which to consider as 'other' or the 'Bright Space' 
  * Which to use as source Text-Mining articles to consider
* Integration of text-mining score into relevance score
* Calculation of interest scores
  * Use triplets not only individual protein interest
* Integration of relevance and interest into one combined score
  * Consider a thresholding strategy
