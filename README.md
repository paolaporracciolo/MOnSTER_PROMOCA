# !!! Please visit https://github.com/Plant-Net/MOnSTER_PROMOCA.git for the final version of the tool !!!


# MOnSTER_PROMOCA
This repository contains the folders of the tools MOnSTER and PRO-MOCA

## About MOnSTER
MOnSTER is a pipeline to find CLUMPs (CLUster of Motifs of Proteins)
that are discriminative between two datasets of sequences
that will be called positive and negative.
       
The CLUMPs are found by gathering together motifs that show
similar physicochemical properties.
These physicochemical properties can show different ranges of
values, for this reason, we propose known 3 methods of scaling
the data for this step.
       
After clustering, the user will visualize the resulting dendogram.
To find the best compromise between number of CLUMPs and having
sufficiently different CLUMPs, MOnSTER employes the Davies Bouldin score.
It is then up to the user to choose number of intervals
of distances to find the best cut of the tree to obtain the most
reasonable number of CLUMPs.
(see doc for more information)
       
To determine which are the most discriminative CLUMPs, MOnSTER
calculates the MOnSTER score for each of them. The score
is between 0 and 2 and takes into account: the physicochemical
properties of the CLUMPs and the two datasets, and the occurrences
of the motifs and the CLUMPs in the two datasets.
The physicochemical properties that are employed for the MOnSTER_score
are those that show a statistically significant difference
in the two datasets.
       
MOnSTER scores the CLUMPs the highest if they show values of these
physicochemical properties that are in line with the positive
dataset and if they are more present in the positive dataset
compared to the negative one.

### Demo
In the folder MOnSTER/demo a demo is available step by step.

## About PRO-MOCA
PRO-MOCA (PROtein-MOtif Characteristics Aligner) is an
aligner of motif sequences.
PRO-MOCA aligns short sequences based on the characteristics
of the AAs rather than the AAs themselves.
       
The algorithm doesn't employ a substitution matrix.
Instead, PRO-MOCA employes a classification scheme of the 
AAs to convert the motifs into patterns of letters. Then
looks for a central position of the patterns, where a 
certain letter is repeated more than in the others.
Based on the central position of the alignment, it then
calculates the number of gaps ('_') to add to the right and 
to the left of all the patterns, after centering them.
At the end it reconverts the patterns into motifs, 
and gives a list of motifs with gaps to the left and 
to the right.
       
No gaps are admitted into the sequences, only at the 
extremities. This is because motifs are very short sequences
hence introducing gaps could mean losing the information
of the motif itself.

### Demo
In the folder PRO_MOCA/notebooks a demo is available step by step.
