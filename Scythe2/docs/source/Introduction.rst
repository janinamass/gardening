============
Introduction
============
Scythe picks best matching transcripts for one-to-one orthologous genes from two or more species.
The goal is to provide the best, i.e. most homologous set of sequences, for a subsequent multiple sequence alignment in order
to minimize sources for misalignment in an automated fashion.
Scythe will perform pairwise global alignments using the Needleman-Wunsch algorithm [NE1970] as implemented in needleall as part of the EMBOSS package [EMB2000].

Please see the :ref:`tutorial` on how to use Scythe.
Important! needleall seems to be broken in EMBOSS 6.6.0.0; this is working with needleall from EMBOSS 6.4.0.0.

.. [NEE1970] Needleman, Saul B.; and Wunsch, Christian D. (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins". Journal of Molecular Biology 48 (3): 443–53. doi:10.1016/0022-2836(70)90057-4. PMID 5420325.

.. [EMB2000] EMBOSS: The European Molecular Biology Open Software Suite (2000) Rice,P. Longden,I. and Bleasby, A.Trends in Genetics 16, (6) pp276--277

Algorithms and Application
==========================
Scythe performs pairwise global alignments [NEED1970] to measure the similarity between transcripts.
Transcripts are only compared between species.

There are three main strategies implemented in Scythe to deal with the sequences after scoring:
A single-linkage approach either starting with reference transcripts as seed (sl_ref_) and adding best matching sequences from non-reference species 
or starting out with the best matching transcript pair for a gene (sl_glob_). 
A `reference species` is defined as a species that has only one transcript (`reference transcript`) for a gene.
Every gene is processed individually, a `reference species` is only local to a gene and is automatically derived from the input data.
Alternatively, the maximum-sum (mx_sum_) approach calculates the score for all transcript pairings between the species and return a maximum-scoring set.
Please note that this approach might not be feasable for large data sets.

.. mx_sum:

mx_sum
------
mx_sum (maximum-sum)
max sum algorithm returns an optimal solution for the problem of score maximization. 
There are, however, scenarios where this would not represent the desired result: In cases of single-transcript outliers,
their pairwise score to all other sequences is always taken into account and may favor sequences that are less dissimilar to the outlier.

.. sl_ref:

sl_ref
------
sl_ref(single-linkage-reference)... 
It works well in cases where the user wants to
focus on finding similar gene models given one reference species, e. g. when the user
has already decided on a splice variant in one species.

.. sl_glob:

sl_glob
-------
sl_glob (single-linkage-global) due to the single linkage selection approach, the algorithm starting with the best match-
ing pair should provide a good starting point in absence of a reference gene model and
may be able circumvent outliers.

.. bf:

bf_msa
------
!todo
brute force all vs all msa comparison using mafft

