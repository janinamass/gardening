.. _converters: 

Converters
==========

Scythe uses simple, human readable  formats to store gene-transcript and ortholog information.
See also :ref:`format`.

loc
---
Scripts to convert the following formats to `loc` format are included:

* gff3
* tab-separated (eg ENSEMBL BioMart)

run ::

    scythe_loc_gff.py -f GFF

or ::

    scythe_loc_tsv.py -f FILE.tsv

See also: :ref:`loc_format`.

grp
---
Please note that the `grp` converters  need a concatenated `loc` file in addition to the orthology information.
Scripts to convert the following formats to `grp` are included:

* orthomcl
* proteinortho
* tab-separated (eg ENSEMBL BioMart)

run ::
    
    scythe_grp_orthomcl.py

::
    
    scythe_grp_proteinortho.py

::

    scythe_grp_tsv.py


See also :ref:`grp_format`.


Downloading from ENSEMBL without the GUI
=====================================
To download sequences (pep and cds fasta files) from ENSEMBL without the graphical user interface, use ... TODO
TODO example
-------------
::
    sthsth.py -s homo_sapiens,mus_musculus -r 75

TODO: -h output

Use `scythe_ensembl_ortho_mysql.py` to download pairwise orthology information. 
scythe_ensembl_ortho_mysql.py example
-------------------------------------
::
     scythe_ensembl_ortho_mysql.py -s A,B,C,D,E -r 75

TODO:rename to scythe_ensembl_ortho_mysql.py
TODO: output of -h

Manual merge of tab-separated files  to one  `.grp` file
======================================================
If you have pairwise (two-species) files ready and want to
merge them into a multi-species `.grp` file you can do so via the 
`scythe_ensembl2grp`  and `scythe_mergeSubsets` scripts.

TODO rename both
scythe_ensembl2grp example
--------------------------
::
    scythe_ensembl2grp.py TODO

scythe_mergeSubsets example
-----------------------
::
    scythe_mergeSubsets.py TODO



.. automodule:: scythe2.convert
    
    :undoc-members:

.. program-output:: ../../scythe2/convert/scythe_loc_gff.py -H

.. program-output:: ../../scythe2/convert/scythe_loc_tsv.py -H


