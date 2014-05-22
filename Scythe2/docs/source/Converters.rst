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


.. automodule:: scythe2.convert
    
    :undoc-members:

.. program-output:: ../../scythe2/convert/scythe_loc_gff.py -H

.. program-output:: ../../scythe2/convert/scythe_loc_tsv.py -H


