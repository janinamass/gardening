========
Tutorial
========


Installation
============
Via pip: ::
    
    pip install scythe
    
From source: ::
    
    wget TODO
    tar xvf TODO
    cd TODO/DIR
    python setup.py install

Requirements
------------

python modules
~~~~~~~~~~~~~~
* configparser

external software
~~~~~~~~~~~~~~~~~
Scythe will perform pairwise global alignments using the Needleman-Wunsch algorithm [NE1970] as implemented in needleall as part of the EMBOSS package [EMB2000]

* needleall `<http://emboss.sourceforge.net/>`_


Summary
~~~~~~~
On debian-based systems: ::

    sudo apt-get install emboss  python-configparser

Using Scythe
=============

Command line interface
----------------------
Please note that you cannot automatically download sequences from ENSEMBL via the command line interface.
Use either the GUI or download tab-separated files via BioMart.

.. program-output:: ../../scythe2/core/scythe_nv.py

Graphical user interface
------------------------
!todo


Converters
----------
Scythe uses a simple formats to store gene-transcript and ortholog information.
See the `Format` section for details. !todo
For the `loc` format, converters from tab-separated files (eg. download from ENSEMBL via BioMart) and gff3 are included.

.. automodule:: scythe2.convert
    
    :undoc-members:

loc
~~~

.. program-output:: ../../scythe2/convert/scythe_loc_gff.py -H

.. program-output:: ../../scythe2/convert/scythe_loc_tsv.py -H

grp
~~~

