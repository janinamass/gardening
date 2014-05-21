.. _tutorial:

========
Tutorial
========


Installation
============

Install Scythe in a virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Install virtualenv via pip: ::
    
    pip install virtualenv

Create a virtual environment that uses Python 3 and activate: ::
 
    virtualenv -p /usr/bin/python3 python3env
    source python3env/bin/activate

Finally, install the scythe package: ::
    
    pip install scythe

Without virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Via pip: ::
    
    pip install scythe
    
From source: ::
    
    wget TODO
    tar xvf TODO
    cd TODO/DIR
    python setup.py install

Requirements
============
* Python 3    

python modules
~~~~~~~~~~~~~~
* configparser
* mysql-connector-python
* httplib2

::

 pip install --allow-external mysql-connector-python mysql-connector-python
 pip install httplib2
 pip install configparser

external software
~~~~~~~~~~~~~~~~~
* needleall `<http://emboss.sourceforge.net/>`_

Scythe will perform pairwise global alignments using the Needleman-Wunsch algorithm [NE1970] as implemented in needleall as part of the EMBOSS package [EMB2000].
Please make sure `needleall` is in your PATH.

Install emboss on debian-based systems: ::

    sudo apt-get install emboss  


Using Scythe
=============

Graphical user interface
~~~~~~~~~~~~~~~~~~~~~~~~

:: 
 
 scythe-gui.py

or with configuration file `conf.py` ::

 scythe-gui.py conf.scy

(Read more about configuration files under :ref:`ConfigurationFiles`)

Command line interface
~~~~~~~~~~~~~~~~~~~~~~

.. program-output:: ../../scythe2/scythe2.py -h

Please note that you cannot automatically download sequences from ENSEMBL via the command line interface.
Use either the GUI or download tab-separated files via BioMart.



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

