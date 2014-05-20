========
Tutorial
========


Installation
------------
Via pip: ::
    
    pip install scythe
    
From source: ::
    
    wget TODO
    tar xvf TODO
    cd TODO/DIR
    python setup.py install

Requirements
~~~~~~~~~~~~
Scythe will pairform pairwise global alignments using the Needleman-Wunsch algorithm [NE1970] as implemented in needleall as part of the EMBOSS package [EMB2000]

* needleall `<http://emboss.sourceforge.net/>`_

On debian-based systems: ::

    sudo  apt-get install emboss

.. [NEE1970] Needleman, Saul B.; and Wunsch, Christian D. (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins". Journal of Molecular Biology 48 (3): 443–53. doi:10.1016/0022-2836(70)90057-4. PMID 5420325.

.. [EMB2000] EMBOSS: The European Molecular Biology Open Software Suite (2000) Rice,P. Longden,I. and Bleasby, A.Trends in Genetics 16, (6) pp276--277
