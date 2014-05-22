.. _ConfigurationFiles:

Configuration Files
===================
The configuration files allow you save your settings/restore your settings.
Most entries should be self-explanatory.

Example file
-------------

.. confexample:

:: 

    [Mode]
    use_ensembl=yes
    use_local_files=no

    [Paths]
    fasta_directory=/home/nin/py3env/fa/
    loc_directory=/home/nin/py3env/loc/
    grp_file=/home/nin/py3env/magohopopafe.shared_by_all.grp
    output_directory=/home/nin/py3env
    
    [Cleanup]
    clean_up_directories=yes
    
    [Run_options]
    max_threads=1
    split_input=1
    
    [Penalties]
    gap_open_cost=10
    gap_extend_cost=0.5
    substitution_matrix=EBLOSUM62
    
    [Algorithm]
    use_sl_glob=unset
    use_sl_ref=unset
    use_mx_sum=yes
    
    [Fasta_header]
    fasta_header_delimiter=" "
    fasta_header_part=0


[Fasta_header]
--------------
In case the ids in your `loc` file are truncated or your fasta headers carry information in addition to the sequence identifier, you can define where to split the fasta header.

Example: ::
    
    >other:xyz_sequenceID_length:123_species:x
    MYQVLQAYDWKYLHDNHDCNFQVGGADQLGNI...

and: ::
    
    [Fasta_header]
    fasta_header_delimiter = "_"
    fasta_header_part = 1

Would result in taking only `sequenceID` into account.

[Run_options]
--------------

Multi core support is not yet implemented.

[Algorithm]
------------

See :ref:`algo`

