=======
Scythe
=======

*Scythe*, Selection of Conserved Transcript by Homology Evaluation, tries to find the best matching set of transcripts for each one-to-one orthologous gene from two or more species.


Usage (GUI): ::
	
	scythe-gui.py


Usage (command line): ::
    
    ######################################
    # scythe.py v0.1                     #
    ######################################
  usage:
     scythe.py -i DIR -g .grpFILE --cleanup

  usage with configuration file:
     scythe.py --config configuration.scy

  general options:
    -C, --config                     use configuration file instead of
                                     command line parameters
    -c, --cleanup                    remove temporary files when done
    -h, --help                       prints this
    -i, --in_dir=DIR                 folder w/ subfolders "fa" and "loc"
    -m, --global_max                 always start with best pair
    -M, --global_sum                 optimize sum of transcript pairs
    -o, --out_dir=DIR                output directory [default:./]
    -v, --verbose                    be wordy

  alignment options:
     -O, --gap_open=FLOAT           needleall gap opening cost [default 10]
     -E, --gap_extend=FLOAT         needleall gap extension cost

  fasta options:
    -d, --delim=STRING               split fasta headers at STRING
    -a, --asID=INT                   use INTth part of fasta header as transcript-ID
                                     (default:0)
  debug options:
    -s, --stop_after=NUM             stop after NUM groups

  further help:
    Please see documentation in the 'docs'-directory.


Installation
============
Scythe needs Python 3. You might want to set up a virtual environment with Python 3 for easier package management. 
Please see the documentation in the 'docs'-directory for more details.
(In your virtual environment) install Scythe via ::
 pip install scythe

Requirements
============

Python modules
~~~~~~~~~~~~~~~
* configparser
* mysql-connector-python
* httplib2 

 ::
    pip install configparser
    pip install --allow-external mysql-connector-python mysql-connector-python 
    pip install httplib2

External programs
------------------
* needleall (EMBOSS package, http://emboss.sourceforge.net/)


Todo
====
* document this


Contributors
============


Thanks also to
==============

