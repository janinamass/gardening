from setuptools import setup, find_packages
setup(
        name='scythe2',
        version='0.1dev',
        author='Janina Mass',
        author_email='janina.mass@hhu.de',
        packages=find_packages(),
        scripts=['scythe2/convert/scythe_loc_gff.py', 'scythe2/convert/scythe_loc_gff.py'],
        license='GPLv3',
        description='Find best set of transcripts for one-to-one orthologous genes from two or more species',
        long_description=open('README.txt').read(),
        classifiers=[
            'Topic :: Scientific/Engineering :: Bio-Informatics'
            ],
        )
