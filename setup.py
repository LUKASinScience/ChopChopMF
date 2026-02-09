from setuptools import setup, find_packages

setup(
    name="ChimeraX-ChopChopMF",
    version="1.1",
    description="Visualization of a sequence alignment in 3D structure",
    long_description=(
        "This tool conducts an alignment of two amino acid sequences. "
        "The result is visible in the sequence viewer and the 3D structure. "
        "The model is color-coded according to the conservation score of the alignment."
    ),
    author="Lukas W. Bauer, ISTA Bioinformatics, Philipp Bauer",
    author_email="Lukas.Bauer@ist.ac.at, it@ist.ac.at",
    url="https://ist.ac.at/en/research/schur-group/",
    packages=find_packages(where="src"),  # Locate all Python packages in 'src'
    package_dir={"": "src"},  # Map root package to 'src'
    data_files=[
        # Include bundle_info.xml from the root directory
        (".", ["bundle_info.xml"]),
        # Include Python scripts
        ("chimerax/chopchop_mf", [
            "src/__init__.py",
            "src/chopalignment.py",
            "src/chopgetmissense.py",
            "src/chopmissense.py",
            "src/cropstructure.py",
            "src/duplicate.py",
            "src/alphafoldinfo.py",
            "src/pdbepisa.py",
            "src/foldseekanalysis.py",
            "src/pae_analysis.py",
             
        ]),
        # Include icon files
        ("chimerax/chopchop_mf/icons", [
            "src/icons/ChopAlignIcon.png",
            "src/icons/ChopGetMissense.png",
            "src/icons/ChopMissenseIcon.png",
            "src/icons/ChopUndo.png",
            "src/icons/AlphaFoldIconChop.png",
            "src/icons/crop.png",
            "src/icons/copyicon.png",
            "src/icons/pisa.png",
            "src/icons/foldseeklogo.png",
            "src/icons/pae_icon.png",
             
        ]),
    ],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: Freeware",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Framework :: ChimeraX",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    install_requires=[
        "ChimeraX-Core~=1.9",
        "ChimeraX-UI~=1.0",
    ],
    python_requires=">=3.7",  # Specify the minimum Python version
)

