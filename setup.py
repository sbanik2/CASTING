from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / "README.md").read_text(encoding="utf-8")


setup(
    name="CASTING",
    version="1.0.0",  
    description="A continuous action spce tree search for inverse design",
    author="Suvo Banik",
    author_email="sbanik2@uic.edu", 
    keywords="Crystal structure prediction, Monte Carlo Tree Search, Nano Cluster, Optimization",
    install_requires=[
            "ase>=3.21.1",
            "networkx>=2.0",
            "numpy>=1.23.1",
            "pandas>=1.4.3",
            "pymatgen>=2022.7.25"
            "scipy>=1.8.0",
    ],
    
    
    scripts =[
        "CASTING/clusterfun.py",
        "CASTING/lammpsEvaluate.py",
        "CASTING/MCTS.py",
        "CASTING/perturb.py",
        "CASTING/utilis.py",
    ],
        
    
    
    classifiers=[

    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.8',
    
    ],

    
    long_description=long_description,  
    long_description_content_type="text/markdown",
    url="https://github.com/sbanik2/CASTING",
   
    packages=find_packages(),  
    python_requires=">=3.7",

)
