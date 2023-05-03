from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / "README.md").read_text(encoding="utf-8")


setup(
    name="CASTING",
    version="0.1.3",  
    description="A continuous action spce tree search for inverse design",
    author="Suvo Banik",
    author_email="sbanik2@uic.edu", 
    install_requires=[
                    "ase==3.22.1",
                    "mendeleev==0.13.1",
                    "networkx==2.8.4",
                    "numpy==1.24.3",
                    "pandas==1.4.4",
                    "pymatgen==2023.2.22",
                    "scipy==1.9.1",
                    "tqdm==4.64.1",

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
