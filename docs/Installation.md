
### Manual Installation
[Install the anaconda package](https://docs.anaconda.com/anaconda/install/). Then, 
```
conda env create --name CASTING
conda activate CASTING
git clone https://github.com/sbanik2/CASTING.git
pip install -r requirements.txt
python setup.py install
```

### Installation with pypi

```
pip install CASTING

```
***The package requires python lammps binding to run. First, lammps package needs to be downloaded from [LAMMPS download](https://www.lammps.org/download.html) and compiled. The instructions on python integration can be found here [LAMMPS-Python](https://docs.lammps.org/Python_install.html).

