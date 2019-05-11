This project aims to demonstrate the the relation between proteome allocation and rate-yield tradeoff of E.coli.
doi: https://doi.org/10.1101/414912

We focus on the k_eff parameters in the ME-model

In the ME-model, the metabolic pathways generate the precursors for protein synthesis. And the protein expressions, which are part of the biomass, back constraint to the metabolic fluxes. In this case, solving the ME-model is a non-linear optimization problem. Here, we demonstrate the basic features of ME-model and 

## Setting up the ME-model and solver
### Create a virtual environment
mkdir rate_yield
virtualenv ./rate_yield/
cd rate_yield/
source ./bin/activate

### Install cobrame package
git clone https://github.com/SBRG/cobrame.git
cd cobrame
python setup.py develop

git clone https://github.com/SBRG/ecolime.git
cd ecolime/
python setup.py develop
