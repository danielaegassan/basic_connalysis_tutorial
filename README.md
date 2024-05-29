# Basic connalysis tutorial

A small tutorial with the basic usage of connalysis for biological neural networks (BNNs) 

## Setup 

***In your own computer***

Install Connectome Utilities and Connectome analysis.  Warning: After this test, we will migrate connalysis to https://github.com/BlueBrain
   ```
   pip install git+https://github.com/danielaegassan/connectome_analysis.git
   pip install Connectome-Utilities
   ```

For more information on installation check: 

https://github.com/danielaegassan/connectome_analysis

https://github.com/BlueBrain/ConnectomeUtilities

***In BB5***

Flagser need specific gcc and cmake versions.  To smoothly run this on BB5 for this tutorials follow these steps: 

1. Load the correct modules
   
```
module purge
module load archive/2023-07
module load python/3.10.8
module load gcc/11.3.0
module load cmake/3.24.3
```

2. Create a virtual environment.  Run these upgrades in your virtual environment.

```
pip install --upgrade pip
pip install --upgrade wheel setuptools
```

3. Install necessary packages

```
pip install git+https://github.com/danielaegassan/connectome_analysis.git
pip install Connectome-Utilities
```

4. If you want to run things in jupyter notebooks, install that too.  Though you can run everything using scripts instead.  See below some notes on one way to run notebooks on BB5. 

```
pip install jupyterlab
pip install notebook
```

## Getting started 

1. Follow LINK TO NOTEBOOK, for instructions on how to quickly get the Celegans connectomes and possibly also the MICrONS connectome. This will also show you how to load them easily with Connectome-Utilities.
2. Follow notebooks 2, 3 and 4 to get an broad overview on using connalysis.
3. If time allows follow notebook 5 to get a quick example on generating distance dependent control models.

## Extra information 

#### Using jupyter notebooks in BB5.

This is how I run notebooks on BB5.  For an alternative (and possibly official way) on how to do this check Confluence and/or ask Christoph.

1. Get an interactive node.  Use the correct project, memory and and time requirements
```salloc --mem=100G  --time=3:00:00 -p interactive --account=proj83 -c4  srun --pty --preserve-env --mpi=none $SHELL```
2. Call Jupyter notebook, choose an available port
```jupyter notebook --no-browser --port=8896```
3. In a different terminal  connect to ssh for port forwarding.  Check out the ``port``, the ``node`` you were allocated and your ``username``.
```ssh -N -f -L localhost:8896:localhost:8896 <your_user_name>@<your_node>.bbp.epfl.ch```
4. Open the link provided to you in the first terminal









