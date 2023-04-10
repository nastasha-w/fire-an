# fire_an
code and scripts for finding main halos, making histograms and 2D
projected maps, and analysing those products in FIRE simulations

'installation':
---------------
This is an in-development set of scripts, not a polished package. Many
of the scripts are in flux, especially those in `queuerun` and 
`makeplots`, and I do have a tendency to do development on the main 
branch (especially if it's unlikely to break something important). 
Because I'm updating this a lot, it's simpler to just use git pulls 
instead of running an installation every time. However, I have provided
a `condaenv_firep1.yml` file that includes the python packages required 
to run these scripts. It can be used to set up a conda environment 
where this code should run.

Note that a few functions requires additional packages, such as Andrew
Wetzel's FIRE analysis code. I haven't included these here, since these
are no longer the defaults I use.

Note that the main projection code and ion balance calculations call 
some compiled C code. I don't have that publicly available right now, 
because I did not write the SPH projection code myself (I only made 
small modifications), and I have not gotten around to asking the authors
for permission.

setup:
------
The `utils` directory should contain a file called `opts_locs.py`. This
file simply lists a number of files and directories needed to run the
code, such as where to search for simulation outputs, and where to store
halo data. Some example files `opts_locs_<system>.py` are provided.
Note that the `kernel_list`, `desngb`, `solar_abunds_ea`, and `Zsun_ea`
variables do not need to be changed. They're basically just stored here
because this started off as a parameter file.


running the scripts:
--------------------
To run, for example, `run.py`, options are:
- run from the directory above `fire_an` as
  ```
  python -m fire_an/queuerun/run -3
  ```
  (This will raise an error because `run.py` only takes postive integer
  arguments, but it also won't accidentally start an hour-long process
  on a login node.)
- or add the directory above fire_an (e.g. `foo`) to `PYTHONPATH`:
  ```
  export PYTHONPATH="/path/to/foo/:${PYTHONPATH}"
  ```
  and then, from anywhere, run: 
  ```
  python -m fire_an.queuerun.run -3
  ```
  or specify an absolute or relative path to `run.py` from your 
  current directory:
  ```
  python /path/to/fire_an/queuerun/run.py -3
  ```

This is basically an issue of the different scripts in different 
directories needing to be able to find each other. Running from the
directory above `fire_an` adds that directory to `PYTHONPATH`, like 
adding it to `PYTHONPATH` explicitly does. 


main functions:
---------------
TODO

warning:
--------
This repo is mainly for me to sync up in-progress work. Much of this
code is untested, half-finished projects I abandoned, or work in 
progress. Also, I have a tendency to do work and testing on the main
branch.

Note that tests are usually of the 'plot and see if it makes sense'
form, not functions that will return a boolean value.

