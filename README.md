# fire_an
code and scripts for finding main halos, making histograms and 2D
projected maps, and analysing those products in FIRE simulations

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
  and then either run 
  ```
  python -m fire_an.queuerun.run -3
  ```
  from whereever, or run
  ```
  python /path/to/fire_an/queuerun/run.py -3
  ```
  specifying a an absolute or relative path to `run.py` from your 
  current directory.

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

