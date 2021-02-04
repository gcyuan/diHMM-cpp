Installation
Go into the build dir and run

cmake ..
make
Then you can open a Python shell in the same dir and do

>>> import greet_ext
>>> greet_ext.greet()
'Hello world'
Training a model
Training a diHMM model can be done by using the script train.py, after making necessary changes to input data path and other parameters.

Applying a diHMM model for chromatin state annotation
Annotation can be done with the script annotation.py, after making necessary changes to input data path and other parameters.
