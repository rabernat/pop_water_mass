# Water Mass Analysis in Hi-Res POP Model #

## Instructions for Running on Yellowstone ##

### The Git Repository ###

First you need to clone the repository onto your home directory on yellowstone

    $ git clone https://ryanaberanthey@bitbucket.org/ryanaberanthey/pop_water_mass.git

This will create a directory called pop_water_mass. In order to update the code from git, you do

    $ git pull

If you make any changes to the repository, you will need to think about how they will be merged with the master repository. Git is a powerful but complicated tool. I recommend [reading some tutorials](https://www.atlassian.com/git/).

### Modules ###

In order for python to work properly, you need to load the proper modules.

    $ module load python all-python-libs

### Installing _watermasstools_ Package ###

In order for the notebooks and scripts to work, you need to "install" the watermasstools package. The best way to do this is the following.

    $ cd pop_water_mass
    $ python setup.py develop --user

Now you will be able to import the watermasstools modules.

### Using the IPython Notebooks ###

(Note: I figured out how to do this by reading the [CISL documentation](https://www2.cisl.ucar.edu/resources/yellowstone/software/ipython#notebook))

The IPython notebook runs in your web browser. In order to use it on Yellowstone, you need to find a way to establish an http connection from your local machine to a Yellowstone compute node where the python kernel is actually running. The steps are as follows

1.    Get a an [interactive node](https://www2.cisl.ucar.edu/resources/geyser_caldera/applications#startapp) on geyser. The easiest way to do this is to just type

         $ execgy

    This will start an interactive session and tell you the name of the host your session is running on. For example:

        <<Waiting for dispatch ...>>
        <<Starting on geyser08-ib>>
    Here the host name is geyser08-ib. You need this in a future step.

2.    Load the python modules in your interactive shell:

        $ module load python all-python-libs

3.   Start the ipython notebook server

        $ ipython notebook --no-browser --ip='*'

    You should see a message like

        [NotebookApp] The IPython Notebook is running at: http://[all ip addresses on your system]:8888/

    This means that the server is running on the default port of 8888.

4.    Now your notebook server is now running, but unfortunately you have no way to connect to it from your local machine. To do this, you need to establish an [SSH tunnel](http://www.revsys.com/writings/quicktips/ssh-tunnel.html). To do this, from your *local machine*, ssh again to yellowstone, with the following special flags:

        $ ssh -L 8890:geyser08-ib:8888 yellowstone.ucar.edu

    This syntax will establish a tunnel connection from port 8890 of your local machine to port 8888 on geyser08-ib via yellowtone.ucar.edu. (Obviously geyser08-ib should be replaced with whatever host your interactive job got placed on in step 1.) The local port 8890 is also your choice.

5.   Finally, connect to the notebook server by entering the address http://localhost:8890/ in your browser. You should see an IPython notebook server. Navigate to pop_water_mass/notebooks.

At this point you can try to run any of the notebooks, although they won't all work on yellowstone. (Some of the file paths refer to my local machine.) One notebook that definitely _does_ work with your model is transformation_3Dpop.ipynb.

### Using the Scripts ###

There are launch scripts in the launch directory.