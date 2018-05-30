Installation
============

There are three ways to install ANNOgesic. Please refer to the following 
sections. ANNOgesic can only work when the requirements are installed properly. If
you install ANNOgesic through source code or ``pip3``, please install the pre-required 
tools by yourself.


Github
----------

All the source code including a run script (contains all the commands which are presented in tutorial) 
of ANNOgesic can be retrieve from our Git repository. Using the following commands can clone the 
source code easily.

::

    $ git clone https://github.com/Sung-Huan/ANNOgesic.git

or

::

    $ git clone git@github.com:Sung-Huan/ANNOgesic.git

In order to make ANNOgesic runnable, we should create a soft link of ``annogesiclib`` in ``bin``.

::

    $ cd ANNOgesic/bin
    $ ln -s ../annogesiclib .

Docker
----------

Some modules of ANNOgesic need third-party tools. In order to avoid all the possible issue caused by the dependencies, 
a Docker image is provided. For the details of Docker image, please check `here <https://www.docker.com/>`_.

For using Docker image, please use one of the following commands:

1. You can simply pull the Docker image as following

::

    $ docker pull silasysh/annogesic

2. Alternatively, you can build the image via Dockerfile.
Please Download the `Dockerfile <https://github.com/Sung-Huan/ANNOgesic>`_ from our Git repository.
Then switch to the folder which Dockerfile are located. For the following commands, please 
refer to `here <https://github.com/Sung-Huan/ANNOgesic/blob/master/docs/source/docker.rst>`_.

If you want to check other commands of Docker, please refer to  `here <https://docs.docker.com/>`_.

Singularity
-----------

`Singularity <https://singularity.lbl.gov/index.html>`_ is another way to install ANNOgesic via 
Docker image without root permission.

::

    $ singularity build \
        annogesic.img \
        docker://silasysh/annogesic:latest

After building Singularity image of ANNOgesic, the user just needs to put the following line before 
the command that needs to be executed.

::

    singularity exec -B $STORAGE_PATH annogesic.img

Please put the storage path of your home directory to ``$STORAGE_PATH``. ``df`` can be used to check the 
storage system. 

pip3
----------

ANNOgesic is also hosted in PyPI server. Thus, it can be simply installed via ``pip3``.

::

    $ pip3 install ANNOgesic
    $ pip3 install ANNOgesic --upgrade

You can also install ANNOgesic without root permission.

::

    $ pip3 install --user ANNOgesic
    $ pip3 install ANNOgesic --user --upgrade

Install Dependencies
====================

If the user want to install ANNOgesic via source code, ``get_package_database.sh`` can 
provide a way to install tools and download database automatically. The required versions 
of the tools will be shown on the screen as well.
