Installation
============

There are three ways to install ANNOgesic. Please refer to the following 
sections. Only via the requirements are also installed. If
you install ANNOgesic through on of the other ways, 
please install the pre-required 
tools by yourself.


Github
----------

::

    $ git clone https://github.com/Sung-Huan/ANNOgesic.git

or

::

    $ git clone git@github.com:Sung-Huan/ANNOgesic.git

Then create a soft link of ``annogesiclib`` in ``bin``.

::

    $ cd ANNOgesic/bin
    $ ln -s ../annogesiclib .

Docker
----------

You can simply pull the Docker image as following

::

    $ docker pull silasysh/annogesic

Alternatively, you can build the image via Dockerfile.
Please Download the `Dockerfile <https://github.com/Sung-Huan/ANNOgesic>`_ from our Git repository.
Then switch to the folder which Dockerfile are located. For the following commands, please 
refer to `here <https://github.com/Sung-Huan/ANNOgesic/blob/master/docs/source/docker.rst>`_.

If you want to check other commands of Docker, please refer to  `here <https://docs.docker.com/>`_.

pip3
----------

::

    $ pip3 install ANNOgesic
    $ pip3 install ANNOgesic --upgrade
