Installation
============

There are three ways to install ANNOgesic. Please refer to the following 
sections. Only Dockerfiles can install the requirments automatically. If 
you install ANNOgesic through other ways, please install the pre-required 
tools by yourself.


Github
----------

::

    git clone https://github.com/Sung-Huan/ANNOgesic.git

or

::

    git clone git@github.com:Sung-Huan/ANNOgesic.git

Then create a soft link of ``annogesiclib`` in ``bin``.

::

    cd ANNOgesic/bin
    ln -s ../annogesiclib .

Dockerfile
----------

Please Download the `Dockerfile <https://github.com/Sung-Huan/ANNOgesic>`_ and 
`RATT  <https://github.com/Sung-Huan/ANNOgesic>_` in our Github first.
Then switch to the folder which Dockerfile and RATT folder are located. For the following commands, please 
refer to `here <https://github.com/Sung-Huan/ANNOgesic/blob/master/docs/source/docker.rst>`_.

If you want to check Other commands of Docker, please refer to  `here <https://docs.docker.com/>`_.

pip3
----------

::

    pip3 install ANNOgesic
