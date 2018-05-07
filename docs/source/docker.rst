Docker image
==============

`Docker <https://www.docker.com>`_ is a platform for distributing package. 
It is light and easy to manage. ``ANNOgesic`` includes a ``Dockerfile`` which 
is for build up a environment and install all required tools for running ``ANNOgesic``.

You can simply pull the Docker image by running

::

    $ docker pull silasysh/annogesic

Alternatively, you can build the image by ``Dockerfile``.
Please go to the folder where ``Dockerfile`` are located. Then type

::

    $ sudo docker build -t="annogesic" .

It will build up an image called annogesic. You can see the images by typing ``docker images``

::

   REPOSITORY          TAG                 IMAGE ID            CREATED             VIRTUAL SIZE
   annogesic           latest              d35f555694ad        3 days ago          2.782 GB
   ubuntu              14.04               d2a0ecffe6fa        11 days ago         188.4 MB

Then we can use the image to create a container for running ``ANNOgesic``. Please type 

::

    $ docker run -t -i annogesic bash

Then you will jump into the container.

::

    root@c9de31fcd7e3:~# ls
    ANNOgesic

If you want to mount the files from your host to the container, just add ``-v`` to the command.

::

    $ docker run -t -i -v /host/path/target:/file/path/within/container annogesic bash

The paths should be absolute path. If we go to ``root`` in container. We can see the file.


If you want to copy the files from container to host, you can use ``cp``.

::

    $ docker cp <containerId>:/file/path/within/container /host/path/target

If you have no root permission for running Docker, Singularity is another way to 
build up the image without root permission.

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
