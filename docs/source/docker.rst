Dockfiles
==============

`Docker <https://www.docker.com>`_ is a platform for distributing package. 
It is light and easy to manage. ``ANNOgesic`` includes a ``Dockfile`` which 
is for build up a environment and install all required tools for running ``ANNOgesic``.

Please go to the folder where ``Dockfile`` are located. Then type

::

    sudo docker build -t="annogesic" .

It will build up an image called annogesic. You can see the images by typing ``docker images``

::

   REPOSITORY          TAG                 IMAGE ID            CREATED             VIRTUAL SIZE
   annogesic           latest              d35f555694ad        3 days ago          2.782 GB
   ubuntu              14.04               d2a0ecffe6fa        11 days ago         188.4 MB

Then we can use the image to create a container for running ``ANNOgesic``. Please type 

::

    docker run -t -i annogesic bash

Then you will jump into the container.

::

    root@c9de31fcd7e3:~# ls
    ANNOgesic

If you want to mount the files from your host to the container, just add ``-v`` to the command.

::

    docker run -t -i -v '/home/silas/dockfiles/test2:/root/test' ubuntu:14.04 bash

The file in host is ``/home/silas/dockfiles/test2``. The mount file in container is 
``/root/test``. If we go to ``root`` in container. We can see the file.

::

    root@9a50d77ef14f:~# ls
    test

If you want to copy the files from container to host, you can use ``cp``.

::

    docker cp <containerId>:/file/path/within/container /host/path/target
