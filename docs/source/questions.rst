Q and A
=======

- **Q1:**

I got the error about complaining no module can be found as following:

::

    Traceback (most recent call last):
      File "ANNOgesic.py", line 6, in <module>
        from annogesiclib.controller import Controller
    ImportError: No module named annogesiclib.controller

- **A1:**

Switch to Python 3.3 or higher. The newer versions can handle this.

- **Q2:**

I try to install ANNOgesic in my own directory. However, when I execute ANNOgesic, an error comes out as following:

::

    Traceback (most recent call last):
      File "bin/annogesic", line 7, in <module>
        from annogesiclib.controller import Controller
    ImportError: No module named 'annogesiclib'

- **A2:**

Please generate a soft link of ``annogesiclib`` to ``bin``.

::

    $ cd ANNOgesic/bin
    $ ln -s ../annogesiclib .

- **Q3:**

When I finished sRNA detection, there are some files in ``sRNAs/gffs/for_classes`` and ``sRNAs/tables/for_classes``. 
What do the classes mean?

- **A3:**

Please check ``stat_sRNA_class_$GENOME.csv`` in ``stat`` folder. The information of classes can be found in it. The 
classes are based on the input information.

- **Q4:**

I got a sRNA gff file which only contain ``##gff-version 3``. What does it mean?

- **A4:**

If you only got ``##gff-version 3`` in sRNA gff file, it means that ANNOgesic can not find any sRNA in this genome. 
This case can be applied to other genomic feature prediction.

- **Q5:**

How long does sRNA target prediction take?

- **A5**

It depends on the input data size. However, most of the sRNA target prediction tools need to spend around one hour for detecting 
targets of one sRNA. RNAplex is a fast one which can detect the targets of one sRNA in several minutes. If you want to reduce the 
running time, you can specify ``--program RNAplex``. However, the accuracy may be influenced.

- **Q6:**

When I run sRNA target prediction, I found many files in ``RNAplfold``. Are they necessary?

- **A6:**

These files are intermediate input for RNAplex. They will be deleted when RNAplex is done.
