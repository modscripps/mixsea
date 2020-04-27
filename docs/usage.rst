=====
Usage
=====

Import mixsea as::

    import mixsea

or::

    import mixsea as mx

Test running code example:

.. ipython:: python

    x = 2
    x**2

Works? Yes!

Let's try to generate a plot. First load modules in hidden mode. You shouldn't
be able to see this!

.. ipython:: python
    :suppress:

    import numpy as np
    import matplotlib.pyplot as plt

Then plot:

.. ipython:: python

    fig, ax = plt.subplots(
        nrows=1, ncols=1,
        figsize=(7, 4),
        constrained_layout=True
    )
    ax.plot([1, 2, 3], marker='o')
    # add directive
    @savefig plotting_example.png width=7in
    ax.set(xlabel='x')

Should show a simple figure.