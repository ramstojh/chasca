chasca
======

.. image:: https://img.shields.io/badge/powered%20by-astroplan-orange
    :target: https://astroplan.readthedocs.io/en/latest/

    
*chasca* means star in the Inca–Andean–Quechua cosmovision. *chasca* is a simple tool based on `astroplan <https://astroplan.readthedocs.io/en/latest/>`_  to plan astronomical observations. To exploid the full `astroplan <https://astroplan.readthedocs.io/en/latest/>`_ features plase see its tutorials.


Installation
------------
Installing *chasca* requires only one step. Please run the following pip command::

    pip install chasca

Note that you will need Python (>=3.7-3.9) installed.
If you already have *chasca* installed, you should consider upgrading to the latest version via::

    pip install chasca --upgrade

Dependencies
------------
The main dependencies of *kanchay* is  `astroplan <https://astroplan.readthedocs.io/en/latest/>`_. However, there are other dependences such as `NumPy <https://numpy.org/>`_, `Astropy <https://www.astropy.org/>`_, `matplotlib <https://matplotlib.org/>`_, `regions <https://pypi.org/project/regions/>`_, `SciPy <https://scipy.org/>`_, `pandas <https://pandas.pydata.org/>`_, `PyPDF2 <https://pypi.org/project/PyPDF2/>`_, and `os <https://docs.python.org/3/library/os.html>`_. Please, make sure that all these dependences are properly installed.

    
Example usage
-------------

.. code-block:: python

    from kanchay import kan
    
    # Search TESS light curves of a given star
    # The tool downloads and plots all the SPOC light curves (LC) observed by TESS in several sectors
    # The tool also normalizes and applies sigma clipping to the LCs. The resultant LC is stored in arrays in x (time), y (flux) and yerr (flux error).
    starid='Gaia DR2 2395082091638750848'
    x, y, yerr = kan.tess_lk(starid, exptime=120, author='SPOC')
    
    # You can also plot a region of interest of the LC
    import matplotlib.pyplot as plt
    plt.scatter(x[0][2000:8500], y[0][2000:8500])
    
    # You can determine rotational periods with only one command (thanks to the starspot code)
    # The rotational period is estimated using three methods (LS, ACF, PDM)
    kan.rotation_calc(x[0][2000:8500], y[0][2000:8500], yerr[0][2000:8500])
    
    # You can also determine rotational periods using GP (thanks to the exoplanet code)
    kan.rotation_calc(x[0][2000:8500], y[0][2000:8500], yerr[0][2000:8500], gp='yes')
    
    #For more details please see the kanchay's tutorial
    

Contributing
------------
*kanchay* is a simple tool created to help undergraduate students measure stellar rotational periods in an easy and simple way. Therefore, the tool needs input to improve. Please contact me (ramstojh@alumni.usp.br) if you have questions. Users are welcome to propose new features or report bugs by opening an issue on GitHub.


Authors
-------
- `Jhon Yana Galarza <https://github.com/ramstojh>`_


License & attribution
---------------------

Copyright 2022, Jhon Yana Galarza.

The source code is made available under the terms of the MIT license.

If you make use of this code, please cite this package and its dependencies.
