.. _user-intro:

###################################################
The first thing to know: Spheral is a Python Module
###################################################

Spheral is a Python extension rather than a standalone code.  This is an unusual choice in the world of physics modeling software, which are generally provided as standalone software packages with their own interface languages and conventions.  Spheral was developed as a Python extension for a few reasons:

 - Python provides a complete, well-developed, and mature language which many users are already familiar with.

 - It is very convenient to develop sophisticated algorithms in Python for tasks such as laying down your initial conditions in the manner most convenient for your problem at hand.  This provides experienced users with the freedom to develop and experiment with the most effective ways to set up their problems in a fully featured programming environment.

 - Python provides a familiar environment for users to analyze, explore, and interpret their simulation outcomes, with direct access to all the physics variables and state internal to Spheral from the Python interface.

 - Users can take advantage of the many other wonderful Python packages in their Spheral scripts, such as NumPy, SciPy, Matplotlib, and so on.

The choice of Spheral existing as a Python module has consequences that, while natural for users accustomed to Python, might be surprising for folks coming from a more traditional standalone physics modeling software background.  First up, the Spheral executable ``spheral`` is simply an alias for the Python executable, so in order to access any of Spheral's compiled assets you must first import the Spheral Python module.  Importing Spheral will print a banner with some useful information, and changes the ordinary Python prompt to make it clear Spheral has been loaded; here is an example starting from an ordinary ``bash`` prompt starting up an interactive Spheral session::

    bash{owen}41: spheral
    Python 3.9.10 (main, Jun  7 2023, 11:19:53) 
    [GCC 10.3.1 20210422 (Red Hat 10.3.1-1)] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import Spheral
    /------------------------------------------------------------------------------\
    |  Spheral version             : v2024.01.0 017f9e506 HEAD                     |
    |    number of MPI tasks       : 1                                             |
    |    number of threads per rank: 112                                           |
    \------------------------------------------------------------------------------/
    Spheral> dir()
    ['Spheral', '__annotations__', '__builtins__', '__doc__', '__loader__', '__name__', '__package__', '__spec__']
    Spheral> quit()
    bash{owen}42: 

Another point to keep in mind is that Spheral is not watching over your shoulder when you execute ordinary Python commands, and therefore in order to affect the outcome of a Spheral simulation you need to directly manipulate Spheral objects.  For instance, in a standalone code very likely typing a statement like ``gamma = 1.4`` means you have set some physical quantity ``gamma`` in the problem, which might be the ratio of specific heats in a gamma-law gaseous equation of state .  At the Spheral Python prompt however this simply declares an ordinary Python variable ``gamma`` to be a floating value of ``1.4``.  If the intent is to set :math:`\gamma = 1.4` in an equation of state, we must explicitly manipulate a Spheral object::

  # Create a gamma-law gas equation of state with gamma = 5.0/3.0 and mean molecular weight of 1.8 in CGS units
  Spheral> eos = GammaLawGas(gamma = 5.0/3.0, mu = 1.8, constants = CGS())
  Spheral> print(eos.gamma)
  1.6666666666666667

  # Now for some reason we changed our mind about gamma
  Spheral> eos.gamma = 1.4
  Spheral> print(eos.gamma)
  1.4
  Spheral> 

As mentioned above a major advantage of choosing to use Python as our interface is that knowing Python immediately gives a potential Spheral user a huge advantage in using the code.  However, for those not familiar with Python probably the first order of business is to gain a working knowledge of Python.  There is a wealth of information available in order to learn Python, from on-line references, video tutorials, or good old fashioned books.  If you are not familiar with Python you should seek out some of these resources and at least learn the basics: a good place to start is the documentation section of the main `Python site <https://www.python.org/doc/>`_.
