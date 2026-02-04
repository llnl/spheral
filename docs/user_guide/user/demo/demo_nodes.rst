.. _nodelist_construction:

=====================
NodeList construction
=====================

The next section of the script constructs the NodeList we will need to hold the SPH nodes that make up our calculation:

.. code-block::
   :linenos:
   :lineno-start: 85

   #-------------------------------------------------------------------------------
   # Create a NodeList and associated Neighbor object.
   #-------------------------------------------------------------------------------
   nodes = makeSolidNodeList(name = "gamma-law gas points",
                             eos = eos, 
                             kernelExtent = WT.kernelExtent,
                             nPerh = nPerh)


Note at this point we are building a Spheral ``NodeList`` object (a ``SolidNodeList`` to be precise), but we have actually not created any nodes yet.  This NodeList object is the object we will populate with points to represent the geometry and initial conditions of our problem to evolve from in a later section of the script (:ref:`point_generation`), but at this point we are simply building our choice of NodeList (which could be one of several types).

The ``SolidNodeList`` can be used to represent either a fluid (as is the case here) or a solid material. In general you will specify a NodeList such as this for each separate material or part in your model. For instance if we had two different gasses in our calculation (such as in a fluid mixing calculation) with different equations of state we would want at least two NodeLists, each of which would then be given one of distinct equations of state.  Even in a case where we have multiple regions consisting material with the same equation of state we may want different NodeLists for each region, since the interface between materials/regions can be treated distinctly if Spheral knows they are different parts (depending on the choice of hydrodyanmics algorithm).

In this case however we only have a single material consisting of our simple gamma-law gas, so we create one NodeList to hold all the points in the calculation. Each NodeList must have a unique name (the first argument) -- if you try to give the same string as the name of two different NodeLists Spheral will raise an error condition.  The succeeding arguments give the equation of state for the material (required for ``SolidNodeList`` or ``FluidNodeList``), the ``kernelExtent`` or dimensionless distance over which the interpolation kernel is non-zero, and target nodes per smoothing scale (``nperh``) describing how the smoothing scale for each point in the NodeList should be adjusted.  Note the ``SolidNodeList`` also can take a strength model as an argument for a material that has material strength; if none is provided the ``SolidNodeList`` defaults to a strengthless fluid that behaves identically to a ``FluidNodeList``.

All these concepts are discussed in more detail later in this document.  See for instance :ref:`NodeLists` for more discussion of the different types of NodeLists supported in Spheral.
