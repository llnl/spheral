.. _Sedov_boundaries:

===================
Boundary conditions
===================

Most problems have boundary conditions of one sort or another. There are quite a few boundary options in Spheral, most of which involve establishing planes over which those boundaries apply. In this example we are running one quadrant of the full physical situation, so we want to create two reflecting boundaries along the :math:`x=0` and :math:`y=0` planes:

.. code-block::
   :linenos:
   :lineno-start: 164

   #-------------------------------------------------------------------------------
   # Create boundary conditions.
   #-------------------------------------------------------------------------------
   xPlane0 = Plane(Vector(0.0, 0.0), Vector(1.0, 0.0))
   yPlane0 = Plane(Vector(0.0, 0.0), Vector(0.0, 1.0))
   xbc0 = ReflectingBoundary(xPlane0)
   ybc0 = ReflectingBoundary(yPlane0)
   for p in packages:
       p.appendBoundary(xbc0)
       p.appendBoundary(ybc0)

Lines 167--168 create our planes in point-normal form: for instance, ``xPlane0`` uses a point at coordinates ``(0.0, 0.0)`` (the first argument) and inward-pointing normal ``(1.0, 0.0)`` (the second argument).  Boundary conditions are set on a Physics package by Physics package basis in case there are some unique to a particular package.  In general however you most likely will be applying the same boundary conditions to all packages such as shown here in the loop on lines 171--173.
