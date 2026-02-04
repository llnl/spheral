.. _Sedov_database:

=========================
Constructing the DataBase
=========================

The next section of code constructs a Spheral ``DataBase`` and adds our NodeList ``nodes`` to it.

.. code-block::
   :linenos:
   :lineno-start: 142

   #-------------------------------------------------------------------------------
   # Construct a DataBase to hold our node list
   #-------------------------------------------------------------------------------
   db = DataBase()
   db.appendNodeList(nodes)
   output("db")
   output("  db.numNodeLists")
   output("  db.numFluidNodeLists")

In Spheral a DataBase is a container of the NodeLists you want to run a calculation on, and is what you hand to Physics packages, the time integrator, etc. to tell them what NodeLists (or materials) are in the problem and should be modeled. In general you might have any number of NodeLists, each of which should be added to the DataBase in order for them to be ''active'' in the sense that Spheral will run physics problems on them. If you want to have a NodeList in your Python script for some other purpose but don't want to actively run any physics model on that NodeList, don't add it to the DataBase: only NodeLists in a DataBase are used in the actual physics models.
