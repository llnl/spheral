Python and C++ container bindings
=================================

Spheral uses PYB11Generator to describe which C++ classes, functions, and
containers Python can use. C++ containers, such as ``std::vector<T>`` and
``std::map<Key, Value>``, store collections of C++ values. PYB11Generator
creates C++ binding code from Python descriptions under ``src/PYB11``;
pybind11 implements the connection to Python.

Three objects can have the same contents but different behavior:

* A Python list, such as ``[1, 2]``, owns Python elements.
* A C++ ``std::vector<int>`` owns C++ integer elements.
* A Python ``Spheral.vector_of_int`` object represents a C++
  ``std::vector<int>``. Depending on how it was created, it owns a vector or
  refers to one owned by another object. It supports operations such as
  indexing, iteration, and ``append``; it is a different Python type from
  ``list``.

A **binding** exposes a C++ type or operation to Python. A pybind11
``type_caster`` is the Python-C++ conversion mechanism: it converts an argument
from Python to a C++ value or reference, or a result in the other direction.

``Spheral.vector_of_int`` above is a **bound Python class**: a class supplied by
the compiled Spheral extension to represent the C++ ``std::vector<int>`` type.
It is not a Python built-in. The developer requests it in
``src/PYB11/CXXTypes/CXXTypes_PYB11.py``:

.. code-block:: python

   vector_of_int = PYB11_bind_vector("int", opaque=True, local=False)

The vector class and its conversion definition have separate sources:

* ``opaque=True`` marks ``std::vector<int>`` as **opaque** in the binding
  description. An opaque declaration disables pybind11's automatic conversion
  of that C++ type to and from Python lists. For Spheral, the developer must
  also add the corresponding ``PYBIND11_MAKE_OPAQUE`` declaration to
  ``src/PYB11/SpheralTypeCasters.hh``. The compiler reads that shared header;
  module-specific opaque declarations emitted by the generator are suppressed.
* ``PYB11_bind_vector`` requests a ``py::bind_vector`` registration call.
  PYB11Generator writes that C++ call during the build; pybind11 executes it
  when the compiled extension is imported, creating the Python class named
  ``vector_of_int``. The description itself is not that Python class.

``local=False`` becomes ``py::module_local(false)`` in that registration call.
The association between the Python class and the C++ type is shared across
extension modules using the same pybind11 type registry. For example, another
Spheral module accepting ``std::vector<int>`` can use the registered
``vector_of_int`` class. ``local=True`` instead permits separate module-local
class registrations for the same C++ type.

Together, the opaque declaration and class registration make a C++ vector
result appear in Python as a ``Spheral.vector_of_int`` object, rather than a
list. Spheral separately enables implicit vector input conversion: when
permitted by the binding, a list argument creates a temporary C++ vector for
the call without an explicit ``vector_of_int(...)`` constructor call in Python.

Adding a class or data type
---------------------------

The following recipe adds an illustrative, dimension-independent
``Spheral::Widget`` to the existing Neighbor package. All paths are relative
to the Spheral repository root; use the corresponding package paths for a
class belonging elsewhere.

1. **Implement and register the C++ type.** Add
   ``src/Neighbor/Widget.hh`` and, if needed, ``src/Neighbor/Widget.cc``.
   In ``src/Neighbor/CMakeLists.txt``, add ``Widget.hh`` to
   ``Neighbor_headers`` and ``Widget.cc`` to ``Neighbor_sources``.
   Header-only types need no new ``.cc`` entry.
2. **Describe the Python class.** Add ``src/PYB11/Neighbor/Widget.py``.
   For a C++ ``Widget(double)`` constructor and ``double weight() const``
   method, the description is:

   .. code-block:: python

      from PYB11Generator import *

      @PYB11namespace("Spheral")
      class Widget:
          def pyinit(self, weight="double"):
              return

          @PYB11const
          def weight(self):
              return "double"

3. **Connect it to the package module.** In
   ``src/PYB11/Neighbor/Neighbor_PYB11.py``, add
   ``'"Neighbor/Widget.hh"'`` to ``PYB11includes`` and
   ``from Widget import Widget`` alongside the other class imports.
   If the Python interface exposes vectors of this type, also add:

   .. code-block:: python

      vector_of_Widget = PYB11_bind_vector(
          "Spheral::Widget", opaque=True, local=False)
      vector_of_vector_of_Widget = PYB11_bind_vector(
          "std::vector<Spheral::Widget>", opaque=True, local=False)

   The nested-vector binding is needed only when that container is exposed.
   Reuse an existing global container registration instead of registering the
   same C++ container class in a second module.
4. **Add the shared container declarations.** In
   ``src/PYB11/SpheralTypeCasters.hh``, add ``class Widget;`` inside the
   existing ``namespace Spheral`` forward-declaration block. Add the
   corresponding declarations alongside the Neighbor containers:

   .. code-block:: cpp

      PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Widget>);
      PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Widget>>);

   Add only the container forms used by the interface. A pointer vector such
   as ``std::vector<Spheral::Widget*>`` or a map is a separate C++ type and
   needs its own declaration and, if Python should construct it, its own
   vector or map registration. Binding ``Widget`` alone does not require
   an opaque container declaration.
5. **Regenerate and check the interface.** Reconfigure when the new class
   changes the generated source-file list, then regenerate, build, and install
   the bindings. The existing
   ``src/PYB11/Neighbor/CMakeLists.txt`` already registers this module through
   ``spheral_add_pybind11_library``. Check construction and methods, the
   requested vector classes, compatible and incompatible list inputs, and the
   actual return, mutation, and ownership behavior.

For a dimension-dependent class, follow the template description in
``src/PYB11/Neighbor/GridCellIndex.py`` and the
``PYB11TemplateClass`` instantiations in
``src/PYB11/Neighbor/Neighbor_PYB11.py``. Write the fully qualified
container declarations in ``src/PYB11/SpheralTypeCasters.hh`` for the relevant
1D, 2D, and 3D types, following its existing families. The shared header covers
all three dimensions because common and Geometry declarations already do;
do not trim it to the dimensions enabled in one build.

For an enum or another type described directly in the package module, edit
``src/PYB11/Neighbor/Neighbor_PYB11.py`` rather than adding a separate
class-description file. General-purpose containers of primitive types belong
in ``src/PYB11/CXXTypes/CXXTypes_PYB11.py``, with matching declarations in
``src/PYB11/SpheralTypeCasters.hh``. A C++ alias of an existing type, such as
an alias of ``double``, reuses its conversion and any container registration
for the underlying type; it does not create another distinct container type.

Shared-header maintenance
~~~~~~~~~~~~~~~~~~~~~~~~~

Keep the package's container descriptions and
``src/PYB11/SpheralTypeCasters.hh`` in the same change. This checked-in header
is maintained explicitly; configuration and builds do not regenerate its
declarations. Adding ``opaque=True`` or a ``PYB11opaque`` entry only in Python
does not update Spheral's shared conversion definitions.

Use fully qualified C++ names in the shared header. Its declarations appear
before the package's ``PYB11includes`` and ``using namespace`` statements.
Introduce a new class name with a forward declaration before the opaque
declarations. An alias or nested type can require a supporting C++ include
there instead; adding that include only to the package module can be too late.
The code that binds a vector also needs the complete element definition.
Avoid including Python conversion code before the shared declarations.

``PYB11opaque`` lists remain available for containers without a corresponding
Python class registration. Existing shared entries live in
``src/PYB11/SpheralCommon.py``; package-specific lists also exist. To add such
a type, update the appropriate list and add its fully qualified declaration
to ``src/PYB11/SpheralTypeCasters.hh``. The list retains the generator's
module-specific declaration path when a shared header is not selected.
An opaque declaration alone creates no Python container class and enables
no implicit list input. Registered vectors and maps already describe their
module-specific opaque declarations, so they need no extra list entry.

When removing or changing a container declaration, check all package
descriptions and handwritten conversion consumers before removing it from
the shared header. Different aliases of the same C++ type must not introduce
duplicate specializations.

Configuration, build, and runtime
---------------------------------

Configuration: CMake records the work
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``spheral_add_pybind11_library`` in ``cmake/spheral/SpheralAddLibs.cmake``
connects ``<Package>_PYB11.py`` to ``PYB11Generator_add_module`` and sets the
Python import path and module options. The Generator's CMake integration is
in ``extern/PYB11Generator/cmake/PYB11Generator.cmake``.

For multiple-file output, CMake runs
``extern/PYB11Generator/cmake/generate_cpp.py`` with ``dry_run=True``.
It imports the descriptions to discover the generated
source filenames without emitting the binding implementations. CMake also
records imported Python files as configuration dependencies. Adding a class
can change the source list and require reconfiguration. For single-file
output, the source filename is already known.

CMake creates the ``Spheral<Package>_src`` source-generation target.
Spheral's target helper defines ``PYB11_SHARED_TYPE_CASTER_HEADER`` as
``"PYB11/SpheralTypeCasters.hh"`` for each binding target. This compiler
definition selects the checked-in shared header; it does not generate it.

Build: the generator writes code, then the compiler compiles it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The source-generation target runs
``extern/PYB11Generator/cmake/generate_cpp.py`` with ``dry_run=False``.
It imports the configured package description, calls ``PYB11generateModule``,
and replaces generated sources whose contents changed.
``extern/PYB11Generator/PYB11Generator/PYB11class.py`` turns the class
description above into a registration function equivalent to:

.. code-block:: cpp

   void bindWidget(py::module_& m) {
     py::class_<Spheral::Widget> obj(m, "Widget");
     obj.def(py::init<double>(), py::arg("weight"));
     obj.def("weight", &Spheral::Widget::weight);
   }

The example omits holder types, explicit method-pointer casts, and optional
registration attributes.
``extern/PYB11Generator/PYB11Generator/PYB11STLmethods.py`` turns the vector
description into code equivalent to:

.. code-block:: cpp

   py::bind_vector<std::vector<Spheral::Widget>>(
       m, "vector_of_Widget", py::module_local(false));
   #ifdef PYB11_VECTOR_BINDING_HOOK
   PYB11_VECTOR_BINDING_HOOK<std::vector<Spheral::Widget>>();
   #endif

Each generated source includes its generated module header. That header
conditionally includes ``PYB11/SpheralTypeCasters.hh``, so the compiler sees
the explicit opaque definitions before compiling package-specific conversion
and binding code. The macro selects Spheral's vector hook, which calls
``Spheral::registerPythonVectorConversion<std::vector<Spheral::Widget>>()``.

The compiler compiles these sources, and the linker combines the binding code
with the C++ implementation into the extension. Installation copies the
checked-in shared header under ``include/PYB11``, as specified in
``src/PYB11/CMakeLists.txt``.

Runtime: importing the compiled extension executes registration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When Python imports the compiled extension, its module initialization executes
the generated registration code. A standalone module uses
``PYBIND11_MODULE``; the aggregated Spheral extension calls generated
``bindSubmodule_<name>`` functions. Those functions call ``bindWidget`` and
the vector registration code. ``py::class_`` registers the Python class and
its C++ type; ``obj.def`` registers constructors and methods.
``py::bind_vector`` registers the Python vector class, after which Spheral's
helper registers implicit vector input conversion.

A later Python call to ``Widget(...)`` invokes the registered C++
constructor; ``widget.weight()`` invokes the registered C++ method.
PYB11Generator does not inspect the descriptions again or generate bindings
during application runtime.

Existing projects without shared definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The shared-header include and post-vector-registration hook are independently
selected interfaces. PYB11Generator does not enable Spheral's policy for other
projects automatically:

* Without ``PYB11_SHARED_TYPE_CASTER_HEADER``, the shared include is inactive.
  The generated module header's ``PYB11opaque`` entries and ``opaque=True``
  container declarations remain active under their preprocessor guard.
* Without ``PYB11_VECTOR_BINDING_HOOK``, the conditional call after
  ``py::bind_vector`` is inactive. Importing the extension still creates the
  vector class, but this hook adds no implicit input conversion. A project
  can supply its own hook.

.. list-table:: Vector conversion with and without the shared Spheral helper
   :header-rows: 1
   :widths: 35 33 32

   * - Arrangement
     - Ordinary vector argument
     - Vector result
   * - Automatic pybind11 conversion, without an opaque declaration for the type
     - An accepted Python input is converted into separate C++ vector storage.
     - A Python list populated from the C++ vector.
   * - A module-specific opaque vector declaration and bound vector class,
       without an additional implicit conversion
     - A matching bound vector. A list can explicitly construct that vector,
       but is not granted implicit function-argument conversion.
     - A bound vector, with ownership determined by the return value policy.
   * - Spheral's shared opaque declaration, bound vector class, and input helper
     - A matching bound vector or a permitted implicitly converted input,
       such as a list. Converting a list creates temporary vector storage.
     - A bound vector, with ownership determined by the return value policy.

For example, an opaque-vector API without an input helper can accept
``consume(vector_of_int([1, 2]))`` even though it does not accept
``consume([1, 2])`` directly. Explicit construction and implicit argument
conversion are different operations.

Why every source file needs the same definitions
------------------------------------------------

A **translation unit** is a C++ source file together with the headers included
while compiling it. C++'s **One Definition Rule (ODR)** requires definitions of
the same class or template specialization to agree across translation units.
For a vector declared opaque in one module and automatically converted to a
list in another, pybind11 would otherwise define two different
``type_caster<std::vector<T>>`` specializations.

Every generated C++ source first includes its generated module header. That
module header includes the shared header after pybind11 headers, including
``stl.h``, and before module-specific C++ headers, declarations, or binding
code. The shared header's relevant order is:

.. code-block:: cpp

   #include <pybind11/pybind11.h>
   #include <pybind11/stl.h>
   // Supporting C++ includes and element forward declarations come next.
   PYBIND11_MAKE_OPAQUE(std::vector<int>);
   // Other explicit container declarations follow.
   // The implicit vector input helper follows the opaque declarations.

When the macro selects the shared header, generated module-specific opaque
declarations are excluded by ``#ifndef PYB11_SHARED_TYPE_CASTER_HEADER``.
Vector and map class registrations are still emitted. Handwritten conversion
headers ``src/FileIO/PyFileIO.hh`` and
``src/PYB11/DataOutput/RestartableObject.hh`` also include the shared header
before using container conversion code. New handwritten conversion code must
see the same definitions before instantiating a relevant ``type_caster``.

**Link-time optimization (LTO)** allows the compiler to optimize across compiled
source files. It exposed the inconsistent conversion definitions; it did not
create permission to use different definitions without LTO. Binding targets
retain LTO options inherited from their dependencies. In the Ubuntu 26.04
environment, MPICH supplies such options through dependency targets.

How a Python list becomes a vector argument
-------------------------------------------

``SpheralTypeCasters.hh`` defines ``PYB11_VECTOR_BINDING_HOOK`` to select
``Spheral::registerPythonVectorConversion<Vector>()``. This helper registers
``pybind11::implicitly_convertible<PythonVectorInput, Vector>()``.
``PythonVectorInput`` is a C++ marker type used to recognize permitted input
categories; it holds no converted data and has no Python-visible class.

For an ordinary Python call to a bound C++ function, argument conversion works
as follows:

* A matching bound vector already identifies the required C++ vector.
  pybind11 obtains that object directly.
* For another input, implicit conversion must be allowed for the argument.
  The marker's ``load`` method uses
  ``pybind11::detail::object_is_convertible_to_std_vector`` to check the input
  category. pybind11 then invokes the target vector class's iterable
  constructor. That constructor converts each element and inserts it into a
  new C++ vector.
* pybind11 keeps the temporary vector object alive during the C++ call and
  passes the vector according to the C++ parameter type. Its storage is
  separate from the original Python container.

Input behavior for value and reference parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following table applies to registered vector types and ordinary bound
function calls that permit implicit conversion. ``T`` must be convertible to
the vector's element type.

.. list-table:: Python inputs and the vector received by C++
   :header-rows: 1
   :widths: 22 38 40

   * - C++ parameter
     - Python list or tuple
     - Matching bound vector
   * - ``const std::vector<T>&``
     - C++ reads a temporary vector populated from the Python input.
     - C++ reads the vector represented by the Python object.
   * - ``std::vector<T>&``
     - C++ can modify the temporary. Those container changes do not update the
       Python list or tuple.
     - C++ modifies the represented vector. Python sees those changes.
   * - ``std::vector<T>``
     - C++ receives a vector value populated from the input. Container changes
       inside the function do not update the input.
     - The C++ value parameter is a copy of the represented vector. Container
       changes inside the function do not update the caller's vector.

Copying a container is not necessarily a deep copy of its elements; pointer
and object-element behavior is discussed below.

Here are illustrative C++ signatures:

.. code-block:: cpp

   int sumValues(const std::vector<int>& values);
   void appendOne(std::vector<int>& values);
   std::vector<int> appendedCopy(std::vector<int> values);
   std::vector<int> makeValues();

Assume ``sumValues`` sums its input, ``appendOne`` appends ``3``,
``appendedCopy`` appends ``3`` and returns its local vector, and ``makeValues``
returns ``{1, 2}``. ``api`` below denotes the module exposing those illustrative
functions; these names are not additional Spheral functions.

.. code-block:: python

   from Spheral import vector_of_int

   items = [1, 2]
   assert api.sumValues(items) == 3
   api.appendOne(items)
   assert items == [1, 2]             # The temporary received the appended 3.

   values = vector_of_int([1, 2])
   api.appendOne(values)
   assert list(values) == [1, 2, 3]   # C++ changed the represented vector.

   original = vector_of_int([1, 2])
   result = api.appendedCopy(original)
   assert list(original) == [1, 2]   # The C++ parameter was passed by value.
   assert list(result) == [1, 2, 3]
   assert isinstance(result, vector_of_int)

A function that stores a pointer or reference to its argument after returning
needs an explicit lifetime design. It cannot safely retain a vector created
temporarily from a Python list. If that API requires persistent vector storage,
require a bound vector and ensure its owner stays alive. Use
``py::arg("values").noconvert()`` in an explicit binding, or
``@PYB11noconvert`` in a PYB11 function or method description, to disable
implicit argument conversion. The PYB11 decorator applies to that description's
arguments. This restriction alone does not keep the supplied vector alive.

Other inputs and overload selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With the inspected pybind11 3.0.1 implementation, the category check accepts
lists, tuples, bound vectors, other Python sequence objects, generators, sets,
frozensets, dictionary views, ``map``, and ``zip`` objects. Strings, bytes, and
ordinary dictionaries are rejected. Not every iterable passes this category
check; for example, a list iterator is different from a generator. The
vector's explicit iterable constructor and implicit function-argument
conversion therefore need not accept precisely the same objects.

Accepted elements must still convert to the required C++ type. An empty vector
of another element type can pass conversion to an empty target vector, but a
populated vector fails when an element cannot convert. Generator conversion
consumes its elements, and converting a set does not promise a meaningful
order. Overloaded functions can choose an overload without conversion before
trying implicit conversions; an accepted input category alone does not
determine which overload is called. Review the internal category predicate
when upgrading pybind11.

How a C++ vector becomes a Python result
----------------------------------------

For an opaque vector, the C++ return type and the binding's **return value
policy** determine whether Python owns independent storage or refers to
existing storage. The policy specifies how pybind11 copies, borrows, or owns
the returned object; accepting a list argument does not change the result type.

.. list-table:: C++ results and the Python object
   :header-rows: 1
   :widths: 30 38 32

   * - C++ return and binding policy
     - Python receives
     - Effect of Python container mutation
   * - ``std::vector<T>`` returned by value, with normal automatic handling
     - A bound vector owning the returned value, normally transferred by move.
     - Changes affect that returned vector, not a source container copied by
       the C++ function.
   * - ``std::vector<T>&`` or ``const std::vector<T>&``, with ``copy``
     - A bound vector with copied container storage.
     - Changes affect the copy, not the original C++ vector.
   * - A member vector reference, with ``reference_internal``
     - A bound vector referring to the member. Python keeps the owning
       instance alive while the result is alive.
     - Changes through mutable vector operations affect the C++ member.
   * - A vector reference or pointer, with ``reference``
     - A bound vector borrowing existing storage.
     - Changes through mutable operations affect the existing vector. The
       developer must arrange for its owner to outlive every use.
   * - An explicit adapter returning ``py::list``
     - A Python list populated by the adapter.
     - List insertion, removal, and replacement do not modify the vector.

For ordinary bound methods, the default automatic policy copies an lvalue
reference result and can take ownership of a raw pointer result. A pointer to
a member or externally owned vector must not accidentally be treated as newly
owned storage. Properties have their own default reference policy, so method
and property getters should be considered separately. An existing Python
wrapper for the same C++ object can also be reused; a policy is not a command
to create a new wrapper every time.

Keeping the owning instance alive does not protect references to elements
when the vector reallocates, or storage the owner destroys independently.
C++ ``const`` on a vector reference does not create a separate immutable Python
vector class. If Python must not mutate the original C++ vector, return a copy.
If the Python result itself must be immutable, provide an explicit read-only
interface rather than relying on ``const`` alone.

For the illustrative ``makeValues`` function:

.. code-block:: python

   result = api.makeValues()
   assert isinstance(result, vector_of_int)
   result.append(3)
   assert list(result) == [1, 2, 3]

   snapshot = list(result)
   snapshot.append(4)
   assert list(result) == [1, 2, 3]   # The Python list has separate storage.

Property assignment and property reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider a class with these C++ methods, where the setter copies its argument
into a member vector:

.. code-block:: cpp

   const std::vector<int>& values() const;
   void values(const std::vector<int>& replacement);

A corresponding excerpt from its PYB11 class description is:

.. code-block:: python

   values = PYB11property(
       "const std::vector<int>&", "values", "values",
       returnpolicy="reference_internal")

For an instance ``obj`` of that illustrative class:

.. code-block:: python

   items = [1, 2]
   obj.values = items                 # Convert, then let the setter copy.
   items[0] = 99
   assert list(obj.values) == [1, 2]   # The member does not alias the list.

   view = obj.values                  # The getter exposes the member vector.
   view[0] = 7
   assert list(obj.values) == [7, 2]

Assigning a matching bound vector also makes this setter copy into the
member; the member does not share the supplied vector's storage.

If the public Python API requires list-valued reads, use an explicit adapter.
For integer vectors, the helpers defined in the ``SpheralCommon.py`` C++
preamble can support a getter that returns a list:

.. code-block:: python

   @PYB11const
   @PYB11implementation(
       "[](const Spheral::VectorExample& self) { "
       "return Spheral::PYB11utils::to_list(self.values()); }")
   def valuesList(self):
       return "py::list"

``VectorExample`` is an illustrative class with the ``values`` getter above.
Its adapter returns an integer list snapshot. A list-taking setter can use
``Spheral::PYB11utils::from_list<int>`` and then copy into the member.
Setting ``opaque=False`` in another module cannot override the shared
declaration for the C++ vector type.

Nested vectors, object elements, and maps
-----------------------------------------

An outer and an inner vector are different C++ types and need their own
registrations. ``CXXTypes_PYB11.py`` registers, for example, both
``std::vector<int>`` and ``std::vector<std::vector<int>>``. A way to construct
a matching nested bound vector without relying on conversion of inner lists is:

.. code-block:: python

   from Spheral import vector_of_int, vector_of_vector_of_int

   rows = vector_of_vector_of_int([
       vector_of_int([1, 2]),
       vector_of_int([3, 4]),
   ])

Each outer-vector element is a C++ vector value, so inserting a bound row
copies that row's container storage. Nested bound-vector results remain bound
vectors at the corresponding levels. ``list(rows)`` creates only an outer
Python list; it does not recursively convert its rows to Python lists.
For integer rows, ``[list(row) for row in rows]`` creates independent lists at
both levels. Whether a list of lists is accepted on a particular constructor,
argument, or callback path needs its own check: inner element conversion and
temporary lifetime handling introduce additional steps.

For ``std::vector<Widget>``, element construction follows ``Widget``'s C++ copy
and Python conversion behavior. For ``std::vector<Widget*>``, copying the
container copies pointers, not the pointed-to objects. A C++ function can change
a pointed-to object even when it received its pointer vector through a
temporary list conversion. Such changes can be visible through the original
Python element objects although list membership stays unchanged. Vector
ownership does not by itself preserve the lifetime of pointed-to objects.

The existing Integrator bindings illustrate the distinction between container
and element ownership: ``physicsPackages`` returns a
``const std::vector<Physics<Dimension>*>&`` with ``reference_internal``, while
``uniqueBoundaryConditions`` returns a vector of Boundary pointers by value.
The first refers to an Integrator member; the second returns independent
pointer-vector storage. Neither vector policy alone establishes ownership of
the Physics or Boundary objects.

For a map registered through
``PYB11_bind_map(..., opaque=True, local=False)``, add its matching opaque
declaration to ``src/PYB11/SpheralTypeCasters.hh`` too. The vector conversion
hook is not called after map registration, so it does not grant ordinary
dictionaries implicit conversion to map arguments. Use the matching bound map
or an explicit dictionary adapter.

Calls from C++ back into Python
-------------------------------

For a Python override of a C++ virtual method, PYB11Generator emits a C++
forwarding class, commonly called a **trampoline**. It finds the Python method,
converts the C++ arguments to Python objects, calls that method, and converts
its result back to the declared C++ return type.

Registered vector arguments exposed to the Python override use bound-vector
representations. Copying and borrowing depend on the forwarding call's argument
conversion; Python changes reach the caller only if the forwarded object
refers to that caller's vector. For a virtual method returning a vector by
value, an explicit matching bound vector gives the return conversion a C++
vector to copy:

.. code-block:: cpp

   virtual std::vector<int> chooseIndices() const;

.. code-block:: python

   def chooseIndices(self):
       return vector_of_int([1, 2])

Returning ``[1, 2]`` uses pybind11's override-return and temporary-object
lifetime handling. Successful ordinary list arguments do not establish this
path; implicit list returns have not been separately established for all
affected Spheral overrides on pybind11 3.0.1.

If the C++ virtual method instead returns a vector reference or pointer, the
Python implementation must preserve a suitable existing vector beyond the
call, for example as an instance attribute. A local list or local temporary
vector does not provide storage that a C++ reference return can safely retain.

Reviewing a binding change
--------------------------

Check the Python description, the shared declaration, generated code, and
runtime behavior as separate parts of the connection. In particular:

* Compare each new or changed opaque container description with
  ``src/PYB11/SpheralTypeCasters.hh``. Confirm that the shared header contains
  one declaration for each distinct C++ type and that element names are known
  before it. Check the generated include placement, class registration, and
  vector hook.
* Exercise matching bound vectors and list inputs with compatible and
  incompatible elements. For mutable references, check which object's
  contents changed. For by-value parameters, check that the caller's container
  did not change.
* Inspect return types, policies, property assignment, and property reads.
  Check nested container behavior and ownership of pointer elements separately.
  Treat callback and virtual-override returns as distinct conversion paths.
* When validating LTO behavior, inspect compilation objects and the final
  extension link to establish that LTO is active, and inspect diagnostics for
  One Definition Rule and LTO type-mismatch warnings. Passing runtime tests
  alone does not establish consistent conversion definitions.

The examples describe behavior for their stated signatures and policies;
they do not establish test coverage for every Spheral API using those types.

Implementation and library references
-------------------------------------

The principal Spheral source files are:

* ``src/PYB11/SpheralTypeCasters.hh``: supporting type declarations, explicit
  shared opaque declarations, and runtime vector input registration.
* ``cmake/spheral/SpheralAddLibs.cmake``: binding module integration and the
  compiler definition selecting the shared header.
* ``src/PYB11/CMakeLists.txt``: package selection, aggregated extension
  registration, and installation of the shared header.
* ``src/PYB11/SpheralCommon.py``: common includes, opaque lists, and explicit
  list-conversion helpers.
* ``src/FileIO/PyFileIO.hh`` and
  ``src/PYB11/DataOutput/RestartableObject.hh``: handwritten conversion
  consumers that also include the shared definitions.

Within ``extern/PYB11Generator``, the implementation files are:

* ``cmake/PYB11Generator.cmake``: configuration-time output discovery and
  source-generation targets.
* ``cmake/generate_cpp.py``: configured module import and source generation.
* ``PYB11Generator/PYB11class.py``: class, constructor, and method registrations.
* ``PYB11Generator/PYB11STLmethods.py``: vector and map registration and the
  conditional vector hook.
* ``PYB11Generator/__init__.py``: module initialization code, shared include
  placement, and guards around module-specific opaque declarations.

``extern/PYB11Generator/docs/stl.rst`` documents the generator's general
container-binding interfaces.

For pybind11's underlying behavior, see its documentation for
`container conversions and opaque types <https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html>`_
and
`return value policies and implicit conversions <https://pybind11.readthedocs.io/en/stable/advanced/functions.html>`_.
The inspected 3.0.1 implementations are in ``include/pybind11/stl.h`` for input
categories, ``stl_bind.h`` for vector construction and access, ``pybind11.h``
for implicit conversion registration and override forwarding, and
``detail/type_caster_base.h`` for bound-object conversion and temporary lifetime
handling. Check the version actually used by the affected Spheral build when
reviewing a library upgrade or an API-specific failure.
