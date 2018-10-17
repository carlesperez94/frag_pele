.. _sec-conditionsByExample:

*********************************
Conditions in Pele++ control file
*********************************

The conditions are comprised of Sensor tags, values, arithmetic
operators, logical operators, a special frequency condition description,
and parenthesis.

Examples
========

Let's see some examples.

Example 1
---------

This is the most basic logical condition possible.

.. code-block:: text

   rmsd1 < 5.5

It is comprised of a Sensor tag ("rmsd1"), an arithmetic operator ("<")
and a value ("5.5").

The corresponding Sensor must have been defined and tagged somewhere
else in the configuration file, (see
:ref:`Sensors <sec-sensors>`).

This logical condition will evaluate to true whenever the current value
of the Sensor tagged as "rmsd1" is less than 5.5.

The defined arithmetic operators are: "<", "<=", ">", ">=" and "=", but
"=" is only defined for integer values.

Example 2
---------

.. code-block:: text

   rmsd1 < 5 or rmsd2 > 1

This condition uses the logical operator "or", to compose the two more
simple logical conditions: "rmsd1 < 5" and "rmsd2 > 1", so that, it will
evaluate to true if either of the two more simple conditions evaluate to
true.

There are three logical operators in Pele++: "not", "and" and "or".

Their preference is: "not" >> "and" >> "or".

If you want to change the order in which the operators are evaluated,
you must use parenthesis, ().

Example 3
---------

Here we use "not" and "or" logical operators.

.. code-block:: text

   not rmsd1 < 5 or rmsd2 > 1

According to the given operator preference, that is the same as:

.. code-block:: text

   (not rmsd1 < 5) or (rmsd2 > 1)

Example 4
---------

In this example, we are using parenthesis, to change the order in which
the operators are evaluated.

.. code-block:: text

   not (rmsd1 < 5 or rmsd2 > 1)

Here, the expression between parentheses is evaluated first and then its
result is negated using the not operator. Compare this with the default
behavior shown in the previous example.

Example 5
---------

More of the same

.. code-block:: text

   not (rmsd1 < 5 or not rmsd2 > 1)

which is the same as:

.. code-block:: text

   not (rmsd1 < 5 or not (rmsd2 > 1))

Example 6
---------

In this example we start using the "and" operator

.. code-block:: text

   rmsd1 < 5 and rmsd2 > 2 or rmsd3 > 1

According to the logical operator preference this is the same as:

.. code-block:: text

   (rmsd1 < 5 and rmsd2 > 2) or rmsd3 > 1

Notice how the precedence of the "and" operator is higher that the
preference of the "or" operator. That's why the "and" operator is
evaluated before.

Example 7
---------

Here we use the "=" arithmetic operator.

.. code-block:: text

   steps = 4

Remember that this operator can only be used with Sensors that measure
integer values. For instance, "rmsd = 4.1" would be incorrect.

Example 8
---------

Using the available logical and arithmetic operators, we can build more
complex conditions like the following:

.. code-block:: text

   not (rmsd1 < 5.0 or not energy >= 3) and bindingEnergy <= 0.1

.. _sec-conditionsByExample-example9:

Example 9
---------

You can express a condition based on a frequency. You need a sensor that
represents a counter (for example, the "currentStep" tracker). The
condition will be true with the provided frequency (expressed by the
period). The counter is assumed to start at 1, so the condition will be
true at 1, 1+period, 1+2\*period, etc. You can add an offset to the
start counter value (which is by default 0), so that instead of starting
at 1 it starts at 1+offset.

This condition uses the special keywords *"_every_"* and *"_offset_"*.
Notice the underscores at the beginning and end of the keywords: they
are needed so that the keyword is recognized.

The following condition uses the default offset and a period of 5:

.. code-block:: text

   currentStep _every_ 5

This condition will be true for currentStep 1, 6, 11, etc.

If you want an offset of 2 (that is, starting at the counter value 3),
you would say:

.. code-block:: text

   currentStep _every_ 5 _offset_ 2

This condition will be true for currentStep 3,8,13, etc.

Notice that for the specific case of currentStep, the metrics are
evaluated just before starting the next step, and currentStep is already
updated to reflect that. Therefore, before starting the first step,
currentStep is 1, so the condition holds (unless you specify an offset).

Configuration blocks where logical conditions are used
======================================================

-  :ref:`JumpController's conditions to jump <sec-jumpController>`
-  :ref:`Dynamic changes in simulation parameters <sec-PELEDynamicChangesInSimulationParameters>`
-  :ref:`Exit conditions <sec-exitConditions>`

