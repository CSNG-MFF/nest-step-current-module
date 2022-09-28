NEST Extension Module with Neuron Models with Built-In Step Currents
====================================================================

.. attention::

   The code in this repository has been tested with NEST 3.3.

This repository contains neuron models with a built-in step current
generator. Behavior of these models is the same as if the underlying
neuron model is driven by a `step_current_generator` with the same
current profile, provided the current switching times are adjusted
according to the delay between `step_current_generator` and neuron.
