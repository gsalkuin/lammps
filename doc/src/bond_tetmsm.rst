.. index:: bond_style tetmsm

bond_style tetmsm command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style tetmsm

Examples
""""""""

.. code-block:: LAMMPS

   bond_style tetmsm
   bond_coeff 1 2.5e4

   bond_style tetmsm
   bond_coeff * 1.8

Description
"""""""""""

The *tetmsm* bond style computes a harmonic spring interaction designed
for tetrahedral mass-spring models (MSM) (:ref:`Lloyd2007 <Lloyd2007>`).
The energy and force are given by:

.. math::

   U      &= \frac{1}{2}\,\kappa_E\,r_0 \left(1 - \frac{r}{r_0}\right)^2 \\
   F      &= -\kappa_E\,\frac{r - r_0}{r_0}

where :math:`r` is the current bond length, :math:`r_0` is
the reference length, and :math:`\kappa_E` is the spring
stiffness (energy/length). The conventional spring constant is
:math:`k_E = \kappa_E / r_0` (force/length).

This command should be used with :doc:`improper_style tetmsm <improper_tetmsm>`
to model an isotropic linear-elastic solid. The coefficient :math:`\kappa_E`
can be derived from the target shear modulus :math:`G` using
(:ref:`Lloyd2007 <Lloyd2007>`, :ref:`Golec2020 <Golec2020>`):

.. math::

   k_E \approx \frac{\sqrt{2}}{5} G \sum_\Delta^\text{adj. tets} 
   \left(\frac{12 V_\Delta}{\sqrt{2}} \right)^{1/3}

where the summand is an effective edge length of an irregular tetrahedron
:math:`\Delta` with volume :math:`V_\Delta` (:ref:`Lloyd2007 <Lloyd2007>`).

The following coefficient must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`\kappa_E` (energy/length)

The equilibrium bond length :math:`r_0` is **not** specified as a coefficient;
instead, it is computed from the initial atom positions on the first timestep
and cached internally (keyed by global atom tags) similar to
:doc:`bond_style bpm/spring <bond_bpm_spring>`.

.. note::

   The equilibrium distance :math:`r_0` for each bond is determined from
   the atom positions at the beginning of the first
   :doc:`run <run>` or :doc:`minimize <minimize>` command after the bond
   style has been defined.  Subsequent ``run`` commands in the same input
   script reuse the cached values.  To reset :math:`r_0`, re-issue the
   ``bond_style tetmsm`` command before the next ``run``.

----------

Restart info
""""""""""""

This bond style supports the :doc:`write_restart <write_restart>` and
:doc:`read_restart <read_restart>` commands. The :math:`\kappa_E` 
coefficient for each bond type and the per-bond reference lengths
:math:`r_0` are stored.

Restrictions
""""""""""""

This bond style maintains internal data to determine the original bond
lengths :math:`r_0`.  This information will be written to binary restart
files but not to data files.  Thus, continuing a simulation from a
deformed state is only possible with :doc:`read_restart <read_restart>`.
When using :doc:`read_data <read_data>`, the reference lengths will be
re-initialized from the current geometry.

This bond style requires that atoms have tags (``atom_modify id yes``,
which is the default).

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`,
:doc:`bond_style harmonic/restrain <bond_harmonic_restrain>`,
:doc:`improper_style tetmsm <improper_tetmsm>`

Default
"""""""

none

----------

.. _Lloyd2007:

**(Lloyd2007)** Lloyd, Szekely, Harders, IEEE Trans Vis Comput Graph, 13(5),
1081-1094 (2007).

.. _Golec2020:

**(Golec2020)** Golec, Palierne, Zara, Nicolle, Damiand, Vis Comput, 36,
809-825 (2020).
