=============
scipion-pkpd
=============

This plugin allows to use cryoSPARC2 programs within the Scipion framework

You will need to use `Scipion 3.0.0 <https://scipion-em.github.io/docs/release-3.0.0/docs/scipion-modes/how-to-install.html>`_ to run the plugin protocols.


**Plugin Versions**
===================

**v1.0.2**
----------
January 20, 2021

**new**: Compatibility with Scipion v3.0.0

**new**: Migrated to python3


**Installing the plugin**
=========================

In order to install the plugin follow these instructions:

1. **Install the plugin**

.. code-block::

     scipion3 installp -p scipion-pkpd

or through the **plugin manager** by launching Scipion and following **Configuration** >> **Plugins**


**To install in development mode**

- Clone or download the plugin repository

.. code-block::

          git clone https://github.com/cossorzano/scipion-pkpd.git

- Install the plugin in developer mode.

.. code-block::

  scipion3 installp -p local/path/to/scipion-pkpd --devel


===============
Buildbot status
===============

Status devel version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/pkpd_devel.svg

Status production version:

.. image:: http://scipion-test.cnb.csic.es:9980/badges/pkpd_prod.svg

==========
SonarCloud
==========
To analize a cloud-based code quality click `here <https://sonarcloud.io/dashboard?id=cossorzano_scipion-pkpd>`_
