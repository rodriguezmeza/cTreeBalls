
Takahashi Catalog Tutorial
==========================

This tutorial shows how to prepare a Takahashi realization for a **cBalls**
correlation-function run.  The download is approximately 3 GB, so use the
compact :doc:`tutorials/minimal-cli` workflow first when validating a build.



Let us test one of the `Takahashi <https://arxiv.org/pdf/1706.01472>`_ realizations. We download the realization using unix ``wget`` command in the terminal::

	wget http://cosmo.phys.hirosaki-u.ac.jp/takahasi/allsky_raytracing/sub1/nres12/allskymap_nres12r000.zs9.mag.dat

Then we will have in our working directory the file ``allskymap_nres12r000.zs9.mag.dat``, 3 Gb in size with ~200 million points distributed on the surface of a unit sphere.
