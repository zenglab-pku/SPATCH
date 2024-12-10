Key functions for the data processing used in SPATCH

1_load_data.py: load raw expression or transcripts data of different platforms using scanpy at 8-μm bin or cell resolution and save in h5ad format

2_8um_bin.py: process transcriptomic data from Xenium or CosMx into 8 µm resolution and remove signals outside tissue region.

3_registeration.py: perform the registration for paired images and apply the transformation to other channels of CODEX or spatial coordinates of ST data

4_diffusion.py: calculate the relative diffusion effects and minimum distance between spots inside and outside tissue for sST platforms.

5_correlation_with_codex.py: calculate the correlation between ST data and CODEX data over the spatial grids.

6_scrna.r: perform preprocessing for scRNA-seq data.

7_cluster.py: perform clustering for ST data.

8_cell_shape.py: calculate the statistics for describing the cell shape.

9_st_annotation.py: transfer the annotations from scRNA-seq to ST data.

9_st_annotation_consistency.r: assess the consistency across different annotation tools.

10_spatial_cluster.py: perform spatial clustering for ST and CODEX data, and assess their consistency.
