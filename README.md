Two-Step Watershed Segmentation of Epithelial Cells
===================================================

USAGE
-----

Place your images into HECD1, ACTIN, DAPI subfolders in the folder
input_images. Corresponding images should have the same name.

Modify parameters.R and double click on 'Run_segmentation.bat' to run the
segmentation. Resulting segmentation should then appear in output_dir. Then
double click on 'Run_analyse.bat' to obtain a file output_dir/data.txt which you
can open in Excel.


TROUBLESHOOTING
---------------

Error: Access denied

Solution: If you have a firewall, it might be blocking the execution of the
programs. Try disabling it. If the firewall is indeed the issue, a proper
install for R and its dependencies would be a safer long-term solution.

REFERENCES
----------

Keraudren, K., Spitaler, M., Braga, V. M., Rueckert, D., & Pizarro, L.:
<i>Segmenting Epithelial Cells in High-Throughput RNAi Screens</i>, MIAAB 2011.               
<a href="http://www.doc.ic.ac.uk/~kpk09/publications/MIAAB-2011.pdf">PDF</a>
<a href="http://www.doc.ic.ac.uk/~kpk09/cell_segmentation_tool.zip">code</a> 
<a href="http://www.doc.ic.ac.uk/~kpk09/publications/MIAAB-2011_slides.pdf">slides</a>

    

