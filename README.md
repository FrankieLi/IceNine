IceNine
=======

-----------
Authors:      Frankie Li (li31@llnl.gov)
-----------

-----------
Requirements:  
-----------    XDM++
	       Boost 1.43
               gcc > 4.4

-----------
Descriptions:
-----------
-This is a forward model reconstruction code that implements [1].  Detailed reconstruction
characteristics can be found in [1].  

-Note that IceNine requires XDM++, which is a library developed in Carnegie Mellon University
during S. F. Li's PhD Thesis.  This library hosted by Carnegie Mellon University on GitHub.

-Because of lagacy reasons, file formats are currently limited to *.mic for output and *.bin and
*.txt for input.  This will hopefully change in the forthcoming releases.



-----------
Content:
-----------
Both XDM++ and Boost must be acquired separately and are not part
of this distribution.

Src/           Source code for VPFFT
               -- Driver.cpp         --  main driver for the IceNine application
               -- SearchTrait.h      --  definition of an orientation search
               -- DiscreteAdaptive.h --  implementation of the adpative orientation
                                         search method. 

ConfigFiles/   Example configuration files and detector files.


-----------
References
-----------

[1]  S. F. Li, and R. M. Suter, Adaptive reconstruction method for three-dimensional
orientation imaging, Journal of Applied Crystallography, 2013.


