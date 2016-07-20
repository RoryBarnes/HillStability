This code calculates the relative proximity of a three-body system
to "Hill stability." A Hill stable system is one for which the    
the ordering of the planets remains constant, i.e. the most       
distent body may escape to infinity and the system would still be 
Hill stable. Hill stability may be calculated for any system, but 
this code is optimized for planetary systems. This code was used  
in R. Barnes & R. Greenberg, 2006, ApJ, 674, L163-L166 to         
demonstrate that for a system of a solar-type star and two        
Jupiter-ish mass planets that Hill stability approximates         
"Lagrange stability," which means no swapping && no ejections.    

The user inputs the central mass in solar units and the orbiters  
parameters are entered into the next two lines with the format:   
ArgumentPericenter MeanAnomaly. The units are Jupiter masses, AU, 
and degrees. The final line must state either "bodycentric" or    
"barycentric" to indicate the coordinate system of the orbital    
elements. There are no command line options.                      

This code will output 2 numbers: Exact and Approx. Both numbers   
represent relative proximity to the boundary, with unity on the   
boundary, values < 1 are Hill unstable, and > 1 are Hill stable.   
Exact is the value of beta/beta_{crit} from BG06, which is        
computed by calculating the energy and angular momentum including 
the primary. Approx is the value of delta/delta_{crit} from BG06, 
which is computed assuming the central mass dominates, see Gladman
(1993). As this option assumes the central body's center is very  
close to the system's center-of-mass, Approx uses the input       
elements.
