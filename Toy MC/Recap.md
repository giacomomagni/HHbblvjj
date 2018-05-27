Toy MC simulation to estime which variable gives the best relative uncertainty on the number of Signal evts predicted.

The pdf are fitted separately using a polN. 
Random generation of N exepected events for background and signal
Fit of B+S distribution and fit with f(x,N_s, N_b) = N_s*pdf_s(x)+N_b*pdf_b(x) (N_s and N_b are now parameters to be estimated)
Evaluation of relative uncertainty on N_s at different values of luminosity. 

Calulation of epected N_s and N_b 
  N = luminosity * cross section * branching ratios * cuts selection 
  	
	1. four different values of integrated luminosity:
    For ex:
			-- 100 fb^-1  actual luminosity 
			-- 150 fb^-1  end of 2018
			-- 300 fb^-1  end of 2nd LHC run
 			-- 3000 fb^-1 high-luminosity LHC
	2.cross section  at E = 13 TeV:
			t tbar		831.76*10e^3 fb
 			H H 			33.5 fb
	3.Branching ratios:
			The Branching ratio for each decay in the signal is: 
 				H >> b bbar							58.24+-0.38	%	(B1)		 		
				H >> W W							21.37+-0.21	%	(B2)
 				W >> q1 q2							67.60+-0.27	%	(B3)
 				W >> l nu	(only e- e+ mu- mu+)	21.32+-0.19	%	(B4)
			The total branching ratio for the signal event is given by the product 2*(B1*B2)*2*(B3*B4):
				H H >> b bbar q1 q2 l+ nu OR b bbar q3 q4 l- nu				7.17+-0.11	%

			The Branching ratio for each decay in the backgroung is: 
 				t >> b W								91+-4	%	(B1)			
 				W >> q1 q2							67.60+-0.27	%	(B2)
 				W >> l nu	(only e- e+ mu- mu+)	21.32+-0.19	%	(B3)
			The total branching ratio for the background is given by the product B1^2*2*(B2*B3):
 				t tbar >> b bbar q1 q2 l+ nu OR b bbar q3 q4 l- nu		23.86+-1.7	%	
     
     4.Cut selection factors are computed when cut on masses are applied ( cut_s around 100%, cut_b around 0.6%)  
