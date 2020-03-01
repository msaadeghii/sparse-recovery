IFD-AMPM version 1
November 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This MATLAB collection includes the implementation of IFD-AMPM (Incoherent Frame Design via an Alternating 
Minimization Penalty Method) [1] for designing frames of very low mutual coherences.

IFD_AMPM.m			                        : main function

ProjectOntoL1Ball.m                     : computes projection onto the l1 norm ball [2]
UNTF_iterative_projection.m     : creates a unit-norm tight-frame using iterative projection
UNTF_frame_potential.m            : creates a unit-norm tight-frame based on [3]
mucoh.m                                         : computes the mutual coherence of a frame	
demo.m			                                 : a demo of how to use the main function


References:

   [1] M. Sadeghi and M. Babaie-Zadeh, “Incoherent Unit-norm Frame Design via an
         Alternating Minimization Penalty Method”,
         IEEE Signal Processing Letters, 2016.

   [2] J. Duchi, S. Shalev-Shwartz, Y. Singer, and T. Chandra, “Efficient
         projections onto the l1-ball for learning in high dimensions,” in
         International Conference on Machine Learning (ICML), 2008.

   [3] J. J. Benedetto and M. Fickus, “Finite normalized tight frames,” 
        Adv. Comput. Math., vol. 18, pp. 357–385, 2003.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Mostafa Sadeghi
  EE Department, Sharif University of Technology, Tehran, Iran
  Email: m.saadeghii@gmail.com
