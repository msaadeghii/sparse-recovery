%�¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢�
% Estimate the sparse signal x using generalized OMP
%
% y		: observation
% Phi	: sensing matrix
% K		: sparsity
% S		: selection length
%
%	Output parameters
% x_omp		: estimated signal
% iter_count: iteration count during estimating
%
% Written by Suhyuk (Seokbeop) Kwon
% Information System Lab., Korea Univ.
% http://isl.korea.ac.kr
%�¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢�
function [x_ommp iter_count] = islsp_EstgOMP(y, Phi, K, S, zero_threshold)

% Check the parameters
	if nargin < 3
		error('islsp_EstgOMP : Input arguments y ,Phi and K must be specified.');
	end
		
	if nargin < 4
		S = max(K/4, 1);
	end
	
	if nargin < 5
		zero_threshold  = 1e-6;
	end

% Initialize the variables
	[nRows nCols]	= size(Phi);
	x_ommp			= zeros(size(Phi,2), 1);
	residual_prev	= y;
	supp			= [];
	iter_count		= 0;

	while (norm(residual_prev) > zero_threshold && iter_count < K)
		iter_count	= iter_count+1;
		[supp_mag supp_idx]	= sort(abs(Phi'*residual_prev), 'descend');
		supp_n				= union(supp, supp_idx(1:S));

		if (length(supp_n) ~= length(supp)) && (length(supp_n) < nRows )
			x_hat			= Phi(:,supp_n)\y;
			residual_prev	= y - Phi(:,supp_n)*x_hat;

			supp	= supp_n;
		else
			break;
		end
	end

	x_ommp(supp)	= Phi(:,supp)\y;
	
	if nargout < 2
		clear('iter_count');
	end
end
