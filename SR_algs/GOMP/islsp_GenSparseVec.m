%�¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢�
% Generate K-sparse vector
%
% N			: original signal size.
% K			: sparsity level
%
%	Output parameters
% x_omp		: estimated signal
% iter_count: iteration count during estimating
%
% Written by Suhyuk (Seokbeop) Kwon
% Information System Lab., Korea Univ.
% http://isl.korea.ac.kr
%�¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢�
function [x x_pos] = islsp_GenSparseVec(N, K)

	KPos	= K;

	if N/2 < K
		KPos	= N-K;
	end

	randPos		= ceil(N*rand( KPos, 1 ));
	randPos		= union(randPos,randPos);
	leftPOsLen	= KPos-length(randPos);

	while leftPOsLen > 0
		tmpPos	= ceil(N*rand( leftPOsLen, 1 ));

		randPos	= union(tmpPos,randPos);
		leftPOsLen = KPos-length(randPos);
	end

	if KPos < K
		randPos	= setxor((1:N), randPos);
	end

	x				= zeros( N, 1 );
	x(randPos)		= randn( K, 1 );
	x_pos			= randPos;

end