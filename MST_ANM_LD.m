function [freq, pow] = MST_ANM_LD(Y, beam_w, K)

% function [freq, pow] = MST_ANM_LD(Y, beam_w, K)
%
% MST_GL_ANM implements the MS Atomic Norm Minimization for
% DOA estimation
%
% Input
% Y: observation
%beam_w; EADF matrix
%
% Output:
% freq: frequency estimate
% pow: power estimate
% sigma: noise variance estimate
%
% Written by Pan Jie, 2019   E-mail:panjie@yzu.edu.cn



[M, N] = size(Y);


if N > 1
    Rhat = Y * Y' / N;
end

cvx_quiet true
cvx_precision default


[freq, pow] = ANM_MS(Rhat, beam_w, N,K);




end



function [freq, pow] = ANM_MS(Rhat,beam_w, N,K)
% SPA when Rhat is nonsingular
    M = size(beam_w, 2);
    B = size(beam_w,1);
    Rhat =Rhat/sum(real(diag(Rhat)));
    delta = svd(Rhat);
    delta = sqrt(N/3600)*B*delta(end);
    eta = (delta);
    degrad = pi/180;
    
    cvx_solver sdpt3
    cvx_begin sdp
    variable x(B,B) hermitian,
    variable u(M) complex,
    variable R(B,B) hermitian
    
    [x R; R beam_w*toeplitz(u)*beam_w'] >= 0,
    norm(R-Rhat,'fro')<=eta
    minimize trace(x) + real(trace(beam_w*toeplitz(u)*beam_w'));
    cvx_end
    sig = 0;
    R_cov = beam_w*toeplitz(u)*beam_w';
    [U SS V] = svd(R_cov);
    En = U(:,K+1:end);
    GG=beam_w'*En*En'*beam_w;
    MM=M+1;
    a = zeros(2*MM-1,1);
    for j=-(MM-1):(MM+1)
        a(j+MM) = sum( diag(GG,j) );
    end
    
    ra=roots([a]);
    rb=ra(abs(ra)<1);
    [~,I]=sort(abs(abs(rb)-1));
    w=angle(rb(I(1:K)));
    freq=sort((w/degrad-1));
    pow = 0;
end





