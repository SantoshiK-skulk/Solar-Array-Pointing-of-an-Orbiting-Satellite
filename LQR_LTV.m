function [K, P] = lqr_LTV(AFun,BFun,Q,R,tSpan)
    nSteps = length(tSpan);

    P{nSteps} = 1/2 * Q;
    K{nSteps} = zeros(length(R),length(Q));
    
    for i = nSteps-1:-1:1
        A_ = AFun(i);
        B_ = BFun(i);
        P_ = P{i+1};
        
        K{i} = ( 1/2 * R + B_' * P_ * B_ )^(-1) * B_' * P_ * A_;
        P{i} = A_' * P_ * ( A_ - B_ * K{ i } ) + Q * 1/2;
    end
end