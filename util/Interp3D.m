function Vq = Interp3D(X, Y, Z, V, Xq, Yq, Zq, method, ExtrapolationMethod)

arguments
    X % X values of original data set
    Y % Y values of original data set
    Z % Z values of original data set
    V % Value of function at points in original data set
    Xq % X query points
    Yq % Y query points
    Zq % Z query points
    method = 'linear';
    ExtrapolationMethod = 'none';
end

if all(size(X)==size(Y)) && all(size(Y) == size(Z))
    %
else
    error('Dimensions of input domain are not the same')
end
if all(size(Xq)==size(Yq)) && all(size(Yq) == size(Zq))
    %
else
    error('Dimensions of output domain are not the same')
end

F = griddedInterpolant(X, Y, Z, V, method, ExtrapolationMethod);

Vq = F(Xq, Yq, Zq);

% if any(isnan(Vq), 'all')
%     warning('A point outside of the original domain was queried. The corresponding value is NaN')
% end
Vq( isnan(Vq) ) = 0;

end