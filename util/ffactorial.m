function out = ffactorial(n)
%FFACTORIAL FFactorial (double factorial) function.

N = n(:);
if any(fix(N) ~= N) || any(N < 0) || ~isa(N,'double') || ~isreal(N)
  error('MATLAB:factorial:NNegativeInt', ...
        'N must be a matrix of non-negative integers.')
end

if n==0 || n==1
    out=1;
else
    out=n*ffactorial(n-2);
end

%n(N>170) = 171;
%m = max([1; n(:)]);
%N = [1 1 cumprod(2:m)];
%n(:) = N(n+1);


