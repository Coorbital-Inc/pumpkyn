function y = chebyspace(a,b,n)
%% Purpose:
%
% chebyspace:  create a chebyshev node spaced vector of n points between
% a and b.  
%
% Source References:
%
%  https://en.wikipedia.org/wiki/Chebyshev_nodes
%
%% Inputs:
%
%  a                    double              Lower bound of interval
%
%  b                    double              Upper bound of interval
%
%  n                    integer             Number of points used for
%                                           chebyshev nodes
%
%% Outputs:
%
%  y                    [1 x n]             n equally spaced points between
%                                           a and b.
%  
%
%% Revision History:
%  Darin C. Koblick                                         (c) 11-26-2024
%% ------------------------------- Begin Code Sequence --------------------
if nargin == 0
    a = -10;
    b = +10;
    n = 10;
    y = chebyspace(a,b,n);
    figure('color',[1 1 1]);
    plot(1:n,y); hold on;
    plot(1:n,y,'.r','markersize',20);
    grid on;
    legend('','Cheyshev Nodes','location','NW');
    xlabel('Node Number');
    ylabel('Value');
    return;
end
%Scale inputs:
k = n-1:-1:0;
y = (a+b)./2 + ((b-a)./2)*cos(pi.*k./(n-1));
end