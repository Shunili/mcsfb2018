function peakloc = zipeaks(y)
%Finds locations of local maxima
%yD=1 at maxima, yD=0 otherwise, end point maxima excluded
% http://www.mathworks.com/matlabcentral/fileexchange/24789-fast-peak-locator
    N=length(y)-2;
    yD=[0 (sign(sign(y(2:N+1)-y(3:N+2))-sign(y(1:N)-y(2:N+1))-.1)+1)/2 0];
%Locate indices of maxima
    Y=logical(yD);
    I=1:length(Y);
    peakloc=I(Y);
