function [dst,DEN,T] = solveLinearEquation( NUM, wy, lambda ,method)

if (~exist('lambda','var'))	% alpha -- parameter for shape
    lambda = 1;
end
if (~exist('method','var'))	% alpha -- parameter for shape
    method = 'direct';
end

%%
[ r, c] = size(NUM);
hw = c;
uvy = lambda * wy( : );
uy = padarray(uvy, 1, 'pre'); uy = uy(1:end-1);
D = uvy+uy+1;
T = spdiags([-uvy],[-1],hw,hw);
DEN = T + T' + spdiags(D, 0, hw, hw);

switch method
    case 'pcg'
        L = ichol(DEN,struct('michol','on'));
        [dst,~] = pcg(DEN, NUM(:), 0.01, 40, L, L');
    case 'minres'
        [dst,~] = minres(DEN,NUM(:), 0.01, 40);
    case 'bicg'
        [L,U] = ilu(DEN,struct('type','ilutp','droptol',0.01));
        [dst,~] = bicg(DEN,NUM(:), 0.01, 40, L, U);
    case 'direct'
        dst = DEN\NUM(:); %#ok<RHSFN>
end

end