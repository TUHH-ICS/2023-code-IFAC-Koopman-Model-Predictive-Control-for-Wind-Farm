lx = 3;
ly = 4;
nx = 3;
ny = 1;
kk = (1:7:34);
nz = length(kk);

n1 = nan(lx,ly,kk);
aVec = (1:12)';


n = aVec .*kk ;

 n2 = reshape(n,lx,ly,5);

testVec = squeeze(n1(ny,nx,:));
disp(testVec)