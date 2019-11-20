% test of gain in speed by moving backslash operation out of for loop

N = 10000;
dta = randn(3,N);
ba = randn(3,1);
D = randn(3,3);

tic
a1 = D\(dta - repmat(ba, 1, N)); 
toc

tic
for k=1:N-1
a2 = D\(dta(:,k) - ba); 
end
toc
