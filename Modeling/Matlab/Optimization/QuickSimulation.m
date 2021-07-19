global tx

FileName = '../Data/ParameterSetAll.csv';
EntryNumber = 1;
p = GetParameters(EntryNumber, FileName);

mu= p(1);
z = p(2);
g = p(3);
alpha = p(4);
c = p(5);
epsilon = p(6); 
b_R = p(7);
beta =p(8);
a = p(9);
b_T = p(10);
e_T = p(11);
e_R = p(12);
kA = p(13);
j = p(14);
kB = p(15);
n = p(16);
d = p(17);

%{
mu= p(1);
z = p(2);
g = p(3);
alpha = p(4);
c = p(5);
epsilon = p(6); 
b_R = p(7);
beta =p(8);
a = p(9);
b_T = p(10);
e_T = p(11);
e_R = p(12);
kA = p(13);
j = p(14);
kB = p(15);
n = p(16);
d = p(17);
%}

%Do not change this order
p0 = [alpha, a, kA, e_T, e_R, g, b_T, b_R, epsilon, mu, beta, c, kB, j, z, n, d];

tx = 0:432; %Maximum amount of time - 18 days

Genotype = [1, 2];

for i = Genotype
    PlottingResults(p0, i)
end

