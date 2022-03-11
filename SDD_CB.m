clear all
yalmip('clear')
sdpvar x y
%%
%load('sdp_data.mat')
mSize = 10;
E = randn(mSize,mSize); E = E + E';
F = randn(mSize,mSize); F = F + F';
%%
matA = eye(mSize)+x*E+y*F;

%% psd constrain
cons_psd = [matA>=0];
var = [x y];
cost = x+y;
optimize(cons_psd,-cost);
cost_psd = value(cost);
figure
plot(cons_psd)

%% SDD constrain
% define SDD Q
tstart = tic;
n_2by2 = nchoosek(mSize,2);
M = sdpvar(2,2,n_2by2)
cons_M = [];
for i = 1:n_2by2;
cons_M = [cons_M, M(1,1,i)>=0];
cons_M = [cons_M, norm([2*M(1,2,i); M(1,1,i)-M(2,2,i)])<=(M(1,1,i)+M(2,2,i))];
end

In = eye(mSize);
Eij_temp = [];
mat_SDD = zeros(mSize,mSize);
iter_n = 0;
for i = 1:mSize-1
    for j = 1:mSize-i
        iter_n = iter_n+1;
        E_ij = [In(i,:);In(i+j,:)];
        mat_SDD = mat_SDD+E_ij'*M(:,:,iter_n)*E_ij;
    end
end

%%
% define initial U
U_0 = eye(mSize);
% define cost
cost = x+y;
% plot psd constrain
%%
figure
plot(cons_psd,var,[0.5 0.2 0.25])
hold on
title('SDD')
%title('SDD iterative change of basis')
grid on
%%
U = U_0;
% save optimal cost each time
cost_temp = [];
% start iterative change of basis
for i = 1:3
% constrain the result
M_SDD = U'*mat_SDD*U;
cons_SDD = [cons_M, M_SDD == matA];
% plot constrain each time
plot(cons_SDD,var,[0.5 (0.2+i*0.1) 0.25])
% optimize and get result
optimize(cons_SDD,-cost);
cost_temp = [cost_temp, value(cost)];
% update x y optimal value and new U
x_1 = value(x);
y_1 = value(y);
matA_1 = eye(mSize)+x_1*E+y_1*F;
U = chol(matA_1);
end

tend = toc(tstart)
cost_psd - cost_temp(end)
%%
figure
plot(cost_temp,'-o')
yline(cost_psd)
grid on
title('SDD optimize cost after each iteration')
xlabel('iteration')
ylabel('optimal cost')
legend('SDD optimal','psd optimal','Location','southeast')

