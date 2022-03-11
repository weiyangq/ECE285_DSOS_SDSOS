clear all
yalmip('clear')
sdpvar x y
%%
mSize = 10;
E = randn(mSize,mSize); E = E + E';
F = randn(mSize,mSize); F = F + F';
%save('sdp_data.mat','E','F','mSize')
%%
sdpvar x y
matA = eye(mSize)+x*E+y*F;

%% psd constrain
cons_psd = [matA>=0];
var = [x y];
cost = x+y;
optimize(cons_psd,-cost);
cost_psd = value(cost);
figure
plot(cons_psd)

%% DD constrain
% define DD Q
%tstart = tic;
n_Vi = mSize^2;
% formulate the unique vi vector basis list
I1 = eye(mSize);
Im1 = -eye(mSize);
vi_temp = [];
for i = 1:mSize-1
    for j = 1:mSize-i
        vi_temp = [vi_temp;I1(i,:)+I1(i+j,:)];
    end
end
for i = 1:mSize-1
    for j = 1:mSize-i
        vi_temp = [vi_temp;I1(i,:)+Im1(i+j,:)];
    end
end
vi_temp = [vi_temp;I1];
%
eta = sdpvar(n_Vi,1);
mat_DD = zeros(mSize,mSize);
for i = 1:n_Vi
    mat_DD = mat_DD + eta(i)*vi_temp(i,:)'*vi_temp(i,:);
end
cons_M = [eta>=0];

%%
% define initial U
U_0 = eye(mSize);
% define cost
cost = x+y;
%%
%plot psd constrain
figure
plot(cons_psd,var,[0 0.1 0.3])
hold on
title('DD')
%title('DD iterative change of basis')
grid on
%%
U = U_0;
% save optimal cost each time
cost_temp = [];
% start iterative change of basis
for i = 1:5
% constrain the result
M_DD = U'*mat_DD*U;
cons_DD = [cons_M, M_DD == matA];
% plot constrain each time
plot(cons_DD,var,[0 (0.1+i*0.1) 0.3])
% optimize and get result
optimize(cons_DD,-cost);
cost_temp = [cost_temp, value(cost)];
% update x y optimal value and new U
x_1 = value(x);
y_1 = value(y);
matA_1 = eye(mSize)+x_1*E+y_1*F;
U = chol(matA_1);
end

%tend = toc(tstart)
%cost_psd - cost_temp(end)

%%
figure
plot(cost_temp,'-o')
yline(cost_psd)
grid on
title('DD optimize cost after each iteration')
xlabel('iteration')
ylabel('optimal cost')
legend('DD optimal','psd optimal','Location','southeast')
