clear all
yalmip('clear')
sdpvar x y
%%
%load('sdp_data.mat')
mSize = 10;
E = randn(mSize,mSize); E = E + E';
F = randn(mSize,mSize); F = F + F';
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
% start iterative column generation
C = eye(mSize);
X = sdpvar(mSize,mSize,1);
cost = sum(dot(C,X));
A_x = -E;
A_y = -F;
b_x = 1;
b_y = 1;
cons_dual_B0 = [];
for i = 1:n_Vi
    Vi = vi_temp(i,:)'*vi_temp(i,:);
    cons_dual_B0 = [cons_dual_B0, sum(dot(Vi,X))>=0];
end
cons_dual_0 = [sum(dot(A_x,X)) == b_x, sum(dot(A_y,X)) == b_y];
cons_dual_0 = [cons_dual_0, cons_dual_B0];
cons_dual = cons_dual_0;

b_j_sav = [];
for i = 1:10
optimize(cons_dual,cost)
X_opt = value(X);
[V_eig,D_eig] = eig(X_opt);
b_j = V_eig(:,1);
b_j_sav = [b_j_sav, b_j];
B_j = b_j*b_j';
cons_dual = [cons_dual, sum(dot(B_j,X))>=0];
end

%%

%%
% plot psd constrain
% figure
% plot(cons_psd,var,[0 0.1 0.3])
% hold on
% title('SDD iterative change of basis')
% grid on

%%
% save optimal cost each time
cost_temp = [];
cost = x+y;
alpha = sdpvar(11,1);
M_DD = mat_DD;
for i = 1:7
cons_DD = [cons_M, alpha>=0];
cons_DD = [cons_DD, M_DD == matA];
% optimize and get result
optimize(cons_DD,-cost);
cost_temp = [cost_temp, value(cost)];
B_j = b_j_sav(:,i)*b_j_sav(:,i)';
M_DD = M_DD + alpha(i)*B_j;
end
% optimize(cons_DD,-cost);
% cost_temp = [cost_temp, value(cost)];
%tend = toc(tstart)
%cost_psd - cost_temp(end)

%%
% alpha = sdpvar(21,1);
% M_DD = mat_DD;
% for i = 1:10
% B_j = b_j_sav(:,i)*b_j_sav(:,i)';
% M_DD = M_DD + alpha(i)*B_j;
% end
% cons_DD = [cons_M, alpha>=0];
% cons_DD = [cons_DD, M_DD == matA];
% optimize(cons_DD,-cost);
% value(cost)


%%
% plot psd constrain
figure
plot(cons_psd,var,[0 0.1 0.3])
hold on
title('DD column generation')
grid on
M_DD_back = M_DD;
for i = 1:7
B_j = b_j_sav(:,8-i)*b_j_sav(:,8-i)';
M_DD_back = M_DD_back - alpha(8-i)*B_j;
cons_DD = [cons_M, alpha>=0];
cons_DD = [cons_DD, M_DD_back == matA];
plot(cons_DD,var,[0 (0.9-i*0.1) 0.3])
end

%%
figure
plot(cost_temp,'-o')
yline(cost_psd)
grid on
title('DD optimize cost after each iteration')
xlabel('iteration')
ylabel('optimal cost')
legend('DD optimal','psd optimal','Location','southeast')


