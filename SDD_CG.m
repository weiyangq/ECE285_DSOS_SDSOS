clear all
yalmip('clear')
sdpvar x y
%%
load('sdp_data.mat')
% mSize = 10;
% E = randn(mSize,mSize); E = E + E';
% F = randn(mSize,mSize); F = F + F';
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
In = eye(mSize);
iter_n = 0;
for i = 1:mSize-1
    for j = 1:mSize-i
        iter_n = iter_n+1;
        V_temp = [In(i,:);In(i+j,:)];
        E_ij(:,:,iter_n) = V_temp';
    end
end
eta = sdpvar(n_2by2,3);

mat_SDD = zeros(mSize,mSize);
for i = 1:n_2by2
    Eta = [eta(i,1) eta(i,2);eta(i,2) eta(i,3)];
    mat_SDD = mat_SDD + E_ij(:,:,i)*Eta*E_ij(:,:,i)';
end

cons_M = [];
for i =1:n_2by2
Eta = [eta(i,1) eta(i,2);eta(i,2) eta(i,3)];
cons_M = [cons_M, Eta>=0];
% cons_M = [cons_M, eta(i)>=0];
% cons_M = [cons_M, norm([2*eta(i)^2; eta(i)-eta(i)^3])<=(eta(i)+eta(i)^3)];
end


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
for i = 1:n_2by2
VXV_i = E_ij(:,:,i)'*X*E_ij(:,:,i);
cons_dual_B0 = [cons_dual_B0, VXV_i>=0];
% cons_dual_B0 = [cons_dual_B0, VXV_i(1,1)>=0];
% cons_dual_B0 = [cons_dual_B0, norm([2*VXV_i(1,2)^2; ...
%     VXV_i(1,1)-VXV_i(2,2)^3])<=(VXV_i(1,1)+VXV_i(2,2)^3)];
end
cons_dual_0 = [sum(dot(A_x,X)) == b_x, sum(dot(A_y,X)) == b_y];
cons_dual_0 = [cons_dual_0, cons_dual_B0];
cons_dual = cons_dual_0;

%%
cost = sum(dot(C,X));
for i = 1:8
optimize(cons_dual,cost)
X_opt = value(X);
[V_eig,D_eig] = eig(X_opt);
D_eig;
b_j = V_eig(:,1:2);
b_j_sav(:,:,i) = b_j;
B_j = b_j*b_j';
VXV_new = b_j'*X*b_j;
cons_dual = [cons_dual, VXV_new>=0];
% cons_dual = [cons_dual, VXV_new(1,1)>=0];
% cons_dual = [cons_dual, norm([2*VXV_new(1,2)^2; ...
%     VXV_new(1,1)-VXV_new(2,2)^3])<=(VXV_new(1,1)+VXV_new(2,2)^3)];
end

%%
% plot psd constrain
% figure
% plot(cons_psd,var,[0 0.1 0.3])
% hold on
% title('SDD iterative change of basis')
% grid on

%%
% save optimal cost each time
cost = x+y;
cost_temp = [];
alpha = sdpvar(7,3);
cons_Alpha = [];
for i = 1:7
Alpha = [alpha(i,1) alpha(i,2);alpha(i,2) alpha(i,3)];
cons_Alpha = [cons_Alpha, Alpha>=0];
end

M_SDD = mat_SDD;
for i = 1:5
cons_SDD = [cons_M, cons_Alpha];
cons_SDD = [cons_SDD, M_SDD == matA];
% optimize and get result
%option = sdpsettings('solver','mosek');
%optimize(cons_DD,-cost,option);
%plot(cons_SDD,var,[0.5 (0.9-i*0.1) 0.25])
%%
% optimize(cons_SDD,-cost);
% cost_temp = [cost_temp, value(cost)];
%%
Alpha = [alpha(i,1) alpha(i,2);alpha(i,2) alpha(i,3)];
VAV_i = b_j_sav(:,:,i)*Alpha*b_j_sav(:,:,i)';
M_SDD = M_SDD + VAV_i;
end
optimize(cons_SDD,-cost);
cost_temp = [cost_temp, value(cost)];

tend = toc(tstart)
cost_psd - cost_temp(end)

%%
% plot psd constrain
figure
plot(cons_psd,var,[0.5 0.1 0.3])
hold on
title('SDD column generation')
grid on
M_SDD_back = M_SDD;

for i = 1:7
Alpha = [alpha(8-i,1) alpha(8-i,2);alpha(8-i,2) alpha(8-i,3)];
VAV_i = b_j_sav(:,:,8-i)*Alpha*b_j_sav(:,:,8-i)';
M_SDD_back = M_SDD_back - VAV_i;
cons_SDD = [cons_M, cons_Alpha];
cons_SDD = [cons_SDD, M_SDD_back == matA];
plot(cons_SDD,var,[0.5 (0.9-i*0.1) 0.25])
end



%%
figure
plot(cost_temp,'-o')
yline(cost_psd)
grid on
title('SDD optimize cost after each iteration')
xlabel('iteration')
ylabel('optimal cost')
legend('DD optimal','psd optimal','Location','southeast')

