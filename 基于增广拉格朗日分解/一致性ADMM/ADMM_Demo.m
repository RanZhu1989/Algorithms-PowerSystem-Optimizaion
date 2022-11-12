% 普通ADMM测试算例
% min (x-1)^2 + (y-2)^2
% st 0<=x<=3
%    1<=y<=4
%    2x+3y==5

%% 直接用求解器的对比结果
clear all
ops = sdpsettings('solver','cplex');
x =sdpvar;
y = sdpvar;
%% 直接用求解器的对比结果
optimize([0<=x<=3, 1<=y<=4, 2*x+3*y==5], (x-1)^2+(y-2)^2, ops);
CX = value(x);
CY = value(y);

% x=0.5384, y = 1.3077

%% Yalmip 代码
k_max = 30; 
xk = zeros(k_max,1); % 迭代初值
xk(1)=0;
yk =zeros(k_max,1);
yk(1)=0;
lambda = 1; % 橙子
rho = 0.5;
k = 1;
gk=ones(k_max,1);
Cons_x = [0<=x<=3];
Cons_y = [1<=y<=4];
while k<=k_max && abs(gk(k))>=1e-4
    % 解y-fixed问题，得x
    Obj_x = (x-1)^2 + (yk(k)-2)^2 + lambda*(2*x+3*yk(k)-5) + 0.5*rho*(2*x+3*yk(k)-5)^2;
    optimize(Cons_x, Obj_x, ops);
    xk(k+1) = value(x);
    % 解x-fixed问题，得y
    Obj_y = (xk(k+1)-1)^2 + (y-2)^2 + lambda*(2*xk(k+1)+3*y-5) + 0.5*rho*(2*xk(k+1)+3*y-5)^2;
    optimize(Cons_y, Obj_y, ops);
    yk(k+1) = value(y);
    % 更新橙子
    lambda = lambda + rho*(2*xk(k+1)+3*yk(k+1)-5);
    gk(k+1) = 2*xk(k+1)+3*yk(k+1)-5;
    k = k+1;
end
figure;
yyaxis left
plot(1:k,gk(1:k),'-b','LineWidth',2);
hold on
yyaxis right
plot(1:k,xk(1:k),'--r','LineWidth',1);
hold on
plot(1:k,yk(1:k),'--g','LineWidth',1);
xlabel('迭代次数')
hh = legend('耦合约束不平衡度','x','y');%
hh.Orientation = 'horizontal';


