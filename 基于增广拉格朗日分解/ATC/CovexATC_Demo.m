% 两个主体的ATC演示算例

%  ------- 待求解问题 --------
% min (x-1)^2 + (y-2)^2
% st 0<=x<=3
%    1<=y<=4
%    2x+3y==5
clear all
ops = sdpsettings('solver','cplex');
x =sdpvar;
x1 = sdpvar;
x2 = sdpvar;
xc = sdpvar;
y = sdpvar;
%% 直接用求解器的对比结果
optimize([0<=x<=3, 1<=y<=4, 2*x+3*y==5], (x-1)^2+(y-2)^2, ops);
CX = value(x);
CY = value(y);

% x=0.5384, y = 1.3077
%% 上层视为协调层的ATC
% 由于只有两个主体，所以上层视为协调层
% 以x作为耦合变量
% 上层问题： min (x1-1)^2 , st. 0<=x1<=3
% 下层问题:   min (y-2)^2, st. 1<=y<=4,   2*x2+3y==5
k_max = 200; 
x1k = zeros(k_max,1); % 保存变量的结果
x2k = zeros(k_max,1); % 保存变量的结果
yk = zeros(k_max,1); % 保存变量的结果
x1k(1) = 0; % 初始化
x2k(1) = 0;

lambda = 0.1; % 乘子初始值
alpha = 1.1; % 学习率
beta = 0.1; 
PR = ones(k_max,1); % 原始残差
Cons_x = [0<=x1<=3];
Cons_y = [1<=y<=4, 2*x2+3*y==5];

k = 1;

while k<=k_max && PR(k)>=1e-6
    % 解下层问题 用上层的x1*作为耦合量, 因此第一项为上层问题返回的x1*
    Obj_y = (y-2)^2 + lambda*(x1k(k)-x2) + (beta^2)*(x2-x1k(k))^2;
    optimize(Cons_y, Obj_y, ops);
    x2k(k+1) = value(x2);
    yk(k+1) = value(y);
    % 解上层问题 得x1
    Obj_x = (x1-1)^2 + lambda*(x1-x2k(k+1)) + (beta^2)*(x1-x2k(k+1))^2;
    optimize(Cons_x, Obj_x, ops);
    x1k(k+1) = value(x1);
    % 更新橙子
    lambda = lambda + 2*(beta^2)*(x1k(k+1)-x2k(k+1));
    beta = alpha*beta;
    PR(k+1) = abs(x1k(k+1)-x2k(k+1));
    k = k+1;
end

figure;
yyaxis left
plot(2:k,PR(2:k),'-o','LineWidth',0.5);
hold on
yyaxis right
plot(2:k,x1k(2:k),'--r','LineWidth',1);
hold on
plot(2:k,x2k(2:k),'--g','LineWidth',1);
xlabel('迭代次数')
hold on
plot(2:k,yk(2:k),'--b','LineWidth',1);
xlabel('迭代次数')
hh = legend('原始残差','x1','x2', 'y');%
hh.Orientation = 'horizontal';
DX1 = value(x1k(k));
DX2 = value(x2k(k));
DY1 = value(yk(k));


%% 含协调层的ATC
% 
% 以x作为耦合变量
% 上层问题： min (x1-1)^2 , st. 0<=x1<=3
% 下层问题:   min (y-2)^2, st. 1<=y<=4,   2*x2+3y==5
k_max = 200; 
x1k = zeros(k_max,1); % 保存变量的结果
x2k = zeros(k_max,1); % 保存变量的结果
xck = zeros(k_max,1);
yk = zeros(k_max,1); % 保存变量的结果


lambda1 = 0.1; % 乘子初始值
lambda2 = 0.1; 
alpha = 1.1; % 学习率
beta1 = 0.1; 
beta2 = 0.1; 
PR = ones(k_max,1); % 原始残差
Cons_x = [0<=x1<=3];
Cons_y = [1<=y<=4, 2*x2+3*y==5];

k = 1;

while k<=k_max && PR(k)>=1e-6
    % 用协调层xc  解上层问题 得x1
    Obj_x = (x1-1)^2 + lambda1*(xck(k)-x1) + (beta1^2)*(x1-xck(k))^2;
    optimize(Cons_x, Obj_x, ops);
    x1k(k+1) = value(x1);
    % 用协调层xc  解下层问题 得x2*
    Obj_y = (y-2)^2 + lambda2*(xck(k)-x2) + (beta2^2)*(x2-xck(k))^2;
    optimize(Cons_y, Obj_y, ops);
    x2k(k+1) = value(x2);
    yk(k+1) = value(y);
    % 解协调层, 得xc*
    Obj_c = lambda1*(xc-x1k(k+1)) + lambda2*(xc-x2k(k+1)) + (beta1^2)*(xc-x1k(k+1))^2 + (beta2^2)*(xc-x2k(k+1))^2;
    optimize([], Obj_c, ops);
    xck(k+1) = value(xc);
    % 更新橙子
    lambda1 = lambda1 + 2*(beta1^2)*(xck(k+1)-x1k(k+1));
    lambda2 = lambda2 + 2*(beta2^2)*(xck(k+1)-x2k(k+1));
    beta1 = alpha*beta1;
    beta2 = alpha*beta2;
    PR(k+1) = abs(x1k(k+1)-x2k(k+1));
    k = k+1;
end

figure;
yyaxis left
plot(2:k,PR(2:k),'-o','LineWidth',0.5);
hold on
yyaxis right
plot(2:k,x1k(2:k),'--r','LineWidth',1);
hold on
plot(2:k,x2k(2:k),'--g','LineWidth',1);
xlabel('迭代次数')
hold on
plot(2:k,yk(2:k),'--b','LineWidth',1);
xlabel('迭代次数')
hh = legend('原始残差','x1','x2', 'y');%
hh.Orientation = 'horizontal';

DX1 = value(x1k(k));
DX2 = value(x2k(k));
DY1 = value(yk(k));



