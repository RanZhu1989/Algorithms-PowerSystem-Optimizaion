% 一致性凸ADMM演示算例

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
y = sdpvar;
%% 直接用求解器的对比结果
optimize([0<=x<=3, 1<=y<=4, 2*x+3*y==5], (x-1)^2+(y-2)^2, ops);
CX = value(x);
CY = value(y);

% x=0.5384, y = 1.3077


%% Consensus-based ADMM
% 以x作为一致变量
% 问题1： min (x1-1)^2 , st. 0<=x1<=3
% 问题2:   min (y-2)^2, st. 1<=y<=4,   2*x2+3y==5
k_max = 200; 
x1k = zeros(k_max,1); % 保存变量的结果
x2k = zeros(k_max,1); % 保存变量的结果
yk =zeros(k_max,1);
xaux_k = zeros(k_max,1); % 一致变量的更新
xaux_k(1) = 3;
lambda1 = 0; % 两个橙子
lambda2 = 0;
rho = 0.5*ones(k_max,1);

% 自适应rho更新  Boyd. 2010
% 可以看出，采用自适应更新后只要8步就收敛，否则需要28步收敛
tau = 2;
mu = 10;


PR = ones(k_max,1); % 原始残差
DR = ones(k_max,1); % 对偶残差
Cons_x = [0<=x1<=3];
Cons_y = [1<=y<=4, 2*x2+3*y==5];

k = 1;
while k<=k_max && (PR(k)>=1e-4 || DR(k)>=1e-4)
    % 解只关于x的问题1
    Obj_x = (x1-1)^2 + lambda1*(x1-xaux_k(k)) + 0.5*rho(k)*(x1-xaux_k(k))^2;
    optimize(Cons_x, Obj_x, ops);
    x1k(k+1) = value(x1);
    % 解只关于y的问题2
    Obj_y = (y-2)^2 + lambda2*(x2-xaux_k(k)) + 0.5*rho(k)*(x2-xaux_k(k))^2;
    optimize(Cons_y, Obj_y, ops);
    x2k(k+1) = value(x2);
    yk(k+1) = value(y);
    % 更新一致量
    xaux_k(k+1) = 0.5*(x1k(k+1)+x2k(k+1));
    % 更新橙子
    lambda1 = lambda1 + rho(k)*(x1k(k+1)-xaux_k(k+1));
    lambda2 = lambda2 + rho(k)*(x2k(k+1)-xaux_k(k+1));
    PR(k+1) = max( abs(x1k(k+1)-xaux_k(k+1)), abs(x2k(k+1)-xaux_k(k+1)));
    DR(k+1) = rho(k)*(xaux_k(k+1)-xaux_k(k)); % 这里只有一个一致量，故不用max了
    % 更新rho 可选
    rho(k+1) = rho(k);
    if PR(k+1)>mu*DR(k+1)
        rho(k+1) = tau*rho(k);
    end
    if DR(k+1)>mu*PR(k+1)
        rho(k+1) = rho(k)/tau;
    end
    k = k+1;
end

figure;
yyaxis left
plot(1:k,PR(1:k),'-o','LineWidth',0.5);
hold on
plot(1:k,DR(1:k),'-*','LineWidth',0.5);
hold on
yyaxis right
plot(1:k,x1k(1:k),'--r','LineWidth',1);
hold on
plot(1:k,x2k(1:k),'--g','LineWidth',1);
xlabel('迭代次数')
hh = legend('原始残差','对偶残差','x1','x2');%
hh.Orientation = 'horizontal';
DX = value(x1k(k-1));
DY = value(yk(k-1));








