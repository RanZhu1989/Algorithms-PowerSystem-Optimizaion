% PH 算法测试的简单算例
% Jean-Paul Watson, Progressive Hedging Innovations for a Class of Stochastic Resource Allocation Problems 
clear all
x1s = sdpvar(2,1,'full'); % 第一阶段变量
x1 = sdpvar;
x2s = sdpvar(2,1,'full');
x2 = sdpvar;
y1 = sdpvar(2,1,'full'); % 第二阶段变量
y2 = sdpvar(2,1,'full');
W_X1 = sdpvar(2,1,'full'); % PHA x1 动态系数
W_X2 = sdpvar(2,1,'full'); % PHA x2 动态系数
q1 = [24;28];
q2 = [28;32];
d1 = [500;300];
d2 = [100;300];
Pi_s = [0.4;0.6]; % 场景概率
rho =0.5; % PHA静态系数
rho_x1 = 0; % SEP-PHA的X1变量的静态系数
rho_x2 = 0; % SEP-PHA的X1变量的静态系数
k = 1; % PHA迭代计数器
kmax = 200; % PHA最大迭代次数
X1k = zeros(2,kmax);
X2k = zeros(2,kmax);
X1_bar = zeros(kmax,1);
X2_bar = zeros(kmax,1);
Gk = 999.*ones(kmax,1); % PHA阈值，这里设置为X1和X2的总偏差
epsilon = 1e-3;
ops = sdpsettings('solver','cplex');
% 第一阶段
% Obj_S1 = 100*x1 + 150*x2;
Cons_S1 = [x1s>=40, x2s>=20, x1s+x2s<=120];

% 第二阶段
% Obj_S2 = -sum( rho.*(q1.*y1+q2.*y2) ,'all');
Cons_S2 = [0<=y1<=d1, 0<=y2<=d2, 6.*y1+10.*y2<=60.*x1s, 8.*y1+5.*y2<=80.*x2s];

% 对比求解结果
Obj_c=100*x1+150*x2-sum( Pi_s.*(q1.*y1+q2.*y2) ,'all');
Cons_c = [ x1>=40, x2>=20, x1+x2<=120,  0<=y1<=d1, 0<=y2<=d2, 6.*y1+10.*y2<=60*x1, 8.*y1+5.*y2<=80*x2];
optimize(Cons_c, Obj_c,ops);
RES_C_X1 = value(x1);
RES_C_X2 = value(x2);
RES_C_Y1 = value(y1);
RES_C_Y2 = value(y2);

% PH算法
Obj01 = 100*x1s(1) + 150*x2s(1) - (q1(1)*y1(1)+q2(1)*y2(1)); % 注意这里不需要带场景概率！
Cons_PHA1 = []; 
Cons_PHA1 = [Cons_PHA1, x1s(1)>=40, x2s(1)>=20, x1s(1)+x2s(1)<=120];
Cons_PHA1 = [Cons_PHA1, 0<=y1(1)<=d1(1), 0<=y2(1)<=d2(1), 6*y1(1)+10*y2(1)<=60*x1s(1), 8*y1(1)+5*y2(1)<=80*x2s(1)];
Obj02 = 100*x1s(2) + 150*x2s(2) - (q1(2)*y1(2)+q2(2)*y2(2));
Cons_PHA2 = []; 
Cons_PHA2 = [Cons_PHA2, x1s(2)>=40, x2s(2)>=20, x1s(2)+x2s(2)<=120];
Cons_PHA2 = [Cons_PHA2, 0<=y1(2)<=d1(2), 0<=y2(2)<=d2(2), 6*y1(2)+10*y2(2)<=60*x1s(2), 8*y1(2)+5*y2(2)<=80*x2s(2)];

% 标准PH算法 
while Gk(k)>=epsilon && k<=kmax-1
    if k==1
        % k=1, 直接分别解两个场景问题
        optimize(Cons_PHA1, Obj01, ops); % 解场景1问题
        X1k(1,k) = value(x1s(1)); 
        X2k(1,k) = value(x2s(1)); 
        optimize(Cons_PHA2, Obj02, ops); % 解场景2问题
        X1k(2,k) = value(x1s(2)); 
        X2k(2,k) = value(x2s(2)); 
        X1_bar(k) = sum(Pi_s.*X1k(:,k),'all'); % 计算X_Bar
        X2_bar(k) = sum(Pi_s.*X2k(:,k),'all');
        W_X1 = rho.*(X1k(:,k)-X1_bar(k)); % 设定初始W
        W_X2 = rho.*(X2k(:,k)-X2_bar(k));
        k = k+1;
        Gk(k) = sum(Pi_s.*abs(X1k(:,k-1)-X1_bar(k-1)),'all')+sum(Pi_s.*abs(X2k(:,k-1)-X2_bar(k-1)),'all'); 
    else
        Obj1 = Obj01+W_X1(1)*x1s(1)+W_X2(1)*x2s(1)+0.5*rho*((x1s(1)-X1_bar(k-1))^2+(x2s(1)-X2_bar(k-1))^2); % 更新两个场景的增广拉格朗日目标函数
        Obj2 = Obj02+W_X1(2)*x1s(2)+W_X2(2)*x2s(2)+0.5*rho*((x1s(2)-X1_bar(k-1))^2+(x2s(2)-X2_bar(k-1))^2);
        optimize(Cons_PHA1, Obj1, ops); % 解场景1问题
        X1k(1,k) = value(x1s(1)); 
        X2k(1,k) = value(x2s(1)); 
        optimize(Cons_PHA2, Obj2, ops); % 解场景2问题
        X1k(2,k) = value(x1s(2)); 
        X2k(2,k) = value(x2s(2)); 
        X1_bar(k) = sum(Pi_s.*X1k(:,k),'all'); % 计算X_Bar
        X2_bar(k) = sum(Pi_s.*X2k(:,k),'all');
        W_X1 = W_X1 + rho.*(X1k(:,k)-X1_bar(k)); % 更新W
        W_X2 = W_X2 + rho.*(X2k(:,k)-X2_bar(k));
        k = k+1;
        Gk(k) = sum(Pi_s.*abs(X1k(:,k-1)-X1_bar(k-1)),'all')+sum(Pi_s.*abs(X2k(:,k-1)-X2_bar(k-1)),'all'); 
    end
end
Gk(1)=[];
figure;
title("标准PHA");
yyaxis left
plot(1:k-1, Gk(1:k-1),'r-','LineWidth',2);
hold on
yyaxis right
plot(1:k-1, X1k(:,1:k-1),'--g','LineWidth',1);
hold on
plot(1:k-1, X2k(:,1:k-1),'--b','LineWidth',1);
xlabel('迭代次数')


k = 1; % PHA迭代计数器
% 基于改进rho值的PH算法( rho值固定)
while Gk(k)>=epsilon && k<=kmax-1
    if k==1
        % k=1, 直接分别解两个场景问题
        optimize(Cons_PHA1, Obj01, ops); % 解场景1问题
        X1k(1,k) = value(x1s(1)); 
        X2k(1,k) = value(x2s(1)); 
        optimize(Cons_PHA2, Obj02, ops); % 解场景2问题
        X1k(2,k) = value(x1s(2)); 
        X2k(2,k) = value(x2s(2)); 
        X1_bar(k) = sum(Pi_s.*X1k(:,k),'all'); % 计算X_Bar
        X2_bar(k) = sum(Pi_s.*X2k(:,k),'all');
        rho_x1 = 100/max(1, sum(Pi_s.*abs(X1k(:,k)-X1_bar(k)),'all') ); % 确定"个性化" Rho
        rho_x2 = 100/max(1, sum(Pi_s.*abs(X2k(:,k)-X2_bar(k)),'all') );
        W_X1 = rho_x1.*(X1k(:,k)-X1_bar(k)); % 设定初始W
        W_X2 = rho_x2.*(X2k(:,k)-X2_bar(k));
        k = k+1;
        Gk(k) = sum(Pi_s.*abs(X1k(:,k-1)-X1_bar(k-1)),'all')+sum(Pi_s.*abs(X2k(:,k-1)-X2_bar(k-1)),'all'); 
    else
        Obj1 = Obj01+W_X1(1)*x1s(1)+W_X2(1)*x2s(1)+0.5*rho_x1*(x1s(1)-X1_bar(k-1))^2+0.5*rho_x2*(x2s(1)-X2_bar(k-1))^2; % 更新两个场景的增广拉格朗日目标函数
        Obj2 = Obj02+W_X1(2)*x1s(2)+W_X2(2)*x2s(2)+0.5*rho_x1*(x1s(2)-X1_bar(k-1))^2+0.5*rho_x2*(x2s(2)-X2_bar(k-1))^2;
        optimize(Cons_PHA1, Obj1, ops); % 解场景1问题
        X1k(1,k) = value(x1s(1)); 
        X2k(1,k) = value(x2s(1)); 
        optimize(Cons_PHA2, Obj2, ops); % 解场景2问题
        X1k(2,k) = value(x1s(2)); 
        X2k(2,k) = value(x2s(2)); 
        X1_bar(k) = sum(Pi_s.*X1k(:,k),'all'); % 计算X_Bar
        X2_bar(k) = sum(Pi_s.*X2k(:,k),'all');
        W_X1 = W_X1 + rho_x1.*(X1k(:,k)-X1_bar(k)); % 更新W
        W_X2 = W_X2 + rho_x2.*(X2k(:,k)-X2_bar(k));
        k = k+1;
        Gk(k) = sum(Pi_s.*abs(X1k(:,k-1)-X1_bar(k-1)),'all')+sum(Pi_s.*abs(X2k(:,k-1)-X2_bar(k-1)),'all'); 
    end
end
Gk(1)=[];
figure;
title("改进Rho的PHA");
yyaxis left
plot(1:k-1, Gk(1:k-1),'r-','LineWidth',2);
hold on
yyaxis right
plot(1:k-1, X1k(:,1:k-1),'--g','LineWidth',1);
hold on
plot(1:k-1, X2k(:,1:k-1),'--b','LineWidth',1);
xlabel('迭代次数')
    
    
    
    
    