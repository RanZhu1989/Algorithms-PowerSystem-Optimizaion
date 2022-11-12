%% 经典Benders分解算法Demo

% 问题
% min x1 + x2 + x3 + x4 + x5 + 7*(y1+y2+y3+y4+y5)
% st.              x1 + x4 + x5 = 8            (mu1)
%                  x2 + x5 = 3                     (mu2)
%                  x3 + x4 = 5                     (mu3)
%     x1<=8y1 (mu4) ,   x2<=3y2   (mu5),    x3<=5y3  (mu6)   ,  x4<=5y4    (mu7),    x5<=3y5    (mu8)
%    x1~x5>=0,   y1~y5 = binary
clear all
ops = sdpsettings('solver','cplex');
x=sdpvar(5,1,'full');
y=binvar(5,1, 'full');
mu=sdpvar(8,1,'full');
rho=sdpvar(1);
%% 求解器求解
optimize([x(1) + x(4) + x(5) == 8, x(2) + x(5) == 3, x(3) + x(4) == 5, x(1)<=8*y(1), x(2)<=3*y(2), x(3)<=5*y(3), x(4)<=5*y(4), x(5)<=3*y(5), x>=0],sum(x,1) + 7*sum(y,1), ops);
C_Obj = value( sum(x,1) + 7*sum(y,1));

%% 经典BD分解
% MP : min  7*(y1+y2+y3+y4+y5) + rho ,   st.  y1~y5 = binary, BD cuts
% SP: max 8*mu1+3*mu2+5*mu3+8*y1_star*mu4+3*y2_star*mu5+5*y3_star*mu6+5*y4_star*mu7+3*y5_star*mu8
%      st.       mu1+mu4<=1,  mu2+mu5<=1, mu3+mu5<=1, mu1+mu3+mu7<=1,  mu1+mu2+mu8<=1
%                 mu4~mu8<=0
k_max = 1000;
LB = ones(k_max,1);
LB(1) = -100;
UB = ones(k_max,1);
UB(1) = 100;

k = 1;
e = 1e-4;
BD_Cuts = [];
MP_Obj = 7*sum(y,1) + rho;
MP_Cons = [rho>=0];
y_star = zeros(5,1);
mu_star = zeros(8,1);
while UB(k)-LB(k)>=e && k<=k_max
    MP = optimize([MP_Cons, BD_Cuts], MP_Obj, ops);
    if MP.problem ~= 0
        fprintf('MP infeasiable in iteration %k ! \n', k);
     error('MP Infeasiable');
    end
    y_star = value(y);
    LB(k+1) = value(MP_Obj);
    SP_Obj = 8*mu(1)+3*mu(2)+5*mu(3)+8*y_star(1)*mu(4)+3*y_star(2)*mu(5)+5*y_star(3)*mu(6)+5*y_star(4)*mu(7)+3*y_star(5)*mu(8);
    SP_Cons = [mu(1)+mu(4)<=1,  mu(2)+mu(5)<=1, mu(3)+mu(5)<=1, mu(1)+mu(3)+mu(7)<=1,  mu(1)+mu(2)+mu(8)<=1, mu(4:8,1)<=0];
    SP = optimize(SP_Cons, -SP_Obj, ops);
    if SP.problem ~= 0 % 判断SP是否有解
        % 如果无解，则需要返回feasibility cuts
       %  求解一个feasibility subproblem的对偶问题， 返回feasibility cuts
        %  ------- feasibility subproblem (FS) ---------
        % min sum(S,1)
        % st. x1 + x4 + x5 +s1+ -s1- = 8 (mu1),  x2 + x5 +s2+ - s2-= 3 (mu2),  x3 + x4 +s3+ - s3- = 5 mu(3), 
        %     x1+s4+ - s4- <=8y1 (mu4),   x2+s5+ - s5- <=3y2 (mu5),  x3+s6+ - s6- <=5y3 (mu6),  x4+s7+ - s7- <=5y4  (mu7),    x5+s8+ - s8- <=3y5 mu(8)
        %    S>=0.
        % ---------------- Dual of FS is ------------------------
        % max 8*mu1+3*mu2+5*mu3+8*y1_star*mu4+3*y2_star*mu5+5*y3_star*mu6+5*y4_star*mu7+3*y5_star*mu8
        % st.    mu1+mu4<=0,  mu2+mu5<=0, mu3+mu5<=0, mu1+mu3+mu7<=0,  mu1+mu2+mu8<=0
        %          mu1~mu8<=1       mu4~mu8<=0,   
        FSP_Cons = [mu(1)+mu(4)<=0,  mu(2)+mu(5)<=0, mu(3)+mu(5)<=0, mu(1)+mu(3)+mu(7)<=0,  mu(1)+mu(2)+mu(8)<=0, mu<=1, mu(4:8,1)<=0];      
        optimize(FSP_Cons, -SP_Obj, ops);
        mu_star = value(mu);
        BD_Cuts = [BD_Cuts, 8*mu_star(1)+3*mu_star(2)+5*mu_star(3)+8*y(1)*mu_star(4)+3*y(2)*mu_star(5)+5*y(3)*mu_star(6)+5*y(4)*mu_star(7)+3*y(5)*mu_star(8)<=0]; 
        UB(k+1) = UB(k);
    else % 若有解，则返回optimality cut
        mu_star = value(mu);
        BD_Cuts = [BD_Cuts, 8*mu_star(1)+3*mu_star(2)+5*mu_star(3)+8*y(1)*mu_star(4)+3*y(2)*mu_star(5)+5*y(3)*mu_star(6)+5*y(4)*mu_star(7)+3*y(5)*mu_star(8)<=rho];
        UB(k+1) = 7*sum(y_star,1)+value(SP_Obj);
    end
    k=k+1;
end
plot(1:k, LB(1:k),'-o','LineWidth',1);
hold on
plot(1:k, UB(1:k),'-x','LineWidth',1);
xlabel('迭代次数')
hh = legend('LB','UB');%
hh.Orientation = 'horizontal';