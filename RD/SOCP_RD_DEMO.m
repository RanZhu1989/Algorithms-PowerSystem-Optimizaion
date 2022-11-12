% SOCP-RD算法收敛性简单算例

close all ; clear all
ops = sdpsettings('solver','cplex');
imax = 10;
NT=10;
x=binvar(NT,1,'full');
y1=sdpvar(NT,1,'full');
y2=sdpvar(NT,1,'full');
z =binvar(NT,1,'full');
lamb1=sdpvar(NT,imax,'full');
lamb2=sdpvar(NT,imax,'full');
lamb3=sdpvar(NT,imax,'full');
lamb4=sdpvar(NT,imax,'full');
lamb5=sdpvar(NT,imax,'full');
u=sdpvar(NT,imax,'full');
v=sdpvar(NT,imax,'full');
AUX = sdpvar(NT,imax,'full'); %双线性项
BigM = 1e3;
cons_DV = [lamb1<=0, lamb2>=0, lamb3<=0, lamb4>=0, lamb5<=0];
UB = 1e3;
LB = 0;
UB_save = [];
LB_save = [];
UB_save(1)=UB;
LB_save(1)=LB;
RD_Error = 0.01;
K = 1; % 当前RD迭代次数
RD_Cuts = [cons_DV];
while UB-LB>=RD_Error && K<imax
    if K==1
        MP_cons =[x+y1+y2+z<=3.5, y1>=-1, y1<=2, y2>=-2, y2<=5]; % lamb1, lamb2, lamb3,  lamb4, lamb5, u, v
        for i = 1:NT
             MP_cons =[MP_cons,  cone(y1(i),y1(i)+y2(i))];
        end
        CostMP=sum(2*x-y1+3*y2, 'all');
        MP = optimize(MP_cons, CostMP,ops); % 求解MP
    end
    LB = value(CostMP);
    LB_save(K)=LB;
    x_star = value(x);
    CostMP0 = sum(2*x_star,'all');
    CostSP1=sum(y1-2*y2,'all');
    SP1_cons=[x_star+y1+y2+z<=3.5, y1>=-1, y1<=2, y2>=-2, y2<=5];
    for i=1:NT
        SP1_cons=[SP1_cons, cone(y1(i),y1(i)+y2(i))];
    end
    optimize(SP1_cons, CostSP1,ops);
    CostSP1_star = value( CostSP1);
    SP2_cons=[SP1_cons, CostSP1<=CostSP1_star];
    CostSP2 = sum(-y1+3*y2,'all');
    optimize(SP2_cons, CostSP2, ops);
    z_star = value(z);
    CostSP2_star = value(CostSP2);
    UB = CostMP0 + CostSP2_star;
    UB_save(K)=UB;
    Dual_cons = [ lamb1(:,K)+lamb2(:,K)+lamb3(:,K)+u(:,K)+v(:,K)==1,  lamb1(:,K)+lamb4(:,K)+lamb5(:,K)+v(:,K)==-2];
    for i =1:NT
         Dual_cons = [Dual_cons, cone(u(i,K),v(i,K))];
    end
    Dual_cons = [Dual_cons, -BigM.*(1-x)<=AUX(:,K)-lamb1(:,K), AUX(:,K)-lamb1(:,K)<=BigM.*(1-x), -BigM.*x<=AUX(:,K), AUX(:,K)<=BigM.*x];
    Dual_cost = sum(  (3.5-z_star).*lamb1(:,K)-AUX(:,K)-lamb2(:,K)+2.*lamb3(:,K)-2.*lamb4(:,K)+5.*lamb5(:,K), 'all'); % 包含双线性项
    RD_Cuts = [RD_Cuts, Dual_cons, CostSP1<=Dual_cost];
    MP_cons = [MP_cons, RD_Cuts];
    optimize(MP_cons,CostMP,ops);
    K=K+1;
end
figure
plot(1:K-1, UB_save);
hold on
plot(1:K-1, LB_save);
