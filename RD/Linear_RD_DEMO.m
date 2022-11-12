close all; clear all;
ops = sdpsettings('solver','cplex');
imax = 10;
ND =32;
NT = 24;
x = binvar(ND,NT,'full');
y = sdpvar(ND,NT,'full');
z = binvar(ND,NT,'full');
lamb1=sdpvar(ND,NT,imax,'full');
lamb2=sdpvar(ND,NT,imax,'full');
lamb3=sdpvar(ND,NT,imax,'full');
AUX = sdpvar(ND,NT,imax,'full'); %双线性项
BigM = 1e5;
cons_DV = [lamb1>=0, lamb2<=0, lamb3>=0];
UB = 1e5;
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
        MP_cons =[y>=4-2.*x+3.*z, y<=10-z, y>=3-z]; % lamb1,
        CostMP=sum(x-y,'all');
        MP = optimize(MP_cons, CostMP,ops); % 求解MP
    end
    LB = value(CostMP);
    LB_save(K)=LB;
    x_star = value(x);
    CostMP0 = value(sum(x,'all'));
    CostSP1=sum(y,'all');
    SP1_cons=[y>=4-2.*x_star+3.*z, y<=10-z, y>=3-z];
    optimize(SP1_cons, CostSP1,ops);
    CostSP1_star = value( CostSP1);
    SP2_cons=[SP1_cons, CostSP1<=CostSP1_star];
    CostSP2 = sum(-y,'all');
    optimize(SP2_cons, CostSP2,ops);
    z_star = value(z);
    CostSP2_star = value(CostSP2);
    UB = CostMP0 + CostSP2_star;
    UB_save(K)=UB;
    Dual_cons = [lamb1(:,:,K)+lamb2(:,:,K)+lamb3(:,:,K)==1];
    Dual_cons = [Dual_cons, -BigM.*(1-x)<=AUX(:,:,K)-lamb1(:,:,K), AUX(:,:,K)-lamb1(:,:,K)<=BigM.*(1-x), -BigM.*x<=AUX(:,:,K), AUX(:,:,K)<=BigM.*x];
    Dual_cost = sum((4*ones(ND,NT)+3.*z_star).*lamb1(:,:,K) + (10*ones(ND,NT)-z_star).*lamb2(:,:,K) + (3*ones(ND,NT)-z_star).*lamb3(:,:,K) - 2.*AUX(:,:,K),'all') ; % 包含双线性项
    RD_Cuts = [RD_Cuts, Dual_cons, CostSP1<=Dual_cost];
    MP_cons = [MP_cons, RD_Cuts];
    optimize(MP_cons,CostMP,ops);
    K=K+1;
end
figure
plot(1:K-1, UB_save);
hold on
plot(1:K-1, LB_save);
