% This is the matlab code for the multi-objective bonobo optimizer 
% with grid-index approach(MOBO1).
% This is written for solving unconstrained optimization problems. 
% However, it can also solve constrained optimization 
% problems with penalty function approaches.
% For details of the MOBO1 algorithm, kindly refer and cite as mentioned below:
% Das, A.K., Nikum, A.K., Krishnan, S.V. et al. Multi-objective Bonobo Optimizer 
%(MOBO): an intelligent heuristic for multi-criteria optimization. 
% Knowl Inf Syst (2020). https://doi.org/10.1007/s10115-020-01503-x
% For any query, please email to: amit.besus@gmail.com

% I acknowledge that this version of MOBO1 has been written using
% a large portion of the following code:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for                                                  %
%                                                                   %
%  Multi-Objective Particle Swarm Optimization (MOPSO)              %
%  Version 1.0 - Feb. 2011                                          %
%                                                                   %
%  According to:                                                    %
%  Carlos A. Coello Coello et al.,                                  %
%  "Handling Multiple Objectives with Particle Swarm Optimization," %
%  IEEE Transactions on Evolutionary Computation, Vol. 8, No. 3,    %
%  pp. 256-279, June 2004.                                          %
%                                                                   %
%  Developed Using MATLAB R2009b (Version 7.9)                      %
%                                                                   %
%  Programmed By: S. Mostapha Kalami Heris                          %
%                                                                   %
%         e-Mail: sm.kalami@gmail.com                               %
%                 kalami@ee.kntu.ac.ir                              %
%                                                                   %
%       Homepage: http://www.kalami.ir                              %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;close all;clear all;
tic;  % CPU time measure

M=6;  % No. of Objectives
I = 2; %Dimesionality of Pareto-Optimal Surface
k = 20;
dim = M-1+k;  % No. of variables - dimensionality
fobj=@(x)dtlz7(x,M);% Objective function
Var_min=zeros(1,dim);  % Lower variable Boundaries 
Var_max=ones(1,dim);   % Upper variable Boundaries
%% Common parameters of  MOBO similar to  other optimization algorithms
N=156; % No. of bonobos in the population, i.e. population size
max_it=250;  % Maximum number of iterations
%% Algorithm-specific Parameters for MOBO (user should set suitable values of the parameters for their problem)
p_xgm_initial=1/dim; % Initial probability for extra-group mating (generally 1/d for higher dimensions)
scab=1.55;  % Sharing cofficient for alpha bonobo (Generally 1-2)
scsb=1.4;   % Sharing coefficient for selected bonobo(Generally 1-2)
rcpp=0.004; % Rate of change in  phase probability (Generally 1e-3 to 1e-2)
tsgs_factor_max=0.07;% Max. value of temporary sub-group size factor
nGrid=10; % Number of grids
nRep=N;% Repository size
alpha=0.1; % Inflation Rate
beta=2; % Leader Selection Pressure
gamma=2;% Deletion Selection Pressure
%% There is no need to change anything below this %%
%% Initialization
empty_particle.bonobo=[];
empty_particle.Cost=[];
empty_particle.isdominated=[];
empty_particle.GridIndex=[];
empty_particle.GridSubIndex=[];
empty_particle.CostOriginal=[];

pop=repmat(empty_particle,N,1);
for i=1:N
    pop(i).bonobo=unifrnd(Var_min,Var_max,[1 dim]);
    pop(i).Cost = fobj(pop(i).bonobo');
    pop(i).CostOriginal = pop(i).Cost;
end

% N-D to 2-D Conversion
temp=som2d(pop);
for i=1:N
    pop(i).Cost = temp(i).cost2d';
end

pop=Determine_Domonation(pop);
rep=pop(~[pop.IsDominated]);
Grid=GenerateGrid(rep,nGrid,alpha);
num=numel(rep);
for i=1:num
    rep(i)=GridIndexFinding(rep(i),Grid);
end
for i=1:N
    pop(i)=GridIndexFinding(pop(i),Grid);
end
%% Initialization of other parameters
npc=0; % Negative phase count
ppc=0; % Positive phase count
p_xgm=p_xgm_initial; % Probability for extra-group mating
tsgs_factor_initial=0.5*tsgs_factor_max; % Initial value for temporary sub-group size factor
tsgs_factor=tsgs_factor_initial; % Temporary sub-group size factor
p_p=0.5; % Phase probability
p_d=0.5; % Directional probability
nObj=numel(rep(1).Cost);
it=1;
nem=zeros(max_it+1,1);
tempdv=zeros(max_it+1,1);
nem(it)=min(nRep,size(rep,1));
rep_costs=[rep.Cost]';
ss=minmax(rep_costs');
tempdv(it)=mean(std(rep_costs)./(ss(:,2)-ss(:,1))');
%% Main Loop of BO
while(it<=max_it)
    %% loop for sorting
    %%select alpha
    tsgs_max=max(2,ceil(N*tsgs_factor)); % Maximum size of the temporary sub-group
    for i=1:N
        [~,alphabonobo]=LeaderSelection(rep,beta);
        newbonobo=zeros(1,dim);
        B = 1:N;
        %% Determining the actual size of the temporary sub-group
        tsg=randi([2 tsgs_max]);
        %% Selection of pth Bonobo using fission-fusion social strategy & flag value determination
        q=randsample(B,tsg);
        temp=pop(q);
        temp_q2=zeros(tsg,1);
        for j=1:tsg
            if (~(temp(j).IsDominated))
                temp_q2(j)=q(j);
            end
        end
        temp_q2=nonzeros(temp_q2);
        p1=numel(temp_q2);
        if (p1==0)
            if(tsg==2)
                p=q(randi(2));
            else
                p=selectp(temp,tsg);
                p=q(p);
            end
        else
            if (p1==1)
                p=temp_q2;
            elseif (p1==2)
                p=temp_q2(randi(2));
            else
                temp_rep2=pop(temp_q2);
                p=selectp(temp_rep2,p1);
                p=temp_q2(p);
            end
        end
        flag=0;
        %% Creation of newbonobo
        if(rand<=p_p)
            flag=1;
            r1=rand(1,dim); %% Promiscuous or restrictive mating strategy
            newbonobo=pop(i).bonobo+scab*r1.*(alphabonobo.bonobo-pop(i).bonobo)+scsb*(1-r1).*(pop(i).bonobo-pop(p).bonobo);
        else
            for j=1:M
                if(rand<=p_xgm)
                    rand_var=rand; %% Extra group mating strategy
                    if(alphabonobo.bonobo(j)>=pop(i).bonobo(j))
                        if(rand<=(p_d))
                            beta1=exp(((rand_var)^2)+rand_var-(2/rand_var));
                            newbonobo(1,j)=pop(i).bonobo(j)+beta1*(Var_max(j)-pop(i).bonobo(j));
                        else
                            beta2=exp((-((rand_var)^2))+(2*rand_var)-(2/rand_var));
                            newbonobo(1,j)=pop(i).bonobo(j)-beta2*(pop(i).bonobo(j)-Var_min(j));
                        end
                    else
                        if(rand<=(p_d))
                            beta1=exp(((rand_var)^2)+(rand_var)-2/rand_var);
                            newbonobo(1,j)=pop(i).bonobo(j)-beta1*(pop(i).bonobo(j)-Var_min(j));
                        else
                            beta2=exp((-((rand_var)^2))+(2*rand_var)-2/rand_var);
                            newbonobo(1,j)=pop(i).bonobo(j)+beta2*(Var_max(j)-pop(i).bonobo(j));
                        end
                    end
                else
                    if((rand<=p_d)) %% Consortship mating strategy
                        newbonobo(1,j)=pop(i).bonobo(j)+(exp(-rand))*(pop(i).bonobo(j)-pop(p).bonobo(j));
                    else
                        newbonobo(1,j)=pop(p).bonobo(j);
                    end
                end
            end
        end
        %% Clipping
        for j=1:M
            if(newbonobo(1,j)>Var_max(j))
                newbonobo(1,j)=Var_max(j);
            end
            if(newbonobo(1,j)<Var_min(j))
                newbonobo(1,j)=Var_min(j);
            end
        end
        
        newcost.CostOriginal = fobj(newbonobo'); 
        r = norm(newcost.CostOriginal);
%         x1=0; y1=0;
%         if(i == 1)
%             ref = pop(2).CostOriginal;
%             x1 = pop(1).Cost(1);
%             y1 = pop(1).Cost(2);
%         else 
%             ref = pop(i-1).CostOriginal;
%             x1 = pop(i-1).Cost(1);
%             y1 = pop(i-1).Cost(2);
%         end
%         distance = norm(newcost.CostOriginal - ref);
%         syms sol;
%         sol = solve((r*cos(sol) - x1)^2 + (r*sin(sol) - y1)^2 == distance^2, sol, 'Real', true);
%         sol = double(sol);
%         if(numel(sol) && sol(1) >0 && sol(1)<pi/2)
%             alpha = sol(1);
%         elseif (numel(sol) && sol(2) >0 && sol(2)<pi/2)
%                 alpha = sol(2);
%         else
%                 alpha = rand*range([0 pi/2]);
%         end 
        alpha = rand*range([0 pi/2]);
        newcost.Cost = [r*cos(alpha) ; r*sin(alpha)];
           
        %% New bonobo acceptance criteria
        if(flag==1 || Dominate(newcost,pop(i)) || (rand<=(p_xgm)))
            pop(i).CostOriginal=newcost.CostOriginal;
            pop(i).Cost = newcost.Cost;
            pop(i).bonobo=newbonobo;
        end
    end
    pop=Determine_Domonation(pop);
    rep=[rep
        pop(~[pop.IsDominated])];
    [~, idx] = unique([rep.Cost].', 'rows');
    rep = rep(idx);
    rep=Determine_Domonation(rep);
    % Keep only Non-Dminated Memebers in the Repository
    rep=rep(~[rep.IsDominated]);
    Grid=GenerateGrid(rep,nGrid,alpha);
    % Update Grid Indices
    for i=1:numel(rep)
        rep(i)=GridIndexFinding(rep(i),Grid);
    end
    % Check if Repository is Full
    if numel(rep)>nRep
        Extra=numel(rep)-nRep;
        for e=1:Extra
            rep=RepMemberDeletion(rep,gamma);
        end
    end
    nem(it+1)=numel(rep);
    if(nem(it+1)>1)
        rep_costs=[rep.Cost]';
        ss=minmax(rep_costs');
        tempdv(it+1)=mean(std(rep_costs)./(ss(:,2)-ss(:,1))');
        if(nem(it+1)>=nem(it) && tempdv(it+1)>tempdv(it))
            pp=1;
        else
            pp=0;
        end
    else
        pp=0;
    end
    %% Parameters updation
    if(pp==1)
        npc=0; %% Positive phase
        ppc=ppc+1;
        cp=min(0.5,(ppc*(rcpp)));
        p_xgm=p_xgm_initial;
        p_p=0.5+cp;
        tsgs_factor=min(tsgs_factor_max,(tsgs_factor_initial+ppc*(rcpp^2)));
    else
        npc=npc+1; %% Negative phase
        ppc=0;
        cp=-(min(0.5,(npc*(rcpp))));
        p_xgm=min(0.5,p_xgm_initial+npc*(rcpp^2));
        tsgs_factor=max(0,(tsgs_factor_initial-npc*(rcpp^2)));
        p_p=0.5+cp;
    end
    clear temp temp_rep2 temp_q2;
    % Show Iteration Information
    disp(['Iteration = ' num2str(it)]);
    it=it+1;
end

% Plot of the results
PlotObjectiveFunction(pop,rep,nObj);
toc;
