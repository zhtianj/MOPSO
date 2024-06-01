global data1 data2 data3 data4 data5 data6 data7 data8 u a b h kt po xy xx yy rw an bn
%============================= 
%Clogging parameters
%============================= 
tmax=24;
t0=0;
u=7000;
a=20;
b=5;
h=0.1;
kt=30/86400;
po=0.3;
c=0.045;
O2=0.212;
fe2=1.047;
an=1;
bn=50;
rw=0.015;
b=b+rw;
xx=linspace(rw,b+rw,bn+1);
yy=linspace(0,a,an+1);
[x,y] = meshgrid(xx,yy);
xy = [x(:),y(:)]';
%a,b,u,h,k0,poroini,O2,c,fe2,t0,tmax,an,bn
[data1,data2,data3,data4,data5,data6,data7,data8,data9,data10]= Feclogg(a,b,u,h,kt,po,O2,c,fe2,t0,tmax,an,bn,rw);
CostFunction = @(x) evaluate_objective(x);  %function
%============================= 
%MOPSO parameters
%============================= 
nVar = 2;                                     %Number of variables
VarSize = [1 nVar];                            
VarMin = [6 0.1];                                    %Definitional domain of a variable value
VarMax = [24 1];                                  %Note: The function variable cannot have a negative value.
VarMin2 = 0.1;                                    %Definitional domain of a variable value
VarMax2 = 1;  
MaxIt = 50;                                   %iterations
N = 30;                                        %swarm size
nRep =20;                                     %archive size
w = 0.7298;                                       %Inertia weighting factor
wdamp = 0.99;                                  %Inertia weight decay rate
c1 = [1.496 1.496];                                      %Individual Learning Factor
c2 = [1.496 1.496];                                      %Global learning factor
nGrid = 5;                                    
alpha = 0.05;                                  
beta = 2;                                     
gamma = 2;                                     
mu = 0.1;                                      
empty_particle.Position = [];                  
empty_particle.Velocity = [];                  
empty_particle.Cost = [];  
empty_particle.Best.Position = [];  
empty_particle.Best.Cost = [];            
empty_particle.IsDominated = [];              
empty_particle.GridIndex = [];                
empty_particle.GridSubIndex = [];             
pop = repmat(empty_particle,N,1);              
   Position=zeros(MaxIt,nRep,2);
   Value=zeros(MaxIt,nRep,2);
for i = 1:N 
     pop(i).Position = [unifrnd(VarMin(1),VarMax(1),1);unifrnd(VarMin(2),VarMax(2),1)];
     pop(i).Velocity = zeros(VarSize);
     pop(i).Cost = CostFunction(pop(i).Position);
     pop(i).Best.Position = pop(i).Position;
     pop(i).Best.Cost = pop(i).Cost;
end

pop = DetermineDomination(pop);
rep = pop(~[pop.IsDominated]);
Grid = CreateGrid(rep,nGrid,alpha);
for i = 1:numel(rep)
 rep(i) = FindGridIndex(rep(i),Grid);
end

%MOPSO
 for it = 1:MaxIt
     for i = 1:N 
         leader = SelectLeader(rep,beta);   
         rep = [rep;pop(~[pop.IsDominated])];   
         pop(i).Velocity =[w*pop(i).Velocity(1) + ...
             c1(1)*rand(1).*(pop(i).Best.Position(1)-pop(i).Position(1))+ ...
             c1(2)*rand(1).*(leader.Position(1)-pop(i).Position(1));w*pop(i).Velocity(2) + ...
             c2(1)*rand(1)*0.05*(pop(i).Best.Position(2)-pop(i).Position(2))+ ...
             c2(2)*rand(1)*0.05*(leader.Position(2)-pop(i).Position(2))];   
         pop(i).Position = pop(i).Position+pop(i).Velocity;  
%          pop(i).Position = limitToPosition(pop(i).Position,VarMin,VarMax);   
         pop(i).Position(1) = max(min(pop(i).Position(1),24), 6);
         pop(i).Position(2) = max(min(pop(i).Position(2), 1), 0.1);
         pop(i).Cost = CostFunction(pop(i).Position); 
         
         pm = (1-(it-1)/(MaxIt-1)^(1/mu));  
         NewSol.Position = [Mutate(pop(i).Position(1),pm,VarMin(1),VarMax(1));Mutate(pop(i).Position(2),pm,VarMin(2),VarMax(2))];
         NewSol.Cost = CostFunction(NewSol.Position);   
         if Dominates(NewSol,pop(i))
             pop(i).Position = NewSol.Position;
             pop(i).Cost  = NewSol.Cost;
         else
             if rand < 0.5
                 pop(i).Position = NewSol.Position;
                 pop(i).Cost = NewSol.Cost;
             end
         end
         if Dominates(pop(i),pop(i).Best)   
             pop(i).Best.Position = pop(i).Position;
             pop(i).Best.Cost = pop(i).Cost;
         else 
             if rand <0.5
                 pop(i).Best.Position = pop(i).Position;
                 pop(i).Best.Cost = pop(i).Cost;
             end
         end
     end   
     
     rep =  DetermineDomination(rep);
     rep = rep(~[rep.IsDominated]);
     Grid = CreateGrid(rep,nGrid,alpha); 
    for i =1:numel(rep) 
        rep(i) = FindGridIndex(rep(i),Grid); 
    end 
    if numel(rep) > nRep 
        Extra = numel(rep)-nRep; 
        for e = 1:Extra 
            rep = DeleteOneRepMemebr(rep,gamma); 
        end 
    end 
     
    disp(['iterations =',num2str(it)]); 
    w = w*wdamp; 

   for j=1:N
         A1(j)=sqrt((pop(j).Best.Cost(1)-pop(j).Cost(1))^2+(pop(j).Best.Cost(2)-pop(j).Cost(2))^2);
    end
    MID1(it)=mean(A1);
    for j=1:nRep
         A2(j)=sqrt((rep(j).Best.Cost(1)-rep(j).Cost(1))^2+(rep(j).Best.Cost(2)-rep(j).Cost(2))^2);
    end
    MID2(it)=mean(A2);
%     record(it)=sum;
fielPosition=zeros(nRep,2);
fieldValue=zeros(nRep,2);
for i=1:nRep
    fieldValue(i,:) = rep(i).Cost;
end
for i=1:nRep
    fielPosition(i,:) = rep(i).Position;
end
   con1(it) = calculateConvergence(pop, leader);
   con2(it) = calculateConvergence(rep, leader);
   CON1(it)=calculateMID(fieldValue);
   CON2(it)=calculateMID(fielPosition);
Position(it,:,:)=fielPosition;
  Value(it,:,:)=fieldValue;

end 

figure(1); 
location = [rep.Cost];   
scatter(location(1,:),location(2,:),'filled','b'); 
xlabel('f1');ylabel('f2'); %zlabel('f3');
title('Pareto optimal boundary'); 
box on; 



%============================= 
%Calculate the value of the objective function  
%============================= 
function f= evaluate_objective(t)
 
    global data1 data2 data3 data4 data5 data6 data7 data8  a b h kt po c an bn u xx yy xy rw
    t0=0;
    dt=floor((t(1)-0)*5+1);
    coo1=data1(dt,:);
    coo2=data2(dt,:);
    coo3=data8(dt,:)+data7(dt,:);
    coo4=data4(dt,:);
    coo5=data5(dt,:);
    ty1=reshape(coo1,an+1,bn+1);
    ty2=reshape(coo2,an+1,bn+1);
    ty3=reshape(coo3,an+1,bn+1);
    ty4=reshape(coo4,an+1,bn+1);
    ty5=reshape(coo5,an+1,bn+1);
    typ1=ty1(:,2);
    ty1(:,1)=typ1;
    typ2=ty2(:,2);
    ty2(:,1)=typ2;
    typ3=ty3(:,2);
    % typ5=ones(1,an+1)*c2;
    ty3(:,1)=typ3';
    typ4=ty4(:,2);
    ty4(:,1)=typ4;
    typ5=ty5(:,2);
    ty5(:,1)=typ5;
    x=xx(1:bn+1);
    y=yy(1:an+1);
    mphinterpolationfile('notekt.txt','sectionwise',ty2,x,y);
    mphinterpolationfile('notepo.txt','sectionwise',ty1,x,y);
    mphinterpolationfile('noteco2.txt','sectionwise',ty3,x,y);
    mphinterpolationfile('notecs.txt','sectionwise',ty4,x,y);
    mphinterpolationfile('notecst.txt','sectionwise',ty5,x,y);
   y2=commodel2(a,b,u,h,kt,po,c,t0,t(2),xy,rw);
   Q=data7(dt,2);
   for j=1:length(y2) 
       ff(j)=po-y2(j,3);
   end
   fk=sum(ff);
   fq=t(1)-t(2)/t(1);
%    values(i,1)=1/fk*fq;
   f(1)=sum(ff)/101;
  f(2)=t(2)/t(1);
  f=[f(1);f(2)];
end

function flow2=commodel2(a,b,u,h,kt,po,c2,t0,t2,xx,rw)
model = mphopen('dfe3');
model.param.set('rx',[num2str(b),'[m]']);
model.param.set('rh',[num2str(a),'[m]']);
model.param.set('H',[num2str(h),'[m]']);  
model.param.set('rw',[num2str(rw),'[m]']);
model.param.set('k0',[num2str(kt),'[m/s]']);
model.param.set('po', po);
model.param.set('t0',[num2str(t0),'[h]']);
model.param.set('t1', [num2str(t2),'[h]']);
model.func('int1').discardData;
model.func('int1').importData;
model.func('int2').discardData;
model.func('int2').importData;
model.func('int3').discardData;
model.func('int3').importData;
model.func('int4').discardData;
model.func('int4').importData;
model.func('int5').discardData;
model.func('int5').importData;
model.func('int6').discardData;
model.func('int6').importData;
model.study('std1').run();
da1= mphinterp(model,'K','coord',xx);
da2= mphinterp(model,'n','coord',xx);
da3= mphinterp(model,'c','coord',xx);
da4= mphinterp(model,'dl.Hp','coord',xx);
da5= mphinterp(model,'c2','coord',xx);
da6= mphinterp(model,'C_s','coord',xx);
da7= mphinterp(model,'C_st','coord',xx);
dt=floor((t2-t0)*10+1);
coo1=da1(dt,:);
coo2=da2(dt,:);
coo3=da3(dt,:);
coo4=da4(dt,:);
coo5=da5(dt,:);
flow2=[xx',coo1',coo2',coo3',coo4',coo5'];
end
%============================= 
%Return to the dominant solution
%============================= 
function pop =DetermineDomination(pop) 
        nPop = numel(pop); 
        for i =1:nPop 
            pop(i).IsDominated = false;   
        end 
        for i = 1:nPop-1 
            for j = i+1:nPop 
                if Dominates(pop(i),pop(j)) 
                    pop(j).IsDominated = true; 
                end 
                    if Dominates(pop(j),pop(i)) 
                        pop(i).IsDominated = true; 
                    end 
            end 
        end 
end 

%============================= 
%Determine the dominance of two target values x,y 
%============================= 
function b = Dominates(x,y) 
        if isstruct(x) 
            x=x.Cost; 
        end 
        if isstruct(y) 
            y=y.Cost; 
        end 
        b=all(x<=y) && any(x<y); 
end 
%============================= 
%Creating a Raster Matrix 
%============================= 
function Grid = CreateGrid(pop,nGrid,alpha) 
        c = [pop.Cost]; 
        cmin = min(c,[],2); 
        cmax = max(c,[],2); 
        dc = cmax-cmin; 
        cmin = cmin-alpha*dc; 
        cmax = cmax+alpha*dc; 
        nObj = size(c,1); 
        empty_grid.LB = []; 
        empty_grid.UB = []; 
        Grid = repmat(empty_grid,nObj,1); 
         
        for j = 1:nObj 
            cj = linspace(cmin(j),cmax(j),nGrid+1); 
            Grid(j).LB = [-inf cj]; 
            Grid(j).UB = [cj +inf]; 
        end 
end 
%============================= 
%Raster Index Positioning 
%============================= 
function particle = FindGridIndex(particle,Grid) 
        nObj = numel(particle.Cost); 
        nGrid = numel(Grid(1).LB); 
        particle.GridSubIndex = zeros(1,nGrid); 
        for j = 1:nObj 
            particle.GridSubIndex(j) = find(particle.Cost(j)<=Grid(j).UB,1,'first'); 
        end 
        particle.GridIndex = particle.GridSubIndex(1); 
        for j = 2:nObj   
            particle.GridIndex = particle.GridIndex-1; 
            particle.GridIndex = nGrid*particle.GridIndex; 
            particle.GridIndex = particle.GridIndex + particle.GridSubIndex(j); 
        end 
end

%============================= 
%Restricting the scope of change of a variable to be within the domain of definition
%============================= 
function Position = limitToPosition(Position,VarMin,VarMax)     
        for i =1:size(Position,2) 
            if Position(i)<VarMin 
                Position(i) = VarMin; 
            elseif Position(i) > VarMax 
                Position(i) = VarMax; 
            end 
        end 
end 
%============================= 
%Identify an optimal individual from the globally dominant individuals
%============================= 
function leader = SelectLeader(rep,beta) 
        GI = [rep.GridIndex]; 
        OC = unique(GI); 
        N = zeros(size(OC)); 
        for k =1:numel(OC) 
            N(k) = numel(find(GI == OC(k))); 
        end 
        P = exp(-beta*N); 
        P = P/sum(P);  
        sci = RouletteWheelSelection(P); 
        sc = OC(sci);  
        SCM = find(GI==sc); 
        smi = randi([1 numel(SCM)]); 
        sm = SCM(smi); 
        leader = rep(sm);  
%============================= 
%Roulette chooses a better dominant individual  
%============================= 
function i = RouletteWheelSelection(P) 
        r = rand; 
        C = cumsum(P); 
        i = find(r<=C,1,'first'); 
end

%============================= 
%Use of mutation strategies  
%============================= 
function xnew = Mutate(x,pm,VarMin,VarMax) 
        nVar = numel(x); 
        j = randi([1 nVar]); 
        dx = pm*(VarMax-VarMin); 
        lb = x(j)-dx; 
        if lb<VarMin 
            lb=VarMin; 
        end 
        ub = x(j)+dx; 
        if ub > VarMax 
            ub = VarMax; 
        end 
        xnew = x; 
        xnew(j) = unifrnd(lb,ub); 
end
%============================= 
%Delete an individual from the archive  
%============================= 
function rep = DeleteOneRepMemebr(rep,gamma) 
        GI = [rep.GridIndex]; 
        OC = unique(GI); 
        N = zeros(size(OC)); 
        for k = 1:numel(OC) 
            N(k) = numel(find(GI == OC(k))); 
        end 
        P = exp(gamma*N); 
        P = P/sum(P); 
        sci = RouletteWheelSelection(P); 
        sc = OC(sci); 
        SCM = find(GI == sc); 
        smi = randi([1 numel(SCM)]); 
        sm = SCM(smi); 
        rep(sm) = []; 
end 
function [convergence] = calculateConvergence(particles, globalBest)
% Computing Convergence of Particle Swarms
% particles
% globalBest

n = size(particles, 1);
distances = zeros(n, 1);
for i = 1:n
    distances(i) = norm(particles(i).Best.Cost - globalBest.Cost);
end
convergence = mean(distances);
end

function midValue = calculateMID(pop)
    % pop: Solution matrix with each row representing a solution and each column representing an objective




    [numSolutions, numObjectives] = size(pop);

    % Calculate the minimum and maximum values for each target
    minObjectives = min(pop);
    maxObjectives = max(pop);

    % Calculate the Ideal point for each solution
    idealPoint = minObjectives;

    % Calculate the Euclidean distance from each solution to the Ideal point
    distancesToIdeal = sqrt(sum((pop - repmat(idealPoint, numSolutions, 1)).^2, 2));

    % MID
    midValue = mean(distancesToIdeal);
end


