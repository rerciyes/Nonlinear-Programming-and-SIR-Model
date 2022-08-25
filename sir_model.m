tic
IVP_ODE_solver(3,0,40);
toc
%%
function [vector]=rhs(c,C,w)
dimension=size(w,1);
N=size(w,2);
ode=zeros(dimension,N+1);
for i=1:dimension
    f=str2func(C{i});
    for j=1:N
        c1=cell2mat(c{i});
        b=[];
        for k=1:length(c1)
            b(k)=w(c1(k),j);
        end
        cell1=num2cell(b);
        ode(i,j)=f(cell1{:});
    end
end
vector=ode;
end
function IVP_ODE_solver(dimension,t0,tend)
    C = {};c = {};
    for i=1:dimension
        prompt = "Initial Value:";
        x = input(prompt);
        initial_values(i)=x;
        func = input('Enter a function: ','s');
        C{i}=func;
        %Which independent variables are used in equation of ODE
        Cases = input('Enter cases with [ ] around them ');
        c{i}=num2cell(Cases);
    end
    t=0;
    %h=[0.1,0.01,0.005,0.0025,0.0001,0.00005,0.000001];
    h=[0.1];
    for j=1:length(h)
        Nh=(tend-t0)/h(j);
        w=zeros(dimension,Nh+1);
        for i=1:dimension
            w(i,1)=initial_values(i);
        end
        %update the ode's
        ode=rhs(c,C,w);
        %assign the second initial values
        for i=1:dimension
            w(i,2)=w(i,1)+h(j)*(ode(i,1));
        end
        %update the ode's
        ode=rhs(c,C,w);
        t=0;
        for i=2:Nh
            for k=1:dimension
                w(k,i+1)=w(k,i-1)+2*h(j)*(ode(k,i));
                %update the ode's
                ode=rhs(c,C,w);
            end
            t=t+h(j);
        end
        for i=1:dimension
            last_values(i,j)=w(i,Nh+1);
        end
    end
    h=transpose(h);
    y_1=transpose(last_values(2,:));
    x_1=transpose(last_values(1,:));
    %outcomes=table(h,y_1,x_1)
    
    figure
    t=t0:h(1):tend;
    y = w(2,:);
    plot(t,y,':')
    
    hold on 
    y2 = w(1,:);
    plot(t,y2,'--ro')
    y3=w(3,:);
    scatter(t,y3) 
    hold off
end