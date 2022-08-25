%firstly put the general inputs
%func = input('Enter a model option: ','s');
option="Sinusoidal Model";
%p = input('Enter the parameters with [ ] around them ');
%initial_values = input('Enter the initials with [ ] around them ');
initial_values=[0 1];

%generate data for just an trial
tspan= 0:0.01:10;
[t,y] = ode45(@vdp_Sinusoidal,tspan,[0;1]);
y=transpose(y);
index = randi([1 length(t)],1,7);

beta=[1.5 2.3]
i=1;
while i<15
    
    %calculate numerical outcome with beta
    %p=[beta];
    %numerical_values=IVP_ODE_solver(0,100,option,p,initial_values);
    p=[beta];
    numerical_values=IVP_ODE_solver(0,10,option,p,initial_values);
    
    %calculate residue
    for j=1:7
        r(j)=y(1,index(j))-numerical_values(1,index(j));
    end
    rss(i)=sum(r.^2);
    r=transpose(r);
    while (size(r,2)==7)
        r=transpose(r);
    end
    
    %calculate jacobian
    h=0.001;
    p_1_1=[beta(1)+h beta(2)];
    p_1_2=[beta(1)-h beta(2)];
    x1_1=IVP_ODE_solver(0,10,option,p_1_1,initial_values);
    x1_2=IVP_ODE_solver(0,10,option,p_1_2,initial_values);
    for k=1:7
        diff1(k)=(-x1_1(1,index(k))+x1_2(1,index(k)))./(2*h);%symmetric difference
    end
    
    h=0.001;
    p_2_1=[beta(1) beta(2)+h];
    p_2_2=[beta(1) beta(2)-h];
    x2_1=IVP_ODE_solver(0,10,option,p_2_1,initial_values);
    x2_2=IVP_ODE_solver(0,10,option,p_2_2,initial_values);
    for k=1:7
        diff2(k)=(-x2_1(1,index(k))+x2_2(1,index(k)))./(2*h);%symmetric difference
    end
    
    change_of_r_1=transpose(diff1);
    change_of_r_2=transpose(diff2);
    J=[change_of_r_1 change_of_r_2];% speed of change in a sense of derivative
    
    %put iteration here
    beta=transpose(beta);
    beta=beta-pinv(J)*r;
    beta=transpose(beta);
    i=i+1;
end
beta
rss

%%
function value=rhs_new(option,p,dimension,j,w)
if(option=="SIR Model")
    if(dimension==1)
        value=-(p(1)/p(3))*w(1,j)*w(2,j);
    elseif(dimension==2)
        value=(p(1)/p(3))*w(1,j)*w(2,j)-p(2)*w(2,j);
    else
        value=p(2)*w(2,j);
    end
end
if(option=="Sinusoidal Model")
    if(dimension==1)
        value=p(1)*w(2,j);
    else
        value=-p(2)*w(1,j);
    end
end

end

function [numerical_values]=IVP_ODE_solver(t0,tend,option,parameters,initials)
    if(option=="SIR Model")
        dimension=3;
        %p = input('Enter the parameters with [ ] around them ');
        p=parameters;
    end
    if(option=="Sinusoidal Model")
        dimension=2;
        %p = input('Enter the parameters with [ ] around them ');
        p=parameters;
    end
    %{
    for i=1:dimension
        prompt = "Initial Value:";
        x = input(prompt);
        initial_values(i)=x;
    end
    %}
    initial_values=initials;
    t=0;
    %h=[0.1,0.01,0.005,0.0025,0.0001,0.00005,0.000001];
    h=[0.01];
    for j=1:length(h)
        Nh=(tend-t0)/h(j);
        w=zeros(dimension,Nh+1);
        for i=1:dimension
            w(i,1)=initial_values(i);
        end
        %assign the second initial values
        for i=1:dimension
            w(i,2)=w(i,1)+h(j)*(rhs_new(option,p,i,1,w));
        end
        t=0;
        for i=2:Nh
            for k=1:dimension
                w(k,i+1)=w(k,i-1)+2*h(j)*(rhs_new(option,p,k,i,w));
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
    
    %{
    figure
    t=t0:h(1):tend;
    y=w(1,:);
    plot(t,y,':')
    
    hold on 
    y2 = w(1,:);
    plot(t,y2,'--ro')
    y3=w(2,:);
    scatter(t,y3) 
    hold off
    %}
    numerical_values=w;
end

function dydt = vdp_SIR(t,y)
dydt = [(-0.4/1000)*y(1)*y(2); 0.4/1000*y(1)*y(2)-0.04*y(2);0.04*y(2)];
end
function dydt = vdp_Sinusoidal(t,y)
dydt = [2*y(2); -2*y(1)];
end