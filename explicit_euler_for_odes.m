%initial point settings
a=0.0;b=10.0;
%h=0.1;
w_for_x(1)=0.0;w_for_y(1)=1.0;
t=0;
h=[0.1,0.01,0.005,0.0025,0.0001,0.00005,0.000001];
C=cell(2,length(h));
for j=1:length(h)
    % # of times loop iteration
    N=(b-a)/h(j);
    w_for_y=[];w_for_x=[];
    t=0;w_for_x(1)=0.0;w_for_y(1)=1.0;
    for i=1:N
        w_for_x(i+1)=w_for_x(i)+h(j)*(2*w_for_y(i));
        w_for_y(i+1)=w_for_y(i)+h(j)*(-2*w_for_x(i));
        t=t+h(j);
    end
    last_value_of_x(j)=w_for_x(N+1);
    last_value_of_y(j)=w_for_y(N+1);
    
    C(1,j)={[w_for_x]};
    C(2,j)={[w_for_y]};
end
h=transpose(h);
y_1=transpose(last_value_of_y);
x_1=transpose(last_value_of_x);
outcomes=table(h,y_1,x_1)
solution_of_error=log(sqrt(((cos(20))-last_value_of_y).^2+((sin(20))-last_value_of_x).^2));
new_h=log(transpose(h));
%loglog(new_h,solution_of_error)
plot(new_h,solution_of_error)

figure
t=a:h(7):b;
y = w_for_y;
plot(t,y,':')

hold on 
y2 = w_for_x;
plot(t,y2,'--ro')
hold off

tic
%implementing the bisection algorithm with discrete function (estimating pi)
a=find(t==1.5);b=find(t==2);
tol=1;N0=20;
i=0;FA=w_for_x(a);
while(i<N0)
    p=ceil((b+a)/2);
    FP=w_for_x(p);
    if(FP==0 ||(b-a)/2 < tol)
        value=p;
        break
    end
    i=i+1;
    if FA*FP >0
        a=p;FA=FP;
    else
        b=p;
    end 
end
toc
vpa(t(p))
w_for_x(p)