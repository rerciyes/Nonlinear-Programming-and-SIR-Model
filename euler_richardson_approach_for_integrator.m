tic
%initial point settings
a=0.0;b=10.0;
%h=0.1;
w_for_x(1)=0.0;w_for_y(1)=1.0;
t=0;
h=[0.1,0.01,0.005,0.0025,0.0001,0.00005,0.000001];

for j=1:length(h)
    % # of times loop iteration
    Nh=(b-a)/h(j);Nh_2=(b-a)/(h(j)/2);
    w_for_yh=[];w_for_xh=[];w_for_yh_2=[];w_for_xh_2=[];
    
    t=0;w_for_xh(1)=0.0;w_for_yh(1)=1.0;
    for i=1:Nh
        w_for_xh(i+1)=w_for_xh(i)+h(j)*(2*w_for_yh(i));
        w_for_yh(i+1)=w_for_yh(i)+h(j)*(-2*w_for_xh(i));
        t=t+h(j);
    end
    estimated_xh=w_for_xh(Nh+1);
    estimated_yh=w_for_yh(Nh+1);
    
    t=0;w_for_xh_2(1)=0.0;w_for_yh_2(1)=1.0;
    for i=1:Nh_2
        w_for_xh_2(i+1)=w_for_xh_2(i)+(h(j)/2)*(2*w_for_yh_2(i));
        w_for_yh_2(i+1)=w_for_yh_2(i)+(h(j)/2)*(-2*w_for_xh_2(i));
        t=t+(h(j)/2);
    end
    estimated_xh_2=w_for_xh_2(Nh_2+1);
    estimated_yh_2=w_for_yh_2(Nh_2+1);
    
    last_value_of_y(j)=2*estimated_yh_2-estimated_yh;
    last_value_of_x(j)=2*estimated_xh_2-estimated_xh;
end
h=transpose(h);
y_1=transpose(last_value_of_y);
x_1=transpose(last_value_of_x);
%outcomes=table(h,y_1,x_1)
solution_of_error=log(sqrt(((cos(20))-last_value_of_y).^2+((sin(20))-last_value_of_x).^2));
new_h=log(transpose(h));
%loglog(new_h,solution_of_error)
plot(new_h,solution_of_error)
toc