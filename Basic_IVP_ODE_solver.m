function IVP_ODE_solver(dimension,t0,tend)
    for i=1:dimension
        prompt = "Initial Value:";
        x = input(prompt);
        initial_values(i)=x;
        Cases = input('Enter cases with [ ] around them ');
        %c_1 = input(prompt);
        c(i,:)=Cases;
        prompt = "Give the independent variable of ODE:";
        v_1 = input(prompt,'s');
        if(v_1=='x')
            v(i)=1;
        end
        if(v_1=='y')
            v(i)=2;
        end
        if(v_1=='z')
            v(i)=3;
        end
    end
    t=0;
    %h=[0.1,0.01,0.005,0.0025,0.0001,0.00005,0.000001];
    h=[0.1,0.01,0.005,0.0025,0.0001];
    for j=1:length(h)
        Nh=(tend-t0)/h(j);
        w=zeros(dimension,Nh+1);
        for i=1:dimension
            w(i,1)=initial_values(i);
        end
        %assign the ode's
        ode=zeros(dimension,Nh+1);
        for i=1:dimension
            ode(i,:)=polyval(c(i,:),w(v(i),:));
        end
        %assign the second initial values
        for i=1:dimension
            w(i,2)=w(i,1)+h(j)*(ode(i,1));
        end
        for i=1:dimension
            ode(i,:)=polyval(c(i,:),w(v(i),:));
        end
        t=0;
        for i=2:Nh
            for k=1:dimension
                w(k,i+1)=w(k,i-1)+2*h(j)*(ode(k,i));
                for m=1:dimension
                    ode(m,:)=polyval(c(m,:),w(v(m),:));
                end
            end
            t=t+h(j);
        end
        for i=1:dimension
            last_values(i,j)=w(i,Nh+1);
        end
end