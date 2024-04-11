%% calculate degree
function MS_new=calculate_direction(MS)
%first using component
x=MS(:,4);
y=MS(:,5);
len=length(x);
de2=zeros(len,1);
for i = 1 : len
    if x(i)>0
        de2(i)=atand(y(i)/x(i));
    else
        if y(i) > 0 
            de2(i)=atand(y(i)/x(i))+180;
        else
            de2(i)=atand(y(i)/x(i))-180;
        end
    end
    
end

%using amplitude
x=MS(:,6);
y=MS(:,7);
len=length(x);
de=zeros(len,1);
for i = 1 : len
    if x(i)>0
        de(i)=atand(y(i)/x(i));
    else
        if y(i) > 0 
            de(i)=atand(y(i)/x(i))+180;
        else
            de(i)=atand(y(i)/x(i))-180;
        end
    end
end
MS_new=[MS,de,de2];
end