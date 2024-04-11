%% calculate degree
function MS_new=calculate_direction_RE(MS)
%first using component
x=MS(:,4);
y=MS(:,5);
len=length(x);
de2=zeros(len,1);

for i = 1:len
    if x(i) == 0
        if y(i) > 0
            de2(i) = 90;
        elseif y(i) < 0
            de2(i) = 270;
        else
            de2(i) = NaN; % Undefined direction
        end
    else
        angle_rad = atan2d(y(i), x(i)); % atan2d returns angles in degrees
        if angle_rad < 0
            de2(i) = angle_rad + 360; % Convert negative angles to positive
        else
            de2(i) = angle_rad;
        end
    end
end


%using amplitude
x=MS(:,6);
y=MS(:,7);
len=length(x);
de=zeros(len,1);


for i = 1:len
    if x(i) == 0
        if y(i) > 0
            de(i) = 90;
        elseif y(i) < 0
            de(i) = 270;
        else
            de(i) = NaN; % Undefined direction
        end
    else
        angle_rad = atan2d(y(i), x(i)); % atan2d returns angles in degrees
        if angle_rad < 0
            de(i) = angle_rad + 360; % Convert negative angles to positive
        else
            de(i) = angle_rad;
        end
    end
end

MS_new=[MS,de,de2];
end