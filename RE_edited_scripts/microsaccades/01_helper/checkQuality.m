function checkQuality(screenCenter, dvaPerPx, new_trial, ms2plot, i, tab)

durStim = 0;

locationdegrees = {315, 135, 225, 45, 270, 90, 180, 0};
%directiondegrees = {45, 225, 135, 315, 90, 270, 0, 180};

location = locationdegrees{tab(i,9)};
direction = tab(i,12); %directiondegrees(tab(i,10));
correct = tab(i,14);

% colors based on whether it is before, during or after stimulus
colors = {[0 1 0], [1 0 0], [0 1 1]};

        if size(new_trial,1) > 105

            x_fulltrial=dvaPerPx*(new_trial(:,2)-screenCenter(1));
            y_fulltrial=dvaPerPx*(new_trial(:,3)-screenCenter(2));
            xpl = filtfilt(fir1(35,0.05),1,x_fulltrial);
            ypl = filtfilt(fir1(35,0.05),1,y_fulltrial);
            plot(xpl, ypl, 'b')
            hold on
            
            [currentNum,~] = size(ms2plot);
            
            for mm=1:currentNum
                if ms2plot(mm,2)<1300 % end of ms before stimulus
                    currColor = colors{1};
                elseif ms2plot(mm,1)<1800 % start of ms before response
                    currColor = colors{2};
                    durStim = 1;
                else
                    currColor = colors{3};
                end
                tp1 = ms2plot(mm,1); %tt1; % ms start (subtract trial start so index aligns)
                tp2 = ms2plot(mm,2); %-tt1; % ms end
                plot(xpl(tp1:tp2),ypl(tp1:tp2), 'Color', currColor)
                hold on
                plot(xpl(tp2),ypl(tp2), 'k.')

            end

            xlim([-1.5 1.5])
            ylim([-1.5 1.5])
            title(sprintf('Trial %i Location %i Direction %i Correct %i', i, location, direction, correct))
            f1 = gcf;
            f1.Position = [42 286 1268 1051];
            axis square
        end
        hold off



end