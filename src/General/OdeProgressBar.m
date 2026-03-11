function status = OdeProgressBar(t,~,flag,bar_scale,update_frequency)
persistent t0 t_end start_wct wct_last counter
if isempty(counter)
    counter = 1;
end
if isempty(flag)
    if((datetime("now")-wct_last)>seconds(1))
        figure(95);
        if bar_scale == "log"
            perc=(log10(t(end))-log10(t0))/(log10(t_end)-log10(t0));
        elseif bar_scale == "lin"
            perc=(t(end)-t0)/(t_end-t0);
        end
        elapsed_time_seconds = seconds(datetime("now")-start_wct);
        area([t(end) t_end-t(end);t(end) t_end-t(end)]);
        set(findobj('Tag','perc'),'String',num2str(perc*100) + " %");
        set(findobj('Tag','eltime'),'String',SecondsToString(elapsed_time_seconds));
        set(findobj('Tag','esttime'),'String',SecondsToString(elapsed_time_seconds*((1-perc)/perc)));
        wct_last=datetime("now");
    end
    if mod(counter,update_frequency) == 0
        clear ComputeOmegaKozhevnikov
        clear ComputeOmegaTownsend
        clear ComputeOmegaParent
        clear ComputeOmegaParentConst
        clear ComputeOmegaBenilov
    end
else
    if flag == "init"
        t0 = t(1);
        t_end = t(end);
        start_wct=datetime("now");
        wct_last=datetime("now");
        fig = figure(95);
        set(fig,'Position',[4,40,100,545],"Resize","off","ToolBar","none","DockControls","off","MenuBar","none");
        ax = axes(fig);
        ax.Units = "pixels";
        ax.Position = [50, 135, 48, 360];
        axis([1,2,t0,t_end]);
        if bar_scale == "log"
            set(gca,'XTickLabel',[],'YScale','Log','NextPlot','replacechildren');
        elseif bar_scale == "lin"
            set(gca,'XTickLabel',[],'NextPlot','replacechildren');
        end
        ylabel('Simulation Progress - Time (s)');
        area([t(1) t(end);t(1) t(end)]);
        uicontrol('Style', 'pushbutton', 'String', 'Abort','Position', [7 515 103 30], 'Callback', 'close(gcf)')
        uicontrol('Style', 'text', 'String', 'Progress','Position', [7 105 90 15])
        uicontrol('Style', 'text', 'Tag', 'perc', 'String', "0 %",'Position', [7 85 90 15])
        uicontrol('Style', 'text', 'String', 'Elapsed Time','Position', [7 65 90 15])
        uicontrol('Style', 'text', 'Tag', 'eltime', 'String', "0s",'Position', [7 45 90 15])
        uicontrol('Style', 'text', 'String', 'Time Remaining','Position', [7 25 90 15])
        uicontrol('Style', 'text', 'Tag', 'esttime', 'String', num2str(inf),'Position', [7 5 90 15])
        if mod(counter,update_frequency) == 0
            clear ComputeOmegaKozhevnikov
            clear ComputeOmegaTownsend
            clear ComputeOmegaParent
            clear ComputeOmegaParentConst
        end
        pause(0.1);
    elseif flag == "done"
        if(ishandle(95))
            close(95);
        end
    end
end
counter = counter + 1;
status = 0;
drawnow;
end
 